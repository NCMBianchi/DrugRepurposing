"""
NETWORK MODEL MODULE: RUNS UNSUPERVISED LEARNING (XGBoost)
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,platform,datetime,logging,builtins,time,multiprocessing

from drugapp.unique import unique_elements
from drugapp.filepaths import initialise_base_directories
from flask import request, jsonify
from flask import current_app as app

global base_directories
global disease_directories

@app.route('/run', methods=['POST'])
def networkmodel_service_run():
    """
    Flask route to run Network Model service for distributed architecture
    """
    # Extract payload from the request
    payload = request.json

    # Extract parameters from payload with defaults
    input_seed = payload.get('input_seed', 'MONDO:0007739')
    today = payload.get('today', datetime.date.today().strftime('%Y-%m-%d'))
    date_str = payload.get('date_str', today)
    nodes = payload.get('nodes', [])
    edges = payload.get('edges', [])
    drug_nodes = payload.get('drug_nodes', [])
    drug_edges = payload.get('drug_edges', [])
    negs_toggle = payload.get('negs_toggle', False)

    # Use initialise_base_directories to set up paths
    base_directories = initialise_base_directories(date_str)
    today_directory = base_directories['today_directory']
    
    # Helper function to manage global-like variables
    def setup_networkmodel_globals(today_directory, date_str, input_seed):
        global today_directory, date_str, input_seed, input_file_path, disease_name_label
        
        # Set up input_file_path
        input_file_path = os.path.join(today_directory, 'inputs.txt')
        
        # Set disease_name_label 
        disease_name_label = payload.get('disease_name_label', 'Huntington disease')
        
        return today_directory, date_str, input_seed, input_file_path, disease_name_label

    # Setup globals
    today_directory, date_str, input_seed, input_file_path, disease_name_label = setup_networkmodel_globals(
        today_directory, date_str, input_seed
    )

    try:
        # Call the appropriate network model function based on negative samples toggle
        if negs_toggle:
            network_edges, network_nodes, ranked_drugs, plots = run_network_model_with_NS(
                input_seed, today, 
                run_jobs=payload.get('run_jobs', None),
                run_depth=payload.get('run_depth', None),
                run_seed=payload.get('run_seed', None),
                prob_input=payload.get('prob_input', None),
                clust_input=payload.get('clust_input', None),
                emb_toggle=payload.get('emb_toggle', 1),
                param_toggle=payload.get('param_toggle', 1),
                ou_toggle=payload.get('ou_toggle', 1)
            )
        else:
            network_edges, network_nodes, ranked_drugs, plots = run_network_model(
                input_seed, today, 
                run_jobs=payload.get('run_jobs', None),
                run_depth=payload.get('run_depth', None),
                run_seed=payload.get('run_seed', None),
                prob_input=payload.get('prob_input', None),
                clust_input=payload.get('clust_input', None),
                emb_toggle=payload.get('emb_toggle', 1),
                param_toggle=payload.get('param_toggle', 1),
                ou_toggle=payload.get('ou_toggle', 1)
            )

        # Prepare the response
        response = {
            'network_edges': network_edges,
            'network_nodes': network_nodes,
            'ranked_drugs': ranked_drugs,
            'plots': plots,
            # Re-include payload data to pass to next service
            **{k: v for k, v in payload.items() if k not in ['nodes', 'edges', 'drug_nodes', 'drug_edges']}
        }

        return jsonify(response)

    except Exception as e:
        logging.error(f"Error in Network Model service: {str(e)}")
        return jsonify({
            'error': str(e),
            'status': 'failed'
        }), 500

def get_network(input_nodes,input_edges,exclude=None):
    '''
    This function builds a network object from lists of nodes and edges in a
    given format. Examples:
    
    NODES = [{'id': 'chembl:CHEMBL348436', 'label': 'CIRSIMARITIN'},
    {'id': 'chembl:CHEMBL1503190', 'label': 'CHEMBL1503190'},...]
    
    EDGES = [[{'id': 'MONDO:0007739', 'label': 'Huntington disease'},
    {'label': 'biolink:has_mode_of_inheritance'},
    {'id': 'HP:0000006', 'label': 'Autosomal dominant inheritance'},
    {'notes': 'first_layer'}],[{'id': 'MONDO:0007739', 'label': 'Huntington disease'},
    {'label': 'biolink:subclass_of'},
    {'id': 'MONDO:0005395', 'label': 'movement disorder'},
    {'notes': 'first_layer'}]...]
    
    It can remove a given edge type in order to build a partial network: it should be
    in the URI format (e.g. 'dgidb:'), and the formula checks if it's in the entire
    relations label.

    :param input_nodes: list of nodes.
    :param input_edges: list of edges.
    :param exclude: edge type to exclude (default = None).
    :return: networkX object.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    G = nx.DiGraph()

    global node_type_dict

    node_type_colours = {
        'disease': '#00008b',  # darkblue
        'gene': '#add8e6',     # lightblue
        'phenotype': '#0000ff',# blue
        'drug': '#ff0000',     # red
        'else': '#ffa500'      # orange
    }

    edge_type_colours = {
        'biolink': '#0000ff',  # blue
        'dgidb': '#ff0000',    # red
        'SMILES': '#ff4500',   # lightred
        'xgboost': '#00ff00'   # green 
    }
    
    # add nodes
    for node in input_nodes:
        if node is None:
            logging.warning("Encountered NoneType node, skipping...")
            continue
        node_id = node.get('id')
        label = node.get('label')
        if node_id is None or label is None:
            logging.warning(f"Encountered node with missing id or label: {node}, skipping...")
            continue
        node_type = 'else'
        
        for type_key, prefixes in node_type_dict.items():
            if any(node_id.startswith(prefix + ':') for prefix in prefixes):
                node_type = type_key
                break
                
        G.add_node(node_id, label=label, node_type=node_type, colour=node_type_colours[node_type])
    
    # add edges
    for edge in input_edges:
        subj_id = edge[0]['id']
        obj_id = edge[2]['id']
        rel_label = edge[1]['label']
        notes = edge[3]['notes']
    
        if subj_id is None or obj_id is None or rel_label is None:
            logging.warning(f"Encountered edge with missing id or label: {edge}, skipping...")
            continue
        
        if exclude is not None and exclude in rel_label:
            continue  # skip this edge if exclusion criteria match
        
        rel_type = rel_label.split(':')[0]
        colour = edge_type_colours.get(rel_type, '#000000')  # default to black if not found
        G.add_edge(subj_id, obj_id, label=rel_label, notes=notes, colour=colour)

    return G


def plot_network(network):
    pos = nx.spring_layout(network)
    ordered_node_colors = [network.nodes[node].get('colour', '#000000') for node in network.nodes()]
    node_sizes = [30 if node == 'MONDO:0007739' else 10 for node in network.nodes()]
    plt.figure(figsize=(3000/300, 2000/300))  # 300dpi
    nx.draw_networkx_nodes(network, pos, node_size=node_sizes, node_color=ordered_node_colors, alpha=0.6)
    for u, v, data in network.edges(data=True):
        edge_color = data['colour']
        if edge_color == '#00ff00':  # Adjust the properties for the green edges
            nx.draw_networkx_edges(network, pos, edgelist=[(u, v)], width=0.8, alpha=0.5, edge_color=edge_color)
        else:  # Default properties for other edges
            nx.draw_networkx_edges(network, pos, edgelist=[(u, v)], width=0.2, alpha=0.2, edge_color=edge_color)
    plt.axis('off')
    plot_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_full_network_plot.png')
    plt.savefig(plot_path, format="png", dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)


def get_embeddings(input_network,emb_t,node_type_select=None):
    '''
    This function builds embeddings (vectors) for node informations -also considering
    the edges they are involved in- from a network object via Node2Vec.
    It can filter to a selected type of nodes

    :param input_network: networkX object.
    :param emb_t: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check.
    :param node_type_select: node type to focus on (default = None): used within run_network_model(),
        not by the user.
    :return: embedding pickle object.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    global network_directory, node_type_dict
    valid_node_types = set(node_type_dict.keys())
    
    if node_type_select and node_type_select not in valid_node_types:
        raise ValueError(f"Invalid node_type_select: {node_type_select}. Valid options are: {valid_node_types}")

    # PKL file path
    if node_type_select is not None:
        embedding_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_embedding_dict_{node_type_select}.pkl')
    else:
        embedding_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_embedding_dict_full.pkl')

    # check if they already exist, if toggled
    if emb_t == 1 and os.path.exists(embedding_path):
        logging.info(f"Loading embeddings from existing file: {embedding_path}")
        
        with open(embedding_path, "rb") as file:
            output_embeddings = pickle.load(file)

    # otherwise compute
    else:   
        # (eventually) filter nodes
        if node_type_select:
            nodes_to_include = [node for node in input_network.nodes() if input_network.nodes[node].get('node_type') == node_type_select]
            subgraph = input_network.subgraph(nodes_to_include)
        else:
            subgraph = input_network

        # generate node embeddings using Node2Vec
        node2vec = Node2Vec(subgraph, dimensions=64, walk_length=30, num_walks=200, workers=2)
        model = node2vec.fit(window=10, min_count=1, batch_words=4)
        output_embeddings = {node: model.wv[node] for node in subgraph.nodes()}

        # save the embeddings dictionary to a PKL file
        with open(embedding_path, "wb") as file:
            pickle.dump(output_embeddings, file)
        logging.info("PKL files saved in network directory.")

    return output_embeddings


def get_embeddings_with_NS(input_network,emb_t,node_type_select=None):
    '''
    This function builds embeddings (vectors) for node informations -also considering
    the edges they are involved in- from a network object via Node2Vec.
    It can filter to a selected type of nodes

    :param input_network: networkX object.
    :param emb_t: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check.
    :param node_type_select: node type to focus on (default = None): used within run_network_model(),
        not by the user.
    :return: embedding pickle object.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    global network_with_NS_directory, node_type_dict
    valid_node_types = set(node_type_dict.keys())
    
    if node_type_select and node_type_select not in valid_node_types:
        raise ValueError(f"Invalid node_type_select: {node_type_select}. Valid options are: {valid_node_types}")

    # PKL file path
    if node_type_select is not None:
        embedding_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_embedding_dict_NS_{node_type_select}.pkl')
    else:
        embedding_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_embedding_dict_NS_full.pkl')

    # check if they already exist, if toggled
    if emb_t == 1 and os.path.exists(embedding_path):
        logging.info(f"Loading embeddings from existing file: {embedding_path}")
        
        with open(embedding_path, "rb") as file:
            output_embeddings = pickle.load(file)

    # otherwise compute
    else:   
        # (eventually) filter nodes
        if node_type_select:
            nodes_to_include = [node for node in input_network.nodes() if input_network.nodes[node].get('node_type') == node_type_select]
            subgraph = input_network.subgraph(nodes_to_include)
        else:
            subgraph = input_network

        # generate node embeddings using Node2Vec
        node2vec = Node2Vec(subgraph, dimensions=64, walk_length=30, num_walks=200, workers=2)
        model = node2vec.fit(window=10, min_count=1, batch_words=4)
        output_embeddings = {node: model.wv[node] for node in subgraph.nodes()}

        # save the embeddings dictionary to a PKL file
        with open(embedding_path, "wb") as file:
            pickle.dump(output_embeddings, file)
        logging.info("PKL files saved in network directory.")

    return output_embeddings


def is_interaction_present(gene_id, drug_id):
    return (
        ((DGIdb_edges_df['subject_id'] == gene_id) & (DGIdb_edges_df['object_id'] == drug_id))
        | ((DGIdb_edges_df['subject_id'] == drug_id) & (DGIdb_edges_df['object_id'] == gene_id))
    ).any()


def is_interaction_present_with_NS(gene_id, drug_id, edges_df):
    # check if the interaction is present and determine the type of interaction
    positive_interaction = (
        ((edges_df['subject_id'] == gene_id) & (edges_df['object_id'] == drug_id) & (edges_df['relation'] != 'biolink:valid_negative_association'))
        | ((edges_df['subject_id'] == drug_id) & (edges_df['object_id'] == gene_id) & (edges_df['relation'] != 'biolink:valid_negative_association'))
    ).any()
    
    valid_negative_interaction = (
        ((edges_df['subject_id'] == gene_id) & (edges_df['object_id'] == drug_id) & (edges_df['relation'] == 'biolink:valid_negative_association'))
        | ((edges_df['subject_id'] == drug_id) & (edges_df['object_id'] == gene_id) & (edges_df['relation'] == 'biolink:valid_negative_association'))
    ).any()
    
    if positive_interaction:
        return 1
    elif valid_negative_interaction:
        return 0
    else:
        return -1


def fuse_embeddings(gene_embedding_dict, drug_embedding_dict, alldrug_embedding_dict, DGIdb_edges_list, nodes_list, alldrug_nodes_list,emb_t):
    '''
    This function fuses the embeddings for the genes and drugs that have a known link,
    and adds additional drug-gene pairs with random non-known links for prediction.
    Structure of embedding dictionaries: {key: gene/drug, value: embedding}

    :param gene_embedding_dict: all gene embeddings associated to 'input_seed'.
    :param drug_embedding_dict: all drug embeddings associated to 'input_seed'.
    :param alldrug_embedding_dict: all drug embeddings with entire 'DGIdb'.
    :param DGIdb_edges_list: list of gene-to-drug edges from 'run_dgidb()'.
    :param nodes_list: list of nodes associated to 'input_seed'.
    :param alldrug_nodes_list: list of nodes with entire 'DGIdb'.
    :param emb_t: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check.
    :return: dataframes used to train the XGboost model, and subsequently test it.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    embedding_start_time = time.time()

    logging.info('Embeddings are being fused, be patient.')

    # CSV file paths
    train_df_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_training_df.csv')
    predict_df_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_prediction_df.csv')
    
    # check if already present, if toggled:
    if emb_t == 1 and os.path.exists(train_df_path) and os.path.exists(predict_df_path):
        logging.info(f"Loading training and prediction dataframes from existing files: {train_df_path} and {predict_df_path}")
        train_df = pd.read_csv(train_df_path)
        predict_df = pd.read_csv(predict_df_path)

        # Convert fused_embedding column from string to numpy array
        train_df['fused_embedding'] = train_df['fused_embedding'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
        predict_df['fused_embedding'] = predict_df['fused_embedding'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))

    # otherwise fuse
    else:
        array_drugs = np.array(list(drug_embedding_dict.keys()))
        array_alldrugs = np.array(list(alldrug_embedding_dict.keys()))
        array_genes = np.array(list(gene_embedding_dict.keys()))

        DGIdb_edges_df = pd.DataFrame({
            'subject_id': [edge[0]['id'] for edge in DGIdb_edges_list],
            'object_id': [edge[2]['id'] for edge in DGIdb_edges_list],
            'relation': [edge[1]['label'] for edge in DGIdb_edges_list]
        })

        # create gene-to-drug interaction dataframe for training
        train_rows = []

        for gene in array_genes:
            gene_emb = gene_embedding_dict[gene]
            for drug in array_drugs:
                fused_emb = np.multiply(gene_emb, drug_embedding_dict[drug])
                class_label = 1 if is_interaction_present(gene, drug) else 0
                train_rows.append({
                    'gene': gene, 'drug': drug, 'fused_embedding': fused_emb, 'class': class_label
                })
        train_df = pd.DataFrame(train_rows)

        # create interaction dataframe for prediction with additional drugs
        predict_rows = []

        for gene in array_genes:
            gene_emb = gene_embedding_dict[gene]
            for drug in array_alldrugs:
                fused_emb = np.multiply(gene_emb, alldrug_embedding_dict[drug])
                predict_rows.append({
                    'gene': gene, 'drug': drug, 'fused_embedding': fused_emb, 'class': np.nan
                })
        predict_df = pd.DataFrame(predict_rows)

        # save the training and prediction dataframes to CSV files
        train_df.to_csv(train_df_path, index=False)
        predict_df.to_csv(predict_df_path, index=False)

    embedding_end_time = time.time()
    duration = embedding_end_time - embedding_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"Embedding fusion finished in {minutes} minutes and {seconds} seconds.")

    return train_df, predict_df


def fuse_embeddings_with_NS(gene_embedding_dict, drug_embedding_dict, alldrug_embedding_dict, DGIdb_edges_list, nodes_list, alldrug_nodes_list, emb_t):
    '''
    This function fuses the embeddings for the genes and drugs that have a known link,
    and adds additional drug-gene pairs with random non-known links for prediction.
    Structure of embedding dictionaries: {key: gene/drug, value: embedding}

    :param gene_embedding_dict: all gene embeddings associated with 'input_seed'.
    :param drug_embedding_dict: all drug embeddings associated with 'input_seed'.
    :param alldrug_embedding_dict: all drug embeddings with entire 'DGIdb'.
    :param DGIdb_edges_list: list of gene-to-drug edges from 'run_dgidb()'.
    :param nodes_list: list of nodes associated with 'input_seed'.
    :param alldrug_nodes_list: list of nodes with entire 'DGIdb'.
    :param emb_t: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check.
    :return: dataframes used to train the XGboost model, and subsequently test it.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    embedding_start_time = time.time()

    logging.info('Embeddings are being fused, be patient.')

    # CSV file paths
    train_df_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_training_df_with_NS.csv')
    predict_df_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_prediction_df_with_NS.csv')
    
    # check if already present, if toggled:
    if emb_t == 1 and os.path.exists(train_df_path) and os.path.exists(predict_df_path):
        logging.info(f"Loading training and prediction dataframes from existing files: {train_df_path} and {predict_df_path}")
        train_df = pd.read_csv(train_df_path)
        predict_df = pd.read_csv(predict_df_path)

        # Convert fused_embedding column from string to numpy array
        train_df['fused_embedding'] = train_df['fused_embedding'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
        predict_df['fused_embedding'] = predict_df['fused_embedding'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))

    # otherwise fuse
    else:
        array_drugs = np.array(list(drug_embedding_dict.keys()))
        array_alldrugs = np.array(list(alldrug_embedding_dict.keys()))
        array_genes = np.array(list(gene_embedding_dict.keys()))

        DGIdb_edges_df = pd.DataFrame({
            'subject_id': [edge[0]['id'] for edge in DGIdb_edges_list],
            'object_id': [edge[2]['id'] for edge in DGIdb_edges_list],
            'relation': [edge[1]['label'] for edge in DGIdb_edges_list]
        })

        # create gene-to-drug interaction dataframe for training
        train_rows = []

        for gene in array_genes:
            gene_emb = gene_embedding_dict[gene]
            for drug in array_drugs:
                fused_emb = np.multiply(gene_emb, drug_embedding_dict[drug])
                class_label = is_interaction_present_with_NS(gene, drug, DGIdb_edges_df)
                train_rows.append({
                    'gene': gene, 'drug': drug, 'fused_embedding': fused_emb, 'class': class_label
                })
        train_df = pd.DataFrame(train_rows)

        # create interaction dataframe for prediction with additional drugs
        predict_rows = []

        for gene in array_genes:
            gene_emb = gene_embedding_dict[gene]
            for drug in array_alldrugs:
                fused_emb = np.multiply(gene_emb, alldrug_embedding_dict[drug])
                predict_rows.append({
                    'gene': gene, 'drug': drug, 'fused_embedding': fused_emb, 'class': np.nan
                })
        predict_df = pd.DataFrame(predict_rows)

        # save the training and prediction dataframes to CSV files
        train_df.to_csv(train_df_path, index=False)
        predict_df.to_csv(predict_df_path, index=False)

    embedding_end_time = time.time()
    duration = embedding_end_time - embedding_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"Embedding fusion finished in {minutes} minutes and {seconds} seconds.")

    return train_df, predict_df


def ml_prediction(train_df,predict_df,param_t,input_jobs=None,depth=None,input_seed=None):
    '''
    This function builds a supervised learning model using as training data 
    the network of interactions between biological associations (via 'Monarch.py'),
    gene-to-drug associations (via 'DGIdb.py'), and drug-to-drug similarity (via
    'drug_similarity.py') –to then predict new interactions.
    
    :param train_df: dataframe of (gene, drug, embedding, class) used for training.
    :param predict_df: dataframe of (gene, drug, embedding, class) used for prediction.
    :param param_t: 1 to look for existing optimised parameter file to load instead of performing
        a new parameter optimisation step, 0 to ignore the parameter check.
    :param input_jobs: how many CPU cores to use ('-1' uses all available, just for repetition).
    :param depth: based on this, different sets of parameters, splits, repeats and iterations
        are performed (default = 'full', factor of 'full', 'light' and 'ultralight', just for repetition).
    :param random_seed: seed to have consistent results (default = None, random, just for repetition).
    :return: dataframe of newly predicted interactions.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    global nodes
    global drug_nodes

    if input_seed == 'random' or input_seed is None:
        input_seed = np.random.randint(1, 1e6)

        ## LOGGING: random seed
        with open(input_file_path, 'w') as file:
            file.write(f"Random seed generated for the run: {input_seed}\n")
            file.write("\n\n")

    if depth is None:
        depth = 'light'

    if input_jobs is None:
        input_jobs = num_cores // 2

    with open(input_file_path, 'a') as file:
        file.write(f"The XGBoost prediction is run on ({input_jobs}) cores with seed ({input_seed}).\n")
        file.write(f"A ({depth}) parameter optimisation step is performed.\n")

    logging.info('ML model is being trained, be patient.')
    ml_train_start_time = time.time()

    # MODEL TRAINING
    # X = embeddings (turn list of embeddings into columns)
    emb_col = train_df['fused_embedding']
    X = pd.DataFrame(emb_col.tolist())
    X = X.astype(float)

    # y = labels
    y = train_df['class'].astype(int)
    unique_classes = np.unique(y)

    # define parameters to be tuned, and model
    if depth == 'full':
        parameters = {
            'min_child_weight': [2, 3, 5, 8, 13, 20, 30],
            'gamma': [0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.5, 4, 6],
            'reg_alpha': [0, 0.5, 1, 3, 5, 10],
            'reg_lambda': [0, 0.5, 1, 3, 5, 10],
            'subsample': uniform(0.5, 0.5),  # = uniform distribution on (0.5, 1)
            'colsample_bytree': uniform(0.2, 0.8),  # = uniform distribution on (0.2, 1)
            'max_depth': [4, 6, 8, 10, 12, 14, 16],
            'n_estimators': [35, 45, 50, 70, 80, 90, 100],
            'learning_rate': uniform(0, 0.3),
        }
        n_s = 10  # number of splits
        n_r = 5  # number of repeats
        n_i = 20  # number of iterations
    elif depth == 'light':
        parameters = {
            'min_child_weight': [3, 8, 20],
            'gamma': [0.5, 1.2, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20, 30],
            'learning_rate': [0.1],
        }
        n_s = 3  # number of splits
        n_r = 2  # number of repeats
        n_i = 5  # number of iterations
    elif depth == 'ultralight':
        parameters = {
            'min_child_weight': [5, 20],
            'gamma': [0.5, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20],
            'learning_rate': [0.1],
        }
        n_s = 2  # number of splits
        n_r = 1  # number of repeats
        n_i = 5  # number of iterations

    params_file_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_hyperparameters.txt')

    # check if the parameters file is already present: if so and if toggled, run with unique parameters
    if param_t==1 and os.path.exists(params_file_path):
        with open(params_file_path, 'r') as file:
            param_dict = ast.literal_eval(file.read())
        parameters = {k: [v] for k, v in param_dict.items()}
        n_s = 2  # number of splits
        n_r = 1  # number of repeats
        n_i = 1  # number of iterations
    
    xgb_model_hyp = XGBClassifier(objective='multi:softmax', eval_metric='mlogloss', num_class=len(unique_classes), random_state=input_seed, n_jobs=input_jobs)

    rskf = RepeatedStratifiedKFold(n_splits=n_s, n_repeats=n_r, random_state=input_seed)
    randomized_search = RandomizedSearchCV(xgb_model_hyp, param_distributions=parameters,
                                            scoring='f1_weighted', n_iter=n_i, n_jobs=input_jobs,
                                            error_score='raise', cv=rskf.split(X, y), verbose=1,
                                            refit=True, random_state=input_seed)  # refit: train with best hyperparameters found
    # make sure weights for training are added to avoid unbalanced training
    weight = class_weight.compute_sample_weight('balanced', y)
    randomized_search.fit(X, y, sample_weight=weight)

    # save the best parameters
    with open(params_file_path, 'w') as file:
        file.write(str(randomized_search.best_params_))

    best_params = randomized_search.best_params_
    best_score = randomized_search.best_score_
    xgb_model_hyp = randomized_search.best_estimator_
    
    xgb_model_hyp.set_params(**best_params)

    # report best parameters and best model (scores)
    logging.info(f'best found hyperparameters: {best_params}')
    if best_score is not None:
        logging.info(f'score of the model: {best_score}')
        with open(input_file_path, 'a') as file:
            file.write(f"The model with parameters ({best_params}) obtained a ({best_score}) score.\n\n")

    ml_train_end_time = time.time()
    duration = ml_train_end_time - ml_train_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"XGBoost model training finished in {minutes} minutes and {seconds} seconds.")

    logging.info('Prediction is being made, be patient.')
    ml_predict_start_time = time.time()

    # INTERACTION PREDICTION
    emb_col_pred = predict_df['fused_embedding']
    X_pred = pd.DataFrame(emb_col_pred.tolist())
    X_pred = X_pred.astype(float)
    predictions = xgb_model_hyp.predict(X_pred)
    predictions_prob = xgb_model_hyp.predict_proba(X_pred)
    interaction_predictions_df = pd.DataFrame(
        {'drug': predict_df['drug'],
         'gene': predict_df['gene'],
         'predicted_interaction': predictions,
         'prob': np.max(predictions_prob, axis=1)
         })

    # add labels to the drugs and genes for better readability
    all_nodes = nodes + drug_nodes
    id_to_label = {node['id']: node['label'] for node in all_nodes}

    interaction_predictions_df['drug'] = interaction_predictions_df['drug'].apply(
        lambda drug_uri: f'<a href="{drug_uri}">{id_to_label.get(drug_uri, drug_uri)}</a>'
    )
    interaction_predictions_df['gene'] = interaction_predictions_df['gene'].apply(
        lambda gene_uri: f'<a href="{gene_uri}">{id_to_label.get(gene_uri, gene_uri)}</a>'
    )

    # save output files
    interaction_predictions_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_predicted_df.csv')
    interaction_predictions_df.to_csv(interaction_predictions_path, index=False)

    ml_predict_end_time = time.time()
    duration = ml_predict_end_time - ml_predict_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"XGBoost model prediction finished in {minutes} minutes and {seconds} seconds.")

    return interaction_predictions_df


def ml_prediction_with_NS(train_df, predict_df, param_t, input_jobs=None, depth=None, input_seed=None):
    '''
    This function builds a supervised learning model using as training data 
    the network of interactions between biological associations (via 'Monarch.py'),
    gene-to-drug associations (via 'DGIdb.py'), and drug-to-drug similarity (via
    'drug_similarity.py') –to then predict new interactions.
    
    :param train_df: dataframe of (gene, drug, embedding, class) used for training.
    :param predict_df: dataframe of (gene, drug, embedding, class) used for prediction.
    :param param_t: 1 to look for existing optimised parameter file to load instead of performing
        a new parameter optimisation step, 0 to ignore the parameter check.
    :param input_jobs: how many CPU cores to use ('-1' uses all available, just for repetition).
    :param depth: based on this, different sets of parameters, splits, repeats and iterations
        are performed (default = 'full', factor of 'full', 'light' and 'ultralight', just for repetition).
    :param random_seed: seed to have consistent results (default = None, random, just for repetition).
    :return: dataframe of newly predicted interactions.
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    global nodes
    global drug_nodes

    if input_seed == 'random' or input_seed is None:
        input_seed = np.random.randint(1, 1e6)

        ## LOGGING: random seed
        with open(input_file_path, 'w') as file:
            file.write(f"Random seed generated for the run: {input_seed}\n")
            file.write("\n\n")

    if depth is None:
        depth = 'light'

    if input_jobs is None:
        input_jobs = num_cores // 2

    with open(input_file_path, 'a') as file:
        file.write(f"The XGBoost prediction is run on ({input_jobs}) cores with seed ({input_seed}).\n")
        file.write(f"A ({depth}) parameter optimization step is performed.\n")

    logging.info('ML model is being trained, be patient.')
    ml_train_start_time = time.time()

    # MODEL TRAINING
    # X = embeddings (turn list of embeddings into columns)
    emb_col = train_df['fused_embedding']
    X = pd.DataFrame(emb_col.tolist())
    X = X.astype(float)

    # y_1: Convert `0` to `1` and `-1` to `0` in train_df (valid negative as positive)
    logging.info('The training dataframe is being converted: "valid negative" to "positive".')
    y_1_start_time = time.time()
    y_1 = train_df['class'].replace({0: 1, -1: 0}).astype(int)
    y_1_end_time = time.time()
    duration = y_1_end_time - y_1_start_time
    minutes = int(duration // 60)
    seconds = int(duration % 60)
    logging.info(f"Conversion finished in {minutes} minutes and {seconds} seconds.")

    # y_2: Convert only `-1` to `0` in train_df (valid negative as negative)
    logging.info('The training dataframe is being converted: "valid negative" to "negative".')
    y_2_start_time = time.time()
    y_2 = train_df['class'].replace({-1: 0}).astype(int)
    y_2_end_time = time.time()
    duration = y_2_end_time - y_2_start_time
    minutes = int(duration // 60)
    seconds = int(duration % 60)
    logging.info(f"Conversion finished in {minutes} minutes and {seconds} seconds.")

    # define parameters to be tuned, and model
    if depth == 'full':
        parameters = {
            'min_child_weight': [2, 3, 5, 8, 13, 20, 30],
            'gamma': [0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.5, 4, 6],
            'reg_alpha': [0, 0.5, 1, 3, 5, 10],
            'reg_lambda': [0, 0.5, 1, 3, 5, 10],
            'subsample': uniform(0.5, 0.5),
            'colsample_bytree': uniform(0.2, 0.8),
            'max_depth': [4, 6, 8, 10, 12, 14, 16],
            'n_estimators': [35, 45, 50, 70, 80, 90, 100],
            'learning_rate': uniform(0, 0.3),
        }
        n_s = 10
        n_r = 5
        n_i = 20
    elif depth == 'light':
        parameters = {
            'min_child_weight': [3, 8, 20],
            'gamma': [0.5, 1.2, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20, 30],
            'learning_rate': [0.1],
        }
        n_s = 3
        n_r = 2
        n_i = 5
    elif depth == 'ultralight':
        parameters = {
            'min_child_weight': [5, 20],
            'gamma': [0.5, 2.5],
            'reg_alpha': [1, 5],
            'reg_lambda': [1, 5],
            'subsample': [0.5],
            'colsample_bytree': [0.7],
            'max_depth': [3, 5],
            'n_estimators': [10, 20],
            'learning_rate': [0.1],
        }
        n_s = 2
        n_r = 1
        n_i = 5

    params_file_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_hyperparameters.txt')

    # check if the parameters file is already present: if so and if toggled, run with unique parameters
    if param_t == 1 and os.path.exists(params_file_path):
        with open(params_file_path, 'r') as file:
            param_dict = ast.literal_eval(file.read())
        parameters = {k: [v] for k, v in param_dict.items()}
        n_s = 2
        n_r = 1
        n_i = 1

    # define and train model on y_1
    xgb_model_hyp_1 = XGBClassifier(objective='multi:softmax', eval_metric='mlogloss', num_class=2, random_state=input_seed, n_jobs=input_jobs)
    rskf_1 = RepeatedStratifiedKFold(n_splits=n_s, n_repeats=n_r, random_state=input_seed)
    randomized_search_1 = RandomizedSearchCV(xgb_model_hyp_1, param_distributions=parameters,
                                             scoring='f1_weighted', n_iter=n_i, n_jobs=input_jobs,
                                             error_score='raise', cv=rskf_1.split(X, y_1), verbose=1,
                                             refit=True, random_state=input_seed)
    weight_1 = class_weight.compute_sample_weight('balanced', y_1)
    randomized_search_1.fit(X, y_1, sample_weight=weight_1)
    best_model_1 = randomized_search_1.best_estimator_

    # define and train model on y_2
    xgb_model_hyp_2 = XGBClassifier(objective='multi:softmax', eval_metric='mlogloss', num_class=2, random_state=input_seed, n_jobs=input_jobs)
    rskf_2 = RepeatedStratifiedKFold(n_splits=n_s, n_repeats=n_r, random_state=input_seed)
    randomized_search_2 = RandomizedSearchCV(xgb_model_hyp_2, param_distributions=parameters,
                                             scoring='f1_weighted', n_iter=n_i, n_jobs=input_jobs,
                                             error_score='raise', cv=rskf_2.split(X, y_2), verbose=1,
                                             refit=True, random_state=input_seed)
    weight_2 = class_weight.compute_sample_weight('balanced', y_2)
    randomized_search_2.fit(X, y_2, sample_weight=weight_2)
    best_model_2 = randomized_search_2.best_estimator_

    # prediction and averaging the results
    logging.info('Prediction is being made, be patient.')
    ml_predict_start_time = time.time()

    # INTERACTION PREDICTION
    predict_df_filtered = predict_df.copy()
    emb_col_pred = predict_df_filtered['fused_embedding']
    X_pred = pd.DataFrame(emb_col_pred.tolist())
    X_pred = X_pred.astype(float)

    predictions_1 = best_model_1.predict(X_pred)
    predictions_prob_1 = best_model_1.predict_proba(X_pred)

    predictions_2 = best_model_2.predict(X_pred)
    predictions_prob_2 = best_model_2.predict_proba(X_pred)

    # logical 'OR' for predictions: 1 if either model predicts 1
    final_predictions = (predictions_1 | predictions_2)

    # AVERAGE the prediction probabilities
    averaged_predictions_prob = (predictions_prob_1 + predictions_prob_2) / 2

    interaction_predictions_df = pd.DataFrame(
        {'drug': predict_df_filtered['drug'],
         'gene': predict_df_filtered['gene'],
         'predicted_interaction': final_predictions,
         'prob': np.max(averaged_predictions_prob, axis=1)
         })

    # add labels to the drugs and genes for better readability
    all_nodes = nodes + drug_nodes
    id_to_label = {node['id']: node['label'] for node in all_nodes}

    interaction_predictions_df['drug'] = interaction_predictions_df['drug'].apply(
        lambda drug_uri: f'<a href="{drug_uri}">{id_to_label.get(drug_uri, drug_uri)}</a>'
    )
    interaction_predictions_df['gene'] = interaction_predictions_df['gene'].apply(
        lambda gene_uri: f'<a href="{gene_uri}">{id_to_label.get(gene_uri, gene_uri)}</a>'
    )

    # save output files
    interaction_predictions_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_predicted_df.csv')
    interaction_predictions_df.to_csv(interaction_predictions_path, index=False)

    ml_predict_end_time = time.time()
    duration = ml_predict_end_time - ml_predict_start_time
    minutes = int(duration // 60)
    seconds = int(duration % 60)
    logging.info(f"XGBoost model prediction finished in {minutes} minutes and {seconds} seconds.")

    return interaction_predictions_df


def extract_id_label(href):
    '''
    Simple function to extrapolate 'id' and 'label' from href in 'prefict_df'.

    :param href: long string for 'id' and 'label' combined
    :return: 'id', and 'label' short strings
    '''
    
    id_match = re.search(r'<a href="([^"]+)">', href)
    label_match = re.search(r'<a href="[^"]+">([^<]+)</a>', href)
    
    if id_match and label_match:
        return id_match.group(1), label_match.group(1)
    else:
        logging.error(f"Failed to extract id and label from href: {href}")
        return None, None


def filter_rank_drugs(predict_df,prob_threshold=0.65,cluster_threshold=0.8,over_under=1):
    '''
    This function filters predicted drug-gene interactions based on probability
    and ranks drugs based on their interactions with genes.
    
    :param predict_df: dataframe containing predictions with columns 
        ['drug', 'gene', 'predicted_interaction', 'prob']
    :param prob_threshold: probability threshold for filtering predictions
    :param cluster_threshold: cluster threshold for ranking drugs
    :param over_under: toggle to filter above or below cluster_threshold
    :return: ranked list of drugs with counter, and list of predicted edges
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    global nodes
    
    # filter the dataframe
    filt_df = predict_df[(predict_df['predicted_interaction'] == 1) & (predict_df['prob'] >= prob_threshold)]
    
    # remove any drug already present in 'nodes'
    nodes_set = {node['id'] for node in nodes}
    new_df = filt_df[~filt_df['drug'].isin(nodes_set)]
    
    # rank drugs based on how many genes each drug interacts with
    drug_counts = new_df['drug'].value_counts()
    max_count = drug_counts.max()
    threshold_count = max_count * cluster_threshold
    if over_under == 1:
        ranked_drugs = drug_counts[drug_counts >= threshold_count].index.tolist()
    elif over_under == 0:
        ranked_drugs = drug_counts[drug_counts <= threshold_count].index.tolist()
    
    # create the ranked list with counts
    ranked_list = []
    for drug_href in ranked_drugs:
        drug_id, drug_label = extract_id_label(drug_href)
        if drug_id and drug_label:
            ranked_list.append({'id': drug_id, 'label': drug_label, 'count': drug_counts[drug_href]})
    
    # create a ranked dataframe
    ranked_df = new_df[new_df['drug'].isin(ranked_drugs)]
    
    # build the list of new 'edges'
    id_to_label = {node['id']: node['label'] for node in nodes}
    predict_edges = []
    for _, row in ranked_df.iterrows():
        gene_href = row['gene']
        drug_href = row['drug']
        gene_id, gene_label = extract_id_label(gene_href)
        drug_id, drug_label = extract_id_label(drug_href)
        
        edge = [
            {'id': gene_id, 'label': gene_label},
            {'label': 'xgboost:has_association'},
            {'id': drug_id, 'label': drug_label},
            {'notes': row['prob']}
        ]
        predict_edges.append(edge)
    
    return ranked_list, predict_edges


def absolute_value(val, sizes):
    '''
    This function is used in the probability proportion pie chart within 'run_network_model()'.

    :param val: numerical value.
    :param sizes: list of sizes for pie chart segments.
    :return: formatted string.
    '''
    
    total = sum(sizes)
    absolute = int(val / 100. * total)
    
    return f'{absolute}\n({val:.1f}%)'


def run_network_model(monarch_input,date,run_jobs=None,run_depth=None,run_seed=None,prob_input=None,clust_input=None,emb_toggle=1,param_toggle=1,ou_toggle=1):
    '''
    This function runs the entire network_model script and saves nodes and edges files
    after a prediction made via XGBoost.

    :param monarch_input: the input seed from the run_monarch() step.
    :param date: the date of creation of the disease graph.
    :param run_jobs: any value that would override default input_jobs=(num_cores//2) in 
        ml_prediction().
    :param run_depth: any value that would override default depth='full' in ml_prediction().
    :param run_seed: any value that would override default random_seed=123 in ml_prediction().
    :param prob_input: any value that would override default prob_threshold=0.65 in filter_rank_drugs().
    :param clust_input: any value that would override default cluster_threshold=0.8 in filter_rank_drugs().
    :param emb_toggle: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check in get_embeddings() and fuse_embeddings().
    :param param_toggle: 1 to look for existing optimised parameter file to load instead of performing
        a new parameter optimisation step, 0 to ignore the parameter check in ml_prediction().
    :param out_toggle: 1 to consider over clust_input, 0 to consider under clust_input in filter_rank_drugs().
    :return: list of edges after the prediction, list of nodes after the prediction, list of newly
        repurposed drugs found by the prediction.
    '''

    start_time = time.time()
    
    logging.info(f"NOW RUNNING: {current_function_name()} following 'drug_similarity({monarch_input},{today},min_simil={input_min_simil})'.")
    
    global nodes, edges, drug_nodes, drug_edges

    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    
    network_directory = disease_directories['network_directory']

    # generate networks for biological associations and drug database
    embedding_start_time = time.time()
    logging.info('Embeddings are being generated, be patient.')
    
    full_network = get_network(nodes,edges)
    partial_alldrug_network = get_network(drug_nodes,drug_edges,exclude='dgidb:')
    
    drug_embeddings = get_embeddings(full_network,emb_toggle,node_type_select='drug')
    alldrug_embeddings = get_embeddings(partial_alldrug_network,emb_toggle,node_type_select='drug')
    gene_embeddings = get_embeddings(full_network,emb_toggle,node_type_select='gene')

    embedding_end_time = time.time()
    duration = embedding_end_time - embedding_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"Embedding generation finished in {minutes} minutes and {seconds} seconds.")

    DGIdb_edges = [edge for edge in edges if 'dgidb' in edge[1]['label']]
    training_df, prediction_df = fuse_embeddings(gene_embeddings,drug_embeddings,alldrug_embeddings,DGIdb_edges,nodes,drug_nodes,emb_toggle)

    # train the model and make the predictions
    predicted_df = ml_prediction(training_df,prediction_df,param_toggle,input_jobs=run_jobs,depth=run_depth,input_seed=run_seed)

    # generate prediction outcome proportion plot
    num_ones = (predicted_df['predicted_interaction'] == 1).sum()
    num_zeros = (predicted_df['predicted_interaction'] == 0).sum()
    plt.figure(figsize=(8, 8))
    labels = 'Predicted Interactions (1)', 'No Interactions (0)'
    sizes = [num_ones, num_zeros]
    colors = ['#ff9999', '#66b3ff']
    explode = (0.1, 0)
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct=lambda val: absolute_value(val, sizes), shadow=True, startangle=90)
    plt.axis('equal')
    pie_chart_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_interaction_pie_chart.png')
    plt.savefig(pie_chart_path, transparent=True)
    plt.close()
    
    # generate prediction probability distribution plot
    interaction_df = predicted_df[predicted_df['predicted_interaction'] == 1]
    interaction_probs = interaction_df['prob']
    plt.figure(figsize=(10, 6))
    plt.hist(interaction_probs, bins=30, edgecolor='k', alpha=0.7)
    plt.title('Distribution of Probabilities for Predicted Interactions')
    plt.xlabel('Probability')
    plt.ylabel('Frequency')
    plt.grid(True)
    hist_plot_path = os.path.join(network_directory, f'{disease_name_label}_{date_str}_probability_distribution.png')
    plt.savefig(hist_plot_path, transparent=True)
    plt.close()

    plot_paths = [pie_chart_path, hist_plot_path]

    # rank newly found / re-purposed drugs (from prediction)
    if not prob_input:
        prob_input = 0.65
    if not clust_input:
        clust_input = 0.8
    ranked_drug_list, predicted_edges = filter_rank_drugs(predicted_df,prob_threshold=prob_input,cluster_threshold=clust_input, over_under=ou_toggle)
    for ranked_drug_elem in ranked_drug_list:
        logging.info(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")
    
    full_edges = edges + drug_edges + predicted_edges
    full_nodes = [edge[0] for edge in full_edges] + [edge[2] for edge in full_edges]
    full_nodes = unique_elements(full_nodes)

    # save nodes and edges as CSV
    full_edges_df = pd.DataFrame(full_edges)
    full_edges_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_fullnetwork_edges.csv'), index=False)
    full_nodes_df = pd.DataFrame(full_nodes)
    full_nodes_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_fullnetwork_nodes.csv'), index=False)
    
    logging.info("CSV files saved in network directory.")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"'network_model.py' run finished in {minutes} minutes and {seconds} seconds.")
    
    return full_edges, full_nodes, ranked_drug_list, plot_paths


def run_network_model_with_NS(monarch_input,date,run_jobs=None,run_depth=None,run_seed=None,prob_input=None,clust_input=None,emb_toggle=1,param_toggle=1,ou_toggle=1):
    '''
    This function runs the entire network_model script and saves nodes and edges files
    after a prediction made via XGBoost.

    :param monarch_input: the input seed from the run_monarch() step.
    :param date: the date of creation of the disease graph.
    :param run_jobs: any value that would override default input_jobs=(num_cores//2) in 
        ml_prediction().
    :param run_depth: any value that would override default depth='full' in ml_prediction_with_NS().
    :param run_seed: any value that would override default random_seed=123 in ml_prediction_with_NS().
    :param prob_input: any value that would override default prob_threshold=0.65 in filter_rank_drugs().
    :param clust_input: any value that would override default cluster_threshold=0.8 in filter_rank_drugs().
    :param emb_toggle: 1 to look for existing embedding files to load instead of generating new ones,
        0 to ignore the loading check in get_embeddings_with_NS() and fuse_embeddings_with_NS().
    :param param_toggle: 1 to look for existing optimised parameter file to load instead of performing
        a new parameter optimisation step, 0 to ignore the parameter check in ml_prediction_with_NS().
    :param out_toggle: 1 to consider over clust_input, 0 to consider under clust_input in filter_rank_drugs().
    :return: list of edges after the prediction, list of nodes after the prediction, list of newly
        repurposed drugs found by the prediction.
    '''

    start_time = time.time()
    
    logging.info(f"NOW RUNNING: {current_function_name()} following 'drug_similarity({monarch_input},{today},min_simil={input_min_simil})'.")
    
    global nodes, edges, drug_nodes, drug_edges

    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    
    network_with_NS_directory = disease_directories['NSnetwork_directory']

    # generate networks for biological associations and drug database
    embedding_start_time = time.time()
    logging.info('Embeddings are being generated, be patient.')
    
    full_network = get_network(nodes,edges)
    partial_alldrug_network = get_network(drug_nodes,drug_edges,exclude='dgidb:')
    
    drug_embeddings = get_embeddings_with_NS(full_network,emb_toggle,node_type_select='drug')
    alldrug_embeddings = get_embeddings_with_NS(partial_alldrug_network,emb_toggle,node_type_select='drug')
    gene_embeddings = get_embeddings_with_NS(full_network,emb_toggle,node_type_select='gene')

    embedding_end_time = time.time()
    duration = embedding_end_time - embedding_start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"Embedding generation finished in {minutes} minutes and {seconds} seconds.")

    DGIdb_edges = [edge for edge in edges if 'dgidb' in edge[1]['label']]
    training_df, prediction_df = fuse_embeddings_with_NS(gene_embeddings,drug_embeddings,alldrug_embeddings,DGIdb_edges,nodes,drug_nodes,emb_toggle)

    # train the model and make the predictions
    predicted_df = ml_prediction_with_NS(training_df,prediction_df,param_toggle,input_jobs=run_jobs,depth=run_depth,input_seed=run_seed)

    # generate prediction outcome proportion plot
    num_ones = (predicted_df['predicted_interaction'] == 1).sum()
    num_zeros = (predicted_df['predicted_interaction'] == 0).sum()
    plt.figure(figsize=(8, 8))
    labels = 'Predicted Interactions (1)', 'No Interactions (0)'
    sizes = [num_ones, num_zeros]
    colors = ['#ff9999', '#66b3ff']
    explode = (0.1, 0)
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct=lambda val: absolute_value(val, sizes), shadow=True, startangle=90)
    plt.axis('equal')
    pie_chart_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_interaction_pie_chart.png')
    plt.savefig(pie_chart_path, transparent=True)
    plt.close()
    
    # generate prediction probability distribution plot
    interaction_df = predicted_df[predicted_df['predicted_interaction'] == 1]
    interaction_probs = interaction_df['prob']
    plt.figure(figsize=(10, 6))
    plt.hist(interaction_probs, bins=30, edgecolor='k', alpha=0.7)
    plt.title('Distribution of Probabilities for Predicted Interactions')
    plt.xlabel('Probability')
    plt.ylabel('Frequency')
    plt.grid(True)
    hist_plot_path = os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_probability_distribution.png')
    plt.savefig(hist_plot_path, transparent=True)
    plt.close()

    plot_paths = [pie_chart_path, hist_plot_path]

    # rank newly found / re-purposed drugs (from prediction)
    if not prob_input:
        prob_input = 0.65
    if not clust_input:
        clust_input = 0.8
    ranked_drug_list, predicted_edges = filter_rank_drugs(predicted_df,prob_threshold=prob_input,cluster_threshold=clust_input,over_under=ou_toggle)
    for ranked_drug_elem in ranked_drug_list:
        logging.info(f"{ranked_drug_elem['count']} new edge(s) for drug {ranked_drug_elem['id']}.")
    
    full_edges = edges + drug_edges + predicted_edges
    full_nodes = [edge[0] for edge in full_edges] + [edge[2] for edge in full_edges]
    full_nodes = unique_elements(full_nodes)

    # save nodes and edges as CSV
    full_edges_df = pd.DataFrame(full_edges)
    full_edges_df.to_csv(os.path.join(network_with_NS_directory, f'{disease_name_label}_{date_str}_fullnetwork_edges.csv'), index=False)
    full_nodes_df = pd.DataFrame(full_nodes)
    full_nodes_df.to_csv(os.path.join(network_directory, f'{disease_name_label}_{date_str}_fullnetwork_nodes.csv'), index=False)
    
    logging.info("CSV files saved in network directory.")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"'network_model.py' run finished in {minutes} minutes and {seconds} seconds.")
    
    return full_edges, full_nodes, ranked_drug_list, plot_paths