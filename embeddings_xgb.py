"""
EMBEDDINGS XGB MODULE: GENERATES NODES/EDGES EMBEDDINGS FOR XGBOOST
Created on February 17th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,multiprocessing,inspect,requests,json,datetime

import pickle
import pandas as pd
import numpy as np

import networkx as nx
from node2vec import Node2Vec
from joblib import Parallel, delayed

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def format_duration(duration):
    hours = int(duration // 3600)
    minutes = int((duration % 3600) // 60)
    seconds = int(duration % 60)
    
    parts = []
    if hours > 0:
        parts.append(f"{hours} hour{'s' if hours > 1 else ''}")
    if hours > 0 or minutes > 0:
        parts.append(f"{minutes} minute{'s' if minutes > 1 else ''}")
    parts.append(f"{seconds} second{'s' if seconds > 1 else ''}")
    
    return " and ".join(parts)


def current_function_name():
    return inspect.currentframe().f_back.f_code.co_name


def unique_elements(nonUnique_list):
    """
    Short function that remove duplicate elements.
    If the list contains nodes, it will simply convert it into a set{}.
    If the list contains edges, it will remove also edges where subject and object
    are inverted, therefore not being recognised as the same by Python.

    :param nonUnique_list: biomedical entities list, where each entity is either a
        node or an edge in association networks.
    :return: list of the same biomedical entities without duplicates.
    """
    
    # if nonUnique_list is empty
    if not nonUnique_list:
        return []
    
    if isinstance(nonUnique_list[0], dict):
        # Handle list of nodes
        nodes_set = set(tuple(sorted(node.items())) for node in nonUnique_list)
        unique_list = [dict(node) for node in nodes_set]

    elif len(nonUnique_list[0]) == 4 and isinstance(nonUnique_list[0], list):
        # Handle list of edges
        unique_list = []
        seen_edges = set()
        for edge in nonUnique_list:
            subj_id = edge[0]['id']
            obj_id = edge[2]['id']
            norm_edge = tuple(sorted([subj_id, obj_id]))
            if norm_edge not in seen_edges:

                # locally store the simplified/normalised edge for parsing
                seen_edges.add(norm_edge)

                # return the actual full edge
                unique_list.append(edge)
                
    else:
        raise ValueError("Input is not recognised.")
    
    return unique_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def process_drug_block(g, g_emb, d_block, edges_df, ns_t, predict_t=0):
    """
    Processes drugs in a block against a single gene.

    :param gene: gene ID.
    :param gene_emd: gene embedding dictionary pair (vector).
    :param drug_block: dictionary slice of drug embeddings.
    :param edges_df: DataFrame of edges.
    :param ns_t: negative sample toggle.
    :param predict_t: if 1, skip class computation and set to NaN (for prediction).
    :return: list of rows for the pairs in this block.
    """
    block_rows = []
    for drug, drug_emb in d_block.items():
        fused_emb = np.multiply(g_emb, drug_emb)
        if predict_t == 1:
            # skip class computation for prediction dataframe
            class_label = np.nan
        else:
            # compute class for training dataframe
            class_label = is_interaction_present(g, drug, edges_df, ns_t)
        block_rows.append({
            'gene': g,
            'drug': drug,
            'fused_embedding': fused_emb,
            'class': class_label
        })
    return block_rows


def split_dict_into_blocks(d, n_blocks):
    """
    Splits a dictionary into n roughly equal blocks.
    :param d: embedding dictionary.
    :param n_blocks: how many blocks to break it into (same as 'n_jobs', as in
        how many CPU cores to use).
    :return: list of dictionary blocks.
    """
    items = list(d.items())
    avg = len(items) // n_blocks
    remainder = len(items) % n_blocks
    blocks = []
    start = 0
    
    for i in range(n_blocks):
        end = start + avg + (1 if i < remainder else 0)
        blocks.append(dict(items[start:end]))
        start = end
        
    return blocks

    
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
    G = nx.DiGraph()

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

    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    
    # add nodes
    for node in input_nodes:
        if node is None:
            print("Encountered NoneType node, skipping...")
            continue
        node_id = node.get('id')
        label = node.get('label')
        if node_id is None or label is None:
            print(f"Encountered node with missing id or label: {node}, skipping...")
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
            print(f"Encountered edge with missing id or label: {edge}, skipping...")
            continue
        
        if exclude is not None and exclude in rel_label:
            continue  # skip this edge if exclusion criteria match
        
        rel_type = rel_label.split(':')[0]
        colour = edge_type_colours.get(rel_type, '#000000')  # default to black if not found
        G.add_edge(subj_id, obj_id, label=rel_label, notes=notes, colour=colour)

    return G


def is_interaction_present(gene_id, drug_id, edges_df, ns_t):
    """
    Check gene-drug interaction with optional negative sample handling
    
    :param gene_id: gene node ID.
    :param drug_id: drug node ID.
    :param edges_df: DataFrame of edges.
    :param ns_toggle: whether to consider valid and invalid negative samples
        (ns_t==1, labels [0,1]) or not (ns_t==0, labels [-1,0,1]), based on
        wheter run_negsamples() was run or not.
    :return: corresponding interaction label.
    """
    positive_interaction = (
        ((edges_df['subject_id'] == gene_id) & (edges_df['object_id'] == drug_id)) |
        ((edges_df['subject_id'] == drug_id) & (edges_df['object_id'] == gene_id))
    )
    
    # valid negative sample specific check
    if ns_t == 1:
        valid_negative_interaction = (
            ((edges_df['subject_id'] == gene_id) & (edges_df['object_id'] == drug_id) & 
             (edges_df['relation'] == 'biolink:valid_negative_association')) |
            ((edges_df['subject_id'] == drug_id) & (edges_df['object_id'] == gene_id) & 
             (edges_df['relation'] == 'biolink:valid_negative_association'))
        )
        
        if positive_interaction.any():
            return 1
        elif valid_negative_interaction.any():
            return 0
        else:
            return -1
    
    # default binary behavior (ns_t == 0)
    return 1 if positive_interaction.any() else 0
    

def get_embeddings(input_network,dis_dir,emb_l,node_type_select=None,ns_t=0,
                  custom_node_type=None):
    """
    Build embeddings (vectors) for node information from a network object via Node2Vec.

    :param input_network: NetworkX object.
    :param dis_dir: base paths to where data is stored.
    :param emb_l: 1 to load existing embedding files, 0 to generate new ones.
    :param node_type_select: Node type to focus on (default = None).
    :param ns_t: whether to consider valid and invalid negative samples
        (ns_t==1, labels [0,1]) or not (ns_t==0, labels [-1,0,1]), based on
        wheter run_negsamples() was run or not.
    :return: Embedding pickle object.
    """
    # initialise path
    emb_dir = dis_dir['embeddingsXGB_directory']
    date_str = dis_dir['date_string']
    dis_name = dis_dir['disease_name']
    
    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    valid_node_types = set(node_type_dict.keys())
    
    if node_type_select and node_type_select not in valid_node_types:
        raise ValueError(f"Invalid node_type_select: {node_type_select}. Valid options are: {valid_node_types}")

    # PKL file path
    ns_suffix = '_NS' if ns_t else ''
    if node_type_select is not None:
        node_type_for_path = custom_node_type if custom_node_type is not None else node_type_select
        embedding_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_embedding_dict{ns_suffix}_{node_type_for_path}.pkl')
    else:
        embedding_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_embedding_dict{ns_suffix}_full.pkl')

    # check if they already exist, if toggled
    if emb_l == 1 and os.path.exists(embedding_path):
        with open(embedding_path, "rb") as file:
            output_embeddings = pickle.load(file)
        return output_embeddings

    # filter nodes if needed
    if node_type_select:
        nodes_to_include = [node for node in input_network.nodes() 
                             if input_network.nodes[node].get('node_type') == node_type_select]
        subgraph = input_network.subgraph(nodes_to_include)
    else:
        subgraph = input_network

    # generate node embeddings using Node2Vec
    node2vec = Node2Vec(subgraph,dimensions=64,walk_length=30,num_walks=200,workers=2,quiet=True)
    model = node2vec.fit(window=10,min_count=1,batch_words=4)
    output_embeddings = {node: model.wv[node] for node in subgraph.nodes()}

    # save embeddings
    with open(embedding_path, "wb") as file:
        pickle.dump(output_embeddings, file)

    return output_embeddings
    

def fuse_embeddings(gene_embeddings,drug_embeddings,DGIdb_edges,
                    dis_dir,emb_l,mode="train",ns_t=0,n_jobs=None):
    """
    Fuses embeddings for gene-drug interactions with optional negative sample handling.

    :param gene_embeddings: gene node embeddings.
    :param drug_embeddings: drug node embeddings.
    :param DGIdb_edges: list of gene-drug interactions.
    :param dis_dir: base paths to where data is stored.
    :param emb_l: toggle for loading existing files.
    :param mode: "train" or "predict" to specify operation mode.
    :param ns_t: whether to consider valid and invalid negative samples
        (ns_t==1, labels [0,1]) or not (ns_t==0, labels [-1,0,1]), based on
        wheter run_negsamples() was run or not.
    :return: 'training' or 'prediction' DataFrame with fused embeddings.
    """
    # initialise path
    emb_dir = dis_dir['embeddingsXGB_directory']
    date_str = dis_dir['date_string']
    dis_name = dis_dir['disease_name']
    
    if n_jobs is None:
        n_jobs = max(1, multiprocessing.cpu_count() // 2)
    
    # CSV file paths with negative sample suffix
    ns_suffix = '_NS' if ns_t else ''
    if mode == "train":
        df_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_training_df{ns_suffix}.csv')
        report_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_training{ns_suffix}.txt')
        predict_t = 0  # always compute class for 'training'
    elif mode == "predict":
        df_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_prediction_df{ns_suffix}.csv')
        report_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_prediction{ns_suffix}.txt')
        predict_t = 1  # skip class computation for 'prediction'
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'train' or 'predict'")
    
    # check for existing files
    if emb_l == 1 and os.path.exists(df_path):
        df = pd.read_csv(df_path)
        # Convert embeddings from string to numpy array
        df['fused_embedding'] = df['fused_embedding'].apply(
            lambda x: np.fromstring(x.strip('[]'), sep=' ')
        )
        return df

    # prepare edges DataFrame
    DGIdb_edges_df = pd.DataFrame({
        'subject_id': [edge[0]['id'] for edge in DGIdb_edges],
        'object_id': [edge[2]['id'] for edge in DGIdb_edges],
        'relation': [edge[1]['label'] for edge in DGIdb_edges]
    })

    # generate data
    with open(report_path, 'w') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting {mode} fusion process at {current_time}...\n\n")
        
    idx = 0
    rows = []
    for gene, gene_emb in gene_embeddings.items():
        if n_jobs == 1:
            for drug, drug_emb in drug_embeddings.items():
                fused_emb = np.multiply(gene_emb, drug_emb)
                if mode == "train":
                    class_label = is_interaction_present(gene, drug, DGIdb_edges_df, ns_t)
                else:
                    class_label = np.nan
                rows.append({
                    'gene': gene, 
                    'drug': drug, 
                    'fused_embedding': fused_emb, 
                    'class': class_label
                })
        else:
            drug_blocks = split_dict_into_blocks(drug_embeddings, n_jobs)
            rows_blocks = Parallel(n_jobs=n_jobs)(
                delayed(process_drug_block)(
                    gene,
                    gene_emb,
                    drug_block,
                    DGIdb_edges_df,
                    ns_t,
                    predict_t  # 0 = train (compute), 1 = predict (NaN)
                )
                for drug_block in drug_blocks
            )
            rows.extend([row for block in rows_blocks for row in block])
        
        idx += 1
        if idx % max(1, len(gene_embeddings) // 10) == 0:
            with open(report_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"{idx}/{len(gene_embeddings)} genes processed at {current_time}\n")
                
    df = pd.DataFrame(rows)
    df.to_csv(df_path, index=False)
    return df


def run_embeddings_xgb(nodes,edges,drug_nodes,drug_edges,disease_directories,ns_tog,emb_load=0,num_jobs=None):
    """
    This function runs the whole embedding generation scripts for network analysis.
    
    :param nodes: all of the biological network nodes.
    :param edges: all of the biological network edges.
    :param drug_nodes: all the drug nodes from DGIdb.
    :param drug_edges: all the drug edges from DGIdb.
    :param disease_directories: base paths to where data is stored.
    :param emb_load: toggle for loading/generating embeddings.
    :param n_jobs: how many CPU to consider for parallelisation.
    :return: tuple of elements for subsequest run_networkmodel()
        (gene_embeddings, drug_embeddings, alldrug_embeddings,
        training_df, prediction_df).
    """
    start_time = time.time()

    if ns_tog == 1:
        print(f"NOW RUNNING: {current_function_name()} following 'run_negsamples()'.")
    else:
        print(f"NOW RUNNING: {current_function_name()} following 'run_drugsimilarity()'.")
    
    # create network
    full_network = get_network(nodes, edges)
    partial_alldrug_network = get_network(drug_nodes, drug_edges, exclude='dgidb:')
    
    # generate embeddings
    embedding_start_time = time.time()
    drug_embeddings = get_embeddings(full_network,disease_directories,emb_load,
                                     node_type_select='drug',ns_t=ns_tog)
    alldrug_embeddings = get_embeddings(partial_alldrug_network,disease_directories,
                                        emb_load,node_type_select='drug',ns_t=ns_tog,
                                        custom_node_type='alldrug')
    gene_embeddings = get_embeddings(full_network,disease_directories,emb_load,
                                     node_type_select='gene',ns_t=ns_tog)
    embedding_end_time = time.time()
    embedding_duration = embedding_end_time - embedding_start_time  # calculate duration in seconds
    embedding_formatted_duration = format_duration(embedding_duration)  # convert for print
    
    # fuse embeddings
    DGIdb_edges = [edge for edge in edges if 'dgidb' in edge[1]['label']]
    
    train_fusion_start_time = time.time()
    training_df = fuse_embeddings(
        gene_embeddings,
        drug_embeddings,
        DGIdb_edges,
        disease_directories,
        emb_load,
        mode="train",
        ns_t=ns_tog,
        n_jobs=num_jobs
    )
    train_fusion_end_time = time.time()
    train_fusion_duration = train_fusion_end_time - train_fusion_start_time
    train_fusion_formatted_duration = format_duration(train_fusion_duration)

    predict_fusion_start_time = time.time()
    prediction_df = fuse_embeddings(
        gene_embeddings,
        alldrug_embeddings,
        DGIdb_edges,
        disease_directories,
        emb_load,
        mode="predict",
        ns_t=ns_tog,
        n_jobs=num_jobs
    )
    predict_fusion_end_time = time.time()
    predict_fusion_duration = predict_fusion_end_time - predict_fusion_start_time
    predict_fusion_formatted_duration = format_duration(predict_fusion_duration)

    end_time = time.time()
    duration = end_time - start_time
    formatted_duration = format_duration(duration)
    print(f"'embeddings.py' run finished in {formatted_duration} –where the generation of embeddings itself took {embedding_formatted_duration}, and the fusion of embeddings took {train_fusion_formatted_duration} for training and {predict_fusion_formatted_duration} for prediction.")
    
    return (
        gene_embeddings,
        drug_embeddings, 
        alldrug_embeddings, 
        training_df, 
        prediction_df
    )