"""
EMBEDDINGS_PYG MODULE: GENERATES NODE EMBEDDINGS FOR PYTORCH GEOMETRIC
Created on March 3rd 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,multiprocessing,inspect,requests,json,datetime

import pickle
import pandas as pd
import numpy as np
import gc

import torch
import torch.serialization
from torch_geometric.data import Data, HeteroData
from torch_geometric.transforms import ToUndirected, NormalizeFeatures
import networkx as nx
from node2vec import Node2Vec

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

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
        # handle list of nodes
        nodes_set = set(tuple(sorted(node.items())) for node in nonUnique_list)
        unique_list = [dict(node) for node in nodes_set]

    elif len(nonUnique_list[0]) == 4 and isinstance(nonUnique_list[0], list):
        # handle list of edges
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_network(input_nodes, input_edges, exclude=None):
    '''
    This function builds a NetworkX graph object from lists of nodes and edges in a
    given format. Examples:
    
    NODES = [{'id': 'chembl:CHEMBL348436', 'label': 'CIRSIMARITIN'}]
    
    EDGES = [[{'id': 'MONDO:0007739', 'label': 'Huntington disease'},
    {'label': 'biolink:has_mode_of_inheritance'},
    {'id': 'HP:0000006', 'label': 'Autosomal dominant inheritance'},
    {'notes': 'first_layer'}]]
    
    It can remove a given edge type to build a partial network: it should be
    in the URI format (e.g. 'dgidb:'), and the formula checks if it's in the entire
    relations label.

    :param input_nodes: list of nodes.
    :param input_edges: list of edges.
    :param exclude: edge type to exclude (default = None).
    :return: networkX object.
    '''
    G = nx.DiGraph()

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
        
        # determine node type
        node_type = 'else'
        for type_key, prefixes in node_type_dict.items():
            if any(node_id.startswith(prefix + ':') for prefix in prefixes):
                node_type = type_key
                break
                
        G.add_node(node_id, label=label, node_type=node_type)
    
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
        
        G.add_edge(subj_id, obj_id, label=rel_label, notes=notes)

    return G


def determine_node_type(node_id, data):
    """
    Determine node type from ID or from data.node_labels if available.
    
    :param node_id: Node ID
    :param data: PyG data object with node_labels
    :return: Node type (gene, drug, disease, phenotype, or else)
    """
    # Try to get type from node_labels if it's a dictionary
    if node_id in data.node_labels:
        if isinstance(data.node_labels.get(node_id), dict):
            node_type = data.node_labels.get(node_id, {}).get('node_type', '')
            if node_type:
                return node_type
    
    # Infer type from node ID prefix
    if any(prefix in node_id.lower() for prefix in ['chembl:', 'wikidata:', 'drugbank:', 'iuphar.ligand:']):
        return 'drug'
    elif any(prefix in node_id.lower() for prefix in ['hgnc:', 'mgi:', 'go:', 'ncbigene:', 'zfin:', 'xenbase:']):
        return 'gene'
    elif any(prefix in node_id.lower() for prefix in ['mondo:']):
        return 'disease'
    elif any(prefix in node_id.lower() for prefix in ['hp:']):
        return 'phenotype'
    else:
        return 'else'


def is_drug_id(node_id):
    """
    Check if a node ID corresponds to a drug.
    
    :param node_id: Node ID
    :return: True if drug, False otherwise
    """
    return any(prefix in node_id.lower() for prefix in ['chembl:', 'wikidata:', 'drugbank:', 'iuphar.ligand:'])


def create_pyg_data(G, node_type=None, r_path=None):
    """
    Creates a PyTorch Geometric Data object from a NetworkX graph.
    
    :param G: NetworkX graph.
    :param node_type: optional filter for specific node type.
    :param r_path: path to write progress reports.
    :return: PyG Data object and node mapping dictionary.
    """
    # filter nodes by type if specified
    if node_type:
        nodes = [node for node, attr in G.nodes(data=True) if attr.get('node_type') == node_type]
        subgraph = G.subgraph(nodes)
    else:
        subgraph = G
        nodes = list(G.nodes())

    if r_path:
        with open(r_path, 'a') as f:
            f.write(f"Processing {len(nodes)} nodes and {subgraph.number_of_edges()} edges\n")
        
    # create node mapping (string IDs to numeric indices)
    node_mapping = {node: i for i, node in enumerate(nodes)}
    reverse_mapping = {i: node for node, i in node_mapping.items()}
    
    # create edge index
    edge_index = []
    edge_attr = []

    edges_list = list(subgraph.edges(data=True))
    total_edges = len(edges_list)
    report_interval_edges = max(1, total_edges // 10)
    
    for i, (u, v, data) in enumerate(edges_list):
        if u in node_mapping and v in node_mapping:
            edge_index.append([node_mapping[u], node_mapping[v]])
            # store edge type as a numeric feature
            edge_type = 0  # Default
            if 'dgidb:' in data.get('label', ''):
                edge_type = 1
            elif 'biolink:' in data.get('label', ''):
                edge_type = 2
            elif 'smiles:' in data.get('label', ''):
                edge_type = 3
            edge_attr.append([edge_type])

        # report progress every 10%
        if r_path and (i + 1) % report_interval_edges == 0:
            progress_percent = (i + 1) * 100 // total_edges
            with open(r_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"Edge processing: {progress_percent+1}% ({i+1}/{total_edges}) at {current_time}\n")
    
    # if no edges, handle the empty case
    if not edge_index:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr = torch.zeros((0, 1), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)
    
    # generate temporary node features (will be replaced later)
    # use degree as initial feature - a reasonable baseline
    x = []

    total_nodes = len(nodes)
    report_interval_nodes = max(1, total_nodes // 10)
        
    for i, node in enumerate(nodes):
        # one-hot encode node type
        node_type = G.nodes[node].get('node_type', 'else')
        feature = [0, 0, 0, 0, 0]  # one-hot for node types
        if node_type == 'gene':
            feature[0] = 1
        elif node_type == 'drug':
            feature[1] = 1
        elif node_type == 'disease':
            feature[2] = 1
        elif node_type == 'phenotype':
            feature[3] = 1
        else:
            feature[4] = 1
            
        # add degree as additional feature
        feature.append(subgraph.degree(node))
        x.append(feature)
        
        # report progress every 10%
        if r_path and (i + 1) % report_interval_nodes == 0:
            progress_percent = (i + 1) * 100 // total_nodes
            with open(r_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"Node feature processing: {progress_percent+1}% ({i+1}/{total_nodes}) at {current_time}\n")
    
    x = torch.tensor(x, dtype=torch.float)
    
    # create PyG Data object
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    
    # store mapping information
    data.node_mapping = node_mapping
    data.reverse_mapping = reverse_mapping
    data.node_labels = {node: G.nodes[node].get('label', node) for node in nodes}
    
    return data, node_mapping


def create_hetero_data(G):
    """
    Creates a PyTorch Geometric HeteroData object from a NetworkX graph.
    
    :param G: NetworkX graph with node_type attributes.
    :return: PyG HeteroData object and node mappings.
    """
    data = HeteroData()
    
    # group nodes by type
    node_groups = {}
    node_mappings = {}
    reverse_mappings = {}
    
    for node, attr in G.nodes(data=True):
        node_type = attr.get('node_type', 'unknown')
        if node_type not in node_groups:
            node_groups[node_type] = []
            node_mappings[node_type] = {}
            reverse_mappings[node_type] = {}
        
        # map each node to a numeric index within its type
        idx = len(node_groups[node_type])
        node_groups[node_type].append(node)
        node_mappings[node_type][node] = idx
        reverse_mappings[node_type][idx] = node
    
    # add node features for each type
    for node_type, nodes in node_groups.items():
        # Create simple features (just degree for now)
        features = []
        for node in nodes:
            feat = [G.degree(node)]
            features.append(feat)
        
        # add to HeteroData
        data[node_type].x = torch.tensor(features, dtype=torch.float)
        data[node_type].node_ids = nodes
    
    # add edges between different node types
    edge_types = {}
    
    for u, v, attr in G.edges(data=True):
        u_type = G.nodes[u].get('node_type', 'unknown')
        v_type = G.nodes[v].get('node_type', 'unknown')
        
        # determine edge type based on relation label
        rel_label = attr.get('label', 'unknown')
        if 'dgidb:' in rel_label:
            rel_type = 'interacts'
        elif 'biolink:' in rel_label:
            rel_type = 'associated'
        elif 'smiles:' in rel_label:
            rel_type = 'similar'
        else:
            rel_type = 'connected'
        
        edge_type = (u_type, rel_type, v_type)
        
        if edge_type not in edge_types:
            edge_types[edge_type] = []
        
        # get node indices within their respective types
        u_idx = node_mappings[u_type].get(u)
        v_idx = node_mappings[v_type].get(v)
        
        if u_idx is not None and v_idx is not None:
            edge_types[edge_type].append((u_idx, v_idx))
    
    # add edge indices to HeteroData
    for edge_type, edges in edge_types.items():
        if edges:
            src_type, rel_type, dst_type = edge_type
            edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
            data[src_type, rel_type, dst_type].edge_index = edge_index
    
    return data, node_mappings, reverse_mappings


def get_node_features(G, node_list, embeddings_dim=64, method="node2vec"):
    """
    Generate node features/embeddings for PyG.
    
    :param G: NetworkX graph
    :param node_list: List of node IDs
    :param embeddings_dim: Dimension of embeddings
    :param method: Method to use for feature generation
    :return: Tensor of node features
    """
    if method == "node2vec":
        # use Node2Vec to generate embeddings
        node2vec = Node2Vec(G, dimensions=embeddings_dim, walk_length=30, 
                          num_walks=200, workers=2, quiet=True)
        model = node2vec.fit(window=10, min_count=1, batch_words=4)
        
        # extract embeddings for each node
        features = []
        for node in node_list:
            if node in model.wv:
                features.append(model.wv[node])
            else:
                # If node not in model, use zeros
                features.append(np.zeros(embeddings_dim))
        
        return torch.tensor(features, dtype=torch.float)
    
    elif method == "gcn_init":
        # simple feature initialization for GCN
        # for each node, include degree and one-hot node type
        features = []
        for node in node_list:
            feature = [G.degree(node)]
            
            # one-hot encode node type
            node_type = G.nodes[node].get('node_type', 'else')
            if node_type == 'gene':
                feature.extend([1, 0, 0, 0])
            elif node_type == 'drug':
                feature.extend([0, 1, 0, 0])
            elif node_type == 'disease':
                feature.extend([0, 0, 1, 0])
            else:
                feature.extend([0, 0, 0, 1])
                
            features.append(feature)
            
        return torch.tensor(features, dtype=torch.float)
    
    else:
        raise ValueError(f"Unknown feature generation method: {method}")


def process_for_pyg(G, nodes, edges, drug_nodes, drug_edges, disease_dir, ns_toggle=0):
    """
    Process data for PyG training and prediction.
    
    :param G: NetworkX graph of the whole network.
    :param nodes: list of nodes from previous pipeline steps.
    :param edges: list of edges from previous pipeline steps.
    :param drug_nodes: list of all drug nodes.
    :param drug_edges: list of drug-related edges.
    :param disease_dir: directory information for saving files.
    :param ns_toggle: toggle for negative samples (1 = enabled).
    :return: training data, prediction data, training mapping, prediction mapping and run times.
    """
    emb_dir = disease_dir['embeddingsPYG_directory']
    date_str = disease_dir['date_string']
    dis_name = disease_dir['disease_name']
    report_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_embeddings.txt')
    
    # create PyG data objects
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting Tensor indices computation for training at {current_time}...\n")

    train_tensor_start_time = time.time()
        
    homogeneous_data, training_mapping = create_pyg_data(G, r_path=report_path)

    train_tensor_end_time = time.time()
    train_tensor_duration = train_tensor_end_time - train_tensor_start_time
    train_tensor_formatted_duration = format_duration(train_tensor_duration)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending Tensor indices computation for training at {current_time}...\n\n")
    
    # get gene-drug interactions
    with open(report_path, 'w') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting gene-drug interaction computation at {current_time}...\n")
    
    gene_drug_edges = []
    
    for edge in edges:
        if 'dgidb:' in edge[1]['label']:
            gene_id = edge[0]['id']
            drug_id = edge[2]['id']
            
            # skip if nodes not in graph
            if gene_id not in G.nodes or drug_id not in G.nodes:
                continue
                
            # get node indices
            gene_idx = training_mapping.get(gene_id)
            drug_idx = training_mapping.get(drug_id)
            
            if gene_idx is not None and drug_idx is not None:
                gene_drug_edges.append({
                    'gene_id': gene_id,
                    'drug_id': drug_id,
                    'gene_idx': gene_idx,
                    'drug_idx': drug_idx,
                    'label': 1  # positive interaction
                })
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending gene-drug interaction computation at {current_time}...\n\n")
    
    # add negative samples if enabled
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting negative interaction computation at {current_time}...\n")
    
    if ns_toggle == 1:
        negative_edges = []
        for edge in edges:
            if 'biolink:valid_negative_association' in edge[1]['label']:
                gene_id = edge[0]['id']
                drug_id = edge[2]['id']
                
                # skip if nodes not in graph
                if gene_id not in G.nodes or drug_id not in G.nodes:
                    continue
                    
                # get node indices
                gene_idx = training_mapping.get(gene_id)
                drug_idx = training_mapping.get(drug_id)
                
                if gene_idx is not None and drug_idx is not None:
                    negative_edges.append({
                        'gene_id': gene_id,
                        'drug_id': drug_id,
                        'gene_idx': gene_idx,
                        'drug_idx': drug_idx,
                        'label': 0  # valid negative interaction
                    })
        
        # combine positive and negative edges
        gene_drug_edges.extend(negative_edges)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending negative interaction computation at {current_time}...\n\n")
    
    # prepare data for training
    gene_indices = []
    drug_indices = []
    edge_labels = []
    
    for edge in gene_drug_edges:
        gene_indices.append(edge['gene_idx'])
        drug_indices.append(edge['drug_idx'])
        edge_labels.append(edge['label'])
    
    # create PyG edge index for gene-drug interactions
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting edge indices computation at {current_time}...\n")
    
    if gene_indices:
        edge_index = torch.tensor([gene_indices, drug_indices], dtype=torch.long)
        edge_labels = torch.tensor(edge_labels, dtype=torch.float)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_labels = torch.zeros(0, dtype=torch.float)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending edge indices computation at {current_time}...\n\n")
    
    # create training data object
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting training data computation at {current_time}...\n")
    
    train_data_start_time = time.time()
    
    training_data = Data(
        x=homogeneous_data.x,
        edge_index=homogeneous_data.edge_index,
        edge_attr=homogeneous_data.edge_attr,
        train_edge_index=edge_index,
        train_edge_label=edge_labels
    )

    train_data_end_time = time.time()
    train_data_duration = train_data_end_time - train_data_start_time
    train_data_formatted_duration = format_duration(train_data_duration)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending training data computation at {current_time}...\n\n")
    
    # add metadata for reference
    training_data.node_mapping = training_mapping
    training_data.reverse_mapping = homogeneous_data.reverse_mapping
    training_data.node_labels = homogeneous_data.node_labels
    
    # now prepare prediction data
    # add all drugs to the graph for prediction
    all_nodes = nodes.copy()
    all_nodes.extend([drug for drug in drug_nodes if drug not in nodes])

    # create full graph including all drugs
    full_G = get_network(all_nodes, edges + drug_edges)
    
    # create PyG data for prediction
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting prediction data computation at {current_time}...\n")
    
    pred_data_start_time = time.time()
    
    prediction_data, prediction_mapping = create_pyg_data(full_G, r_path=report_path)

    pred_data_end_time = time.time()
    pred_data_duration = pred_data_end_time - pred_data_start_time
    pred_data_formatted_duration = format_duration(pred_data_duration)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending prediction data computation at {current_time}...\n\n")
    
    # extract genes and all drugs for prediction pairs
    genes = [node for node, attr in G.nodes(data=True) if attr.get('node_type') == 'gene']
    all_drugs = [node for node, attr in full_G.nodes(data=True) if attr.get('node_type') == 'drug']
    
    # create all possible gene-drug pairs for prediction
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting gene-drug pairs computation for prediction at {current_time}...\n")
    
    pairs_start_time = time.time()
    
    gene_drug_pairs = []
    total_genes = len(genes)
    report_interval = max(1, total_genes // 10)
    
    for gene_idx, gene in enumerate(genes):
        gene_node_idx = prediction_mapping.get(gene)
        if gene_node_idx is None:
            continue
            
        # report progress every 10% of genes processed
        if (gene_idx + 1) % report_interval == 0:
            progress_percent = (gene_idx + 1) * 100 // total_genes
            with open(report_path, 'a') as f:
                current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                f.write(f"Gene-drug pairs: {progress_percent}% ({gene_idx+1}/{total_genes} genes) processed at {current_time}\n")
        
        for drug in all_drugs:
            drug_node_idx = prediction_mapping.get(drug)
            if drug_node_idx is None:
                continue
                
            # check if this is a known interaction (from training)
            known = False
            for edge in gene_drug_edges:
                if edge['gene_id'] == gene and edge['drug_id'] == drug:
                    known = True
                    break
                    
            if not known:
                gene_drug_pairs.append({
                    'gene_id': gene,
                    'drug_id': drug,
                    'gene_idx': gene_node_idx,
                    'drug_idx': drug_node_idx
                })

    pairs_end_time = time.time()
    pairs_duration = pairs_end_time - pairs_start_time
    pairs_formatted_duration = format_duration(pairs_duration)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending gene-drug pairs computation for prediction at {current_time}...\n\n")
    
    # create prediction edge index
    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Starting Tensor indices computation for prediction at {current_time}...\n")
    
    pred_tensor_start_time = time.time()
    
    pred_gene_indices = []
    pred_drug_indices = []
    
    for pair in gene_drug_pairs:
        pred_gene_indices.append(pair['gene_idx'])
        pred_drug_indices.append(pair['drug_idx'])
    
    if pred_gene_indices:
        pred_edge_index = torch.tensor([pred_gene_indices, pred_drug_indices], dtype=torch.long)
    else:
        pred_edge_index = torch.zeros((2, 0), dtype=torch.long)
    
    # add prediction edge index to prediction data
    prediction_data.pred_edge_index = pred_edge_index
    
    # add metadata for reference
    prediction_data.gene_drug_pairs = gene_drug_pairs

    pred_tensor_end_time = time.time()
    pred_tensor_duration = pred_tensor_end_time - pred_tensor_start_time
    pred_tensor_formatted_duration = format_duration(pred_tensor_duration)

    with open(report_path, 'a') as f:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"Ending Tensor indices computation for prediction at {current_time}...\n\n")

    run_times = [train_data_formatted_duration, pred_data_formatted_duration, pairs_formatted_duration, train_tensor_formatted_duration, pred_tensor_formatted_duration]
    
    return training_data, prediction_data, training_mapping, prediction_mapping, run_times


def run_embeddings_pyg(nodes, edges, drug_nodes, drug_edges, disease_directories, 
                   ns_tog, emb_load=0, num_jobs=None):
    """
    Generates embeddings and prepares data for PyG training and prediction.
    
    :param nodes: all of the biological network nodes.
    :param edges: all of the biological network edges.
    :param drug_nodes: all the drug nodes from DGIdb.
    :param drug_edges: all the drug edges from DGIdb.
    :param disease_directories: base paths to where data is stored.
    :param ns_tog: toggle for negative samples (1 if used, 0 if not).
    :param emb_load: toggle for loading/generating embeddings.
    :param num_jobs: how many CPU to consider for parallelisation.
    :return: PyG training data, prediction data, and node mappings.
    """
    start_time = time.time()

    if ns_tog == 1:
        print(f"NOW RUNNING: {current_function_name()} following 'run_negsamples()'.")
    else:
        print(f"NOW RUNNING: {current_function_name()} following 'run_drugsimilarity()'.")
    
    # initialize paths
    emb_dir = disease_directories['embeddingsPYG_directory']
    date_str = disease_directories['date_string']
    dis_name = disease_directories['disease_name']
    
    # define file paths for PyG data
    ns_suffix = '_NS' if ns_tog else ''
    training_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_pyg_training{ns_suffix}.pt')
    prediction_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_pyg_prediction{ns_suffix}.pt')
    training_mapping_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_pyg_training_mapping{ns_suffix}.pkl')
    prediction_mapping_path = os.path.join(emb_dir, f'{dis_name}_{date_str}_pyg_prediction_mapping{ns_suffix}.pkl')
    
    # check if files already exist and load them if emb_load is enabled
    if emb_load == 1 and os.path.exists(training_path) and os.path.exists(prediction_path) and os.path.exists(training_mapping_path) and os.path.exists(prediction_mapping_path):
        print("Loading existing PyG data files...")
        
        training_data = torch.load(training_path, weights_only=False)
        prediction_data = torch.load(prediction_path, weights_only=False)
        
        with open(training_mapping_path, 'rb') as f:
            training_mapping = pickle.load(f)
            
        with open(prediction_mapping_path, 'rb') as f:
            prediction_mapping = pickle.load(f)
            
        print("PyG data loaded successfully.")
        
        # for compatibility with original function interface, prepare return values
        gene_embeddings = {}
        drug_embeddings = {}
        alldrug_embeddings = {}
        
        # extract gene and drug embeddings from node features
        for node_id, idx in training_mapping.items():
            if idx < training_data.x.size(0):
                node_type = determine_node_type(node_id, training_data)
                if node_type == 'gene':
                    gene_embeddings[node_id] = training_data.x[idx].numpy()
                elif node_type == 'drug':
                    drug_embeddings[node_id] = training_data.x[idx].numpy()

        for node_id, idx in prediction_mapping.items():
            if idx < prediction_data.x.size(0) and is_drug_id(node_id):
                alldrug_embeddings[node_id] = prediction_data.x[idx].numpy()
        
        # create empty DataFrames to maintain interface compatibility
        training_df = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
        prediction_df = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
        
        end_time = time.time()
        duration = end_time - start_time
        formatted_duration = format_duration(duration)
        print(f"'{current_function_name()}' run finished in {formatted_duration}.")
        
        return gene_embeddings, drug_embeddings, alldrug_embeddings, training_df, prediction_df, training_data, prediction_data
    
    # create network
    full_network = get_network(nodes, edges)
    print(f"Network built with {full_network.number_of_nodes()} nodes and {full_network.number_of_edges()} edges.")
    
    # process data for PyG
    training_data, prediction_data, training_mapping, prediction_mapping, emb_run_times = process_for_pyg(
        full_network, nodes, edges, drug_nodes, drug_edges, disease_directories, ns_tog)

    train_dur = emb_run_times[0]
    pred_dur = emb_run_times[1]
    pairs_dur = emb_run_times[2]
    train_tensor_dur = emb_run_times[3]
    pred_tensor_dur = emb_run_times[4]
    
    # save PyG data to files
    torch.save(training_data, training_path)
    torch.save(prediction_data, prediction_path)
    
    # save mappings
    with open(training_mapping_path, 'wb') as f:
        pickle.dump(training_mapping, f)  # training mapping
    
    with open(prediction_mapping_path, 'wb') as f:
        pickle.dump(prediction_mapping, f)  # prediction mapping (full: Monarch + all DGIdb)
    
    print(f"PyG data saved to {emb_dir}.")
    
    # for compatibility with original function interface, prepare return values
    gene_embeddings = {}
    drug_embeddings = {}
    alldrug_embeddings = {}
    
    # extract gene and drug embeddings from TRAINING data using TRAINING mapping
    for node_id, idx in training_mapping.items():
        if node_id in full_network.nodes:
            node_type = full_network.nodes[node_id].get('node_type', '')
            
            if node_type == 'gene':
                gene_embeddings[node_id] = training_data.x[idx].numpy()
            elif node_type == 'drug':
                drug_embeddings[node_id] = training_data.x[idx].numpy()
            
    # extract all drug embeddings from PREDICTION data using PREDICTION mapping
    for node_id, idx in prediction_mapping.items():
        if idx < prediction_data.x.size(0) and is_drug_id(node_id):
            alldrug_embeddings[node_id] = prediction_data.x[idx].numpy()
    
    # create empty DataFrames to maintain interface compatibility
    training_df = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
    prediction_df = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
    
    end_time = time.time()
    duration = end_time - start_time
    formatted_duration = format_duration(duration)
    print(f"'embeddings_pyg.py' run finished in {formatted_duration} –where the generation of training and prediction embeddings themselves respectively took {train_dur} and {pred_dur}, the generation of drug-to-gene pairs took {pairs_dur} and the generation of prediction indices via Torch Tensor took respectively {train_tensor_dur} and {pred_tensor_dur} for training and prediction embeddings.")
    
    return gene_embeddings, drug_embeddings, alldrug_embeddings, training_df, prediction_df, training_data, prediction_data