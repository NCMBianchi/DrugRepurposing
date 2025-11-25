"""
NEGATIVE SAMPLES MODULE: GENERATES NEGATIVE TRIPLES
Updated on February 12th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,inspect,warnings
import pandas as pd
import numpy as np

import networkx as nx

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

def annotate_nodes_with_centrality(nodes, edges, measure, paths):
    """
    Compute the corresponding centrality measure for each node and annotate existing nodes.
    
    :param nodes: list of nodes.
    :param edges: list of edges.
    :param measure: centrality measure (i.e. eigenvector, degree, in-betweenness).
    :param negs_dir: path for negative samples storing. 
    :param dis_label: name of the input disease.
    :return: list of annotated nodes with centrality scores.
    """
    negs_dir = paths['negsamples_directory']
    date_str = paths['date_string']
    dis_label = paths['disease_name']
    
    G = nx.Graph()
    
    for edge in edges:
        G.add_edge(edge[0]['id'], edge[2]['id'])
    
    if measure == "eigen":
        try:
            warnings.filterwarnings('ignore', message='Could not compute weighted eigenvector centrality')
            centrality = nx.eigenvector_centrality_numpy(G)
        except:
            # fallback to unweighted method if numpy method fails
            centrality = nx.eigenvector_centrality(G)
    elif measure == "deg":
        centrality = nx.degree_centrality(G)
    elif measure == "betw":
        centrality = nx.betweenness_centrality(G)
    
    # annotate nodes with centrality scores
    for node in nodes:
        node['centrality_score'] = centrality.get(node['id'], 0)

    # store in a CSV
    centrality_df = pd.DataFrame([
        {'id': node['id'], 'label': node['label'], 'centrality_score': node.get('centrality_score', 0)} 
        for node in nodes])
    centrality_df.to_csv(os.path.join(negs_dir, f'{dis_label}_{date_str}_{measure}_centrality_nodes.csv'),
                         index=False)
    
    return nodes


def generate_negative_samples(positive_edges, dis_dir, similarity_threshold=0.90,
                              centrality_threshold=0.80, method="simple"):
    """
    Generate valid negative samples from the existing edges in the network, based on
    a simple algorithm that considers one degree of distance and similarity between
    drugs.

    :param positive_edges: list of edges.
    :param dis_dir: base paths to where data is stored.
    :param similarity_threshold: minimum similarity score for drug-to-drug subnetwork.
    :param centrality_threshold: minimum centrality score (default = 0.80).
    :param method: method to generate negative samples (i.e. "simple", "eigen", "deg"
        or "betw", deafult = "simple").
    :return: list of valid negative gene-to-drug edges.
    """
    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    
    gene_to_drug_edges = [edge for edge in positive_edges if 'dgidb:' in edge[1]['label']]
    gene_to_gene_edges = [edge for edge in positive_edges if 'biolink:' in edge[1]['label'] and edge[0]['id'].split(':')[0] in node_type_dict['gene'] and edge[2]['id'].split(':')[0] in node_type_dict['gene']]
    drug_to_drug_edges = [edge for edge in positive_edges if 'smiles:' in edge[1]['label']]
    val_neg_edges = []

    all_genes = []
    all_drugs = []
    for edge in gene_to_drug_edges:
        all_genes.append(edge[0])
        all_drugs.append(edge[2])
    for edge in gene_to_gene_edges:
        all_genes.append(edge[0])
        all_genes.append(edge[2])
    for edge in drug_to_drug_edges:
        all_drugs.append(edge[0])
        all_drugs.append(edge[2])
    all_genes = unique_elements(all_genes)
    all_drugs = unique_elements(all_drugs)

    print(f"GENE-to-GENE EDGES: {len(gene_to_gene_edges)}")
    print(f"GENE-to-DRUG EDGES: {len(gene_to_drug_edges)}")
    print(f"DRUG-to-DRUG EDGES: {len(drug_to_drug_edges)}")
    print(f"GENE NODES for negative embeddings: {len(all_genes)}")
    print(f"DRUG NODES for negative embeddings: {len(all_drugs)}")

    if not all_genes or not all_drugs:
        print("ERROR: no genes or drugs found to generate negative samples.")
        return []

    # compute centrality measure (if required by the method)
    if (method == "eigen" or method == "deg" or method == "betw"):
        weighted_nodes = annotate_nodes_with_centrality(all_genes+all_drugs, positive_edges, method,
                                                        dis_dir)
        weighted_nodes_dict = {node['id']: node.get('centrality_score', 0) for node in weighted_nodes}
    elif method == "simple":
        weighted_nodes_dict = {}

    # iterate through gene-to-drug edges to generate valid negative samples
    for edge in gene_to_drug_edges:
        gene = edge[0]
        drug = edge[2]

        # subnetwork of genes
        subnetwork_genes = [gene]
        for g2g_edge in gene_to_gene_edges:
            if g2g_edge[0]['id'] == gene['id']:
                subnetwork_genes.append(g2g_edge[2])
            if g2g_edge[2]['id'] == gene['id']:
                subnetwork_genes.append(g2g_edge[0])
        subnetwork_genes = unique_elements(subnetwork_genes)
        
        # subnetwork of drugs
        subnetwork_drugs = [drug]
        for d2d_edge in drug_to_drug_edges:
            if d2d_edge[0]['id'] == drug['id'] or d2d_edge[2]['id'] == drug['id']:
                similarity_score = float(d2d_edge[3]['notes'].split('similarity score: ')[1])
                if similarity_score >= similarity_threshold:
                    subnetwork_drugs.append(d2d_edge[0])
                    subnetwork_drugs.append(d2d_edge[2])
        subnetwork_drugs = unique_elements(subnetwork_drugs)

        # generate valid negative edges from subnetworks
        for sub_gene in subnetwork_genes:
            for sub_drug in subnetwork_drugs:
                if sub_gene['id'] != gene['id'] and sub_drug['id'] != drug['id']:
                    if (method == "simple" or
                        weighted_nodes_dict.get(sub_gene['id'], 0) >= centrality_threshold or 
                        weighted_nodes_dict.get(sub_drug['id'], 0) >= centrality_threshold or
                        weighted_nodes_dict.get(gene['id'], 0) >= centrality_threshold or
                        weighted_nodes_dict.get(drug['id'], 0) >= centrality_threshold):
                        valid_edge = ([sub_gene,
                                       {'label': 'biolink:valid_negative_association'},
                                       sub_drug,
                                       {'notes': f'related to: {gene["id"]}-{drug["id"]}'}])
                        if valid_edge not in val_neg_edges and valid_edge not in positive_edges:
                            val_neg_edges.append(valid_edge)

    val_neg_edges = unique_elements(val_neg_edges)

    return val_neg_edges


def run_negsamples(positive_edges,d_edges,disease_directories,similarity_t=0.90,
                  centrality_t=0.8,ns_method="simple",ns_load=0):
    """
    This function runs the whole negsamples script and saves valid negative edges files.

    :param positive_edges: all the biolofical network edges resulting from the previous step
        in the pipeline (i.e. 'run_drugsimilarity()').
    :param d_edges: all the drug edges from DGIdb.
    :param disease_directories: base paths to where data is stored.
    :param similarity_t: minimum similarity score (default = 0.90, suggested 0.50 for methods
        other than 'simple').
    :param centrality_t: minimum centrality score (default = 0.80).
    :param ns_method: method to generate negative samples (i.e. "simple", "eigen", "deg"
        or "betw", deafult = "simple").
    :param ns_load: toggle for loading existing files (1) or generating new ones (0).
    :return: semantically valid negative edges between gene and drugs.
    """

    start_time = time.time()

    print(f"NOW RUNNING: {current_function_name()} following 'run_drugsimilarity()'.")
    
    if ns_method == "simple":
        print(f"The minimum similarity threshold for drug-to-drug subnetworks is set to: {similarity_t}.")
    elif (ns_method == "eigen" or ns_method == "deg" or ns_method == "betw"):
        print(f"The minimum similarity threshold for drug-to-drug subnetworks is set to: {similarity_t}.")
        print(f"The minimum centrality threshold for generating negative samples is set to: {centrality_t}.")
    else:
        raise ValueError("Invalid negative sampling method. Choose 'simple', 'eigen', 'deg' or 'betw'.")

    # initialise path
    negsamples_directory = disease_directories['negsamples_directory']
    date_str = disease_directories['date_string']
    disease_name_label = disease_directories['disease_name']

    # define paths for output files
    neg_edges_path = os.path.join(negsamples_directory,f'{disease_name_label}_{date_str}_{ns_method}_negsamples_edges.csv')

    if ns_load == 1 and os.path.exists(neg_edges_path) and os.path.getsize(neg_edges_path) > 0:
        # Load negative edges
        neg_edges_df = pd.read_csv(neg_edges_path)
        valid_negative_edges = []
        
        for _, row in neg_edges_df.iterrows():
            edge = [
                {'id': row['subject_id'], 'label': row['subject_label']},
                {'label': row['relation']},
                {'id': row['object_id'], 'label': row['object_label']},
                {'notes': row['notes']}
            ]
            valid_negative_edges.append(edge)
            
        print(f"Loaded {len(valid_negative_edges)} existing negative edges.")
        
        full_edges = positive_edges + valid_negative_edges
        full_edges = unique_elements(full_edges)
        
        return full_edges 

    try:
        # launch the appropriate method based on the 'ns_method' parameter
        if ns_method == "simple":
            valid_negative_edges = generate_negative_samples(positive_edges, disease_directories,
                                                             similarity_threshold=similarity_t)
        elif (ns_method == "eigen" or ns_method == "deg" or ns_method == "betw"):
            valid_negative_edges = generate_negative_samples(positive_edges, disease_directories,
                                                             similarity_threshold=similarity_t,
                                                             centrality_threshold=centrality_t,
                                                             method=ns_method)
        else:
            raise ValueError("Invalid negative sampling method.")

        # save the valid negative edges as CSV
        edges_df = pd.DataFrame([
            {
                'subject_id': edge[0]['id'],
                'subject_label': edge[0]['label'],
                'relation': edge[1]['label'],
                'object_id': edge[2]['id'],
                'object_label': edge[2]['label'],
                'notes': edge[3]['notes']
            }
            for edge in valid_negative_edges
        ])
        edges_df.to_csv(neg_edges_path,index=False)

        print(f"VALID NEGATIVE EDGES: {len(valid_negative_edges)}")

    except Exception as e:
        print(f"Error in negative samples generation: {str(e)}")
        return []

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    formatted_duration = format_duration(duration)  # convert for print
    print(f"'negsamples.py' run finished in {formatted_duration}.")

    full_edges = positive_edges + valid_negative_edges
    full_edges = unique_elements(full_edges)
    
    return full_edges