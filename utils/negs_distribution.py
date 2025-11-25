"""
NEGS_DISTRIBUTION UTILITY MODULE: ANALYZE NODE DISTRIBUTION IN NEGATIVE SAMPLES
Created on February 17th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os,requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def fetch_disease_name(seed_id):
    """
    Retrieve the name of a specific Monarch Initiative entity.
    
    Args:
        seed_id: A Monarch Initiative ID (e.g., 'MONDO:0007739')
                 https://api-v3.monarchinitiative.org/v3/api/entity/MONDO:0007739
    
    Returns:
        The name of the entity or the seed_id if retrieval fails
    """
    response = requests.get(f'https://api-v3.monarchinitiative.org/v3/api/entity/{seed_id[0]}')
    response.raise_for_status()
    data = response.json()
        
    return data.get('name')
    

def get_centrality_score(node_ids, dis_dir, dis_input):
    """
    Retrieve the centrality measure score for a specific node.
    
    :param node_ids: ID or list of IDs of the nodes to find centrality scores for.
    :param dis_dir: base paths to where data is stored. 
    :param dis_label: name of the input disease.
    :return: none.
    """
    # initialise path
    negs_dir = dis_dir['negsamples_directory']
    date_str = dis_dir['date_string']
    dis_label = dis_dir['disease_name']
    
    method_paths = {
        'degree': os.path.join(negs_dir, f'{dis_label}_{date_str}_deg_centrality_nodes.csv'),
        'eigenvector': os.path.join(negs_dir, f'{dis_label}_{date_str}_eigen_centrality_nodes.csv'),
        'in-betweenness': os.path.join(negs_dir, f'{dis_label}_{date_str}_betw_centrality_nodes.csv')
    }

    
    if isinstance(node_ids, str):
        node_ids = [node_ids]

    node_centrality_data = {}

    for node_id in node_ids:
        node_centrality_data[node_id] = {}
        
        for method, path in method_paths.items():
            try:
                df = pd.read_csv(path)
                max_score = df['centrality_score'].max()
                
                node_row = df[df['id'] == node_id]
                if not node_row.empty:
                    label = node_row['label'].values[0]
                    score = node_row['centrality_score'].values[0]
                    node_centrality_data[node_id][method] = {
                        'label': label,
                        'score': score,
                        'max_score': max_score
                    }
                else:
                    node_centrality_data[node_id][method] = None
            except FileNotFoundError:
                node_centrality_data[node_id][method] = None

    for node_id, methods_data in node_centrality_data.items():
        print(f"\nCentrality Scores for Node: {node_id}")
        print("-" * 80)
        print(f"{'Method':<30} {'Label':<25} {'Centrality Score':<20}")
        print("-" * 80)
        
        # track if node was found in any dataset
        node_found_in_any = False
        
        for method, method_data in methods_data.items():
            if method_data is None:
                print(f"{method.capitalize() + ' Centrality:':<30} {'Measure not computed':<45}")
            else:
                label = method_data['label']
                score = method_data['score']
                max_score = method_data['max_score']
                score_display = f"{score:.4f} (max: {max_score:.4f})"
                print(f"{method.capitalize() + ' Centrality:':<30} {label:<25} {score_display:<20}")
                node_found_in_any = True
        
        # if node was not found in any dataset
        if not node_found_in_any:
            print("\nWARNING: Node not found in any centrality dataset")


def print_node_distribution(counts, node_type):
    print(f"\n{node_type} NODES IN NEGATIVE SAMPLES:")
    print("-" * 45)
    print(f"{'Node ID':<30} {'Count':<10}")
    print("-" * 45)
    
    if len(counts) <= 10:
        # if 10 or fewer nodes, just print all
        for node, count in counts.items():
            print(f"{node:<30} {count:<10}")
    else:
        # otherwise split to show top 5 and bottom 5
        for node, count in counts.head(5).items():
            print(f"{node:<30} {count:<10}")
        print(f"{'(...)':^40}")
        for node, count in counts.tail(5).items():
            print(f"{node:<30} {count:<10}")

            
def report_negs(negs_edges, prnt=True):
    """
    Analyze the distribution of nodes in negative samples.
    
    :param negative_edges: list of negative sample edges from 'negsamples.py'.
    :param prnt: whether to print the report: used to pass the computed stats
        to plot_negs() without printing the table when prnt=False and instead
        returning values.
    :return: none, or dictionaries of gene and drug node counts.
    """ 
    df = pd.DataFrame([
        {'subject_id': edge[0]['id'], 'object_id': edge[2]['id']} 
        for edge in negs_edges
    ])
    
    gene_counts = pd.concat([
        df['subject_id'][df['subject_id'].str.contains('HGNC:|MGI:|GO:|NCBIgene:|ZFIN:|Xenbase:')],
        df['object_id'][df['object_id'].str.contains('HGNC:|MGI:|GO:|NCBIgene:|ZFIN:|Xenbase:')]
    ]).value_counts()
    
    drug_counts = pd.concat([
        df['subject_id'][df['subject_id'].str.contains('chembl:|wikidata:')],
        df['object_id'][df['object_id'].str.contains('chembl:|wikidata:')]
    ]).value_counts()
    
    if prnt:
        # gene node distribution
        print_node_distribution(gene_counts, 'GENE')
        
        # drug node distribution
        print_node_distribution(drug_counts, 'DRUG')
        
    else:
        return gene_counts, drug_counts


def plot_negs(negs_edges):
    """
    Visualize the distribution of nodes in negative samples.
    
    :param negative_edges: list of negative sample edges.
    :return: none.
    """
    # seaborn does not accept palettes, but plt.hist() could
    palette_1 = ['#0D1B2A', '#1B263B', '#415A77', '#778DA9', '#BDC0C6']
    palette_2 = ['#1A0D29', '#281B3B', '#5D4278', '#8D76A8', '#C2BDC7']
    
    gene_counts, drug_counts = report_negs(negs_edges, prnt=False)
    
    plt.figure(figsize=(15, 6))
    
    # gene nodes distribution plot
    plt.subplot(1, 2, 1)
    sns.histplot(gene_counts.values, kde=True, color=palette_1[2])
    plt.title('GENE Nodes Distribution in Negative Samples')
    plt.xlabel('Number of Occurrences')
    plt.ylabel('Frequency')
    plt.yscale('symlog', linthresh=10)
    
    # drug nodes distribution plot
    plt.subplot(1, 2, 2)
    sns.histplot(drug_counts.values, kde=True, color=palette_2[2])
    plt.title('DRUG Nodes Distribution in Negative Samples')
    plt.xlabel('Number of Occurrences')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()