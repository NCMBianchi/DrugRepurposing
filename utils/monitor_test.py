"""
MONITOR_TEST UTILITY MODULE: monitor CPU strain, 
Created on March 11th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import psutil, threading, time
import networkx as nx

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def monitor_cpu(cpu_data, stop_monitoring):
    """
    Helper function to monitor CPU strain while running the `run_test_pipeline()` function in
    the `test-runs.ipynb` notebook.
    
    :param cpu_data: dictionary to store CPU usage data.
    :param stop_monitoring: threading event to signal when to stop monitoring.
    """
    while not stop_monitoring.is_set():
        cpu_percentages = psutil.cpu_percent(interval=0.5, percpu=True)
        cpu_data['usage'].append(cpu_percentages)
        cpu_data['timestamp'].append(time.time())


def get_cpu_stats(cpu_data):
    """
    Process CPU monitoring data and return statistics.
    
    :param cpu_data: dictionary with CPU monitoring data.
    :return: core_stats, cpu_cores.
    """
    cpu_cores = len(cpu_data['usage'][0]) if cpu_data['usage'] else 0
    core_stats = []
    
    for core in range(cpu_cores):
        core_values = [usage[core] for usage in cpu_data['usage']]
        if core_values:
            core_stats.append({
                'core': core,
                'min': min(core_values),
                'max': max(core_values),
                'avg': sum(core_values) / len(core_values)
            })
    return core_stats, cpu_cores


def format_value(param, value, parameters):
    """
    Helper function to format values based on defaults in the table for the `run_test_pipeline()`
    function.

    :param param: parameter name.
    :param value: current parameter value.
    :param parameters: dictionary of all parameters passed to the function.
    :return: formatted value string.
    """
    if param == 'ml_model_type' and parameters.get('ml_method') == 'xgboost':
        return "not applicable for XGBoost"
    
    if param in parameters and parameters.get(param) not in [None, 0]:
        return str(value)
    else:
        if param == 'ds_k':
            default = '10'
        elif param == 'ds_min':
            default = '0.5'
        elif param == 'ds_radius':
            default = '2'
        elif param == 'ds_feat_length':
            default = '4096'
        elif param == 'ns_toggle':
            default = '"with negative samples"'
        elif param == 'ns_method':
            default = 'deg'
        elif param == 'ns_simil_t':
            default = '0.9 or 0.5 depending on method'
        elif param == 'ns_centr_t':
            default = '0.0001'
        elif param == 'emb_jobs':
            default = f'half of cores (max 4)'
        elif param == 'ml_method':
            default = 'xgboost'
        elif param == 'ml_model_type':
            default = 'gcn'
        elif param == 'ml_jobs':
            default = 'half of cores (max 2)'
        elif param == 'ml_depth':
            default = 'ultralight'
        elif param == 'ml_seed':
            default = 'random'
        elif param == 'ml_batch_size':
            default = '2000'
        elif param == 'ml_prob_filter':
            default = '0.65'
        elif param == 'ml_clust_filter':
            default = '0.8'
        elif param == 'ml_iterations':
            default = '1'
        elif param == 'ml_iter_max_edges':
            default = '100'
        else:
            default = 'unknown'
        return f"defaulted to {default}"


def get_high_centrality_gene_edges(nodes, edges, n_top_genes=3, centrality_method="eigen"):
    """
    Identify gene-drug edges involving the top N genes by centrality score.
    
    :param nodes: list of network nodes.
    :param edges: list of network edges.
    :param n_top_genes: number of top genes to select.
    :param centrality_method: centrality method (eigen, deg, betw).
    :return: list of gene-drug edges to remove.
    """
    # create network for centrality calculation
    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0]['id'], edge[2]['id'])
    
    # get only gene nodes
    gene_nodes = [node for node in nodes 
                  if any(prefix in node['id'].lower() 
                         for prefix in ['hgnc:', 'mgi:', 'go:', 'ncbigene:', 'zfin:', 'xenbase:'])]
    
    # calculate centrality for genes
    gene_ids = [node['id'] for node in gene_nodes]
    gene_subgraph = G.subgraph(gene_ids)
    
    # calculate centrality based on the selected method
    if centrality_method == "eigen":
        try:
            centrality = nx.eigenvector_centrality_numpy(gene_subgraph)
        except:
            print("Falling back to standard eigenvector centrality...")
            centrality = nx.eigenvector_centrality(gene_subgraph)
    elif centrality_method == "deg":
        centrality = nx.degree_centrality(gene_subgraph)
    elif centrality_method == "betw":
        centrality = nx.betweenness_centrality(gene_subgraph)
    else:
        print(f"Unknown centrality method: {centrality_method}, using degree centrality")
        centrality = nx.degree_centrality(gene_subgraph)
    
    # sort genes by centrality and get top N
    sorted_genes = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    top_genes = [gene_id for gene_id, score in sorted_genes[:n_top_genes]]
    
    # display top genes and their scores
    print(f"Top {n_top_genes} genes by {centrality_method} centrality:")
    print(f"{'Gene ID':<30} {'Score':<10} {'Label':<30}")
    print("-" * 70)
    for gene_id, score in sorted_genes[:n_top_genes]:
        gene_label = next((node['label'] for node in gene_nodes if node['id'] == gene_id), "Unknown")
        print(f"{gene_id:<30} {score:<10.4f} {gene_label:<30}")
    
    # find all gene-drug edges involving these genes
    gene_drug_edges = []
    for edge in edges:
        # check if edge is a gene-drug interaction
        if 'dgidb:interacts_with' in edge[1]['label']:
            # check if one of our top genes is involved
            if edge[0]['id'] in top_genes:  # gene is subject
                gene_drug_edges.append(edge)
            elif edge[2]['id'] in top_genes:  # gene is object (just in case)
                gene_drug_edges.append(edge)
    
    print(f"Found {len(gene_drug_edges)} gene-drug edges involving top {n_top_genes} genes")
    
    return gene_drug_edges


def print_gene_drug_edge_table(edges):
    """
    Prints a nicely formatted table of gene-drug edges.
    
    :param edges: list of edges in the format [gene_node, relation, drug_node, notes].
    :return: None.
    """
    print(f"{'Gene ID':<15} {'Gene Name':<15} | {'Drug ID':<25} {'Drug Name':<25} | {'Score':<10}")
    print("-" * 95)
    
    for edge in edges:
        gene_id = edge[0]['id']
        gene_name = edge[0]['label']
        drug_id = edge[2]['id']
        drug_name = edge[2]['label']
        
        # extract interaction score from notes
        score_str = edge[3]['notes']
        score = score_str.split(':')[-1] if ':' in score_str else "N/A"
        
        # try to convert score to float and format it
        try:
            score_value = float(score)
            score_formatted = f"{score_value:.4f}"
        except:
            score_formatted = score
            
        print(f"{gene_id:<15} {gene_name:<15} | {drug_id:<25} {drug_name:<25} | {score_formatted:<10}")


def filter_edges_by_pairs(all_edges, gene_drug_pairs):
    """
    Filters edges to keep only those matching the specified gene-drug pairs.
    
    :param all_edges: list of all edges.
    :param gene_drug_pairs: list of tuples (gene_id, drug_id).
    :return: list of edges matching the pairs.
    """
    filtered_edges = []
    
    # convert pairs to a set for faster lookup
    pairs_set = set(gene_drug_pairs)
    
    for edge in all_edges:
        if 'dgidb:interacts_with' in edge[1]['label']:  # snsure it's a gene-drug edge
            gene_id = edge[0]['id']
            drug_id = edge[2]['id']
            
            # check if this pair exists in our set
            if (gene_id, drug_id) in pairs_set:
                filtered_edges.append(edge)
            
            # also check the reverse direction (drug as subject, gene as object)
            elif (drug_id, gene_id) in pairs_set:
                filtered_edges.append(edge)
    
    return filtered_edges


def remove_edges(network_edges, edges_to_remove):
    """
    Removes specific edges from a network.
    
    :param network_edges: original list of network edges.
    :param edges_to_remove: edges to be removed.
    :return: modified list with edges removed.
    """
    # for each edge in the original network
    partial_edges = []
    
    # create identifiers for edges to remove
    remove_ids = set()
    for edge in edges_to_remove:
        # Create a unique identifier based on key elements
        edge_id = (edge[0]['id'], edge[1]['label'], edge[2]['id'])
        remove_ids.add(edge_id)
    
    # keep edges that aren't in the remove set
    for edge in network_edges:
        edge_id = (edge[0]['id'], edge[1]['label'], edge[2]['id'])
        if edge_id not in remove_ids:
            partial_edges.append(edge)
            
    return partial_edges