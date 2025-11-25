"""
NETWORK GRAPH UTILITY MODULE: VISUALIZE NETWORK GRAPHS
Created on March 3rd 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import pandas as pd

import networkx as nx
import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_network(input_nodes, input_edges, exclude=None):
    '''
    This function builds a network object from lists of nodes and edges in a
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
        'smiles': '#ff4500',   # lightred
        'xgboost': '#00ff00',  # green
        'pyg': '#00ff00'       # green (same as xgboost for consistency)
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
        
        rel_type = rel_label.split(':')[0].lower()
        colour = edge_type_colours.get(rel_type, '#000000')  # default to black if not found
        G.add_edge(subj_id, obj_id, label=rel_label, notes=notes, colour=colour)

    return G


def convert_df_to_edges(edges_df):
    """
    Convert edges DataFrame to edge list format.
    
    :param edges_df: DataFrame with 'subject_id', 'subject_label', 'relation',
                    'object_id', 'object_label', 'notes' columns.
    :return: list of edges in the format expected by get_network().
    """
    edges = []
    for _, row in edges_df.iterrows():
        edge = [
            {'id': row['subject_id'], 'label': row['subject_label']},
            {'label': row['relation']},
            {'id': row['object_id'], 'label': row['object_label']},
            {'notes': row['notes']}
        ]
        edges.append(edge)
    return edges


def convert_df_to_nodes(nodes_df):
    """
    Convert nodes DataFrame to node list format.
    
    :param nodes_df: DataFrame with 'id' and 'label' columns.
    :return: list of nodes in the format expected by get_network().
    """
    nodes = []
    for _, row in nodes_df.iterrows():
        node = {'id': row['id'], 'label': row['label']}
        nodes.append(node)
    return nodes


def plot_network(nodes_or_graph, edges=None, output_path=None, title=None, 
                highlight_nodes=None, layout='spring', figsize=(30, 20), 
                node_size=10, edge_width=0.2, alpha=0.6, dpi=300,
                show_labels=False, label_font_size=8,transparent=True,
                special_edge_treatment=True):
    """
    Visualize a network graph.
    
    :param nodes_or_graph: either NetworkX graph or list of nodes.
    :param edges: list of edges (required if nodes_or_graph is a list).
    :param output_path: path to save the image –if None, the plot is displayed.
    :param title: title for the plot.
    :param highlight_nodes: list of node IDs to highlight (larger size).
    :param layout: layout algorithm ('spring', 'circular', 'kamada_kawai', etc.).
    :param figsize: figure size (width, height) in inches.
    :param node_size: base size for nodes.
    :param edge_width: width of edges.
    :param alpha: transparency of nodes and edges.
    :param dpi: resolution of the output image.
    :param show_labels: whether to show node labels.
    :param label_font_size: font size for node labels if shown.
    :param transparent: whether to use transparent background.
    :param special_edge_treatment: highlight prediction edges (xgboost/pyg).
    :return: NetworkX graph object.
    """
    # check input and convert to NetworkX graph if needed
    if isinstance(nodes_or_graph, nx.Graph):
        G = nodes_or_graph
    elif edges is not None:
        # convert DataFrames to lists if needed
        if isinstance(nodes_or_graph, pd.DataFrame):
            nodes_or_graph = convert_df_to_nodes(nodes_or_graph)
        if isinstance(edges, pd.DataFrame):
            edges = convert_df_to_edges(edges)
        
        G = get_network(nodes_or_graph, edges)
    else:
        raise ValueError("If nodes_or_graph is not a NetworkX graph, edges must be provided")
    
    # determine layout
    if layout == 'spring':
        pos = nx.spring_layout(G, k=0.15, iterations=50)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    elif layout == 'spectral':
        pos = nx.spectral_layout(G)
    elif layout == 'shell':
        pos = nx.shell_layout(G)
    else:
        # default to spring layout
        pos = nx.spring_layout(G)
    
    # set up figure
    plt.figure(figsize=figsize)
    
    # determine node sizes (highlight specific nodes if requested)
    if highlight_nodes is None:
        highlight_nodes = []
    
    node_sizes = [node_size * 3 if node in highlight_nodes else node_size for node in G.nodes()]
    
    # get node colors
    node_colors = [G.nodes[node].get('colour', '#000000') for node in G.nodes()]
    
    # draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=alpha)
    
    # special treatment for different edge types
    if special_edge_treatment:
        # prediction edges (xgboost/pyg) get special treatment
        for u, v, data in G.edges(data=True):
            edge_color = data.get('colour', '#000000')
            edge_label = data.get('label', '')
            
            # special treatment for prediction edges (thicker, more visible)
            if 'xgboost:' in edge_label or 'pyg:' in edge_label or edge_color == '#00ff00':
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], 
                                      width=edge_width * 4,  # Thicker
                                      alpha=alpha * 0.8,     # More visible
                                      edge_color=edge_color,
                                      arrowsize=8,
                                      connectionstyle="arc3,rad=0.1")
            # drug similarity edges (slightly thicker)
            elif 'smiles:' in edge_label or edge_color == '#ff4500':
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], 
                                      width=edge_width * 2,
                                      alpha=alpha * 0.6,
                                      edge_color=edge_color,
                                      arrowsize=5,
                                      connectionstyle="arc3,rad=0.1")
            # regular edges (thinner, more transparent)
            else:
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], 
                                      width=edge_width,
                                      alpha=alpha * 0.3,   # More transparent
                                      edge_color=edge_color,
                                      arrowsize=4,
                                      connectionstyle="arc3,rad=0.1")
    else:
        # simple approach: draw all edges the same way
        edge_colors = [G[u][v].get('colour', '#000000') for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, width=edge_width, alpha=alpha*0.5, 
                              edge_color=edge_colors, arrowsize=5, 
                              connectionstyle="arc3,rad=0.1")
    
    # draw labels if requested
    if show_labels:
        # use the 'label' attribute for the labels
        labels = {node: G.nodes[node].get('label', node) for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=label_font_size, 
                                font_color='black', alpha=0.7)
    
    # set title if provided
    if title:
        plt.title(title, fontsize=16)
    
    # turn off axis
    plt.axis('off')
    
    # save or show the plot
    if output_path:
        # create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        plt.savefig(output_path, format="png", dpi=dpi, bbox_inches='tight',
                    pad_inches=0.1, transparent=transparent)
        plt.close()
        print(f"Network graph saved to {output_path}")
    else:
        plt.tight_layout()
        plt.show()
    
    return G


def plot_subnetwork(G, central_node, distance=1, output_path=None, **kwargs):
    """
    Plot a subnetwork centered around a specific node.
    
    :param G: NetworkX graph
    :param central_node: Central node for the subnetwork
    :param distance: Maximum distance from central node to include
    :param output_path: Path to save the image
    :param kwargs: Additional arguments to pass to plot_network()
    :return: Subgraph NetworkX object
    """
    # create subgraph of nodes within distance
    nodes = {central_node}
    curr_distance = 0
    curr_layer = {central_node}
    
    while curr_distance < distance:
        next_layer = set()
        for node in curr_layer:
            next_layer.update(G.neighbors(node))
        nodes.update(next_layer)
        curr_layer = next_layer
        curr_distance += 1
    
    # create the subgraph
    subgraph = G.subgraph(nodes)
    
    # plot the subgraph
    return plot_network(
        subgraph, 
        output_path=output_path, 
        highlight_nodes=[central_node],
        **kwargs
    )

   
def plot_drug_gene_network(edges_df, output_path=None, **kwargs):
    """
    Plot a network focusing on drug-gene interactions.
    
    :param edges_df: DataFrame with edge information
    :param output_path: Path to save the image
    :param kwargs: Additional arguments for plot_network()
    :return: NetworkX graph
    """
    # filter for relevant edges
    filtered_edges = edges_df[
        edges_df['relation'].str.contains('dgidb:|pyg:|xgboost:')
    ]
    
    # extract nodes
    nodes = []
    for _, row in filtered_edges.iterrows():
        nodes.append({'id': row['subject_id'], 'label': row['subject_label']})
        nodes.append({'id': row['object_id'], 'label': row['object_label']})
    
    # remove duplicates
    unique_nodes = []
    seen_ids = set()
    for node in nodes:
        if node['id'] not in seen_ids:
            unique_nodes.append(node)
            seen_ids.add(node['id'])
    
    # convert to edge format
    edges = convert_df_to_edges(filtered_edges)
    
    # create and plot the network
    return plot_network(
        unique_nodes, 
        edges, 
        output_path=output_path,
        title="Drug-Gene Interaction Network",
        **kwargs
    )