"""
NODE_TYPES UTILITY MODULE: PRINT THE NODE TYPES
Created on February 13th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def report_node_types(input_nodes,prnt=True):
    """
    Prints tables for node URI(types) and categories.
    
    :param nodes: list of nodes from Monarch Initiative's database and DGIdb.
    :param prnt: whether to print the report: used to pass the computed stats
        to plot_node_types() without printing the table when prnt=False and instead
        returning values.
    :return: none, or list counts for plot_node_types().
    """
    node_type_dict = {
        'disease': ['MONDO'],
        'gene': ['HGNC', 'MGI', 'GO', 'NCBIgene', 'ZFIN', 'Xenbase'],
        'phenotype': ['HP'],
        'drug': ['chembl', 'wikidata']
    }
    
    uri_counts = {}
    for node in input_nodes:
        if 'id' in node:
            uri = node['id'].split(':')[0]
            uri_counts[uri] = uri_counts.get(uri, 0) + 1
            
    category_counts = {cat: 0 for cat in list(node_type_dict.keys()) + ['-']}
    uri_categories = {}
    
    for uri in uri_counts:
        category_found = False
        for category, prefixes in node_type_dict.items():
            if uri.lower() in [prefix.lower() for prefix in prefixes]:
                category_counts[category] += uri_counts[uri]
                uri_categories[uri] = category
                category_found = True
                break
        if not category_found:
            category_counts['-'] += uri_counts[uri]
            uri_categories[uri] = '-'
            
    if prnt:
        print("\nNODE URI Prefix Counts:")
        print("-" * 45)
        print(f"{'URI':<15} {'Count':<10} {'Category':<15}")
        print("-" * 45)
        for uri, count in sorted(uri_counts.items()):
            print(f"{uri:<15} {count:<10} {uri_categories[uri]:<15}")
            
        print("\nNODE Category Counts:")  
        print("-" * 30)
        print(f"{'Category':<15} {'Count':<10}")
        print("-" * 30)
        for category, count in sorted(category_counts.items()):
            if count > 0:  # only show categories with nodes
                print(f"{category:<15} {count:<10}")
    else:
        return category_counts


def plot_node_types(input_nodes):
    """
    Displays a pie chart for node categories.
    
    :param nodes: list of nodes from Monarch Initiative's database and DGIdb.
    :return: none.
    """
    category_counts = report_node_types(input_nodes,prnt=False)
    colors = ['#cad2c5', '#84a98c', '#52796f', '#354f52', '#2f3e46']
   
    plt.figure(figsize=(10, 8))
    wedges, texts, autotexts = plt.pie(category_counts.values(),
                                       labels=[''] * len(category_counts),
                                       autopct='%1.1f%%',
                                       colors=colors)
    plt.title('Node Type Distribution')
   
    plt.legend(wedges, category_counts.keys(),
               title="Categories",
               loc="center left",
               bbox_to_anchor=(1, 0, 0.5, 1))
   
    plt.axis('equal')
    plt.show()