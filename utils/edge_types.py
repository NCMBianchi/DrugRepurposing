"""
EDGE_TYPES UTILITY MODULE: PRINT THE EDGE TYPES
Created on February 14th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def report_edge_types(input_edges,prnt=True):
    """
    Prints table for edge categories.
    
    :param edges: list of edges where each edge is [subject, relation, object, notes].
        :param prnt: Whether to print the report: used to pass the computed stats
        to plot_edge_types() without printing the table when prnt=False and instead
        returning values.
    :return: none, or list counts for plot_edge_types().
    """
    edge_type_dict = {
        'drug-to-drug': ['smiles: similar to'],
        'gene-to-drug': ['dgidb:interacts_with'],
        'disease-to-gene': ['biolink:gene_associated_with_condition'],
        'gene-to-gene': ['biolink:interacts_with']
    }
    
    category_counts = {cat: 0 for cat in list(edge_type_dict.keys()) + ['-']}
    
    for edge in input_edges:
        if 'label' in edge[1]:
            rel_type = edge[1]['label']
            category_found = False
            
            for category, rel_types in edge_type_dict.items():
                if any(rel_type.startswith(rt) for rt in rel_types):
                    category_counts[category] += 1
                    category_found = True
                    break
                    
            if not category_found:
                category_counts['-'] += 1
    
    if prnt:
        print("\nEDGE Category Counts:")
        print("-" * 30)
        print(f"{'Category':<20} {'Count':<10}")
        print("-" * 30)
        for category, count in sorted(category_counts.items()):
            if count > 0:  #only show categories with edges
                print(f"{category:<20} {count:<10}")
    else:
        return category_counts


def plot_edge_types(input_edges):
    """
    Displays a pie chart for edge categories.
    
    :param nodes: list of edges where each edge is [subject, relation, object, notes].
    :return: none.
    """
    category_counts = report_edge_types(input_edges,prnt=False)
    colors = ['#f9dbbd', '#ff7f51', '#ce4257', '#720026','#450920']
   
    plt.figure(figsize=(10, 8))
    wedges, texts, autotexts = plt.pie(category_counts.values(),
                                       labels=[''] * len(category_counts),
                                       autopct='%1.1f%%',
                                       colors=colors)
    plt.title('Edge Type Distribution')
   
    plt.legend(wedges, category_counts.keys(),
               title="Categories",
               loc="center left",
               bbox_to_anchor=(1, 0, 0.5, 1))
   
    plt.axis('equal')
    plt.show()

