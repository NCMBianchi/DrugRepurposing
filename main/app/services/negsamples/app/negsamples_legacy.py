"""
NEGATIVE SAMPLES MODULE: GENERATES NEGATIVE TRIPLES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys
import os
import platform
import datetime
import logging
import builtins
import time
import multiprocessing
import json
import numpy as np
import pandas as pd
import networkx as nx

from .unique import unique_elements
# from .filepaths import initialise_disease_directories

# # These global variables will be set by the calling function
# global today_directory, date_str, input_seed, input_file_path
# global disease_name_label, disease_directories, base_data_directory
# disease_name_label = 'Huntington disease'  # default value
## comment out since they would not be used, but they are kept for future updates
## and if storing negative samples in their own file is ever necessary

def generate_negative_samples(positive_edges, similarity_threshold=0.90):
    """
    Generate valid negative samples based on the existing edges in the network.

    :param positive_edges: list of edges.
    :param similarity_threshold: minimum similarity score for drug-to-drug subnetwork.
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
    for edge in gene_to_gene_edges:
        all_genes.append(edge[0])
        all_genes.append(edge[2])
    for edge in drug_to_drug_edges:
        all_drugs.append(edge[0])
        all_drugs.append(edge[2])
    all_genes = unique_elements(all_genes)
    all_drugs = unique_elements(all_drugs)

    logging.info(f"Total genes: {len(all_genes)}")  # Debugging print
    logging.info(f"Total drugs: {len(all_drugs)}")  # Debugging print

    if not all_genes or not all_drugs:
        logging.info("No genes or drugs found to generate negative samples.")
        return []

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
                if sub_gene['id'] != gene['id'] or sub_drug['id'] != drug['id']:
                    valid_edge = ([sub_gene,
                                   {'label': 'biolink:valid_negative_association'},
                                   sub_drug,
                                   {'notes': f'related to: {gene["id"]}-{drug["id"]}'}])
                    if valid_edge not in val_neg_edges and valid_edge not in positive_edges:
                        val_neg_edges.append(valid_edge)

    val_neg_edges = unique_elements(val_neg_edges)

    return val_neg_edges
