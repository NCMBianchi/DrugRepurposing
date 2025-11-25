"""
CREATE_SMALL_NETWORK UTILITY MODULE: Generates a smaller representative network
Created on March 10th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os, random, logging, requests, json
import pandas as pd
import numpy as np
import networkx as nx

# import libraries for drug similarity calculation
from biothings_client import get_client
mc = get_client('chem')

from rdkit import Chem
from rdkit.Chem import DataStructs, MolFromSmiles, AllChem, rdFingerprintGenerator

from SPARQLWrapper import SPARQLWrapper, JSON
logging.getLogger('httpx').setLevel(logging.WARNING)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## functions from 'drugsimilarity.py'

def get_smiles(nodes):
    """
    This function performs drug ID conversion from chembl and wikidata ID to SMILES
    chemical structure notation.
    
    :param nodes: list of tuples (id, label) containing all drugs to convert.
    :return: concept dictionary (key = smiles structures: values = old ontology IDs).
    """
    
    chembl = [(id.split(':')[1], label) for id, label in nodes if 'chembl' in id.split(':')[0].lower()]
    wikidata = [(id.split(':')[1], label) for id, label in nodes if 'wikidata' in id.split(':')[0].lower()]

    concept_dct = {}
    duplicate_hits = []
    no_hit_ids = []

    # get SMILES for 'chembl:' drug nodes
    df_chembl = mc.querymany(qterms=[id for id, _ in chembl],
                             scopes=['chembl.molecule_chembl_id'], 
                             fields=['chembl.smiles'], 
                             as_dataframe=True,
                             verbose=False)
    ids_chembl = df_chembl.reset_index().copy()
    
    duplicate_hits = [qterm for qterm, count in df_chembl.index.value_counts().items() if count > 1]
    no_hit_ids = [id for id, _ in chembl if id not in df_chembl.index]

    if duplicate_hits:
        print(f"Duplicate ChEMBL IDs found: {len(duplicate_hits)}.")
    if no_hit_ids:
        print(f"ChEMBL IDs with no matches: {len(no_hit_ids)}.")
    
    if 'chembl.smiles' in ids_chembl.columns:
         for smile, (id, label) in zip(ids_chembl['chembl.smiles'], chembl):
            concept_dct[smile] = {'id': f'chembl:{id}', 'label': label}
    else:
        print('0 chembl IDs could be mapped to smiles.')
    
    # get SMILES for 'wikidata:' drug nodes
    smiles2wikidata = {}
    sparql = SPARQLWrapper("https://query.wikidata.org/sparql")
    for qid, label in wikidata:
        query = f"""
        SELECT ?smiles WHERE {{
        wd:{qid} wdt:P233 ?smiles.
        }}
        """
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        try:
            results = sparql.query().convert()
            if results["results"]["bindings"]:
                smiles = results["results"]["bindings"][0]["smiles"]["value"]
                concept_dct[smiles] = {'id': f'wikidata:{qid}', 'label': label}
            else:
                pass
        except Exception as e:
            #logging.warning(f'no wikidata IDs can be mapped to smiles: {e}')
            pass

    return concept_dct


def compute_similarity(smiles_dict, radius=2, length=4096):
    """
    This function computes the pairwise similarities of the provided list of SMILES
    notations, using Tanimoto coefficient.
    To achieve that, it first converts SMILES into RDKit objects, and then ECFP bitvectors.
    
    :param smiles_dict: the list of smiles.
    :param radius: ECFP fingerprint radius (default = 2, just for repetition).
    :param length: number of ECFP bits (default = 4096, just for repetition).
    :return: symmetric matrix of pairwise similarities. Diagonal is set to zero. 
    """
    
    # extrapolate actual SMILES notations from the concept dictionary
    smiles_list = list(smiles_dict.keys())
    valid_smiles = [
        sm for sm in smiles_list 
        if isinstance(sm, str) and sm.lower() != 'nan' and Chem.MolFromSmiles(sm) is not None
    ]

    # store IDs for the matrix
    id_list = list(smiles_dict.values())
    valid_ids = [
        id_list[smiles_list.index(sm)] 
        for sm in valid_smiles
    ]
    
    if len(smiles_list) > 10000:
        logging.warning(f'Calculating internal similarity on large set of SMILES strings ({len(smiles_list)})')
    
    # define the fingerprint list
    try:
        ## updated to the new RDKit format
        fingerprint_gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,  # Radius for the fingerprint
            countSimulation=False,  # Use counts instead of bits
            fpSize=length,  # Number of bits (replaces old length parameter)
            includeChirality=False,  # Consider stereochemistry
            useBondTypes=True  # Consider bond types
        )
    except (TypeError, ArgumentError):
        fingerprint_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius)

    fingerprint_list = []
    
    i = 0
    while i < len(smiles_list):
        sm = smiles_list[i]
        if pd.notna(sm):
            try:
                mol = Chem.MolFromSmiles(sm)
                if mol:
                    fingerprint = fingerprint_gen.GetFingerprint(mol)
                    fingerprint_list.append(fingerprint)
                    i += 1
                else:
                    smiles_list.pop(i)
                    id_list.pop(i)
            except Exception as e:
                #logging.error(f"Error processing SMILES {sm}: {e}")
                smiles_list.pop(i)
                id_list.pop(i)
        else:
            #logging.warning(f"Skipping NaN SMILES entry: {sm}")
            smiles_list.pop(i)
            id_list.pop(i)
    
    n_fingerprints = len(fingerprint_list)
    similarities = np.ones((n_fingerprints, n_fingerprints))

    for j in range(1, n_fingerprints):
        similarity = DataStructs.BulkTanimotoSimilarity(fingerprint_list[j], fingerprint_list[:j])
        similarities[j, :j] = similarity
        similarities[:j, j] = similarity

    return similarities, id_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def calculate_centrality(graph, method="eigen"):
    """
    Calculate centrality using different methods with fallbacks.
    
    :param graph: NetworkX graph.
    :parm method: centrality method ('eigen', 'deg', or 'betw').
    :return: dictionary of node centrality values.
    """
    if method == "eigen":
        try:
            # try eigenvector centrality with increased iterations
            centrality = nx.eigenvector_centrality_numpy(graph, max_iter=1000)
        except Exception as e:
            print(f"Eigenvector centrality failed with error: {str(e)}")
            # fallback to PageRank (more robust than eigenvector)
            try:
                print("Trying PageRank instead...")
                centrality = nx.pagerank(graph)
            except Exception as e:
                print(f"PageRank failed with error: {str(e)}")
                # final fallback to degree centrality
                print("Using degree centrality as final fallback")
                centrality = nx.degree_centrality(graph)
    
    elif method == "deg":
        centrality = nx.degree_centrality(graph)
    
    elif method == "betw":
        try:
            # use approximation for large graphs
            if len(graph) > 1000:
                centrality = nx.betweenness_centrality(graph, k=100)  # sample 100 nodes
            else:
                centrality = nx.betweenness_centrality(graph)
        except Exception as e:
            print(f"Betweenness centrality failed with error: {str(e)}")
            print("Using degree centrality instead")
            centrality = nx.degree_centrality(graph)
    
    return centrality


def save_small_network(small_network, output_dir, disease_name, date_string):
    """
    Save the smaller network to disk.
    
    :param small_network: dictionary with small network data.
    :param output_dir: directory to save files.
    :param disease_name: name of the disease.
    :param date_string: date string for filenames.
    :return: dictionary with paths to saved files.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # convert to DataFrames for saving
    nodes_df = pd.DataFrame([{
        'id': node['id'],
        'label': node['label']
    } for node in small_network['nodes']])
    
    edges_df = pd.DataFrame([{
        'subject_id': edge[0]['id'],
        'subject_label': edge[0]['label'],
        'relation': edge[1]['label'],
        'object_id': edge[2]['id'],
        'object_label': edge[2]['label'],
        'notes': edge[3]['notes']
    } for edge in small_network['edges']])
    
    drug_nodes_df = pd.DataFrame([{
        'id': node['id'],
        'label': node['label']
    } for node in small_network['all_drugs']])
    
    # save files
    nodes_path = os.path.join(output_dir, f"{disease_name}_{date_string}_small_dgidb_nodes.csv")
    edges_path = os.path.join(output_dir, f"{disease_name}_{date_string}_small_dgidb_edges.csv")
    drugs_path = os.path.join(output_dir, f"dgidb_all_drugs_{date_string}_small_api_response.csv")
    
    nodes_df.to_csv(nodes_path, index=False)
    edges_df.to_csv(edges_path, index=False)
    drug_nodes_df.to_csv(drugs_path, index=False)
    
    print(f"Saved small network files.")
    
    return {
        'nodes_path': nodes_path,
        'edges_path': edges_path,
        'drugs_path': drugs_path
    }


def print_network_statistics(original_nodes, original_edges, small_nodes, small_edges):
    """
    Print comparative statistics for the original and small networks.
    
    :param original_nodes: original network nodes.
    :param original_edges: original network edges.
    :param small_nodes: small network nodes.
    :param small_edges: small network edges.
    """
    # classify nodes by type
    def classify_nodes(nodes):
        genes = []
        drugs = []
        other = []
        
        for node in nodes:
            node_id_prefix = node['id'].split(':')[0].lower()
            
            if node_id_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase']:
                genes.append(node)
            elif node_id_prefix in ['chembl', 'wikidata']:
                drugs.append(node)
            else:
                other.append(node)
                
        return genes, drugs, other
    
    # classify edges by type
    def classify_edges(edges):
        gene_drug = []
        gene_gene = []
        drug_drug = []
        other = []
        
        for edge in edges:
            subj_prefix = edge[0]['id'].split(':')[0].lower()
            obj_prefix = edge[2]['id'].split(':')[0].lower()
            
            if (subj_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase'] and 
                obj_prefix in ['chembl', 'wikidata']):
                gene_drug.append(edge)
            elif (obj_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase'] and 
                  subj_prefix in ['chembl', 'wikidata']):
                gene_drug.append(edge)
            elif (subj_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase'] and 
                  obj_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase']):
                gene_gene.append(edge)
            elif (subj_prefix in ['chembl', 'wikidata'] and 
                  obj_prefix in ['chembl', 'wikidata']):
                drug_drug.append(edge)
            else:
                other.append(edge)
                
        return gene_drug, gene_gene, drug_drug, other
    
    # get statistics for original network
    orig_genes, orig_drugs, orig_other = classify_nodes(original_nodes)
    orig_gene_drug, orig_gene_gene, orig_drug_drug, orig_other_edges = classify_edges(original_edges)
    
    # get statistics for small network
    small_genes, small_drugs, small_other = classify_nodes(small_nodes)
    small_gene_drug, small_gene_gene, small_drug_drug, small_other_edges = classify_edges(small_edges)
    
    # print summary
    print("="*80)
    print(f"NETWORK STATISTICS COMPARISON")
    print("="*80)
    print(f"{'':20} {'ORIGINAL':>15} {'SMALL':>15} {'PERCENTAGE':>15}")
    print("-"*80)
    print(f"{'NODES':}")
    print(f"{'  Genes':20} {len(orig_genes):15d} {len(small_genes):15d} {len(small_genes)/max(1,len(orig_genes))*100:15.2f}%")
    print(f"{'  Drugs':20} {len(orig_drugs):15d} {len(small_drugs):15d} {len(small_drugs)/max(1,len(orig_drugs))*100:15.2f}%")
    print(f"{'  Other':20} {len(orig_other):15d} {len(small_other):15d} {len(small_other)/max(1,len(orig_other))*100:15.2f}%")
    print(f"{'  Total':20} {len(original_nodes):15d} {len(small_nodes):15d} {len(small_nodes)/max(1,len(original_nodes))*100:15.2f}%")
    print("-"*80)
    print(f"{'EDGES':}")
    print(f"{'  Gene-Drug':20} {len(orig_gene_drug):15d} {len(small_gene_drug):15d} {len(small_gene_drug)/max(1,len(orig_gene_drug))*100:15.2f}%")
    print(f"{'  Gene-Gene':20} {len(orig_gene_gene):15d} {len(small_gene_gene):15d} {len(small_gene_gene)/max(1,len(orig_gene_gene))*100:15.2f}%")
    print(f"{'  Drug-Drug':20} {len(orig_drug_drug):15d} {len(small_drug_drug):15d} {len(small_drug_drug)/max(1,len(orig_drug_drug))*100:15.2f}%")
    print(f"{'  Other':20} {len(orig_other_edges):15d} {len(small_other_edges):15d} {len(small_other_edges)/max(1,len(orig_other_edges))*100:15.2f}%")
    print(f"{'  Total':20} {len(original_edges):15d} {len(small_edges):15d} {len(small_edges)/max(1,len(original_edges))*100:15.2f}%")
    print("="*80)


def sample_drugs_by_similarity(selected_drugs, all_drugs, num_to_sample, similarity_threshold=0.5):
    """
    Sample additional drugs based on similarity to the already selected drugs.
    
    :param selected_drugs: list of drug nodes already selected for the small network.
    :param all_drugs: list of all drug nodes to sample from.
    :param num_to_sample: number of additional drugs to sample.
    :param similarity_threshold: minimum similarity score to consider.
    :return: list of additional drug nodes selected based on similarity.
    """
    print("\nSampling additional drugs based on chemical similarity...")
    
    # extract drug IDs and labels for SMILES conversion
    selected_drug_tuples = [(drug['id'], drug['label']) for drug in selected_drugs]
    all_drug_tuples = [(drug['id'], drug['label']) for drug in all_drugs]
    
    # get SMILES notations for selected drugs
    selected_smiles = get_smiles(selected_drug_tuples)
    
    # if we couldn't get SMILES for selected drugs, fall back to random sampling
    if not selected_smiles:
        print("Warning: Could not get SMILES for selected drugs, falling back to random sampling")
        all_drug_ids = set(drug['id'] for drug in all_drugs)
        selected_drug_ids = set(drug['id'] for drug in selected_drugs)
        remaining_drug_ids = list(all_drug_ids - selected_drug_ids)
        
        if len(remaining_drug_ids) <= num_to_sample:
            return [drug for drug in all_drugs if drug['id'] in remaining_drug_ids]
        else:
            sampled_ids = random.sample(remaining_drug_ids, num_to_sample)
            return [drug for drug in all_drugs if drug['id'] in sampled_ids]
    
    # get SMILES notations for all drugs
    batch_size = 2000  # Process in batches to avoid memory issues
    all_smiles = {}
    
    # find all drugs that aren't in selected_drugs
    selected_drug_ids = set(drug['id'] for drug in selected_drugs)
    remaining_drugs = [drug for drug in all_drugs if drug['id'] not in selected_drug_ids]
    
    print(f"Computing similarity for {len(remaining_drugs)} remaining drugs...")
    
    # process in batches
    for i in range(0, len(remaining_drugs), batch_size):
        batch = remaining_drugs[i:i+batch_size]
        batch_tuples = [(drug['id'], drug['label']) for drug in batch]
        batch_smiles = get_smiles(batch_tuples)
        all_smiles.update(batch_smiles)
        print(f"Processed batch {i//batch_size + 1}/{(len(remaining_drugs) + batch_size - 1)//batch_size} ({len(batch_smiles)} SMILES)")
    
    print(f"Got SMILES for {len(all_smiles)} drugs out of {len(remaining_drugs)} remaining drugs")
    
    if not all_smiles:
        print("Warning: Could not get SMILES for any drugs, falling back to random sampling")
        if len(remaining_drugs) <= num_to_sample:
            return remaining_drugs
        else:
            return random.sample(remaining_drugs, num_to_sample)
    
    # compute similarity scores
    max_similarity_scores = {}
    
    # for each drug in all_smiles, compute its maximum similarity to any drug in selected_smiles
    for smiles, drug_info in all_smiles.items():
        max_sim = 0
        for sel_smiles, sel_drug_info in selected_smiles.items():
            try:
                mol1 = Chem.MolFromSmiles(smiles)
                mol2 = Chem.MolFromSmiles(sel_smiles)
                
                if mol1 and mol2:
                    fp1 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096).GetFingerprint(mol1)
                    fp2 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096).GetFingerprint(mol2)
                    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                    max_sim = max(max_sim, sim)
            except Exception as e:
                continue
        
        # store the maximum similarity score for this drug
        max_similarity_scores[drug_info['id']] = max_sim
    
    # sort drugs by similarity score
    sorted_drugs = sorted([(drug_id, score) for drug_id, score in max_similarity_scores.items()], 
                          key=lambda x: x[1], reverse=True)
    
    # apply similarity threshold and select top drugs
    filtered_drugs = [(drug_id, score) for drug_id, score in sorted_drugs if score >= similarity_threshold]
    
    # ff we don't have enough drugs after filtering, relax the threshold
    if len(filtered_drugs) < num_to_sample:
        print(f"Warning: Only {len(filtered_drugs)} drugs meet similarity threshold {similarity_threshold}, using all available similar drugs")
        # Take all available similar drugs, sorted by similarity
        filtered_drugs = sorted_drugs
    
    # take top N drugs
    top_drug_ids = [drug_id for drug_id, _ in filtered_drugs[:num_to_sample]]
    
    # get complete drug objects
    sampled_drugs = [drug for drug in all_drugs if drug['id'] in top_drug_ids]
    
    # if we still don't have enough, add randomly sampled drugs
    if len(sampled_drugs) < num_to_sample:
        remaining_ids = set(drug['id'] for drug in remaining_drugs) - set(top_drug_ids)
        additional_needed = num_to_sample - len(sampled_drugs)
        additional_ids = random.sample(list(remaining_ids), min(additional_needed, len(remaining_ids)))
        additional_drugs = [drug for drug in all_drugs if drug['id'] in additional_ids]
        sampled_drugs.extend(additional_drugs)
    
    print(f"Selected {len(sampled_drugs)} additional drugs based on similarity")
    
    return sampled_drugs


def create_small_network(dgidb_nodes, dgidb_edges, drug_nodes, gene_count=500, drug_count=300, all_drug_count=5000, seed=42):
    """
    Create a smaller representative network while maintaining centrality distributions.
    
    :param dgidb_nodes: list of all nodes from DGIdb.
    :param dgidb_edges: list of all edges from DGIdb.
    :param drug_nodes: list of all drug nodes from DGIdb.
    :param gene_count: number of genes to include.
    :param drug_count: number of drugs from the disease network to include.
    :param all_drug_count: total number of drugs to include.
    :param seed: random seed for reproducibility.
    :return: dictionary with the smaller network data.
    """
    random.seed(seed)
    np.random.seed(seed)
    
    print("Creating smaller representative network...")
    
    # create NetworkX graph for centrality calculations
    G = nx.Graph()
    
    # add all nodes
    for node in dgidb_nodes:
        G.add_node(node['id'], label=node['label'])
    
    # add all edges
    for edge in dgidb_edges:
        G.add_edge(edge[0]['id'], edge[2]['id'])
    
    print(f"Built graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    # classify nodes
    genes = []
    drugs = []
    other_nodes = []
    
    for node in dgidb_nodes:
        node_id_prefix = node['id'].split(':')[0].lower()
        
        if node_id_prefix in ['hgnc', 'mgi', 'go', 'ncbigene', 'zfin', 'xenbase']:
            genes.append(node)
        elif node_id_prefix in ['chembl', 'wikidata']:
            drugs.append(node)
        else:
            other_nodes.append(node)
    
    print(f"Original network: {len(genes)} genes, {len(drugs)} drugs, {len(other_nodes)} other nodes")
    
    # calculate GENES CENTRALITY
    print("Calculating gene centrality measures...")
    gene_ids = [gene['id'] for gene in genes]
    gene_subgraph = G.subgraph(gene_ids)
    connected_components = list(nx.connected_components(gene_subgraph))
    print(f"Gene graph has {len(connected_components)} connected components")
    if len(connected_components) > 1:
        largest_component_size = len(max(connected_components, key=len))
        print(f"Largest component contains {largest_component_size} nodes ({largest_component_size/len(gene_subgraph.nodes())*100:.2f}% of graph)")
    gene_centrality = calculate_centrality(gene_subgraph, method="eigen")
    
    # calculate DRUGS CENTRALITY
    print("Calculating drug centrality measures...")
    drug_ids = [drug['id'] for drug in drugs]
    drug_subgraph = G.subgraph(drug_ids)
    connected_components = list(nx.connected_components(drug_subgraph))
    print(f"Drug graph has {len(connected_components)} connected components")
    if len(connected_components) > 1:
        largest_component_size = len(max(connected_components, key=len))
        print(f"Largest component contains {largest_component_size} nodes ({largest_component_size/len(drug_subgraph.nodes())*100:.2f}% of graph)")
    drug_centrality = calculate_centrality(drug_subgraph, method="eigen")
    
    # create percentile bins for genes
    gene_df = pd.DataFrame({
        'id': gene_ids,
        'centrality': [gene_centrality.get(id, 0) for id in gene_ids]
    })
    try:
        gene_df['percentile'] = pd.qcut(gene_df['centrality'], 10, labels=False, duplicates='drop')
    except:
        # Handle the case where there are too many ties
        print("Warning: Too many ties in gene centrality - using quantile binning instead")
        gene_df['percentile'] = pd.cut(gene_df['centrality'], 5, labels=False)
    
    # create percentile bins for drugs
    drug_df = pd.DataFrame({
        'id': drug_ids,
        'centrality': [drug_centrality.get(id, 0) for id in drug_ids]
    })
    try:
        drug_df['percentile'] = pd.qcut(drug_df['centrality'], 10, labels=False, duplicates='drop')
    except:
        # handle the case where there are too many ties
        print("Warning: Too many ties in drug centrality - using quantile binning instead")
        drug_df['percentile'] = pd.cut(drug_df['centrality'], 5, labels=False)
    
    # sample from each percentile bin to maintain distribution
    print(f"Sampling {gene_count} genes and {drug_count} drugs based on centrality distribution...")
    
    sampled_genes = []
    for percentile in sorted(gene_df['percentile'].unique()):
        percentile_genes = gene_df[gene_df['percentile'] == percentile]['id'].tolist()
        if not percentile_genes:
            continue
            
        # calculate proportional sample size
        orig_count = len(percentile_genes)
        sample_size = max(1, int(round(gene_count * (orig_count / len(genes)))))
        sample_size = min(sample_size, orig_count)
        
        # sample from this percentile
        percentile_sample = random.sample(percentile_genes, sample_size)
        sampled_genes.extend(percentile_sample)
    
    # adjust to match target gene_count
    if len(sampled_genes) > gene_count:
        sampled_genes = random.sample(sampled_genes, gene_count)
    elif len(sampled_genes) < gene_count and len(sampled_genes) < len(genes):
        # add more genes if needed
        additional_needed = gene_count - len(sampled_genes)
        remaining_genes = [g for g in gene_ids if g not in sampled_genes]
        if remaining_genes:
            additional_genes = random.sample(remaining_genes, min(additional_needed, len(remaining_genes)))
            sampled_genes.extend(additional_genes)
    
    # sample drugs using the same approach
    sampled_drugs = []
    for percentile in sorted(drug_df['percentile'].unique()):
        percentile_drugs = drug_df[drug_df['percentile'] == percentile]['id'].tolist()
        if not percentile_drugs:
            continue
            
        # calculate proportional sample size
        orig_count = len(percentile_drugs)
        sample_size = max(1, int(round(drug_count * (orig_count / len(drugs)))))
        sample_size = min(sample_size, orig_count)
        
        # sample from this percentile
        percentile_sample = random.sample(percentile_drugs, sample_size)
        sampled_drugs.extend(percentile_sample)
    
    # adjust to match target drug_count
    if len(sampled_drugs) > drug_count:
        sampled_drugs = random.sample(sampled_drugs, drug_count)
    elif len(sampled_drugs) < drug_count and len(sampled_drugs) < len(drugs):
        # add more drugs if needed
        additional_needed = drug_count - len(sampled_drugs)
        remaining_drugs = [d for d in drug_ids if d not in sampled_drugs]
        if remaining_drugs:
            additional_drugs = random.sample(remaining_drugs, min(additional_needed, len(remaining_drugs)))
            sampled_drugs.extend(additional_drugs)
    
    # sample 10% of other nodes as well
    other_count = max(1, int(round(len(other_nodes) * 0.1)))  # target ~10% of other nodes
    if len(other_nodes) > other_count:
        small_other_nodes = random.sample(other_nodes, other_count)
    else:
        small_other_nodes = other_nodes
        
    print(f"Selected {len(small_other_nodes)} other nodes out of {len(other_nodes)}")
    
    # get the full node objects
    small_genes = [gene for gene in genes if gene['id'] in sampled_genes]
    small_drugs = [drug for drug in drugs if drug['id'] in sampled_drugs]
    
    print(f"Selected {len(small_genes)} genes and {len(small_drugs)} drugs based on centrality")
    
    # sample additional drugs based on similarity to the selected drugs
    sampled_drug_ids = set(sampled_drugs)
    all_drug_ids = set(drug['id'] for drug in drug_nodes)
    
    # only sample from drugs that aren't already in the small network
    remaining_drug_count = all_drug_count - len(sampled_drug_ids)
    
    if remaining_drug_count > 0:
        # use similarity-based sampling for additional drugs
        additional_drugs = sample_drugs_by_similarity(
            small_drugs, 
            drug_nodes, 
            remaining_drug_count,
            similarity_threshold=0.3  # lower threshold for more diverse sampling
        )
        
        all_sampled_drugs = small_drugs + additional_drugs
    else:
        all_sampled_drugs = small_drugs
    
    # create the small network nodes
    small_network_nodes = small_genes + small_drugs + small_other_nodes
    
    # filter edges to only include nodes in the small network
    small_network_node_ids = set(node['id'] for node in small_network_nodes)
    small_network_edges = [
        edge for edge in dgidb_edges 
        if edge[0]['id'] in small_network_node_ids and edge[2]['id'] in small_network_node_ids
    ]
    
    print(f"Small network: {len(small_genes)} genes, {len(small_drugs)} drugs, {len(small_other_nodes)} other nodes, {len(small_network_edges)} edges")
    print(f"Total drugs (including additional): {len(all_sampled_drugs)}")
    
    return {
        'nodes': small_network_nodes,
        'edges': small_network_edges,
        'genes': small_genes,
        'drugs': small_drugs,
        'other_nodes': small_other_nodes,
        'all_drugs': all_sampled_drugs
    }