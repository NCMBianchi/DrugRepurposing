"""
DRUG SIMILARITY MODULE: COMPUTE D-to-D SIMILARITY, AND GENERATE EDGES
Updated on February 10th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,inspect,requests,json,logging
import pandas as pd
import numpy as np

from biothings_client import get_client
mc = get_client('chem')

from rdkit import Chem
from rdkit.Chem import DataStructs, MolFromSmiles, AllChem, rdFingerprintGenerator

from SPARQLWrapper import SPARQLWrapper, JSON
logging.getLogger('httpx').setLevel(logging.WARNING)

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


def compute_similarity(smiles_dict,radius=2,length=4096):
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


def generate_drug_edges(similarity_matrix,sorted_matrix,ids,K,minimum=0.5):
    '''
    This function builds the edges object from the similarity matrix.

    :param similarity_matrix: similarity matrix from compute_similarity().
    :param sorted_matrix: result of KNN on similarity matrix .
    :param ids: list of ids from compute_similarity().
    :param K: number of top scoring drugs based on similarity.
    :param min_simil: minimum similarity score (default = 0.5, just for repetition).
    :return: drug-to-drug edges.
    '''
    
    similarity = []
    added_edges_count = 0
    skipped_edges_count = 0
    chunk_size = 2000
    num_chunks = (similarity_matrix.shape[0] // chunk_size) + 1

    for chunk in range(num_chunks):
        start_idx = chunk * chunk_size
        end_idx = min((chunk + 1) * chunk_size, similarity_matrix.shape[0])
        
        sub_similarity_matrix = similarity_matrix[start_idx:end_idx, :]
        sub_sorted_matrix = sorted_matrix[start_idx:end_idx, :]
        sub_ids = ids[start_idx:end_idx]

        # further update to remove all ==1 and if so, raise K+ how many times ==1
        for i in range(sub_similarity_matrix.shape[0]):
            global_i = start_idx + i  # Adjust index for the original matrix
            for j in sub_sorted_matrix[i, :K]:
                # filter results below the 'min_simil' threshold (default = 0.5)
                # and those at 0 even if 'min_simil == 0.0', as well as any isoform with == 1.0 other than identity
                if sub_similarity_matrix[i, j] >= minimum and sub_similarity_matrix[i, j] != 0 and sub_similarity_matrix[i, j] != 1:
                    new_row = [
                        {'id': ids[global_i]['id'], 'label': ids[global_i]['label']},
                        {'label': 'smiles: similar to'},
                        {'id': ids[j]['id'], 'label': ids[j]['label']},
                        {'notes': f'similarity score: {sub_similarity_matrix[i, j]}'}
                    ]
                    similarity.append(new_row)
    
    return similarity
    

def run_drugsimilarity(nodes,edges,d_nodes,disease_directories,
                       K=10,min_simil=None,input_radius=None,input_length=None, simil_load=0):
    '''
    This function runs the whole drug_similarioty script and saves nodes and edges files.

    :param nodes: all the node resulting from the previous step in the pipeline (i.e. 'run_dgidb()').
    :param edges: all the edges resulting from the previous step in the pipeline (i.e. 'run_dgidb()').
    :param disease_directories: base paths to where data is stored.
    :param K: number of top scoring drugs based on similarity (default = 10).
    :param min_simil: any value that would override default minimum=0.5 in compute_similarity().
    :param input_radius: any value that would override default radius=2 in compute_similarity().
    :param input_length: any value that would override default length=4096 in compute_similarity().
    :param simil_load: toggle for loading existing files (1) or generating new ones (0).
    :return: lists of edges with drug-to-drug associations based on feature similarity.
    '''

    start_time = time.time()

    print(f"NOW RUNNING: {current_function_name()} following 'run_dgidb()'.")
    
    if not input_radius:
        input_radius = 2
    if not input_length:
        input_length = 4096
    if not min_simil:
        min_simil = 0.5

    print(f"The ({K}) top scoring drugs based on similarity are considered.")
    print(f"The minimum similarity threshold is set to: {min_simil}.")
    print(f"The input radius of the ECFP ({input_length} features) is set to: {input_radius}.")

    # initialise path
    drugsimil_directory = disease_directories['drugsimil_directory']
    date_str = disease_directories['date_string']
    disease_name_label = disease_directories['disease_name']

    # define path for output files
    smiles_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_drugsim_smiles.csv')
    alldrugs_smiles_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_alldrugsim_smiles.csv')
    similarities_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_drugsim_similarityMatrix.csv')
    alldrug_simil_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_alldrugsim_similarityMatrix.csv')
    edges_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_drugsim_edges.csv')
    drug_edges_path = os.path.join(drugsimil_directory, f'{disease_name_label}_{date_str}_alldrugsim_edges.csv')

    # check if output files exist and load them if simil_toggle is set to 1
    if simil_load == 1 and all(os.path.exists(f) and os.path.getsize(f) > 0 for f in 
                                [smiles_path, alldrugs_smiles_path, similarities_path, 
                                 alldrug_simil_path, edges_path, drug_edges_path]):
        
        # load SMILES data
        smiles_df = pd.read_csv(smiles_path)
        smiles = {row['SMILES']: {'id': row['ID'], 'label': row['Label']} 
                  for _, row in smiles_df.iterrows() if pd.notna(row['SMILES'])}
        
        alldrugs_smiles_df = pd.read_csv(alldrugs_smiles_path)
        alldrug_smiles = {row['SMILES']: {'id': row['ID'], 'label': row['Label']} 
                         for _, row in alldrugs_smiles_df.iterrows() if pd.notna(row['SMILES'])}
        
        # load similarity matrices
        similarities_df = pd.read_csv(similarities_path, index_col=0)
        similarities = similarities_df.values
        ids = similarities_df.index.tolist()
        id_list = [{'id': id_str, 'label': next((node['label'] for node in nodes if node['id'] == id_str), id_str)} 
                  for id_str in ids]
        
        alldrug_simil_df = pd.read_csv(alldrug_simil_path, index_col=0)
        alldrug_simil = alldrug_simil_df.values
        alldrug_ids = alldrug_simil_df.index.tolist()
        alldrug_id_list = [{'id': id_str, 'label': next((drug['label'] for drug in d_nodes if drug['id'] == id_str), id_str)} 
                          for id_str in alldrug_ids]

    else:
        # convert drug IDs into SMILES notations
        nodes_list = [(node['id'], node['label']) for node in nodes]
        alldrug_list = [(drug_node['id'], drug_node['label']) for drug_node in d_nodes]
        smiles = get_smiles(nodes_list)
        alldrug_smiles = get_smiles(alldrug_list)
        
        # save SMILES notations as CSV
        smiles_data = [{'SMILES': k, 'ID': v['id'], 'Label': v['label']} for k, v in smiles.items()]
        smiles_df = pd.DataFrame(smiles_data)
        smiles_df.to_csv(smiles_path, index=False)
        alldrugs_smiles_data = [{'SMILES': k, 'ID': v['id'], 'Label': v['label']} for k, v in alldrug_smiles.items()]
        alldrugs_smiles_df = pd.DataFrame(alldrugs_smiles_data)
        alldrugs_smiles_df.to_csv(alldrugs_smiles_path, index=False)
    
        # compute similarities
        similarities, id_list = compute_similarity(smiles,radius=input_radius,length=input_length)
        alldrug_simil, alldrug_id_list = compute_similarity(alldrug_smiles,radius=input_radius,length=input_length)
        
        # save the similarity matrix as CSV
        ids = [d['id'] for d in id_list]
        similarities_df = pd.DataFrame(similarities, index=ids, columns=ids)
        similarities_df.to_csv(similarities_path, index=True)
        alldrug_ids = [d['id'] for d in alldrug_id_list]
        alldrug_simil_df = pd.DataFrame(alldrug_simil, index=alldrug_ids, columns=alldrug_ids)
        alldrug_simil_df.to_csv(alldrug_simil_path, index=True)
    
    # sort based on similarity scores (ignore identity, first element for any sorting)
    sortedSimilarities = np.argsort(-similarities, axis=1)[:, 1:K+1]
    alldrug_sorted = np.argsort(-alldrug_simil, axis=1)[:, 1:K+1]
    
    # generate drug-to-drug edges
    drug_similarity = generate_drug_edges(similarities,sortedSimilarities,id_list,K,minimum=min_simil)
    drug_edges = generate_drug_edges(alldrug_simil,alldrug_sorted,alldrug_id_list,K,minimum=min_simil)

    # MERGE new edges with the existing ones
    all_edges = edges + drug_similarity
    unique_edges = unique_elements(all_edges)
    unique_drug_edges = unique_elements(drug_edges)

    # save the unique edges as CSV
    edges_df = pd.DataFrame(unique_edges)
    edges_df.to_csv(edges_path, index=False)
    drug_edges_df = pd.DataFrame(drug_edges)
    drug_edges_df.to_csv(drug_edges_path, index=False)

    print(f"DRUG-TO-DRUG EDGES: {len(drug_similarity)} associated to the disease of interest, {len(drug_edges)} in total in DGIdb.")
    
    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    formatted_duration = format_duration(duration)  # convert for print
    print(f"'drugsimilarity.py' run finished in {formatted_duration}.")

    return unique_edges, unique_drug_edges