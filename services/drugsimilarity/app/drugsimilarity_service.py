"""
DRUG SIMILARITY MODULE: COMPUTE D-to-D SIMILARITY, AND GENERATE EDGES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,platform,datetime,logging,builtins,time,multiprocessing

from drugapp.filepaths import initialise_base_directories, setup_service_globals
from drugapp.unique import unique_elements

def get_smiles(nodes):
    """
    This function performs drug ID conversion from chembl and wikidata ID to SMILES
    chemical structure notation.
    
    :param nodes: list of tuples (id, label) containing all drugs to convert.
    :return: concept dictionary (key = smiles structures: values = old ontology IDs).
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    chembl = [(id.split(':')[1], label) for id, label in nodes if 'chembl' in id.split(':')[0].lower()]
    wikidata = [(id.split(':')[1], label) for id, label in nodes if 'wikidata' in id.split(':')[0].lower()]

    concept_dct = {}

    # get SMILES for 'chembl:' drug nodes
    df_chembl = mc.querymany(qterms=[id for id, _ in chembl],
                             scopes=['chembl.molecule_chembl_id'], size=1,
                             fields=['chembl.smiles'], as_dataframe=True)
    ids_chembl = df_chembl.reset_index().copy()
    if 'chembl.smiles' in ids_chembl.columns:
         for smile, (id, label) in zip(ids_chembl['chembl.smiles'], chembl):
            concept_dct[smile] = {'id': f'chembl:{id}', 'label': label}
    else:
        logging.warning('no chembl IDs can be mapped to smiles')
    
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
                logging.info(f'Wikidata ID {qid} cannot be mapped to smiles')
        except Exception as e:
            logging.warning('no wikidata IDs can be mapped to smiles: {e}')

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

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    # extrapolate actual SMILES notations from the concept dictionary
    smiles_list = list(smiles_dict.keys())

    # store IDs for the matrix
    id_list = list(smiles_dict.values())
    
    if len(smiles_list) > 10000:
        logging.warning(f'Calculating internal similarity on large set of SMILES strings ({len(smiles_list)})')
    
    # define the fingerprint generator parameters
    fingerprint_gen = rdFingerprintGenerator.GetMorganGenerator(fpradius=radius,fpSize=length)

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
                    i += 1  # only increment the index if the item is not removed
                else:
                    # if RDKit fails to create a molecule, but the SMILES is not NaN
                    smiles_list.pop(i)
                    id_list.pop(i)
            except Exception as e:
                logging.error(f"Error processing SMILES {sm}: {e}")
                smiles_list.pop(i)
                id_list.pop(i)
        else:
            # if SMILES is NaN
            logging.warning(f"Skipping NaN SMILES entry: {sm}")
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
    :return: drug-to-drug edges
    '''

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
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
    

def run_drugsimilarity(monarch_input,date,K=10,min_simil=None,input_radius=None,input_length=None):
    '''
    This function runs the whole drug_similarioty script and saves nodes and edges files.

    :param monarch_input: the input seed from the run_monarch() step.
    :param date: the date of creation of the disease graph.
    :param K: number of top scoring drugs based on similarity (default = 10).
    :param min_simil: any value that would override default minimum=0.5 in compute_similarity().
    :param input_radius: any value that would override default radius=2 in compute_similarity().
    :param input_length: any value that would override default length=4096 in compute_similarity().
    :return: list of edges with drug-to-drug associations.
    '''

    start_time = time.time()
    
    logging.info(f"NOW RUNNING: {current_function_name()} following 'run_dgidb({monarch_input},{today},layers={run_layers})'.")

    global nodes
    global edges
    global drug_nodes

    if not input_radius:
        input_radius = 2
    if not input_length:
        input_length = 4096
    if not min_simil:
        min_simil = 0.5

    with open(input_file_path, 'a') as file:
        file.write(f"The ({K}) top scoring drugs based on similarity are considered.\n")
        file.write(f"The minimum similarity threshold is set to: {min_simil}\n")
        file.write(f"The input radius of the ECFP ({input_length} features) is set to: {min_simil}\n\n")

    drugSimil_directory = disease_directories['drugsimil_directory']
    
    # convert drug IDs into SMILES notations
    nodes_list = [(node['id'], node['label']) for node in nodes]
    alldrug_list = [(drug_node['id'], drug_node['label']) for drug_node in drug_nodes]
    smiles = get_smiles(nodes_list)
    alldrug_smiles = get_smiles(alldrug_list)
    
    # save SMILES notations as CSV
    smiles_data = [{'SMILES': k, 'ID': v['id'], 'Label': v['label']} for k, v in smiles.items()]
    smiles_df = pd.DataFrame(smiles_data)
    smiles_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date}_drugSim_smiles.csv'), index=False)
    alldrugs_smiles_data = [{'SMILES': k, 'ID': v['id'], 'Label': v['label']} for k, v in alldrug_smiles.items()]
    alldrugs_smiles_df = pd.DataFrame(alldrugs_smiles_data)
    alldrugs_smiles_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date}_alldrugSim_smiles.csv'), index=False)

    # compute similarities
    similarities, id_list = compute_similarity(smiles,radius=input_radius,length=input_length)
    alldrug_simil, alldrug_id_list = compute_similarity(alldrug_smiles,radius=input_radius,length=input_length)
    
    # save the similarity matrix as CSV
    ids = [d['id'] for d in id_list]
    similarities_df = pd.DataFrame(similarities, index=ids, columns=ids)
    similarities_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date}_drugSim_similarityMatrix.csv'), index=True)
    alldrug_ids = [d['id'] for d in alldrug_id_list]
    alldrug_simil_df = pd.DataFrame(alldrug_simil, index=alldrug_ids, columns=alldrug_ids)
    alldrug_simil_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date}_alldrugSim_similarityMatrix.csv'), index=True)
    
    # sort based on similarity scores (ignore identity, first element for any sorting)
    # later on remove and handle in generate_drug_edges()
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
    edges_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date_str}_drugSim_edges.csv'), index=False)
    drug_edges_df = pd.DataFrame(drug_edges)
    drug_edges_df.to_csv(os.path.join(drugSimil_directory, f'{disease_name_label}_{date_str}_alldrugSim_edges.csv'), index=False)
    
    logging.info("CSV files saved in drug_similarity directory.")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    print(f"'drugsimilarity.py' run finished in {minutes} minutes and {seconds} seconds.")
    logging.info(f"'drugsimilarity.py' run finished in {minutes} minutes and {seconds} seconds.")

    return unique_edges, unique_drug_edges