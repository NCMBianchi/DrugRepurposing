"""
DGIdb MODULE: MAKES API CALL AND ITERATES THROUGH NODES AT 2-3 DOD
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,platform,datetime,logging,builtins,time,multiprocessing

from drugapp.filepaths import initialise_base_directories, initialise_disease_directories, setup_service_globals
from drugapp.unique import unique_elements

import requests
import json
from unittest.mock import Mock
import pandas as pd

def json_reverse(json_dict):
    """
    Short function that converts calls turned into JSON files by the json() function
    into a request.Response object
    """
    json_file = Mock()
    json_file.json.return_value = json_dict
    json_file.status_code = 200
    json_file.ok = True
    return json_file

def hit_dgidb_api(gene_name=None, gene_id=None, max_retries=3):
    """
    This function performs API calls to DGIdb to retrieve interactions for a specific gene
    or all drugs.
    It will perform up to 3 retries if any 'ConnectionError', 'HTTPError' or 'Timeout' occurs.

    :param gene_name: Name of the gene. If None, retrieves all drugs.
    :param gene_id: ID of the gene. If None, retrieves all drugs.
    :param max_retries: The maximum number of retries (integer) (default = 3).
    :return: An API response object.
    """

    if gene_name and gene_id:
        logging.info(f"NOW RUNNING: {current_function_name()} with gene name {gene_name} and gene ID {gene_id}.")
    else:
        logging.info(f"NOW RUNNING: {current_function_name()} to retrieve all drugs in DGIdb.")

    # load strings required for path location (other than in 'filepaths.py')
    global disease_name_label
    #global disease_id_label (can be removed)
    global base_data_directory

    dgidb_directory = disease_directories['dgidb_directory']

    # GraphQL API endpoint and query
    dgidb_link = 'https://dgidb.org/api/graphql'
    
    if gene_name and gene_id:
        query = """
        {
          genes(names: ["%s"]) {
            nodes {
              interactions {
                drug {
                  name
                  conceptId
                }
                interactionScore
                interactionTypes {
                  type
                  directionality
                }
                interactionAttributes {
                  name
                  value
                }
                publications {
                  pmid
                }
                sources {
                  sourceDbName
                }
              }
            }
          }
        }
        """ % (gene_name)
    else:
        query = """
        {
          genes {
            nodes {
              name
              interactions {
                drug {
                  name
                  conceptId
                }
                interactionScore
                interactionTypes {
                  type
                  directionality
                }
                interactionAttributes {
                  name
                  value
                }
                publications {
                  pmid
                }
                sources {
                  sourceDbName
                }
              }
            }
          }
        }
        """

    headers = {'Content-Type': 'application/json'}
    attempt = 0
    timeout_seconds = 10
    while attempt < max_retries:
        try:
            response = requests.post(dgidb_link, json={'query': query}, headers=headers, timeout=timeout_seconds)
            logging.info(f"GraphQL request URL: {dgidb_link} with query {query}")
            logging.info(f"In edges response status: {response.status_code}")

            if response.ok:
                if gene_name and gene_id:
                    gene_id_string = gene_id.replace(':', '_')
                    gene_dir = os.path.join(dgidb_directory, f'{gene_id_string}')
                    drugs_file_path = os.path.join(gene_dir, f'{gene_id_string}_api_response.json')
                    os.makedirs(os.path.dirname(drugs_file_path), exist_ok=True)
                else:
                    drugs_file_path = os.path.join(dgidb_directory, 'all_drugs', f'{disease_name_label}_{date_str}_dgidb_all_drugs.json')
                    os.makedirs(os.path.dirname(drugs_file_path), exist_ok=True)

                with open(drugs_file_path, 'w') as json_file:
                    json.dump(response.json(), json_file)
                logging.info(f"API response saved to {drugs_file_path}")
                return response

            else:
                logging.warning(f"HTTP Error {response.status_code}: {response.reason}")
                if response.status_code == 504:
                    attempt += 1
                    if attempt < max_retries:
                        logging.info("Retrying due to 504 Gateway Timeout...")
                        time.sleep(5)  # Add a delay before retrying
                    else:
                        if not gene_name and not gene_id:  # If the call is for all drugs
                            file_name = "Huntington disease_2024-06-06_dgidb_all_drugs.json"
                            fallback_path = os.path.join(base_data_directory, "reference", "Huntington disease (2024-06-06)", "dgidb", "all_drugs", file_name)
                            logging.info(f"Retrieving data from fallback file: {fallback_path}")
                            with open(fallback_path, 'r') as file:
                                fallback_data = json.load(file)
                            response = json_reverse(fallback_data)
                            # save the fallback response as a JSON file
                            with open(drugs_file_path, 'w') as json_file:
                                json.dump(fallback_data, json_file)
                            logging.info(f"Fallback data saved to {drugs_file_path}")
                            return response
                        else:
                            return None
                else:
                    raise ValueError(f"HTTP Error {response.status_code}: {response.reason}")

        except (requests.ConnectionError, requests.HTTPError, requests.Timeout) as e:
            logging.warning(f"Attempt {attempt + 1} failed: {e}")
            attempt += 1
            if attempt < max_retries:
                time.sleep(5)  # Add a delay before retrying
            if attempt == max_retries:
                logging.error(f"All retries exhausted. Failed to fetch data.")
                if not gene_name and not gene_id:  # If the call is for all drugs
                    file_name = "Huntington disease_2024-06-06_dgidb_all_drugs.json"
                    fallback_path = os.path.join(base_data_directory, "reference", "Huntington disease (2024-06-06)", "dgidb", "all_drugs", file_name)
                    logging.info(f"Retrieving data from fallback file: {fallback_path}")
                    with open(fallback_path, 'r') as file:
                        fallback_data = json.load(file)
                    response = json_reverse(fallback_data)
                    # save the fallback response as a JSON file
                    drugs_file_path = os.path.join(dgidb_directory, 'all_drugs', f'{disease_name_label}_{date_str}_dgidb_all_drugs.json')
                    with open(drugs_file_path, 'w') as json_file:
                        json.dump(fallback_data, json_file)
                    logging.info(f"Fallback data saved to {drugs_file_path}")
                    return response
                else:
                    return None

    logging.warning("Failed to process drug interactions after multiple retries.")
    return None


def run_dgidb(monarch_input,date,layers = 3):
    """
    This function runs the whole DGIdb script and saves nodes and edges files.

    :param monarch_input: the input seed from the run_monarch() step.
    :param date: the date of creation of the disease graph.
    :param layers: how many layers to consider for the API calls (default = 3).
    :return: nodes and edges files in /DGIdb folder.
    """

    start_time = time.time()

    with open(input_file_path, 'a') as file:
        file.write(f"Degrees of distance from the Disease ID considered for DGIdb: {layers}\n\n")
    
    logging.info(f"NOW RUNNING: {current_function_name()} following 'run_monarch({monarch_input})'.")

    global nodes
    global edges

    dgidb_directory = os.path.join(today_directory, f'{disease_name_label} ({date_str})', 'dgidb')
    
    # initialise gene-to-drug list
    drug_edges = []

    query_nodes = nodes
    query_edges = edges

    # adjust edges and nodes based on the layers parameter ('logging.critical()' just to make it explicit)
    if layers == 3:
        logging.info("All 3 layers are considered.")
        print("All 3 layers are considered.")
        pass  # use all layers, no filtering needed
    elif layers == 2:
        query_edges = [edge for edge in query_edges if 'first_layer' in edge[3]['notes'] or 'second_layer' in edge[3]['notes']]
        logging.info("Only layers 1 and 2 are considered.")
        print("Only layers 1 and 2 are considered.")
    elif layers == 1:
        query_edges = [edge for edge in query_edges if 'first_layer' in edge[3]['notes']]
        logging.info("Only layer 1 is considered.")
        print("Only layer 1 is considered.")
    else:
        raise ValueError("Invalid number of layers specified. Choose 1 or 2 –otherwise 3 by default.")
    if layers in [1, 2]:
        node_ids = set(edge[0]['id'] for edge in query_edges) | set(edge[2]['id'] for edge in query_edges)
        query_nodes = [node for node in nodes if node['id'] in node_ids]

    logging.info(f"{len(query_nodes)} nodes and {len(query_edges)} edges from 'run_monarch({monarch_input})'.")
    
    valid_prefixes = ['flybase:', 'wormbase:', 'mgi:', 'hgnc:', 'ensembl:', 'zfin:']

    total_nodes = len(query_nodes)
    error_count_504 = 0
    
    for node in query_nodes:
        node_name = node['label']
        node_id = node['id']
        if any(prefix in node_id.lower() for prefix in valid_prefixes):
            node_id_string = node_id.replace(':', '_')
            gene_dir = os.path.join(dgidb_directory, f'{node_id_string}')
            drug_data = None  # Initialize drug_data to None

            # check if the API call for this gene ID was already done
            if os.path.exists(gene_dir):
                drug_file_path = os.path.join(gene_dir, f'{node_id_string}_api_response.json')
                with open(drug_file_path, 'r') as file:
                    drug_data = json.load(file)
                logging.info(f"Data for {node_id} fetched from local storage.")
            else:
                # make the API call
                try:
                    r_drugs = hit_dgidb_api(node_name, node_id)
                    if r_drugs is not None:
                        drug_data = r_drugs.json()
                        os.makedirs(gene_dir, exist_ok=True)
                        with open(os.path.join(gene_dir, f'{node_id_string}_api_response.json'), 'w') as file:
                            json.dump(drug_data, file)
                    else:
                        logging.info(f"No drug data returned for {node_id}")
                        error_count_504 += 1
                        continue  # Move to the next ID
                except ValueError as e:
                    logging.error(f"HTTP Error encountered for {node_id} - {str(e)}")
                    continue  # Move to the next ID
                except Exception as e:
                    logging.error(f"Failed to process drug interactions for {node_id}: {str(e)}")
                    continue  # Move to the next ID

            # extrapolate the drug-to-gene associations
            if drug_data:
                genes_data = drug_data.get('data', {}).get('genes', {}).get('nodes', [])
                if genes_data:
                    interactions = genes_data[0].get('interactions', [])
                    for interaction in interactions:
                        drug_info = interaction.get('drug', {})
                        drug_node = {'id': drug_info.get('conceptId'), 'label': drug_info.get('name')}
                        relation = {'label': 'dgidb:interacts_with'}
                        notes = {'notes': f"interactionScore:{interaction.get('interactionScore')}"}
                        drug_edges.append((node, relation, drug_node, notes))

    # MERGE new edges with the existing ones
    all_edges = query_edges + drug_edges
    unique_edges = unique_elements(all_edges)  # they should already be unique

    # EXTRACT nodes from edges
    all_nodes = [edge[0] for edge in unique_edges] + [edge[2] for edge in unique_edges]
    unique_nodes = unique_elements(all_nodes)

    # also store the entire list of drugs in DGIdb
    all_drug_dir = os.path.join(dgidb_directory,'all_drugs')
    os.makedirs(all_drug_dir, exist_ok=True)
    all_drug_path = os.path.join(all_drug_dir, f'{disease_name_label}_{date_str}_dgidb_all_drugs.csv')
    if os.path.exists(all_drug_path):
        # fetch from locally stored files
        all_drugs_df = pd.read_csv(all_drug_path)
        all_drugs = all_drugs_df.to_dict('records')
        logging.info(f"Data for all drugs fetched from local storage.")
    else:
        # make the API call to fetch the entire list of drugs in DGIdb
        all_genes_data = hit_dgidb_api()
        if all_genes_data is not None:
            all_genes = all_genes_data.json().get('data', {}).get('genes', {}).get('nodes', [])
            all_drugs = [{'id': interaction['drug']['conceptId'], 'label': interaction['drug']['name']} for gene in all_genes for interaction in gene.get('interactions', []) if 'conceptId' in interaction['drug'] and 'name' in interaction['drug']]
        else:
            all_drugs = []
    
    # save the unique nodes and edges as CSV
    nodes_df = pd.DataFrame(unique_nodes)
    nodes_df.to_csv(os.path.join(dgidb_directory, f'{disease_name_label}_{date_str}_dgidb_nodes.csv'), index=False)
    edges_df = pd.DataFrame(unique_edges)
    edges_df.to_csv(os.path.join(dgidb_directory, f'{disease_name_label}_{date_str}_dgidb_edges.csv'), index=False)

    # save the DGIdb nodes
    all_drugs_df = pd.DataFrame(all_drugs)
    all_drugs_df.to_csv(all_drug_path, index=False)
    
    logging.info("CSV files saved in DGIdb directory.")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    minutes = int(duration // 60)  # convert seconds to whole minutes
    seconds = int(duration % 60)  # get the remaining seconds
    logging.info(f"'DGIdb.py' run finished in {minutes} minutes and {seconds} seconds.")

    logging.info(f"{error_count_504} out of a total of {total_nodes} gene nodes encountered 'HTTP Error 504: Gateway Time-out'.")
    
    return unique_nodes, unique_edges, all_drugs