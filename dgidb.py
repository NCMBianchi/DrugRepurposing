"""
DGIdb MODULE: MAKES API CALL AND ITERATES THROUGH NODES AT 2-3 DOD
Updated on February 5th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,inspect,requests,json
import pandas as pd

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

def hit_dgidb_api(disease_dir, gene_name=None, gene_id=None, max_retries=3):
    """
    This function performs API calls to DGIdb to retrieve interactions for a specific gene
    or all drugs.
    It will perform up to 3 retries if any 'ConnectionError', 'HTTPError' or 'Timeout' occurs.

    :param disease_dir: base paths to where data is stored.
    :param gene_name: Name of the gene. If None, retrieves all drugs.
    :param gene_id: ID of the gene. If None, retrieves all drugs.
    :param max_retries: The maximum number of retries (integer) (default = 3).
    :return: An API response object.
    """
    date_str = disease_dir['date_string']
    dgidb_dir = disease_dir['dgidb_directory']
    
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
    if gene_name and gene_id:
        timeout_seconds = 10
    else:
        timeout_seconds = 90
    while attempt < max_retries:
        try:
            response = requests.post(dgidb_link, json={'query': query}, headers=headers, timeout=timeout_seconds)

            if response.ok:
                if gene_name and gene_id:
                    pass
                else:
                    print(f"ALLDRUGS FETCH: attempt {attempt+1} successful.")
                if gene_name and gene_id:
                    gene_id_string = gene_id.replace(':', '_')
                    gene_dir = os.path.join(dgidb_dir, f'{gene_id_string}')
                    drugs_file_path = os.path.join(gene_dir, f'{gene_id_string}_api_response.json')
                    os.makedirs(os.path.dirname(drugs_file_path), exist_ok=True)
                else:
                    drugs_file_path = os.path.join(dgidb_dir, 'all_drugs', f'dgidb_all_drugs_{date_str}_api_response.json')
                    os.makedirs(os.path.dirname(drugs_file_path), exist_ok=True)

                with open(drugs_file_path, 'w') as json_file:
                    json.dump(response.json(), json_file)
                return response

            else:
                if response.status_code == 504:
                    attempt += 1
                    if gene_name and gene_id:
                        pass
                    else:
                        print(f"ALLDRUGS FETCH: attempt {attempt+1} failed for HTTP error 504: {response.reason}.")
                    if attempt < max_retries:
                        time.sleep(5)  # Add a delay before retrying
                    else:
                        return None
                else:
                    if gene_name and gene_id:
                        pass
                    else:
                        print(f"ALLDRUGS FETCH: attempt {attempt+1} failed for HTTP error {response.status_code}:{response.reason}.")

        except (requests.ConnectionError, requests.HTTPError, requests.Timeout) as e:
            if gene_name and gene_id:
                pass
            else:
                print(f"ALLDRUGS FETCH: attempt {attempt+1} failed: {type(e).__name__}, {str(e)}.")
            attempt += 1
            if attempt < max_retries:
                time.sleep(5)  # Add a delay before retrying
            if attempt == max_retries:
                return None

    return None


def run_dgidb(nodes,edges,disease_directories,layers = 3):
    """
    This function runs the whole DGIdb script and saves nodes and edges files.

    :param nodes: all the node resulting from the previous step in the pipeline (i.e. 'run_monarch()').
    :param edges: all the edges resulting from the previous step in the pipeline (i.e. 'run_monarch()').
    :param disease_directories: base paths to where data is stored.
    :param layers: how many layers to consider for the API calls (default = 3).
    :return: drug nodes and gene-to-drug edges in DGIdb, plus the nodes and edges related to the disease
        of interests in Monarch Initiative's database.
    """

    start_time = time.time()

    # initialise path
    dgidb_directory = disease_directories['dgidb_directory']
    date_str = disease_directories['date_string']
    disease_name_label = disease_directories['disease_name']

    print(f"Degrees of distance from the Disease ID considered for DGIdb: {layers}.")
    
    print(f"NOW RUNNING: {current_function_name()} following 'run_monarch()'.")
    
    # initialise gene-to-drug list
    drug_edges = []

    query_nodes = nodes
    query_edges = edges

    # adjust edges and nodes based on the layers parameter ('logging.critical()' just to make it explicit)
    if layers == 3:
        pass  # use all layers, no filtering needed
    elif layers == 2:
        query_edges = [edge for edge in query_edges if 'first_layer' in edge[3]['notes'] or 'second_layer' in edge[3]['notes']]
    elif layers == 1:
        query_edges = [edge for edge in query_edges if 'first_layer' in edge[3]['notes']]
    else:
        raise ValueError("Invalid number of layers specified. Choose 1 or 2 –otherwise 3 by default.")
    if layers in [1, 2]:
        node_ids = set(edge[0]['id'] for edge in query_edges) | set(edge[2]['id'] for edge in query_edges)
        query_nodes = [node for node in nodes if node['id'] in node_ids]
    
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
            drug_file_path = os.path.join(gene_dir, f'{node_id_string}_api_response.json')
            if os.path.exists(drug_file_path):
                if os.path.getsize(drug_file_path) > 0:
                    with open(drug_file_path, 'r') as file:
                        drug_data = json.load(file)
            else:
                # make the API call
                try:
                    r_drugs = hit_dgidb_api(disease_directories, node_name, node_id)
                    if r_drugs is not None:
                        drug_data = r_drugs.json()
                        os.makedirs(gene_dir, exist_ok=True)
                        with open(os.path.join(gene_dir, f'{node_id_string}_api_response.json'), 'w') as file:
                            json.dump(drug_data, file)
                    else:
                        error_count_504 += 1
                        continue  # Move to the next ID
                except ValueError as e:
                    continue  # Move to the next ID
                except Exception as e:
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
    all_drug_path = os.path.join(all_drug_dir, f'dgidb_all_drugs_{date_str}_api_response.csv')
    if os.path.exists(all_drug_path):
        if os.path.getsize(all_drug_path) > 1:
            # fetch from locally stored files
            all_drugs_df = pd.read_csv(all_drug_path)
            all_drugs = all_drugs_df.to_dict('records')
    else:
        # make the API call to fetch the entire list of drugs in DGIdb
        all_genes_data = hit_dgidb_api(disease_directories)
        if all_genes_data is not None:
            all_genes = all_genes_data.json().get('data', {}).get('genes', {}).get('nodes', [])
            all_drugs = [{'id': interaction['drug']['conceptId'], 'label': interaction['drug']['name']} for gene in all_genes for interaction in gene.get('interactions', []) if 'conceptId' in interaction['drug'] and 'name' in interaction['drug']]
        else:
            all_drugs = []
    
    # save the unique nodes and edges as CSV
    nodes_df = pd.DataFrame(unique_nodes)
    nodes_df.to_csv(os.path.join(dgidb_directory, f'{disease_name_label}_{date_str}_dgidb_nodes.csv'), index=False)
    edges_df = pd.DataFrame([
        {
            'subject_id': edge[0]['id'],
            'subject_label': edge[0]['label'],
            'relation': edge[1]['label'],
            'object_id': edge[2]['id'],
            'object_label': edge[2]['label'],
            'notes': edge[3]['notes']
        }
        for edge in unique_edges
    ])
    edges_df.to_csv(os.path.join(dgidb_directory, f'{disease_name_label}_{date_str}_dgidb_edges.csv'), index=False)

    # save the DGIdb nodes
    all_drugs_df = pd.DataFrame(all_drugs)
    all_drugs_df.to_csv(all_drug_path, index=False)

    print(f"MATCHED DRUGS: {len(drug_edges)} drug-to-gene edges out of a total of {len(all_drugs)} drugs on DGIdb.")
    print(f"{error_count_504} out of a total of {total_nodes} gene nodes encountered 'HTTP Error 504: Gateway Time-out'.")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    formatted_duration = format_duration(duration)  # convert for print
    print(f"'DGIdb.py' run finished in {formatted_duration}.")
    
    return unique_nodes, unique_edges, all_drugs