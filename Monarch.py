"""
MONARCH MODULE: MAKES API CALL AND ITERATES THROUGH DEGREES OF DISTANCE
Updated on February 5th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,time,inspect,requests,json
import pandas as pd
from unittest.mock import Mock

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


def json_reverse(json_dict):
    """
    Short function that converts calls turned into JSON files by the json() function
    into a request.Response object

    :param seed_list: JSON response stored as a dictionary.
    :return: Mock object mimicking requests.Response with the provided JSON data.
    """
    
    json_file = Mock()
    json_file.json.return_value = json_dict
    json_file.status_code = 200
    json_file.ok = True

    return json_file

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def hit_monarch_api(seed, d_path, rows=2000, direct_true=1):
    """
    Originally part of NQR's "bioknowledgeReviewer" tool.
    
    This function performs API calls to Monarch to retrieve 'OUT' and 'IN' edges from
        a query node.
    It retrieves entities plus associations via the Monarch Initiative new V3 API
        (2023 onward).
    
    The original version used the older Biolink API, no longer supported:
    https://api-biolink.monarchinitiative.org/api/association
     
    Moreover:
       - Biolink supported 'OMIM:' URIs, V3 doesn't and 'MONDO:' is used
       - Biolink required 'association/____', V3 requires 'association?____='
       - Biolink supported the ":" character, V3 requires its ASCII version (i.e.%3A)
       
    It hits two endpoints:
       - association?subject= (for 'OUT' edges)
       - association?object= (for 'IN' edges)
    It returns 'OUT' and 'IN' edges.

    :param node: The specific seed URI to query.
    :param d_path: __.
    :param rows: the maximum number of results to return (integer) (default = 2000).
    :param direct_true: if ==1 then 'direct=True' is added to the API call URL to filter
        out all URI aliases
    :return: two API response objects: 'OUT' and 'IN' response objects, in this order.
    """
    input_disease_id = d_path['disease_id']
    input_disease_name = d_path['disease_name']
    disease_path = d_path['disease_directory']
    monarch_path = d_path['monarch_directory']
    
    # API address
    biolink = 'https://api-v3.monarchinitiative.org/v3/api/association'  #V3
    seed_new = seed.split(':')
    seed_call = seed_new[0]+'%3A'+seed_new[1]

    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': rows}

    # OUT edges API call  ('?subject=')
    if direct_true == 1:
        out_url = f'{biolink}?subject={seed_call}&direct=true'
    elif direct_true == 0:
        out_url = f'{biolink}?subject={seed_call}'
    r_out = requests.get(out_url, params=parameters)  # API call

    # IN edges API call (?object=)
    if direct_true == 1:
        in_url = f'{biolink}?object={seed_call}&direct=true'
    elif direct_true == 0:
        in_url = f'{biolink}?object={seed_call}'
    r_in = requests.get(in_url, params=parameters)  # API call

    # Error Handling for API responses
    if not r_out.ok or not r_in.ok:
        error_message = f"Error fetching data for seed {seed}: OUT status {r_out.status_code}, IN status {r_in.status_code}"
        logging.error(error_message)
        # handle error appropriately, could raise an exception or return an error code
        raise Exception(error_message)
    
    # extract disease name and ID from the API response
    r_out_json = r_out.json()
    r_in_json = r_in.json()
    if seed == input_disease_id:
        # store disease ID and name in a file
        disease_id_file_path = os.path.join(disease_path, 'disease_id.txt')
        with open(disease_id_file_path, "w") as text_file:
            # store 'disease_name' and 'disease_id'
            text_file.write(f"{input_disease_id};{input_disease_name}\n")
            # store the list of subject closure for the input_seed
            # (information on the ontology of the disease)
            text_file.write(f"Subject closure:{r_out_json['items'][0]['subject_closure']}")
    
    # create a directory for each seed and store API responses
    seed_new = seed.split(':')
    seed_dir = seed_new[0]+'_'+seed_new[1]
    seed_directory = os.path.join(monarch_path, seed_dir)
    os.makedirs(seed_directory, exist_ok=True)
    
    # store OUT .json file
    out_file_path = os.path.join(seed_directory, f'{seed_dir}_api_response_OUT.json')
    if r_out.ok:
        with open(out_file_path, 'w') as json_file:
            json.dump(r_out_json, json_file)
    else:
        logging.error(f"Failed to fetch OUT edges: {r_out.status_code} - {r_out.reason}")

    # store IN .json file
    in_file_path = os.path.join(seed_directory, f'{seed_dir}_api_response_IN.json')
    if r_in.ok:
        with open(in_file_path, 'w') as json_file:
            json.dump(r_in_json, json_file)
    else:
        logging.error(f"Failed to fetch IN edges: {r_in.status_code} - {r_in.reason}")
    
    return r_out, r_in


def get_neighbours(seed_list,layer,disease_paths):
    """
    Originally part of NQR's "bioknowledgeReviewer" tool.
    
    This function parses through the .json files obtained from the hit_monarch_api()
    function, in order to obtain nodes that are related to any seed provided.
    
    If the API call for a given seed has already been run –and therefore a directory
    named with the same URI is present in
    ~/drugapp/data/YYYY-MM-DD/disease YYYY-MM-DD/monarch– then relevant information
    are fetched from the already locally stored 'URI_api_response_OUT.json' and
    'URI_api_response_IN.json' files. This allows for a reduced overhead.

    Information regarding subject, relation and object for both the 'OUT' and 'IN' calls
    are then stored in a table with a noSQL approach (?).

    :param seed_list: biomedical entities list, where each entity is a URI
        identifier string such as 'MONDO:0007739' or 'HGNC:4851'.
    :param layer: degrees of distance from the input disease, for logging.
    :param disease_name: variable required for correct path generation.
    :param disease_id: variable required for correct path generation.
    :param disease_dir: __.
    :return: neighbours lists, variables required for correct path generation.
    """

    # create the object that will store all the associations found while parsing
    # through the .json file ouputs from hit_monarch_api()
    edges_list = list()
    nodes_list = list()

    # iterate through each seed
    for seed in seed_list:
        seed_new = seed.split(':')
        seed_dir = seed_new[0]+'_'+seed_new[1]
        monarch_path = disease_paths['monarch_directory']
        seed_path = os.path.join(monarch_path, seed_dir)
        
        # run API calls (OUT and IN) for the given seed in the 'FOR' loop
        out_json_path = os.path.join(seed_path, f'{seed_dir}_api_response_OUT.json')
        in_json_path = os.path.join(seed_path, f'{seed_dir}_api_response_IN.json')
        if os.path.exists(out_json_path) and os.path.exists(in_json_path):
            with open(out_json_path, 'r') as file:
                r_out_json = json.load(file)
            with open(in_json_path, 'r') as file:
                r_in_json = json.load(file)
            r_out = json_reverse(r_out_json)
            r_in = json_reverse(r_in_json)
        else:
            r_out, r_in = hit_monarch_api(seed,disease_paths)
        
        # .json file handling
        for associations in [r_out.json()['items'], r_in.json()['items']]:
            for association in associations:
                
                # parse through the files and store the associations in edges_list
                subj = {'id': association['subject'], 'label': association['subject_label']}
                rel = {'label': association['predicate']} # relation ID no longer present
                obj = {'id': association['object'], 'label': association['object_label']}
                notes = {'notes': layer}
                ## reference (link to literature) no longer present
                association_elem = [subj, rel, obj, notes]
                edges_list.append(association_elem)

                # extrapolate list of nodes (both subjects and objects) and store it in nodes_list
                nodes_list.append(subj)
                nodes_list.append(obj)

    # filter for only unique nodes and edges
    nodes_list = unique_elements(nodes_list)
    edges_list = unique_elements(edges_list)
    
    return nodes_list, edges_list


def run_monarch(input_id, disease_directories):
    """
    This function runs the whole Monarch script and saves nodes and edges files.

    :param input_number: The input MONDO: URI  of the disease (default = 'MONDO:0007739' for
        Huntington's disease).
    :param disease_directories: base paths to where data is stored.
    :return: nodes and edges related to the disease of interests in Monarch Initiative's database,
        variables required for correct path generation.
    """
    
    print(f"NOW RUNNING: {current_function_name()} with seed {input_id}.")

    start_time = time.time()
    
    # initialise path
    disease_directory = disease_directories['disease_directory']
    monarch_directory = disease_directories['monarch_directory']
    date_str = disease_directories['date_string']
    disease_name_label = disease_directories['disease_name']
    
    # First layer: use the flaskApp input
    firstLayer_nodes, firstLayer_edges = get_neighbours(input_id,'first_layer',disease_directories)

    # Second layer: use the nodes in the first layer
    firstLayer_seeds = [item['id'] for item in firstLayer_nodes]  # turn into a list of IDs
    firstLayer_seeds = [fl_id for fl_id in firstLayer_seeds if fl_id not in input_id]  # to avoid re-running nodes and counting them in the second layer
    secondLayer_nodes, secondLayer_edges = get_neighbours(firstLayer_seeds,'second_layer',disease_directories)

    # Third layer: use the nodes in the second layer
    secondLayer_seeds = [item['id'] for item in secondLayer_nodes]  # turn into a list of IDs
    secondLayer_seeds = [sl_id for sl_id in secondLayer_seeds if sl_id not in firstLayer_seeds]  # to avoid re-running nodes and counting them in the third layer
    thirdLayer_nodes, thirdLayer_edges = get_neighbours(secondLayer_seeds,'third_layer',disease_directories)
    
    # MERGE the nodes' lists
    nonUnique_nodes = firstLayer_nodes + secondLayer_nodes + thirdLayer_nodes
    unique_nodes = unique_elements(nonUnique_nodes)
    
    # MERGE the edges' lists
    nonUnique_edges = firstLayer_edges + secondLayer_edges + thirdLayer_edges
    unique_edges = unique_elements(nonUnique_edges)
    
    # save the unique nodes and edges as CSV
    nodes_df = pd.DataFrame(unique_nodes)
    nodes_df.to_csv(os.path.join(monarch_directory, f'{disease_name_label}_{date_str}_monarch_nodes.csv'), index=False)
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
    edges_df.to_csv(os.path.join(monarch_directory, f'{disease_name_label}_{date_str}_monarch_edges.csv'), index=False)

    print(f"EDGES: {len(unique_edges)}")
    print(f"NODES: {len(unique_nodes)}")

    end_time = time.time()
    duration = end_time - start_time  # calculate duration in seconds
    formatted_duration = format_duration(duration)  # convert for print
    print(f"'Monarch.py' run finished in {formatted_duration}.")
    
    return unique_nodes, unique_edges