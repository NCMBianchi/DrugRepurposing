# -*- coding: utf-8 -*-
"""
MODULE FOR MONARCH NETWORK PREPARATION AND MANAGEMENT

Last modified on Wed March 27 10:30:00 2024

@authors: Núria Queralt Rosinach, Carmen Reep, Niccolò Bianchi
"""

import requests
import sys,os
import json  # to save the latest .json returned by the API
import datetime
import logging  # to check for running variables
import inspect  # to chec for running functions
import pandas as pd
from biothings_client import get_client
from tqdm import tqdm


# timestamp
today = datetime.date.today()
curr_year = int(str(today)[:4])


# logging configuration
base_data_directory = os.path.join(os.getcwd(), 'drugapp', 'data')
log_directory = os.path.join(base_data_directory, 'logs')
log_filename = datetime.datetime.now().strftime('monarch_app_%Y-%m-%d.log')
log_file_path = os.path.join(log_directory, log_filename)
os.makedirs(log_directory, exist_ok=True)
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(log_file_path),
                        logging.StreamHandler()  # Enables logging to stderr.
                    ])

def current_function_name():
    return inspect.currentframe().f_back.f_code.co_name



def read_connections(filename):  #bioknowledgeReviewer
    """
    This function reads monarch_connections CSV file.
    :param filename: complete path to the monarch connections csv file string
    :return: monarch edges dataframe
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # monarch network
    path = os.getcwd() + '/monarch'
    csv_path = path + '/' + filename
    network_df = pd.read_csv('{}'.format(csv_path))
    #print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    #print('* These are the attributes: {}'.format(network_df.columns))
    #print('* This is the first record:\n{}'.format(network_df.head(1)))
    logging.info(f"Reading monarch connections CSV from: {csv_path}")
    logging.indo(f"DataFrame: {network_df}")
    logging.info(f"DataFrame size: {network_df.shape}")
    logging.info(f"DataFrame attributes: {network_df.columns.tolist()}")
    logging.info(f"First record in the DataFrame:\n{network_df.head(1)}")

    return network_df


## RETRIEVE SUBNETWORK FROM MONARCH KNOWLEDGE GRAPH

def hit_monarch_api(node = 'HGNC:4851', rows = 2000):  #bioknowledgeReviewer
    """
    This function performs api calls to Monarch to retrieve out and in edges from a query node.
    It retrieves entities plus associations via the BioLink API service.
    It hits two endpoints:
        * association/from - for out edges
        * association/to - for in edges
    It returns out and in edges.

    :param node: node id to query (string). Default: 'HGNC:4851' (=HTT gene).
    :param rows: the maximum number of results to return (integer). Default: 2000.
    :return: two api response objects: 'out' and 'in' response objects, in this order.
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    # API address
    #biolink = 'https://api-biolink.monarchinitiative.org/api/association'  #biolink
    biolink = 'https://api-v3.monarchinitiative.org/v3/api/association'  #V3
    node_new = node.split(':')
    node_call = node_new[0]+'%3A'+node_new[1]
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': rows}
    # out edges: from/
    #r_out = requests.get('{}/from/{}'.format(biolink,node_call),params=parameters)
    #out_url = f'{biolink}/from/{node_call}'  #biolink
    out_url = f'{biolink}?subject={node_call}'  #V3
    r_out = requests.get(out_url, params=parameters)
    logging.info(f"Out edges request URL: {out_url} with parameters {parameters}")
    logging.info(f"Out edges response status: {r_out.status_code}")

    # in edges: to/
    #r_in = requests.get('{}/to/{}'.format(biolink,node_call),params=parameters)
    #in_url = f'{biolink}/to/{node_call}'  #biolink
    in_url = f'{biolink}?object={node_call}'  #V3
    r_in = requests.get(in_url, params=parameters)
    logging.info(f"In edges request URL: {in_url} with parameters {parameters}")
    logging.info(f"In edges response status: {r_in.status_code}")

    return r_out, r_in


def get_disease_name_id(disease_input_ID = 'MONDO:0007739'):
    """
    This function finds the disease name and disease ID of the input URI ID.
    It returns two strings, name and ID respectively.

    :param disease_ID: The input ID of the disease
    :return: two strings, name and ID respectively
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # api address
    biolink = 'https://api-v3.monarchinitiative.org/v3/api/association'
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': 1}
    ## out edges: from/
    #r_out_disease = requests.get('{}/from/{}'.format(biolink,disease_input_ID),params=parameters)

    disease_input_ID_new = disease_input_ID.split(':')
    disease_input_ID_call = disease_input_ID_new[0]+'%3A'+disease_input_ID_new[1]
    #request_url = f'{biolink}/from/{disease_input_ID}'  #biolink
    request_url = f'{biolink}?subject={disease_input_ID_call}'  #V3
    logging.info(f"Making API call to: {request_url}")
    r_out_disease = requests.get(request_url, params=parameters)
    #for association in r_out_disease.json()['associations']:
    #    disease_name = association['subject']['label']
    #    disease_id = association['subject']['id']

    response_json = r_out_disease.json()
    #print(response_json)
    #if 'associations' in response_json:
    #    for association in response_json['associations']:
    #        disease_name = association['subject']['label']
    #        disease_id = association['subject']['id']
    #else:
    #    print("Warning: 'associations' key not found in response.")
    if 'items' in response_json and response_json['items']:
        disease_name = response_json['items'][0]['subject_label']
        disease_id = response_json['items'][0]['subject']
    else:
        logging.warning("Warning: 'associations' key not found in response or no associations present.")
        return None, None
        
    # ADDITIONAL STEP TO CHECK THE .json FILE
    #with open('api_response.json', 'w') as json_file:
    #    json.dump(response_json, json_file)
    data_directory_path = os.path.join(os.getcwd(), 'drugapp', 'data')
    json_file_path = os.path.join(data_directory_path, 'api_response.json')
    os.makedirs(data_directory_path, exist_ok=True)
    with open(json_file_path, 'w') as json_file:
        json.dump(response_json, json_file)
    logging.info(f"API response saved to {json_file_path}")

    return disease_name, disease_id


def get_edges_objects(r_out, r_in):  #bioknowledgeReviewer
    """
    This function prepares the api object responses from Monarch.
    It returns four lists, one for subjects, relations, objects, and references.
    Subjects, relations and objects are lists of dictionaries, where each dictionary is a node.
    References list lists strings, where each string is a chain of references for each edge.

    :param r_out: BioLink API 'out' response object
    :param r_in: BioLink API 'in' response object
    :return: subjects, relations, objects and references lists (in this order)
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # variables
    sub_l = list()
    rel_l = list()
    obj_l = list()
    ref_l = list()

    ## ADDITIONAL CODE
    #biolink = 'https://api-biolink.monarchinitiative.org/api/association'
    #response = requests.get(biolink, params=parameters)
    #if response.status_code == 200:
    #    response_json = response.json()
    #    if 'associations' in response_json:
    #
    #        ## out edges: from/
    #        r_out = requests.get('{}/from/{}'.format(biolink,node),params=parameters)
    #    
    #        # in edges: to/
    #        r_in = requests.get('{}/to/{}'.format(biolink,node),params=parameters)
    #
    #    else:
    #        print("Warning: 'associations' key not found in response.")
    #else:
    #    print(f"API call failed with status code {response.status_code}")

    # compose list of dictionaries
    #for associations in [r_out.json()['associations'], r_in.json()['associations']]:
    for associations in [r_out.json()['items'], r_in.json()['items']]:
        for association in associations:
            pub_l = list()
            sub_l.append({'id': association['subject'], 'label': association['subject_label']})
            rel_l.append({'id': association['predicate'], 'label': association['predicate']}) #######relation -> only label given
            obj_l.append({'id': association['object'], 'label': association['object_label']})
            # add references to each association as a list of strings
            if association['publications']:
                for publication in association['publications']:
                    pub_l.append(publication)#['id']) -> only id is given in new version
            else:
                pub_l.append('NA')
            ref_l.append('|'.join(pub_l))

    logging.info("Prepared Monarch API object responses.")
    logging.info(f"Subjects list length: {len(sub_l)}")
    logging.debug(f"Subjects sample: {sub_l[:2]}")
    logging.info(f"Relations list length: {len(rel_l)}")
    logging.debug(f"Relations sample: {rel_l[:2]}")
    logging.info(f"Objects list length: {len(obj_l)}")
    logging.debug(f"Objects sample: {obj_l[:2]}")
    logging.info(f"References list length: {len(ref_l)}")
    logging.debug(f"References sample: {ref_l[:2]}")

    return sub_l, rel_l, obj_l, ref_l


def get_edges(sub_l, rel_l, obj_l, ref_l, attribute='id'):  #bioknowledgeReviewer
    """
    This function builds edges using a user-specified attribute for each node.
    It returns a set of edges, where edges are tuples.

    :param sub_l: subjects (objects) list from the get_edges_objects() function
    :param rel_l: relations (objects) list from the get_edges_objects() function
    :param obj_l: objects (objects) list from the get_edges_objects() function
    :param ref_l: references (strings) list from the get_edges_objects() function
    :param attribute: object attribute, default 'id'
    :return: edges (as tuples) set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    edges = set()
    # compose tuple
    for i in range(len(sub_l)):
        sub = sub_l[i][attribute]
        rel = rel_l[i][attribute]
        obj = obj_l[i][attribute]
        ref = ref_l[i]
        edges.add((sub, rel, obj, ref))

    logging.info(f"Generated edges set with {len(edges)} edges.")
    sample_edges = list(edges)[:5]
    logging.debug(f"Sample edges: {sample_edges}")

    return edges


def keep_edges(keep, new):  #bioknowledgeReviewer
    """
    This function adds edges from a new set to a keep set.
    :param keep: edges set
    :param new: edges set
    :return: updated edges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    logging.info(f"Before adding: keep set size = {len(keep)}, new set size = {len(new)}")

    for edge in new:
        keep.add(edge)

    logging.info(f"After adding: keep set size = {len(keep)}")
    sample_edges = list(keep)[:5]
    logging.debug(f"Sample edges from the updated keep set: {sample_edges}")

    return keep


def keep_nodes(keep, edges, seed):  #bioknowledgeReviewer
    """
    This function collects nodes from the edges that are not in the nodes query list to keep nodes set.
    It filters out: PMID nodes, and nodes related by provenance:
        * None
        * 'dc:source'
        * 'IAO:0000136' or 'is about'
        * 'IAO:0000142' or 'mentions'
    i.e., not biologically related
    It returns a set of nodes.

    :param keep: nodes set
    :param edges: edges set
    :param seed: query nodes list
    :return: updated nodes set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    logging.info(f"Before filtering: keep set size = {len(keep)}")
    initial_size = len(keep)

    for (sub, rel, obj, ref) in edges:
        if 'PMID' in sub or 'PMID' in obj:
            continue
        if rel == None:
            rel = 'None'
        if 'dc:source' in rel:
            continue
        if 'IAO:0000136' in rel:  # is about (?)
            continue
        if 'IAO:0000142' in rel:  # mentions (?)
            continue
        if sub not in seed:
            keep.add(sub)
        if obj not in seed:
            keep.add(obj)
    # https://bioportal.bioontology.org/ontologies/IAO

    logging.info(f"After filtering: keep set size = {len(keep)}; added {len(keep) - initial_size} new nodes")
    sample_nodes = list(keep)[:5]
    logging.debug(f"Sample nodes from the updated keep set: {sample_nodes}")

    return keep


def get_neighbours(seed):  #bioknowledgeReviewer
    """
    This function gets the first layer of neighbours and relations.
    :param seed: query nodes list
    :return: nodes set, edges set (in this order)
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    keepNodes = set()
    keepEdges = set()
    seedNodes = set(seed)
    for node in tqdm(seedNodes):
        try:
            r_out, r_in = hit_monarch_api(node)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            keepEdges = keep_edges(keepEdges, edges)
            keepNodes = keep_nodes(keepNodes, edges, seedNodes)

        #except (ValueError, KeyError):
        #    pass
        #except:
        #    print('error: {}'.format(sys.exc_info()[0]))
        #    print(node)
        except (ValueError, KeyError):  # (!!!) where it seems to break
            logging.warning(f"Skipping node {node} due to ValueError or KeyError.")
        except Exception as e:
            logging.error(f"An error occurred for node {node}: {e}")

    logging.info(f"Final keepNodes set size: {len(keepNodes)}")
    sample_nodes = list(keepNodes)[:5]
    logging.debug(f"Sample nodes from keepNodes set: {sample_nodes}")

    logging.info(f"Final keepEdges set size: {len(keepEdges)}")
    sample_edges = list(keepEdges)[:5]
    logging.debug(f"Sample edges from keepEdges set: {sample_edges}")

    return keepNodes, keepEdges


def filter_edges(nodes, edges):  #bioknowledgeReviewer
    """
    This function filters down edges with both nodes in a nodes list.

    :param nodes: nodes list
    :param edges: edges set
    :return: filtered edges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    nodes = set(nodes)
    keep = set()
    for (start, pred, stop, ref) in edges:
        if {start, stop} <= nodes: # is {..} a subset of {..}
            keep.add((start, pred, stop, ref))

    logging.info(f"Filtered edges set size: {len(keep)}")
    sample_keep = list(keep)[:5]
    logging.debug(f"Sample of filtered edges: {sample_keep}")

    return keep


def add_attributes(sub_l, rel_l, obj_l, edges):  #bioknowledgeReviewer
    """
    This function adds 'label', 'iri', 'category' attribute to each entity in the edge.
    :param sub_l: subjects (object) list
    :param rel_l: relations (object) list
    :param obj_l: objects (object) list
    :param edges: edges set
    :return: metaedges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    metaedges = set()
    for (sub_id, rel_id, obj_id, refs) in edges:
        for i in range(len(sub_l)):
            if sub_l[i]['id'] == sub_id and rel_l[i]['id'] == rel_id and obj_l[i]['id'] == obj_id:
                metaedges.add((sub_l[i]['id'],
                               sub_l[i]['label'],
                               #sub_l[i]['iri'],  # temporary removed for testing
                               sub_l[i]['category'][0], # add [0] because category is in a list
                               rel_l[i]['id'],
                               rel_l[i]['label'],
                               #rel_l[i]['iri'],
                               obj_l[i]['id'],
                               obj_l[i]['label'],
                               #obj_l[i]['iri'],
                               #obj_l[i]['category'][0],
                               refs)
                              )
                break

    logging.info(f"Metaedges set size: {len(metaedges)}")
    sample_metaedges = list(metaedges)[:5]
    logging.debug(f"Sample of metaedges: {sample_metaedges}")

    return metaedges


def keep_node_type(edges, seed, nodeType='ortho'):  #bioknowledgeReviewer
    """
    This function keeps specific node types for objects in edges.

    :param edges: edges set
    :param seed: the query nodes list
    :param nodeType: Introduce node type to keep (string): 'ortho' for orthologs or 'pheno' \
    for phenotypes/diseases, default is 'ortho'
    :return: nodes set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    propertyList = ['RO:HOM0000017', 'RO:HOM0000020']
    if nodeType == 'pheno':
        propertyList = ['RO:0002200', 'RO:0002607', 'RO:0002326', 'GENO:0000840']
    # https://github.com/monarch-initiative/dipper/issues/378 apparently they no longer are in OWL

    keep = set()
    for (sub, rel, obj, ref) in edges:
        if rel == None:
            continue
        if rel in propertyList:
            if sub not in seed:
                keep.add(sub)
            if obj not in seed:
                keep.add(obj)

    logging.info(f"Kept nodes set size: {len(keep)}")
    sample_keep_nodes = list(keep)[:5]
    logging.debug(f"Sample of kept nodes: {sample_keep_nodes}")

    return keep


def get_connections(nodes):  #bioknowledgeReviewer
    """
    This function returns associations retrieved from Monarch among a list of query nodes.
    :param nodes: the query nodes list
    :return: edges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    keep = set()
    for node in tqdm(nodes):
        try:
            r_out, r_in = hit_monarch_api(node, 1000)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            filteredEdges = filter_edges(nodes, edges)
            metaFilteredEdges = add_attributes(sub_l, rel_l, obj_l, filteredEdges)
            keep = keep_edges(keep, metaFilteredEdges)

        except (ValueError, KeyError):
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            #print(node)

    logging.info(f"Final edges set size: {len(keep)}")
    sample_keep_edges = list(keep)[:5]
    logging.debug(f"Sample of final kept edges: {sample_keep_edges}")

    return keep


# NETWORK MANAGEMENT FUNCTIONS

def get_neighbours_list(seed_list):  #bioknowledgeReviewer
    """
    This function returns the first explicit layer of neighbours from a list of query nodes.
    :param seed_list: biomedical entities list, where each entity is the identifier string like 'HGNC:4851'
    :return: neighbours list
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    # print executing function
    print('\nThe function "get_neighbours_list()" is running. Its runtime may take some minutes. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)
    neighbours_list = list(neighbours)
    print('\nFinished get_neighbours_list().\n')

    logging.info(f"Neighbours list size: {len(neighbours_list)}")
    sample_neighbours = neighbours_list[:5]
    logging.debug(f"Sample of neighbours: {sample_neighbours}")

    return neighbours_list


def get_orthopheno_list(seed_list):  #bioknowledgeReviewer
    """
    This function returns orthologs-phenotypes nodes in ortho-pheno relationships for a list of query genes.
    :param seed_list: gene list, where each gene is the identifier string like 'HGNC:4851'
    :return: orthopheno list
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    # print executing function
    print('\nThe function "get_orthopheno_list()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)

    # keep orthologs in the first layer
    orthologs = keep_node_type(relations, seed_list)

    # get second layer from orthologs
    neighbours, relations = get_neighbours(orthologs)

    # keep phenotypes in the second layer
    phenotypes = keep_node_type(relations, orthologs, 'pheno')

    nodes = set()
    nodes.update(orthologs, phenotypes)
    nodes_list = list(nodes)
    print('\nFinished get_orthopheno_list().\n')

    logging.info(f"Orthopheno list size: {len(nodes_list)}")
    sample_nodes = nodes_list[:5]
    logging.debug(f"Sample of orthopheno nodes: {sample_nodes}")

    return nodes_list


def extract_edges(gene_list):  #bioknowledgeReviewer
    """
    This function returns the Monarch network from a list of query nodes. It retrieves connectivity from Monarch, i.e. \
    edges from Monarch between query nodes.
    :param gene_list: gene list
    :return: edges (as tuples) set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    # print executing function
    print('\nThe function "extract_edges()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the edges retrieved '
          'and you should start over the execution of this function.')

    # set network nodes: gene list provided by the user
    nodes = set(gene_list)

    # get connections
    network = get_connections(nodes)
    print('\nFinished extract_edges(). To save the retrieved Monarch edges use the function "print_network()".\n')

    logging.info(f"Extracted network size: {len(network)}")
    sample_network = list(network)[:5]
    logging.debug(f"Sample of extracted network edges: {sample_network}")

    return network


def _print_network2(network, filename):  #bioknowledgeReviewer
    """This function saves the Monarch expanded network into a CSV file. this function save connections file format into
    get-monarch-connections/"""

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path):
        os.makedirs(path)
    #with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
    #    f.write(
    #        'subject_id,subject_label,subject_uri,subject_category,relation_id,relation_label,relation_uri,object_id,object_label,object_uri,object_category,reference_id_list\n')
    #    for edge in network:
    #        edge = ['None' if t is None else '"{}"'.format(str(t)) for t in edge]
    #        f.write('{}\n'.format(','.join(edge)))
    file_path = f"{path}/{filename}_v{today}.csv"
    with open(file_path, 'w') as f:
        f.write('subject_id,subject_label,subject_uri,subject_category,relation_id,relation_label,relation_uri,object_id,object_label,object_uri,object_category,reference_id_list\n')
        for edge in network:
            edge_str = ['None' if t is None else f'"{str(t)}"' for t in edge]
            f.write(f"{','.join(edge_str)}\n")

    logging.info(f"Monarch network saved to CSV file: {file_path}")
    sample_data = [','.join(['None' if t is None else f'"{str(t)}"' for t in edge]) for edge in network[:5]]
    logging.debug("Sample data being saved to CSV:")
    for line in sample_data:
        logging.debug(line)
    logging.debug(f"Number of edges saved: {len(network)}")

    # commented out since the file is saved elsewhere in the function
    #return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def print_network(network, filename):  #bioknowledgeReviewer
    """
    This function saves the Monarch network into a CSV file. It only prints connections file format only.
    :param network: monarch edges dataframe
    :param filename: file name without extension string, e.g. 'monarch_connections'
    :return: None object
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # transform set of tuples to list of dictionaries
    edges = list()
    for tuple in network:
        row = dict()
        row['subject_id'] = tuple[0]
        row['subject_label'] = tuple[1]
        #row['subject_uri'] = tuple[2]  # temporary removed for testing
        row['subject_category'] = tuple[3]
        row['relation_id'] = tuple[4]
        row['relation_label'] = tuple[5]
        #row['relation_uri'] = tuple[6]
        row['object_id'] = tuple[7]
        row['object_label'] = tuple[8]
        #row['object_uri'] = tuple[9]
        #row['object_category'] = tuple[10]
        row['reference_id_list'] = tuple[11]
        edges.append(row)
    
    df=pd.DataFrame(edges).fillna('None')
    
    ## change every 'model' node that does not have 'Coriell' as prefix to genotype
    #df.loc[(df['subject_category'] == 'model') &  (df['subject_id'].str.contains('zfin|mgi|flybase|wormbase', case=False)), 'subject_category'] = 'genotype'
    #df.loc[(df['object_category'] == 'model') &  (df['object_id'].str.contains('zfin|mgi|flybase|wormbase', case=False)), 'object_category'] = 'genotype'
    #
    ## remove every row that has Coriell in either subject or object id
    #df.drop(df[df['object_id'].str.contains('coriell', case=False)].index, inplace=True)
    #df.drop(df[df['subject_id'].str.contains('coriell', case=False)].index, inplace=True)
    # (!!!) Removed this section for testing

    ## print output file
    #path = os.getcwd() + '/monarch'
    #if not os.path.isdir(path): os.makedirs(path)
    #df.to_csv('{}/{}_v{}.csv'.format(path, filename, today), index=False)
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path):
        os.makedirs(path)
    file_path = f'{path}/{filename}_v{today}.csv'
    df.to_csv(file_path, index=False)

    logging.debug("Sample data from DataFrame being saved to CSV:")
    sample_df = df.head(5).to_csv(index=False, line_terminator='\n')
    logging.debug(sample_df)
    logging.info(f"Monarch edges saved at: {file_path}")

    # commented out since the file is saved elsewhere in the function
    #return print("\nSaving Monarch edges at: '{}/{}_v{}.csv'...\n".format(path, filename, today))


def print_nodes(nodes, filename):  #bioknowledgeReviewer
    """
    This function saves Monarch nodes into a CSV file.
    :param nodes: nodes list
    :param filename: file name without path and extension, e.g. 'monarch_nodes'
    :return: None object
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path):
        os.makedirs(path)
    #with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
    #    f.write('{}\n'.format('\n'.join(list(nodes))))
    file_path = f'{path}/{filename}_v{today}.csv'
    with open(file_path, 'w') as f:
        #f.write('\n'.join([str(node) for node in nodes]))
        f.write('{}\n'.format('\n'.join(list(nodes))))  # modified for testing

    logging.info(f"Nodes saved to CSV file: {file_path}")
    sample_nodes = nodes[:5]
    logging.debug("Sample nodes saved to CSV:")
    for node in sample_nodes:
        logging.debug(node)
    logging.info(f"File '{file_path}' saved.")

    # commented out since the file is saved elsewhere in the function
    #return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def expand_edges(seed_list):  #bioknowledgeReviewer
    """
    This function returns the Monarch network expanded with the first layer of neighbors from a list of query nodes.
    This function builds monarch 1shell network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # get 1shell list of nodes or neighbors
    neighbours, relations = get_neighbours(seed_list)

    # network nodes:  seed + 1shell
    nodes = set(seed_list).union(neighbours)

    # get connections for network nodes
    network = get_connections(nodes)

    logging.info(f"Expanded network contains {len(network)} edges.")
    sample_edges = list(network)[:5]
    logging.debug("Sample edges from expanded network:")
    for edge in sample_edges:
        logging.debug(edge)

    return network


def orthopheno_expand_edges(seed_list):  #bioknowledgeReviewer
    """
    This function returns the Monarch network expanded with the orthologs and the ortholog associated phenotypes from
     a list of query nodes.
    This function builds monarch 1shell-animal network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # get ortholog-phenotypes list
    orthophenoList = get_orthopheno_list(seed_list)

    # network nodes:  seed + neighbors + orthologs + phenotypes
    nodes = set(seed_list).union(set(orthophenoList))

    # get connections for network nodes
    network = get_connections(nodes)

    logging.info(f"Orthopheno expanded network contains {len(network)} edges.")
    sample_edges = list(network)[:5]
    logging.debug("Sample edges from orthopheno expanded network:")
    for edge in sample_edges:
        logging.debug(edge)

    return network


# BUILD NETWORK

def build_edges(edges_df, edges_fname):  #bioknowledgeReviewer
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    ## variables
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            #record['subject_uri'] = tuple[2]  # removed for testing
            #record['subject_category'] = tuple[3]
            #record['relation_id'] = tuple[4]  # tuple id modified according to previous steps
            record['relation_id'] = tuple[2]
            #record['relation_label'] = tuple[5]
            record['relation_label'] = tuple[3]
            #record['relation_uri'] = tuple[6]
            #record['object_id'] = tuple[7]
            record['object_id'] = tuple[4]
            #record['object_label'] = tuple[8]
            record['object_label'] = tuple[5]
            #record['object_uri'] = tuple[9]
            #record['object_category'] = tuple[10]
            #record['reference_id_list'] = tuple[11]
            record['reference_id_list'] = tuple[6]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)


    ## generate static variable: uriPrefixes_dct (url references)  # PREVIOUS VERSION
    #uriPrefixes_dct = {
    #    'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  
    #    'react': 'http://reactome.org/content/detail/',  
    #    'zfin': 'http://zfin.org/',
    #    'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  
    #    'mgi': 'http://www.informatics.jax.org/accession/MGI:',  
    #    'flybase': 'http://flybase.org/reports/',
    #    'wormbase': 'http://www.wormbase.org/resources/paper/',
    #    'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
    #    'isbn-10': 'ISBN-10:',
    #    'isbn-13': 'ISBN-13:',
    #    'mondo': 'http://purl.obolibrary.org/obo/MONDO_',  
    #    'rgd': 'https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=', 
    #    'omim': 'http://omim.org/entry/',  
    #    'sgd_ref': 'https://db.yeastgenome.org/reference/',  
    #    'genereviews': 'https://www.ncbi.nlm.nih.gov/books/',  
    #    'omia': 'http://omia.angis.org.au/',  
    #    'hgnc': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:', 
    #}

    # generate static variable: uriPrefixes_dct (url references)
    uriPrefixes_dct = {
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  # 'http://identifiers.org/pubmed/',
        'react': 'http://reactome.org/content/detail/',  # 'http://identifiers.org/reactome/',
        'zfin': 'http://zfin.org/',
        'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  # 'http://identifiers.org/go.ref/GO_REF:',
        'mgi': 'http://www.informatics.jax.org/accession/MGI:',  # 'http://identifiers.org/mgi/MGI:'
        'flybase': 'http://flybase.org/reports/',
        'wormbase': 'http://www.wormbase.org/resources/paper/',
        'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
        'isbn-10': 'ISBN-10:',
        'isbn-13': 'ISBN-13:',
        #'isbn-10': 'https://www.wikidata.org/wiki/Special:BookSources/',
        #'isbn-13': 'https://www.wikidata.org/wiki/Special:BookSources/'
        'mondo': 'http://purl.obolibrary.org/obo/MONDO_',  # http://purl.obolibrary.org/obo/MONDO_0009026
        'rgd': 'https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=', \
        # https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=1600115
        'omim': 'http://omim.org/entry/',  # http://omim.org/entry/61527
        'sgd_ref': 'https://db.yeastgenome.org/reference/',  # https://db.yeastgenome.org/reference/S000124036
        'genereviews': 'https://www.ncbi.nlm.nih.gov/books/',  # https://www.ncbi.nlm.nih.gov/books/NBK1526/
        'omia': 'http://omia.angis.org.au/',  # http://omia.angis.org.au/000214/9913
        'hgnc': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:', \
        # https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:7132
        'orpha': 'Go to ORPHANET web site (https://www.orpha.net/) and in the search field introduce the Orpha number: '
                 'ORPHA:' # no entry in monarch for that edge
    }

    # generate static variable: dbPrefixes_dct (source/database references)
    dbPrefixes_dct = {
        'na': 'NA',
        'nan': 'NA',
        'mgi': 'http://www.informatics.jax.org/',
        'fb': 'http://flybase.org/',
        'rgd': 'http://rgd.mcw.edu/',
        'zfin': 'http://zfin.org/',
        'sgd': 'https://www.yeastgenome.org/',
        'hgnc': 'https://www.genenames.org/',
        'xenbase': 'http://www.xenbase.org/'
    }

    # provenance variables
    ref_text = 'This edge comes from the Monarch Initiative {}.'.format(curr_year) 
    ref_date = today

    ## build graph schema network edges data structure and save edges file
    # prepare dataframe = [ {} ... {} ], where every row = {} = concept
    edges_l = list()
    for edge in edges_df.itertuples():
        # edge or row is a tuple (named and ordered attributes)
        # edge.reference_id_list >> can be 1) np.nan (float type) or 2) str without "|" 3) str with "|"
        ref_s = str(edge.reference_id_list)

        # prepare reference_uri_list attribute
        ref_uri_l = list()
        # expand to uri or NA
        pmid_l = list()
        # reference_id list iteration
        for ref in ref_s.strip().split('|'):
            # NA or database
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("Warning:")
                    print('Detected a new reference database in Monarch not yet implemented in this module. '
                          'The new source should be added to the dictionary of databases.'
                          'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                    print("In the build_edges() method, update 'dbPrefixes_dct' dictionary with '{}'".format(ref))
                    print('The edge that includes this new reference database is {}'.format(edge))
                    print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                          'until the dictionary is updated.')
                ref_uri_l.append(ref_uri)
            # publication uri: pubmed_id or url
            else:
                pref, uriId = ref.split(':')
                # separate pmid from non pmid and detect:
                # pubmed_id
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                # url
                elif ref.lower().startswith('http'):
                    ref_uri_l.append(ref)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("Warning:")
                        print('Detected a new reference source in Monarch not yet implemented in this module. '
                              'The new source should be added to the dictionary of sources.'
                              'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                        print("In the build_edges() method, update 'uriPrefixes_dct' dictionary with '{}'".format(pref))
                        print('The edge that includes this new reference source is {}'.format(edge))
                        print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                              'until the dictionary is updated.')
                    ref_uri_l.append(ref_uri)
        # create multi-term pubmed url
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)

        # prepare edge attributes: sub_id, obj_id, rel_id, rel_label, rel_def, rel_uri
        sub_id = 'NA' if edge.subject_id is None or str(edge.subject_id) == 'nan' or str(edge.subject_id) == 'None' else edge.subject_id
        rel_id = 'NA' if edge.relation_id is None or str(edge.relation_id) == 'nan'  or str(edge.relation_id) == 'None' else edge.relation_id
        obj_id = 'NA' if edge.object_id is None or str(edge.object_id) == 'nan'  or str(edge.object_id) == 'None' else edge.object_id
        rel_label = 'NA' if edge.relation_label is None or str(edge.relation_label) == 'nan'  or str(edge.relation_label) == 'None' else edge.relation_label
        rel_uri = 'NA' if edge.relation_uri is None or str(edge.relation_uri) == 'nan'  or str(edge.relation_uri) == 'None' else edge.relation_uri
        rel_def = 'NA'
        # "or str(edge.subject_id) == 'None'" section missing in the original code, but I doubt it's a problem

        ## make sure there are no 'NA' relation ids
        #if rel_id == 'NA':
        #    if edge.object_category == 'genotype': # if object is genotype, relation is 'has genotype'
        #        rel_id = 'GENO:0000222'
        #        rel_label = 'has_genotype'
        #        rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000222'
        #    else:
        #        rel_id = 'RO:0002610' # else (gene, model objects), relation is 'correlated with'
        #        rel_label = 'correlated with'
        #        rel_uri = 'http://purl.obolibrary.org/obo/RO_0002610'
        #
        ## change every 'contributes to' to 'contributes to condition'
        #if rel_id == 'RO:0002326':
        #    rel_id = 'RO:0003304'
        #    rel_label = 'contributes to condition'
        #    rel_uri = 'http://purl.obolibrary.org/obo/RO_0003304'
        #
        ## change every 'is causal germline mutation in' and 'is causal germline mutation partially giving rise to' and 'pathogenic for condition' to 'causes condition'
        #if rel_id == 'RO:0004013' or rel_id == 'RO:0004016' or rel_id == 'GENO:0000840':
        #    rel_id = 'RO:0003303'
        #    rel_label = 'causes condition'
        #    rel_uri = 'http://purl.obolibrary.org/obo/RO_0003303'
        #
        ## change every 'in orthology relationship with' to 'in 1 to 1 orthology relationship with'
        #if rel_id == 'RO:HOM0000017':
        #    rel_id = 'RO:HOM0000020'
        #    rel_label = 'in 1 to 1 orthology relationship with'
        #    rel_uri = 'http://purl.obolibrary.org/obo/RO_HOM0000020'  
        #
        ## if 'genotype' -->'has phenotype' --> 'disease', change 'has phenotype' to 'has role in modelling'
        #if edge.subject_category == 'genotype' and edge.object_category == 'disease' and rel_id == 'RO:0002200':
        #    rel_id = 'RO:0003301'
        #    rel_label = 'has role in modeling'
        #    rel_uri = 'http://purl.obolibrary.org/obo/RO_0003301'
        #
        ## change every 'is reference allele of' to 'is allele of'
        #if rel_id == 'GENO:0000610':
        #    rel_id = 'GENO:0000408'
        #    rel_label = 'is_allele_of'
        #    rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000408'
        # (!!!) whole section removed for testing

        # build the data structure = list of edges as list of dict, where a dict is an edge
        edge = dict()
        edge['subject_id'] = sub_id
        edge['object_id'] = obj_id
        edge['property_id'] = rel_id
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] = rel_uri
        edge['reference_uri'] = ref_uri_list
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edges_l.append(edge)

    # save edges file
    df = pd.DataFrame(edges_l)
    #print('df',df.shape)
    #path = os.getcwd() + '/monarch'
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', \
             'property_label', 'property_description', 'property_uri']]
    df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,edges_fname,today), index=False)
    #print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    #print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    #print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    #print('\nThe Monarch network edges are built and saved at: {}/monarch_edges_v{}.csv\n'.format(path,today))
    #print('\nFinished build_edges().\n')
    logging.info("Edges DataFrame created for building network.")
    logging.info(f"DataFrame size: {df.shape}")
    logging.info(f"DataFrame columns: {df.columns.tolist()}")
    sample_df = df.head(1)
    logging.debug(f"Sample data from the DataFrame:\n{sample_df.to_string(index=False)}")
    file_path = f'{path}/{edges_fname}_v{today}.csv'
    df.fillna('NA').to_csv(file_path, index=False)
    logging.info(f"Monarch network edges file has been built and saved.")
    logging.info(f"File saved at: {file_path}")
    logging.info(f"DataFrame size (rows, columns): {df.shape}")
    logging.info(f"DataFrame columns: {df.columns.tolist()}")
    logging.debug(f"First record in the DataFrame:\n{sample_df.to_string(index=False)}")

    #return df   # check if required down the line to be in df format
    return edges_l  # in the original bioknowledge reviewer


def build_nodes(edges_df, nodes_fname):  #bioknowledgeReviewer
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    # build semantic groups dictionary
    # collide concepts in a concept dict
    concept_dct = dict()

    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            #record['subject_uri'] = tuple[2]
            #record['subject_category'] = tuple[3]
            #record['relation_id'] = tuple[4]
            record['relation_id'] = tuple[2]
            record['relation_label'] = tuple[3]
            #record['relation_label'] = tuple[5]
            #record['relation_uri'] = tuple[6]
            #record['object_id'] = tuple[7]
            record['object_id'] = tuple[4]
            record['object_label'] = tuple[5]
            #record['object_label'] = tuple[8]
            #record['object_uri'] = tuple[9]
            #record['object_category'] = tuple[10]
            #record['reference_id_list'] = tuple[11]
            record['reference_id_list'] = tuple[6]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        sid = edge.subject_id #fields[0]
        oid = edge.object_id #fields[4]
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))

    ## build concept attributes dict: id integration of sub and obj IDs in a common data structure
    #concept_dct = dict()
    #
    ## read edges from variable
    #for edge in edges_df.itertuples():
    #    #fields = list(edge_tuple)
    #    # id: integration of sub and obj IDs in a unique data structure
    #    sid = edge.subject_id 
    #    slab = edge.subject_label 
    #    sgroup = edge.subject_category
    #    suri = edge.subject_uri
    #    oid = edge.object_id 
    #    olab = edge.object_label 
    #    ogroup = edge.object_category
    #    ouri = edge.object_uri
    #    
    #    # if 'genotype' -->'is allele of' --> 'gene', change 'genotype' to 'variant'
    #    if sgroup == 'genotype' and ogroup == 'gene' and edge.relation_id == 'GENO:0000408':
    #        sgroup = 'variant'
    #    
    #    # build the concept data structure
    #    concept_dct[sid] = {'preflabel': slab,
    #                        #'name': 'NA',
    #                        'semantic_groups': sgroup,
    #                        'uri': suri,
    #                        'synonyms': 'NA', 
    #                        'description': 'NA'}
    #    concept_dct[oid] = {'preflabel': olab,
    #                        #'name': 'NA',
    #                        'semantic_groups': ogroup,
    #                        'uri': ouri,
    #                        'synonyms': 'NA', 
    #                        'description': 'NA'}
    #
    ## biothings: annotate name,synonyms,description to genes
    #print('\nAdding BioThings annotation: name, synonyms, description...')
    ## input: (preflabel) symbol,alias
    #symbols = list()
    #for concept in concept_dct:
    #    if isinstance(concept_dct[concept]['semantic_groups'], list):
    #        for label in concept_dct[concept]['semantic_groups']:
    #            if 'gene' in label:
    #                symbols.append(concept_dct[concept]['preflabel'])
    #    else:
    #        if 'gene' in concept_dct[concept]['semantic_groups']:
    #            symbols.append(concept_dct[concept]['preflabel'])
    #
    #print('symbols:', len(symbols))
    # (!!!) removed for testing, substituted with original below

        # list of concept prefixes with dict
    conceptPrefix_dct = dict()
    for concept in concept_dct:
        conceptPrefix_dct[concept.split(':')[0]] = 1
    print('Number of nodes CURIEs: {}'.format(len(conceptPrefix_dct.keys())))
    print('List of nodes CURIEs: {}'.format(conceptPrefix_dct.keys()))

    # build conceptPrefix2semantic dict
    conceptPrefix2semantic_dct = dict()
    for prefix in conceptPrefix_dct:
        prefix = prefix.lower()
        if 'variant' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'VARI'
        elif 'phenotype' in prefix or 'mondo' in prefix or 'omim' in prefix or 'doid' in prefix or 'mesh' in prefix or 'hp' in prefix or 'mp' in prefix or 'fbcv' in prefix or 'fbbt' in prefix or 'zp' in prefix or 'apo' in prefix or 'trait' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'DISO'
        elif 'gene' in prefix or 'hgnc' in prefix or 'ensembl' in prefix or 'mgi' in prefix or 'flybase' in prefix or 'wormbase' in prefix or 'xenbase' in prefix or 'zfin' in prefix or 'rgd' in prefix or 'sgd' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENE'
        elif 'react' in prefix or 'kegg-path' in prefix or 'go' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'PHYS'
        elif 'uberon' in prefix or 'cl' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'ANAT'
        elif 'geno' in prefix or 'coriell' in prefix or 'monarch' in prefix or 'mmrrc' in prefix or '' in prefix \
                or 'bnode' in prefix:
            conceptPrefix2semantic_dct[prefix] = 'GENO'
        else:
            conceptPrefix2semantic_dct[prefix] = 'CONC'

    # build concept attributes dict: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()
    # read edges from file
    #header = 1
    ##for row in open('../monarch/1shell-animal-hgnc/get-monarch-connections/monarch_connections.tsv').readlines():
    #for row in open('{}'.format(csv_path)).readlines():
    #    if header:
    #        header = 0
    #        continue
    #    fields = row.strip('\n').split('\t')
    # read edges from variable
    for edge in edges_df.itertuples():
        #fields = list(edge_tuple)
        # id: integration of sub and obj IDs in a unique data structure
        sid = edge.subject_id #fields[0]
        slab = edge.subject_label #fields[1]
        oid = edge.object_id #fields[4]
        olab = edge.object_label #fields[5]
        # build the concept data structure
        concept_dct[sid] = {'preflabel': slab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(sid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            'semantic_groups': conceptPrefix2semantic_dct.get(oid.split(':')[0].lower(), 'CONC'),
                            'synonyms': 'NA', 'description': 'NA'}

    # build graph schema network nodes data structure and save nodes file
    # biothings: annotate name,synonyms,description to genes
    print('\nAdding BioThings annotation: gene name, synonyms, description...')
    # input: (preflabel) symbol,alias
    symbols = list()
    for concept in concept_dct:
        if isinstance(concept_dct[concept]['semantic_groups'], list):
            for label in concept_dct[concept]['semantic_groups']:
                if 'GENE' in label:
                    symbols.append(concept_dct[concept]['preflabel'])
        else:
            if 'GENE' in concept_dct[concept]['semantic_groups']:
                symbols.append(concept_dct[concept]['preflabel'])
    print('symbols:', len(symbols))
    # original code above since previous (!!!) comment

    # query biothings
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name,alias,summary', size=1, as_dataframe=True)

    # dictionary: {symbol:name}
    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
    ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
    monarch_s2n = dict(zip(ids.symbol, ids.name))
    monarch_s2s = dict(zip(ids.symbol, ids.synonyms))
    monarch_s2d = dict(zip(ids.symbol, ids.description))

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        preflabel = concept_dct[concept]['preflabel']
        concept_dct[concept]['synonyms'] = monarch_s2s[preflabel] if preflabel in monarch_s2s.keys() else 'NA'
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['uri'] = concept_dct[concept]['uri']
        node['preflabel'] = preflabel
        node['name'] = monarch_s2n[preflabel] if preflabel in monarch_s2n.keys() else preflabel
        node['synonyms'] = '|'.join(list(concept_dct[concept]['synonyms'])) if isinstance(
            concept_dct[concept]['synonyms'], list) else concept_dct[concept]['synonyms']
        node['description'] = monarch_s2d[preflabel] if preflabel in monarch_s2d.keys() else 'NA'
        nodes_l.append(node)

    # save nodes file
    df = pd.DataFrame(nodes_l)
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/monarch'
    #df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,nodes_fname, today), index=False)
    file_path = f'{path}/{nodes_fname}_v{today}.csv'
    df.fillna('NA').to_csv(file_path, index=False)
    #print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    #print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    #print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    #print('\nThe Monarch network nodes are built and saved at: {}/monarch_nodes_v{}.csv\n'.format(path,today))
    #print('\nFinished build_nodes().\n')
    logging.info("Nodes DataFrame created for building network.")
    logging.info(f"DataFrame size: {df.shape}")
    logging.info(f"DataFrame columns: {df.columns.tolist()}")
    sample_df = df.head(1)
    logging.debug(f"Sample data from the DataFrame:\n{sample_df.to_string(index=False)}")
    logging.info(f"Monarch network nodes file has been built and saved.")
    logging.info(f"File saved at: {file_path}")

    #return df   # check if required down the line to be in df format
    return nodes_l  # in the original bioknowledge reviewer


def get_symptoms_disease(disease_id, monarch_edges_csv, monarch_nodes_csv):
    """
    This function finds all symptoms of a disease.
    :param disease_id: id of the disease of which symptoms and drugs should be found (e.g. 'MONDO:0007739' for HD)  
    :param monarch_edges_csv: name of the edges csv file from Monarch
    :param monarch_nodes_csv: name of the nodes csv file from Monarch
    :return: a list of symptom names
    """

    logging.info(f"NOW RUNNING: {current_function_name()}. This might take a while.")

    # open csv file
    fname_edges = monarch_edges_csv
    fname_nodes = monarch_nodes_csv
    edges_df = pd.read_csv(fname_edges)
    nodes_df = pd.read_csv(fname_nodes)

    # find all nodes that are one step away from disease and have as relation 'has phenotype' (='RO:0002200')
    df_symptoms = edges_df.loc[(edges_df['subject_id'] == disease_id) & (edges_df['property_id'] == 'RO:0002200')]
    symptoms_id_lst = df_symptoms['object_id'].tolist()

    #get names of symptoms
    symptoms_name_lst = nodes_df.loc[nodes_df['id'].isin(symptoms_id_lst), 'preflabel'].to_list()

    logging.info(f"Total number of symptoms found for disease ID {disease_id}: {len(symptoms_name_lst)}")
    sample_symptoms = symptoms_name_lst[:5]
    logging.info("Sample symptoms (or complete list if fewer than 5):")
    for symptom in sample_symptoms:
        logging.info(symptom)

    return symptoms_name_lst


def symptom_list_specified_folder(input_folder = 'Huntington disease (2022-06-23)'):
    """
    This function finds the symptom list of the disease specified by user 
    clicking a folder and changes the cwd to the folder specified.
    :param input_folder: The input folder of the disease
    :return: symptom name list
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    date_str = input_folder[-11:-1] # get the date from the disease name
    date = datetime.datetime.strptime(date_str, '%Y-%m-%d').date()
    logging.info(f"Processing folder: {input_folder} with date: {date}")

    #change directory to the folder
    #os.chdir(os.getcwd() +'/drugapp/data/' + input_folder)
    new_path = os.getcwd() +'/drugapp/data/' + input_folder
    os.chdir(new_path)
    logging.info(f"Changed directory to: {new_path}")

    # find the disease ID that was previously written in a txt file
    #with open(os.getcwd()+"/disease_id.txt", 'r') as file:
    #    file_text = file.read()
    #    input_number, disease_id, disease_name = file_text.split(';')
    disease_id_file_path = os.path.join(new_path, "disease_id.txt")
    with open(disease_id_file_path, 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')
    logging.info(f"Read disease information: ID={disease_id}, Name={disease_name}, Input Number={input_number}")

    monarch_fname_edges = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    
    # get list of symptoms for disease to ask user
    symptoms_name_lst=get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)

    logging.info(f"Symptoms list size for {disease_name}: {len(symptoms_name_lst)}")
    sample_symptoms = symptoms_name_lst[:5]
    logging.info("Sample symptoms (or complete list if fewer than 5):")
    for symptom in sample_symptoms:
        logging.info(symptom)

    return symptoms_name_lst, date, input_number


def symptom_list_today():
    """
    This function finds the symptom list of the disease specified by user on 
    that same day.
    :return: symptom name list, date of creation of files required, disease name
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")

    date=today
    logging.info(f"Generating symptom list for today: {date}")
    
    # find the disease ID that was previously written in a txt file
    #with open(os.getcwd()+"/disease_id.txt", 'r') as file:
    #    file_text = file.read()
    #    input_number, disease_id, disease_name = file_text.split(';')
    disease_id_file_path = os.path.join(os.getcwd(), "disease_id.txt")
    with open(disease_id_file_path, 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')
    logging.info(f"Read disease information for today: ID={disease_id}, Name={disease_name}, Input Number={input_number}")

    #disease_name_date = disease_name+' ({})'.format(date)
    disease_name_date = f"{disease_name} ({date})"

    monarch_fname_edges = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    
    # get list of symptoms for disease to ask user
    symptoms_name_lst=get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)

    logging.info(f"Symptoms list size for {disease_name_date}: {len(symptoms_name_lst)}")
    sample_symptoms = symptoms_name_lst[:5]
    logging.info("Sample symptoms (or complete list if fewer than 5):")
    for symptom in sample_symptoms:
        logging.info(symptom)
    
    return symptoms_name_lst, date, disease_name_date


#def run_monarch(input_number = '143100'):
def run_monarch(input_id = 'MONDO:0007739'):
    """
    This function runs the whole Monarch script and saves nodes and edges files.

    :param input_number: The input phenotype MIM number of the disease
    :return: nodes and edges files in /monarch folder
    """

    logging.info(f"NOW RUNNING: {current_function_name()}.")
    
    ## turn input number into input ID
    #input_id = 'OMIM:'+input_number
    # NOTE: changed so that it can take different kinds of URI
        
    # get the disease name (str)
    disease_name, disease_id = get_disease_name_id(disease_input_ID = input_id)
    
    # make folder for HD and change cwd to this folder
    new_path = os.getcwd() + '/drugapp/data/'+disease_name+' ({})'.format(today)

    # check if directory exists, if not create it, then change directory to it
    #os.makedirs(new_path)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    os.chdir(new_path)

    # save disease_ID as text file
    text_file = open(os.getcwd()+"/disease_id.txt", "w")
    #text_file.write(input_number + ';' + disease_id + ';'+ disease_name)
    text_file.write(disease_id + ';'+ disease_name)
    text_file.close()
    
    #build monarch network
    seedList = [input_id]

    neighbourList = get_neighbours_list(seedList)
    orthophenoList = get_orthopheno_list(seedList) 
    geneList = sum([seedList,neighbourList,orthophenoList], [])
    network = extract_edges(geneList) 

    # save network
    print_network(network, 'monarch_orthopeno_network_disease')
    
    file = 'monarch_orthopeno_network_disease_v{}.csv'.format(today)
    monarch_connections = read_connections(file)
    build_edges(monarch_connections, edges_fname = 'monarch_edges_disease')
    build_nodes(monarch_connections, nodes_fname = 'monarch_nodes_disease')


def run_monarch_symptom(input_symptom, date):
    """
    This function runs the whole Monarch script using the disease phenotype MIM number
    and symptom ID as seeds and saves nodes and edges files.

    :param input_symptom: The input symptom
    :return: nodes and edges files in /monarch folder
    """

    logging.info(f"NOW RUNNING: {current_function_name()}")

    # turn input symptom name into ID
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    nodes_df = pd.read_csv(monarch_fname_nodes)
    symptom_id = nodes_df.loc[nodes_df['preflabel'] == input_symptom, 'id'].iloc[0]
    
    #build monarch network
    seedList = [symptom_id]

    neighbourList = get_neighbours_list(seedList)
    orthophenoList = get_orthopheno_list(seedList) 
    geneList = sum([seedList,neighbourList,orthophenoList], [])
    network = extract_edges(geneList) 

    # save network
    print_network(network, 'monarch_orthopeno_network_symptom')
    
    file = 'monarch_orthopeno_network_symptom_v{}.csv'.format(today)
    monarch_connections = read_connections(file)
    build_edges(monarch_connections, edges_fname='monarch_edges_symptom')
    build_nodes(monarch_connections, nodes_fname = 'monarch_nodes_symptom')



if __name__ == '__main__':
    
    input_number = '143100' # input of user
    run_monarch(input_number)



