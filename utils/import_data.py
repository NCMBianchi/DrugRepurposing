"""
IMPORT DATA UTILITY MODULE: Loads and formats data from CSV files
Created on March 10th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os, datetime, requests
import pandas as pd

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def df_to_nodes(nodes_df):
    """Convert DataFrame nodes to pipeline format (list of dicts)"""
    return [{'id': row['id'], 'label': row['label']} for _, row in nodes_df.iterrows()]


def df_to_edges(edges_df):
    """Convert DataFrame edges to pipeline format (list of lists of dicts)"""
    edges = []
    for _, row in edges_df.iterrows():
        edge = [
            {'id': row['subject_id'], 'label': row['subject_label']},
            {'label': row['relation']},
            {'id': row['object_id'], 'label': row['object_label']},
            {'notes': row['notes']}
        ]
        edges.append(edge)
    return edges


def nodes_to_df(nodes_list):
    """Convert pipeline format nodes (list of dicts) to DataFrame"""
    return pd.DataFrame(nodes_list)


def edges_to_df(edges_list):
    """Convert pipeline format edges (list of lists of dicts) to DataFrame"""
    edges_df = pd.DataFrame([
        {
            'subject_id': edge[0]['id'],
            'subject_label': edge[0]['label'],
            'relation': edge[1]['label'],
            'object_id': edge[2]['id'],
            'object_label': edge[2]['label'],
            'notes': edge[3]['notes']
        }
        for edge in edges_list
    ])
    return edges_df


def fetch_disease_name(seed_id):
    """
    Retrieve the name of a specific Monarch Initiative entity.
    
    :param seed_id: A Monarch Initiative ID (e.g., 'MONDO:0007739')
                 https://api-v3.monarchinitiative.org/v3/api/entity/MONDO:0007739
    :return: the name of the entity or the seed_id if retrieval fails.
    """
    response = requests.get(f'https://api-v3.monarchinitiative.org/v3/api/entity/{seed_id[0]}')
    response.raise_for_status()
    data = response.json()
        
    return data.get('name')


def load_and_convert_data(custom_paths=None):
    """
    Load data from CSV files and convert to pipeline format.
    
    :param custom_paths: custom paths to use instead of generating from date and disease.
    :return: dictionary containing all loaded data in both DataFrame and list formats.
    """
    # initialize data dictionary
    data = {}
    
    # check which provided files exist
    path_exists = {}
    for path_type, path in custom_paths.items():
        path_exists[path_type] = os.path.exists(path)
    
    # count how many files exist
    files_found = sum(path_exists.values())
    total_files = len(custom_paths)
    
    # print summary
    print(f"File existence status: {files_found}/{total_files} found")
    for path_type, exists in path_exists.items():
        print(f"  {path_type}: {'✓' if exists else '✗'} - {custom_paths[path_type]}")
    
    # load DataFrames for paths that exist
    try:
        for path_type, path in custom_paths.items():
            if path_exists[path_type]:
                # load file as DataFrame
                df_key = f"{path_type.replace('_path', '')}_df"
                data[df_key] = pd.read_csv(path)
                
                # convert to pipeline format
                # assume path_type ends with '_nodes_path' or '_edges_path'
                if 'nodes_path' in path_type:
                    data_key = path_type.replace('_path', '')
                    data[data_key] = df_to_nodes(data[df_key])
                elif 'edges_path' in path_type:
                    data_key = path_type.replace('_path', '')
                    data[data_key] = df_to_edges(data[df_key])
        
        # print summary
        print("\nLoaded data summary:")
        for key, value in data.items():
            if isinstance(value, pd.DataFrame):
                print(f"  {key}: {len(value)} rows")
            else:
                print(f"  {key}: {len(value)} items")
        
    except Exception as e:
        print(f"Error loading or converting data: {e}")
    
    return data