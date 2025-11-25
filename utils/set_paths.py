"""
SET_PATHS UTILITY MODULE: set the path dictionary
Created on February 21st 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import sys,os,requests

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def initialise_disease_directories(today_date,input_disease_id):
    """
    Initialises and returns the paths for the disease-specific directory and Monarch directory.
    
    :param today_directory: the directory for today's date.
    :param disease_name: the name of the disease.
    :return: a dictionary containing paths for the disease directory and sub directories.
    """
    sys.path.append(os.getcwd())

    input_disease_name = fetch_disease_name(input_disease_id)
    date_string = today_date.strftime('%Y-%m-%d')
    
    base_data_directory = os.path.join(os.getcwd(), 'data')
    today_directory = os.path.join(base_data_directory, date_string)
    
    disease_directory = os.path.join(today_directory, f'{input_disease_name} ({date_string})')
    monarch_directory = os.path.join(disease_directory, 'monarch')
    dgidb_directory = os.path.join(disease_directory, 'dgidb')
    drugsimil_directory = os.path.join(disease_directory, 'drug_similarity')
    negsamples_directory = os.path.join(disease_directory, 'negative_samples')
    embeddingsXGB_directory = os.path.join(disease_directory, 'embeddings_xgb')
    embeddingsPYG_directory = os.path.join(disease_directory, 'embeddings_pyg')
    networkXGB_directory = os.path.join(disease_directory, 'network_xgb')
    networkPYG_directory = os.path.join(disease_directory, 'network_pyg')

    # create the disease directory and its subdirectories if they don't exist
    os.makedirs(base_data_directory, exist_ok=True)
    os.makedirs(today_directory, exist_ok=True)
    os.makedirs(disease_directory, exist_ok=True)
    os.makedirs(monarch_directory, exist_ok=True)
    os.makedirs(dgidb_directory, exist_ok=True)
    os.makedirs(drugsimil_directory, exist_ok=True)
    os.makedirs(negsamples_directory, exist_ok=True)
    os.makedirs(embeddingsXGB_directory, exist_ok=True)
    os.makedirs(embeddingsPYG_directory, exist_ok=True)
    os.makedirs(networkXGB_directory, exist_ok=True)
    os.makedirs(networkPYG_directory, exist_ok=True)

    return {
        'disease_id': input_disease_id,
        'disease_name': input_disease_name,
        'today_date': today_date,
        'date_string': date_string,
        'base_data_directory': base_data_directory,
        'today_directory': today_directory,
        'disease_directory': disease_directory,
        'monarch_directory': monarch_directory,
        'dgidb_directory': dgidb_directory,
        'drugsimil_directory': drugsimil_directory,
        'negsamples_directory': negsamples_directory,
        'embeddingsXGB_directory': embeddingsXGB_directory,
        'embeddingsPYG_directory': embeddingsPYG_directory,
        'networkXGB_directory': networkXGB_directory,
        'networkPYG_directory': networkPYG_directory,
    }


def initialise_test_directories(today_date, input_disease_id):
    """
    Initialises and returns the paths for test runs performed in 'test-runs.ipynb'.
    
    :param today_directory: the directory for today's date.
    :param disease_name: the name of the disease.
    :return: a dictionary containing paths for the test directory and sub directories.
    """
    sys.path.append(os.getcwd())

    input_disease_name = fetch_disease_name(input_disease_id)
    date_string = today_date.strftime('%Y-%m-%d')
    
    base_data_directory = os.path.join(os.getcwd(), 'data')
    today_directory = os.path.join(base_data_directory, date_string)
    
    test_directory = os.path.join(today_directory, f'TEST ({date_string})')
    monarch_directory = os.path.join(test_directory, 'monarch')
    dgidb_directory = os.path.join(test_directory, 'dgidb')
    drugsimil_directory = os.path.join(test_directory, 'drug_similarity')
    negsamples_directory = os.path.join(test_directory, 'negative_samples')
    embeddingsXGB_directory = os.path.join(test_directory, 'embeddings_xgb')
    embeddingsPYG_directory = os.path.join(test_directory, 'embeddings_pyg')
    networkXGB_directory = os.path.join(test_directory, 'network_xgb')
    networkPYG_directory = os.path.join(test_directory, 'network_pyg')

    # create the test directory and its subdirectories if they don't exist
    os.makedirs(base_data_directory, exist_ok=True)
    os.makedirs(today_directory, exist_ok=True)
    os.makedirs(test_directory, exist_ok=True)
    os.makedirs(monarch_directory, exist_ok=True)
    os.makedirs(dgidb_directory, exist_ok=True)
    os.makedirs(drugsimil_directory, exist_ok=True)
    os.makedirs(negsamples_directory, exist_ok=True)
    os.makedirs(embeddingsXGB_directory, exist_ok=True)
    os.makedirs(embeddingsPYG_directory, exist_ok=True)
    os.makedirs(networkXGB_directory, exist_ok=True)
    os.makedirs(networkPYG_directory, exist_ok=True)

    return {
        'disease_id': input_disease_id,
        'disease_name': input_disease_name,
        'today_date': today_date,
        'date_string': date_string,
        'base_data_directory': base_data_directory,
        'today_directory': today_directory,
        'test_directory': test_directory,
        'monarch_directory': monarch_directory,
        'dgidb_directory': dgidb_directory,
        'drugsimil_directory': drugsimil_directory,
        'negsamples_directory': negsamples_directory,
        'embeddingsXGB_directory': embeddingsXGB_directory,
        'embeddingsPYG_directory': embeddingsPYG_directory,
        'networkXGB_directory': networkXGB_directory,
        'networkPYG_directory': networkPYG_directory,
    }


def fetch_disease_name(seed_id):
    """
    Retrieve the name of a specific Monarch Initiative entity.
    
    :param seed_id: a Monarch Initiative ID (e.g., 'MONDO:0007739')
                 https://api-v3.monarchinitiative.org/v3/api/entity/MONDO:0007739
    :return: the name of the entity or the seed_id if retrieval fails.
    """
    response = requests.get(f'https://api-v3.monarchinitiative.org/v3/api/entity/{seed_id[0]}')
    response.raise_for_status()
    data = response.json()
        
    return data.get('name')