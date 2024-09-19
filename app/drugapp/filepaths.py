"""
MAIN MODULE: FILE PATHS GENERATION 
Created on September 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os

def initialise_base_directories(date_string):
    """
    Initialises and returns the paths for the base directory and today's directory.
    
    :param date_string: A string representing the date (e.g., '2024-09-03').
    :return: A dictionary containing paths for the base directory and today's directory.
    """
    base_data_directory = os.path.join(os.getcwd(), 'drugapp', 'data')
    today_directory = os.path.join(base_data_directory, date_string)
    
    # Create today's directory if it doesn't exist
    os.makedirs(today_directory, exist_ok=True)
    
    return {
        'base_data_directory': base_data_directory,
        'today_directory': today_directory
    }

def initialise_disease_directories(today_directory, disease_name, date_string):
    """
    Initialises and returns the paths for the disease-specific directory and Monarch directory.
    
    :param today_directory: The directory for today's date.
    :param disease_name: The name of the disease.
    :param date_string: The date string (for archival purposes).
    :return: A dictionary containing paths for the disease directory and sub directories.
    """
    disease_directory = os.path.join(today_directory, f'{disease_name} ({date_string})')
    monarch_directory = os.path.join(disease_directory, 'monarch')
    dgidb_directory = os.path.join(disease_directory, 'dgidb')
    drugsimil_directory = os.path.join(disease_directory, 'drug_similarity')
    network_directory = os.path.join(disease_directory, 'network')
    NSnetwork_directory = os.path.join(disease_directory, 'network_with_NS')

    # Create the disease directory and its subdirectories if they don't exist
    os.makedirs(disease_directory, exist_ok=True)
    os.makedirs(monarch_directory, exist_ok=True)
    os.makedirs(dgidb_directory, exist_ok=True)
    os.makedirs(drugsimil_directory, exist_ok=True)
    os.makedirs(network_directory, exist_ok=True)
    os.makedirs(NSnetwork_directory, exist_ok=True)

    return {
        'disease_directory': disease_directory,
        'monarch_directory': monarch_directory,
        'dgidb_directory': dgidb_directory,
        'drugsimil_directory': drugsimil_directory,
        'network_directory': network_directory,
        'NSnetwork_directory': NSnetwork_directory
    }

