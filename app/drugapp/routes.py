"""
MAIN MODULE: FLASKAPP FUNCTIONALITIES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

# for setup, logging and date-override
import sys,os,platform,datetime,logging,builtins,time,multiprocessing
from drugapp.platform import *

# for service communication in Docker Swarm
import requests
import json
from urllib.parse import urljoin

# for Flask App operations
from flask import Flask, render_template, url_for, request, session, redirect
from drugapp import app 
from drugapp.forms import *
from drugapp.Monarch import *
from drugapp.dgidb import *
from drugapp.drugsimilarity import *
from drugapp.negsamples import *
from drugapp.networkmodel import *
from drugapp.filepaths import *
from drugapp.unique import *

## PLATFORM INFO
python_executable = sys.executable
python_version = platform.python_version()
num_cores = multiprocessing.cpu_count()

## SERVICE CONFIGURATION
SERVICE_URLS = {
    'monarch': os.environ.get('MONARCH_SERVICE_URL', 'http://monarch:5000'),
    'dgidb': os.environ.get('DGIDB_SERVICE_URL', 'http://dgidb:5000'),
    'drugsim': os.environ.get('DRUGSIM_SERVICE_URL', 'http://drugsim:5000'),
    'negsample': os.environ.get('NEGSAMPLE_SERVICE_URL', 'http://negsample:5000'),
    'networkmodel': os.environ.get('NETWORKMODEL_SERVICE_URL', 'http://networkmodel:5000')
}

def call_service(service_name, endpoint, payload):
    """
    Generic method to call microservices in Docker Swarm
    """
    try:
        service_url = urljoin(SERVICE_URLS[service_name], endpoint)
        response = requests.post(service_url, json=payload)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        logging.error(f"Service {service_name} call failed: {e}")
        raise

@app.route("/", methods = ['GET', 'POST'])
@app.route("/home", methods = ['GET', 'POST'])
def config():
    
    form = user_input()

    if form.validate_on_submit():

        payload = {
            'd_toggle': form.date_t.data,  # Date override toggle
            'date_OR_day': form.date_OR_day.data,
            'date_OR_month': form.date_OR_month.data,
            'date_OR_year': form.date_OR_year.data,
            'input_seed': form.disease_URI.data,  # disease URI
            'run_layers': form.deg_of_dist.data,  # degree of distance
            'input_min_simil': float(form.inp_minimum_sim.data),  # minimum drug similarity
            'negs_toggle': form.ns_toggle.data,  # negative samples toggle
            'sim_threshold': float(form.sim_t.data),  # drug similarity threshold
            'num_jobs': form.n_cores.data,  # number of CPU cores
            'depth_input': form.po_mode.data,  # mode of operation
            'seed_input': form.ML_seed.data  # ML seed
        }

        # generate the override date if toggle is active
        if payload['d_toggle']:
            d_override = datetime.date(
                payload['date_OR_year'], 
                payload['date_OR_month'], 
                payload['date_OR_day']
            )
        else:
            d_override = None

        ## DATES
        today, actual_today, date_str, curr_year, overall_start_time, formatted_start_time = set_date(payload['d_toggle'], d_override)
        logging.info(f"'DrugRepurposing' pipeline started at {formatted_start_time}.\n")

        ## DIRECTORY NAMES
        base_directories = initialise_base_directories(date_str)

        ## LOGGING: set-up
        platform_filename = "running_platform.txt"
        platform_file_path = os.path.join(base_directories['today_directory'], platform_filename)
        log_filename = "logfile.log"
        log_file_path = os.path.join(base_directories['today_directory'], log_filename)
        input_filename = "inputs.txt"
        input_file_path = os.path.join(base_directories['today_directory'], input_filename)
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)  # highest level to capture all logs
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.ERROR)  # only log errors to the console
        console_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        if not logger.handlers:
            logger.addHandler(file_handler)
            logger.addHandler(console_handler)

        ## LOGGING: platform info
        with open(platform_file_path, 'w') as file:
            file.write(f"Run date (real): {actual_today}\n\n")
            file.write(f"Python kernel: {python_executable}\n")
            file.write(f"Python version: {python_version}\n")
            file.write(f"The current machine has {num_cores} CPU cores.\n")
            file.write("\n\n")

        ## LOGGING: inputs info
        with open(input_file_path, 'w') as file:
            file.write("\n\n")

        ## PACKAGES: python-included packages
        import tqdm
        import re
        import json  # to save .json files returned by the API
        import logging  # to check for running variables
        import inspect  # to check for running functions
        import shutil  # to properly copy and rename files in run_monarch_mock()
        import ast  # for syntax trees in 'DGIdb.py' steps
        from unittest.mock import Mock  # to reverse requests.Response()
        import pickle

        ## LOGGING: only non-core packages are loaded with the 'import_and_log()' function:
        import_and_log('pandas')  # for general object handling
        import_and_log('numpy')  # for general operations
        import_and_log('requests')  # for several API calls
        import_and_log('networkx')  # for ' Monarch.py' to plot the network
        import_and_log('matplotlib')  # for ' Monarch.py' to plot the network
        import_and_log('biothings_client','get_client') # for several API calls
        import_and_log('rdkit','Chem')  # for 'drugsimilarity.py'
        import_and_log('rdkit','RDLogger')  # for 'drugsimilarity.py'
        import_and_log('rdkit','DataStructs')  # for 'drugsimilarity.py'
        import_and_log('SPARQLWrapper','SPARQLWrapper') # for 'drugsimilarity.py'
        import_and_log('SPARQLWrapper','JSON') # for 'drugsimilarity.py'
        import_and_log('node2vec','Node2Vec')  # for 'network_model.py'
        import_and_log('sklearn')  # for 'network_model.py'
        import_and_log('xgboost','XGBClassifier')  # for 'network_model.py'
        import_and_log('scipy','stats')  # for 'network_model.py'
        import_and_log('IPython','display')  # for 'network_model.py'

        ## PACKAGES: for ease of use, some packages are reloaded with an alias
        import pandas as pd
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt
        mc = get_client('chem')
        # these submodules are imported normally as the main versions are already logged,
        # and importlib.import_module() doesn't handle well submodules
        from rdkit.Chem import AllChem  # older version, do not place in REQUIREMENTS.txt
        from rdkit.Chem import Draw
        from rdkit.Chem import rdFingerprintGenerator # instead of AllChem
        from sklearn.model_selection import RepeatedStratifiedKFold
        from sklearn.model_selection import RandomizedSearchCV
        from sklearn.utils import class_weight
        from scipy.stats import uniform
        from IPython.display import HTML,Image,display
        ## REMEMBER that 'sklearn' is used to call the package otherwise installed as
        ## 'scikit-learn': the old PyPi 'sklearn' is deprecated 

        ## DISTRIBUTED SERVICE CALLS
        try:
            ## NETWORK CONSTRUCTION
            monarch_result = call_service('monarch', '/run', {
                **payload, 
                'today': today.strftime('%Y-%m-%d'),
                'actual_today': actual_today.strftime('%Y-%m-%d'),
                'date_str': date_str,
                'curr_year': curr_year,
                'overall_start_time': overall_start_time
            })

            # Extract variables from monarch_result
            nodes = monarch_result.get('nodes', [])
            edges = monarch_result.get('edges', [])
            disease_name_id = monarch_result.get('disease_name_id', None)
            disease_name_label = monarch_result.get('disease_name_label', None)
            disease_directories = monarch_result.get('disease_directories', {})

            dgidb_result = call_service('dgidb', '/run', {
                **payload, 
                **monarch_result,  # Pass all previous results
                'today': today.strftime('%Y-%m-%d'),
                'run_layers': payload['run_layers']
            })

            # Extract variables from dgidb_result
            nodes = dgidb_result.get('nodes', nodes)
            edges = dgidb_result.get('edges', edges)
            drug_nodes = dgidb_result.get('drug_nodes', [])

            drugsim_result = call_service('drugsim', '/run', {
                **payload, 
                **dgidb_result,  # Pass all previous results
                'today': today.strftime('%Y-%m-%d'),
                'input_min_simil': payload['input_min_simil']
            })

            # Extract variables from drugsim_result
            edges = drugsim_result.get('edges', edges)
            drug_edges = drugsim_result.get('drug_edges', [])
            
            ## DRUG PREDICTION
            if payload['negs_toggle'] == 1:
                negsample_result = call_service('negsample', '/generate', {
                    **payload, 
                    **drugsim_result,
                    'edges': edges,
                    'similarity_threshold': payload['sim_threshold']
                })

                # Extract valid negative edges
                valid_negative_edges = negsample_result.get('valid_negative_edges', [])
                edges = unique_elements(edges + valid_negative_edges)

                network_result = call_service('networkmodel', '/run-with-ns', {
                    **payload,
                    **negsample_result,
                    'today': today.strftime('%Y-%m-%d'),
                    'run_jobs': payload['num_jobs'],
                    'run_depth': payload['depth_input'],
                    'run_seed': payload['seed_input']
                })
            else:
                network_result = call_service('networkmodel', '/run', {
                    **payload,
                    **drugsim_result,
                    'today': today.strftime('%Y-%m-%d'),
                    'run_jobs': payload['num_jobs'],
                    'run_depth': payload['depth_input'],
                    'run_seed': payload['seed_input']
                })

            # RESULTS EXTRACTIONS
            network_edges = network_result.get('network_edges', [])
            network_nodes = network_result.get('network_nodes', [])
            ranked_drugs = network_result.get('ranked_drugs', [])
            plots = network_result.get('plots', [])

        except Exception as e:
            logging.error(f"Pipeline execution failed: {e}")
            ## FUTURE: make error.html or a custom response
            return render_template('error.html', error=str(e))

        # ## NETWORK CONSTRUCTION
        # nodes, edges, disease_name_id, disease_name_label, disease_directories = run_monarch(input_seed)
        # nodes, edges, drug_nodes = run_dgidb(input_seed,today, layers = run_layers)
        # edges, drug_edges = run_drugsimilarity(input_seed,today,min_simil=input_min_simil)
        
        # ## DRUG PREDICTION
        # if negs_toggle == 1:
        #     valid_negative_edges = generate_negative_samples(edges,similarity_threshold=sim_threshold)
        #     edges = edges + valid_negative_edges
        #     edges = unique_elements(edges)
        #     network_edges, network_nodes, ranked_drugs, plots = run_network_model_with_NS(input_seed,today,run_jobs=num_jobs,run_depth=depth_input, run_seed=seed_input)
        # elif negs_toggle == 0:
        #     network_edges, network_nodes, ranked_drugs, plots = run_network_model(input_seed,today,run_jobs=num_jobs,run_depth=depth_input, run_seed=seed_input)

        ## LOGGING: report run duration
        overall_end_time = time.time()
        overall_duration = overall_end_time - overall_start_time  # calculate duration in seconds
        minutes = int(overall_duration // 60)  # convert seconds to whole minutes
        seconds = int(overall_duration % 60)  # get the remaining seconds
        logging.info(f"PIPELINE run finished in {minutes} minutes and {seconds} seconds.")

    return render_template('home.html', form=form)

@app.route("/about")
def about():
    return render_template('about.html', title='About')

