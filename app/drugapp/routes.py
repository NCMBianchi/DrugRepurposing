"""
MAIN MODULE: FLASKAPP FUNCTIONALITIES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

# for setup, logging and date-override
import sys,os,platform,datetime,logging,builtins,time,multiprocessing
from drugapp.platform import *

# for Flask App operations
from flask import Flask, render_template, url_for, request, session, redirect
from drugapp import app 
from drugapp.forms import *
from drugapp.monarch import *
from drugapp.dgidb import *
from drugapp.drugsimilarity import *
from drugapp.negsamples import *
from drugapp.networkmodel import *

## PLATFORM INFO
python_executable = sys.executable
python_version = platform.python_version()
num_cores = multiprocessing.cpu_count()

@app.route("/", methods = ['GET', 'POST'])
@app.route("/home", methods = ['GET', 'POST'])
def config():

    form = user_input()
    if form.validate_on_submit():
        d_toggle = form.user_input.data #date_t
        d_override = form.user_input.data #date_OR
        input_seed = form.user_input.data #disease_URI
        run_layers = form.user_input.data #deg_of_dist
        input_min_simil = form.user_input.data #inp_minimum_sim
        negs_toggle = form.user_input.data #ns_toggle
        sim_threshold = form.user_input.data #sim_t
        num_jobs = form.user_input.data #n_cores
        depth_input = form.user_input.data #po_mode
        seed_input = form.user_input.data #ML_seed

        ## DATES
        today, actual_today, date_str, curr_year, overall_start_time, formatted_start_time = set_date(d_toggle, d_override)
        logging.info(f"'DrugRepurposing' pipeline started at {formatted_start_time}.\n")

        ## DIRECTORY NAMES
        base_data_directory = os.path.join(os.getcwd(), 'drugapp', 'data')
        today_directory = os.path.join(base_data_directory, date_str)
        os.makedirs(today_directory, exist_ok=True)
        disease_name_label = None #placeholder

        ## LOGGING: set-up
        platform_filename = "running_platform.txt"
        platform_file_path = os.path.join(today_directory, platform_filename)
        log_filename = "logfile.log"
        log_file_path = os.path.join(today_directory, log_filename)
        input_filename = "inputs.txt"
        input_file_path = os.path.join(today_directory, input_filename)
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

        ## PACKAGES: python-included packages
        import tqdm
        import re
        import json  # to save .json files returned by the API
        import logging  # to check for running variables
        import inspect  # to chec for running functions
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

        ## NETWORK CONSTRUCTION
        nodes, edges, disease_name_id, disease_name_label = run_monarch(input_seed)
        nodes, edges, drug_nodes = run_dgidb(input_seed,today, layers = run_layers)
        edges, drug_edges = run_drugsimilarity(input_seed,today,min_simil=input_min_simil)
        
        ## DRUG PREDICTION
        if negs_toggle == 1:
            valid_negative_edges = generate_negative_samples(edges,similarity_threshold=sim_threshold)
            edges = edges + valid_negative_edges
            edges = unique_elements(edges)
            network_edges, network_nodes, ranked_drugs, plots = run_network_model_with_NS(input_seed,today,run_jobs=num_jobs,run_depth=depth_input, run_seed=seed_input)
        elif negs_toggle == 0:
            network_edges, network_nodes, ranked_drugs, plots = run_network_model(input_seed,today,run_jobs=num_jobs,run_depth=depth_input, run_seed=seed_input)

        ## LOGGING: report run duration
        overall_end_time = time.time()
        overall_duration = overall_end_time - overall_start_time  # calculate duration in seconds
        minutes = int(overall_duration // 60)  # convert seconds to whole minutes
        seconds = int(overall_duration % 60)  # get the remaining seconds
        logging.info(f"PIPELINE run finished in {minutes} minutes and {seconds} seconds.")

    return render_template('home.html', folders = folderlist, form=form, folders_output=folders_output)

@app.route("/about")
def about():
    return render_template('about.html', title='About', folders_output=folders_output)

