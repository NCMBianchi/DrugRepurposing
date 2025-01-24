"""
MAIN MODULE: FLASKAPP FUNCTIONALITIES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

# for setup, logging and date-override
import sys,os,platform,datetime,logging,builtins,time,multiprocessing
from logger_utils import *

# for service activation and communication in Docker Swarm
import requests
import json
from urllib.parse import urljoin
import docker
from docker import DockerClient

# for Flask App operations
from flask import Flask, render_template, url_for, request, session, redirect
from __init__ import app
from forms import *
from filepaths import *
from unique import *
from queue_manager import *

## PLATFORM INFO
python_executable = sys.executable
python_version = platform.python_version()
num_cores = multiprocessing.cpu_count()

@app.route("/", methods = ['GET', 'POST'])
@app.route("/home", methods = ['GET', 'POST'])
def config():
    
    form = user_input()

    if form.validate_on_submit():

        d_toggle = form.date_t.data  # Date override toggle
        date_OR_day = form.date_OR_day.data
        date_OR_month = form.date_OR_month.data
        date_OR_year = form.date_OR_year.data

        # generate the override date if toggle is active
        if d_toggle:
            d_override = datetime.date(date_OR_year, date_OR_month, date_OR_day)
        else:
            d_override = None
        
        input_seed = form.disease_URI.data  # disease URI
        run_layers = form.deg_of_dist.data  # degree of distance
        input_min_simil = float(form.inp_minimum_sim.data)  # minimum drug similarity
        negs_toggle = form.ns_toggle.data  # negative samples toggle
        sim_threshold = float(form.sim_t.data)  # drug similarity threshold
        num_jobs = form.n_cores.data  # number of CPU cores
        depth_input = form.po_mode.data  # mode of operation
        seed_input = form.ML_seed.data  # ML seed

        ## DATES
        today, actual_today, date_str, curr_year, overall_start_time, formatted_start_time = set_date(d_toggle, d_override)
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
            ## (...) 
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

        ## placeholder

    return render_template('home.html', form=form)

@app.route("/about")
def about():
    return render_template('about.html', title='About')