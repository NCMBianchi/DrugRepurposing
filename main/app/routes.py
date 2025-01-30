"""
MAIN MODULE: FLASKAPP FUNCTIONALITIES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

# for setup, logging and date-override
import sys,os,platform,datetime,logging,builtins,time,multiprocessing
from logger_utils import *

# for Flask App operations
from flask import Flask, render_template, url_for, request, session, redirect
from __init__ import app
from forms import *
from filepaths import *
from unique import *
from queue_manager import *

# for Docker micro services implementation
from servicerun import ServiceOrchestrator

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

        ## NETWORK CONSTRUCTION
        try:
            # Monarch Service
            monarch_launch_status = service_queue_manager.request_service_instance(
                'monarch',
                {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory']
                }
            )
            monarch_instance_id = monarch_launch_status['instance_id']

            if not monarch_launch_status['can_launch']:
                flash(f"Monarch service is currently busy. You are in queue position {monarch_launch_status.get('queue_position', 'unknown')}. "
                      f"Active instances: {monarch_launch_status['active_instances']}/{monarch_launch_status['max_instances']}")
                return render_template('queue_wait.html', service='Monarch', status=monarch_launch_status)

            monarch_orchestrator = ServiceOrchestrator()
            try:
                result = monarch_orchestrator.process_service_request('monarch', {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory']
                })

                # Unpack the result
                nodes = result['nodes']
                edges = result['edges']
                disease_name_id = result.get('disease_name_id')
                disease_name_label = result.get('disease_name_label')
                disease_directories = result.get('disease_directories')
                service_queue_manager.release_service_instance('monarch', monarch_instance_id)

            except Exception as e:
                service_queue_manager.release_service_instance('monarch', monarch_instance_id)
                logging.error(f"Monarch service failed: {e}")
                flash("Monarch service encountered an error. Please try again.")
                return render_template('home.html', form=form)

            # DGIdb Service
            dgidb_launch_status = service_queue_manager.request_service_instance(
                'dgidb',
                {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'nodes': nodes,
                    'run_layers': run_layers
                }
            )
            dgidb_instance_id = dgidb_launch_status['instance_id']

            if not dgidb_launch_status['can_launch']:
                flash(f"DGIdb service is currently busy. You are in queue position {dgidb_launch_status.get('queue_position', 'unknown')}. "
                      f"Active instances: {dgidb_launch_status['active_instances']}/{dgidb_launch_status['max_instances']}")
                return render_template('queue_wait.html', service='DGIdb', status=dgidb_launch_status)

            dgidb_orchestrator = ServiceOrchestrator()
            try:
                result = dgidb_orchestrator.process_service_request('dgidb', {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'nodes': nodes,
                    'run_layers': run_layers
                })

                # Unpack the result
                nodes = result['nodes']
                edges = result['edges']
                drug_nodes = result['drug_nodes']
                service_queue_manager.release_service_instance('dgidb', dgidb_instance_id)

            except Exception as e:
                service_queue_manager.release_service_instance('dgidb', dgidb_instance_id)

                logging.error(f"DGIdb service failed: {e}")
                flash("DGIdb service encountered an error. Please try again.")
                return render_template('home.html', form=form)

            # Drug Similarity Service
            drugsim_launch_status = service_queue_manager.request_service_instance(
                'drugsimilarity',
                {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'nodes': nodes,
                    'input_min_simil': input_min_simil
                }
            )
            drugsim_instance_id = drugsim_launch_status['instance_id']

            if not drugsim_launch_status['can_launch']:
                flash(f"Drug Similarity service is currently busy. You are in queue position {drugsim_launch_status.get('queue_position', 'unknown')}. "
                      f"Active instances: {drugsim_launch_status['active_instances']}/{drugsim_launch_status['max_instances']}")
                return render_template('queue_wait.html', service='Drug Similarity', status=drugsim_launch_status)

            drugsim_orchestrator = ServiceOrchestrator()
            try:
                result = drugsim_orchestrator.process_service_request('drugsimilarity', {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'nodes': nodes,
                    'input_min_simil': input_min_simil
                })

                # Unpack the result
                edges = result['edges']
                drug_edges = result['drug_edges']
                service_queue_manager.release_service_instance('drugsimilarity', drugsim_instance_id)

            except Exception as e:
                service_queue_manager.release_service_instance('drugsimilarity', drugsim_instance_id)

                logging.error(f"Drug Similarity service failed: {e}")
                flash("Drug Similarity service encountered an error. Please try again.")
                return render_template('home.html', form=form)

            ## DRUG PREDICTION
            network_model_launch_status = service_queue_manager.request_service_instance(
                'networkmodel',
                {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'nodes': nodes,
                    'edges': edges,
                    'drug_nodes': drug_nodes,
                    'drug_edges': drug_edges,
                    'negs_toggle': negs_toggle
                }
            )
            network_model_instance_id = network_model_launch_status['instance_id']

            if not network_model_launch_status['can_launch']:
                flash(f"Network Model service is currently busy. You are in queue position {network_model_launch_status.get('queue_position', 'unknown')}. "
                      f"Active instances: {network_model_launch_status['active_instances']}/{network_model_launch_status['max_instances']}")
                return render_template('queue_wait.html', service='Network Model', status=network_model_launch_status)

            network_model_orchestrator = ServiceOrchestrator()
            try:
                if negs_toggle == 1:
                    # Negative Samples Service
                    negsamples_launch_status = service_queue_manager.request_service_instance(
                        'negsample',
                        {
                            'input_seed': input_seed,
                            'date': date_str,
                            'base_directory': base_directories['today_directory'],
                            'edges': edges,
                            'sim_threshold': sim_threshold
                        }
                    )
                    negsamples_instance_id = negsamples_launch_status['instance_id']

                    if not negsamples_launch_status['can_launch']:
                        flash(f"Negative Samples service is currently busy. You are in queue position {negsamples_launch_status.get('queue_position', 'unknown')}. "
                              f"Active instances: {negsamples_launch_status['active_instances']}/{negsamples_launch_status['max_instances']}")
                        return render_template('queue_wait.html', service='Negative Samples', status=negsamples_launch_status)

                    negsamples_orchestrator = ServiceOrchestrator()
                    try:
                        result = negsamples_orchestrator.process_service_request('negsamples', {
                            'input_seed': input_seed,
                            'date': date_str,
                            'base_directory': base_directories['today_directory'],
                            'edges': edges,
                            'sim_threshold': sim_threshold
                        })

                        # Unpack the result
                        valid_negative_edges = result['valid_negative_edges']
                        service_queue_manager.release_service_instance('negsample', negsamples_instance_id)

                    except Exception as e:
                        # Release all previously acquired services
                        service_queue_manager.release_service_instance('negsample', negsamples_instance_id)

                        logging.error(f"Negative Samples service failed: {e}")
                        flash("Negative Samples service encountered an error. Please try again.")
                        return render_template('home.html', form=form)

                # Run Network Model
                network_edges, network_nodes, ranked_drugs, plots = network_model_orchestrator.run(
                    nodes=nodes,
                    edges=edges,
                    drug_nodes=drug_nodes,
                    drug_edges=drug_edges,
                    negs_toggle=negs_toggle,
                    run_depth=depth_input,
                    num_jobs=num_jobs,
                    seed_input=seed_input,
                    input_min_simil=input_min_simil,
                    sim_threshold=sim_threshold
                )

                # Successfully completed - release the 'networkmodel' service
                service_queue_manager.release_service_instance('networkmodel', network_model_instance_id)

            except Exception as e:
                # Release the 'networkmodel' service
                service_queue_manager.release_service_instance('networkmodel', network_model_instance_id)

                service_queue_manager.release_service_instance('monarch')
                service_queue_manager.release_service_instance('dgidb')
                service_queue_manager.release_service_instance('drugsimilarity')

                logging.error(f"Network Model service failed: {e}")
                flash("Network Model service encountered an error. Please try again.")
                return render_template('home.html', form=form)

        except Exception as e:
            logging.error(f"Overall pipeline execution failed: {e}")
            flash("An unexpected error occurred during pipeline execution.")
            return render_template('home.html', form=form)

        ## LOGGING: report run duration
        overall_end_time = time.time()
        overall_duration = overall_end_time - overall_start_time
        minutes = int(overall_duration // 60)
        seconds = int(overall_duration % 60)
        logging.info(f"PIPELINE run finished in {minutes} minutes and {seconds} seconds.")

    return render_template('home.html', form=form)

@app.route("/about")
def about():
    return render_template('about.html', title='About')
