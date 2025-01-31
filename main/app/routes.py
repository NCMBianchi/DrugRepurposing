"""
MAIN MODULE: FLASKAPP FUNCTIONALITIES
Created on August 3rd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

# for setup, logging and date-override
import sys,os,platform,datetime,logging,builtins,time,multiprocessing
from logger_utils import *

# for run tracking
import uuid
import json
import redis

# for Flask App operations
from flask import Flask, render_template, url_for, request, session, redirect, jsonify
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

def save_run_metadata(run_id, metadata):
    """Save run metadata to Redis"""
    if redis_client:
        redis_key = f"run:{run_id}"
        try:
            redis_client.hmset(redis_key, metadata)
            ## redis_client.expire(redis_key, 86400)
            # currently 24h, can be modified if an expiration threshold is even necessary
        except Exception as e:
            logging.error(f"Failed to save run metadata: {e}")

def update_run_stage(run_id, stage, status):
    """Update specific stage status in Redis"""
    if redis_client:
        redis_key = f"run:{run_id}"
        try:
            stage_metadata = {
                "status": status,
                "timestamp": datetime.now().isoformat(),
                "start_time": datetime.now().isoformat() if status == "running" else None,
                "end_time": datetime.now().isoformat() if status == "completed" else None
            }
            redis_client.hset(redis_key, f"stage:{stage}", json.dumps(stage_metadata))
        except Exception as e:
            logging.error(f"Failed to update run stage: {e}")

def stop_run(run_id):
    """Forcibly stop a running pipeline"""
    try:
        # Update run status in Redis
        if redis_client:
            redis_client.hset(f"run:{run_id}", "status", "stopped")

        # Terminate associated services (you'll need to implement service termination logic)
        # This might involve calling service_queue_manager to release instances
        service_queue_manager.release_service_instance('monarch')
        service_queue_manager.release_service_instance('dgidb')
        # Add other services as needed

        logging.info(f"Run {run_id} forcibly stopped")
        return True
    except Exception as e:
        logging.error(f"Error stopping run {run_id}: {e}")
        return False

def delete_run(run_id):
    """Delete run metadata and associated data"""
    try:
        # Get run details from Redis to locate data directory
        if redis_client:
            run_data = redis_client.hgetall(f"run:{run_id}")

            # Extract necessary information (e.g., date, disease name)
            date_str = run_data.get('date', datetime.now().strftime('%Y-%m-%d'))
            disease_name = run_data.get('disease_name', 'unknown_disease')

            # Construct path to run data
            base_data_directory = os.path.join('/app', 'data')
            run_directory = os.path.join(base_data_directory, date_str, f"{disease_name} ({date_str})")

            # Recursively delete directory
            if os.path.exists(run_directory):
                shutil.rmtree(run_directory)

            # Remove Redis entry
            redis_client.delete(f"run:{run_id}")

            logging.info(f"Run {run_id} deleted. Data removed from {run_directory}")
            return True
    except Exception as e:
        logging.error(f"Error deleting run {run_id}: {e}")
        return False

@app.route("/", methods = ['GET', 'POST'])
@app.route("/home", methods = ['GET', 'POST'])
def config():

    form = user_input()

    if form.validate_on_submit():
        # Generate unique run ID
        run_id = str(uuid.uuid4())

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
        save_run_metadata(run_id, {
            "input_seed": input_seed,
            "start_time": datetime.now().isoformat(),
            "disease_name": input_seed,
            "stages": json.dumps(["monarch", "dgidb", "drugsimilarity", "negsamples", "networkmodel"])
        })

        try:
            # Monarch Service
            update_run_stage(run_id, "monarch", "pending")
            monarch_launch_status = service_queue_manager.request_service_instance(
                'monarch',
                {
                    'input_seed': input_seed,
                    'date': date_str,
                    'base_directory': base_directories['today_directory'],
                    'run_id': run_id
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

            update_run_stage(run_id, "monarch", "completed")

            # DGIdb Service
            update_run_stage(run_id, "dgidb", "pending")
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

            update_run_stage(run_id, "dgidb", "completed")

            # Drug Similarity Service
            update_run_stage(run_id, "drugsimilarity", "pending")
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

            update_run_stage(run_id, "drugsimilarity", "completed")

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
                    update_run_stage(run_id, "negsample", "pending")
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

                    update_run_stage(run_id, "negsample", "completed")

                if negs_toggle == 0:
                    update_run_stage(run_id, "negsample", "skipped")

                # Run Network Model
                update_run_stage(run_id, "networkmodel", "pending")
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

            update_run_stage(run_id, "networkmodel", "completed")

        except Exception as e:
            if redis_client:
                redis_client.hset(f"run:{run_id}", "status", "failed")
                redis_client.hset(f"run:{run_id}", "error", str(e))
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

@app.route("/stop_run/<run_id>", methods=['POST'])
def stop_run_route(run_id):
    """Stop a running run and show its status"""
    try:
        # Update run status in Redis
        if redis_client:
            # Retrieve existing run data
            run_key = f"run:{run_id}"
            run_data = redis_client.hgetall(run_key)

            # Mark run as stopped
            redis_client.hset(run_key, "status", "stopped")

            # Terminate associated services
            service_queue_manager.release_service_instance('monarch')
            service_queue_manager.release_service_instance('dgidb')
            # Add other services as needed

            # Prepare context for run_status template
            stages = []
            stage_keys = [key.decode() for key in redis_client.keys(f"{run_key}:stage:*")]

            for stage_key in stage_keys:
                stage_name = stage_key.split(':')[-1]
                stage_info = json.loads(redis_client.hget(run_key, stage_key))
                stages.append({
                    'name': stage_name,
                    'status': stage_info['status'],
                    'timestamp': stage_info['timestamp']
                })

            context = {
                'run_id': run_id,
                'disease_id': run_data.get('input_seed', 'Unknown'),
                'start_time': run_data.get('start_time', 'N/A'),
                'stages': stages,
                'overall_status': 'stopped'
            }

            logging.info(f"Run {run_id} forcibly stopped")
            return render_template('run_status.html', **context)

        flash("Unable to stop run: Redis unavailable")
        return redirect(url_for('current_runs'))

    except Exception as e:
        logging.error(f"Error stopping run {run_id}: {e}")
        flash(f"Error stopping run: {e}")
        return redirect(url_for('current_runs'))

@app.route("/delete_run/<run_id>", methods=['POST'])
def delete_run_route(run_id):
    """Delete run metadata and associated data"""
    try:
        if redis_client:
            # Get run details from Redis to locate data directory
            run_data = redis_client.hgetall(f"run:{run_id}")

            # Extract necessary information
            date_str = run_data.get('date', datetime.now().strftime('%Y-%m-%d'))
            disease_name = run_data.get('input_seed', 'unknown_disease')

            # Construct path to run data
            base_data_directory = os.path.join('/app', 'data')
            run_directory = os.path.join(base_data_directory, date_str, f"{disease_name} ({date_str})")

            # Recursively delete directory
            if os.path.exists(run_directory):
                shutil.rmtree(run_directory)

            # Remove Redis entry
            redis_client.delete(f"run:{run_id}")

            logging.info(f"Run {run_id} deleted. Data removed from {run_directory}")
            return render_template('run_deleted.html', run_id=run_id)

    except Exception as e:
        logging.error(f"Error deleting run {run_id}: {e}")
        flash(f"Error deleting run: {e}")
        return redirect(url_for('current_runs'))

@app.route("/current_runs")
def current_runs():
    """List all current/recent runs"""
    try:
        # Try to connect to Redis if not already connected
        if 'redis_client' not in globals() or redis_client is None:
            try:
                redis_client = redis.Redis(
                    host=os.environ.get('REDIS_HOST', 'localhost'),
                    port=int(os.environ.get('REDIS_PORT', 6379)),
                    db=int(os.environ.get('REDIS_DB', 0)),
                    decode_responses=True
                )
                redis_client.ping()  # Test connection
            except Exception as e:
                # If Redis connection fails, log the error and proceed with an empty list
                logging.warning(f"Redis connection failed: {e}")
                return render_template('current_runs.html', runs=[])

        # Find all run keys in Redis
        run_keys = redis_client.keys("run:*")

        runs = []
        for key in run_keys:
            run_id = key.decode().split(':')[1]
            run_data = redis_client.hgetall(key)

            # Basic processing to prepare run data for template
            run_info = {
                'id': run_id,
                'input_seed': run_data.get('input_seed', 'Unknown'),
                'start_time': run_data.get('start_time', 'N/A'),
                'status': run_data.get('status', 'running')
            }
            runs.append(run_info)

        return render_template('current_runs.html', runs=runs)

    except Exception as e:
        # Fallback to empty runs list if any unexpected error occurs
        logging.error(f"Error in current_runs: {e}")
        return render_template('current_runs.html', runs=[])

@app.route("/run_status/<run_id>")
def run_status(run_id):
    """Detailed status for a specific run"""
    if not redis_client:
        flash("Redis service unavailable")
        return render_template('home.html')

    # Retrieve run data from Redis
    run_key = f"run:{run_id}"
    run_data = redis_client.hgetall(run_key)

    if not run_data:
        flash("Run not found")
        return redirect(url_for('current_runs'))

    # Process stage information
    stages = []
    stage_keys = [key.decode() for key in redis_client.keys(f"{run_key}:stage:*")]

    for stage_key in stage_keys:
        stage_name = stage_key.split(':')[-1]
        stage_info = json.loads(redis_client.hget(run_key, stage_key))
        runtime = None
        if stage_info['start_time'] and stage_info['end_time']:
            start = datetime.fromisoformat(stage_info['start_time'])
            end = datetime.fromisoformat(stage_info['end_time'])
            runtime = f"{(end - start).total_seconds():.2f}s"

        stages.append({
            'name': stage_name,
            'status': stage_info['status'],
            'timestamp': stage_info['timestamp'],
            'runtime': runtime
        })

    # Compute partial runtime
    total_run_time = (datetime.now() - overall_start_time).total_seconds()
    run_time_display = f"{total_run_time:.2f}s" if total_run_time < 600 else f"{int(total_run_time // 60)} min"

    # Template context
    context = {
        'run_id': run_id,
        'disease_id': run_data.get('input_seed', 'Unknown'),
        'start_time': run_data.get('start_time', 'N/A'),
        'run_time': run_time_display,
        'stages': stages,
        'overall_status': run_data.get('status', 'running')
    }

    return render_template('run_status.html', **context)

@app.route("/run_status_update/<run_id>")
def run_status_update(run_id):
    """
    Provide dynamic updates for a specific run status
    Returns JSON with current stage and overall status
    """
    if not redis_client:
        return jsonify({
            'error': 'Redis service unavailable',
            'overall_status': 'error'
        }), 500

    # Retrieve run data from Redis
    run_key = f"run:{run_id}"
    run_data = redis_client.hgetall(run_key)

    if not run_data:
        return jsonify({
            'error': 'Run not found',
            'overall_status': 'error'
        }), 404

    # Process stage information
    stages = []
    stage_keys = [key.decode() for key in redis_client.keys(f"{run_key}:stage:*")]

    for stage_key in stage_keys:
        stage_name = stage_key.split(':')[-1]
        try:
            stage_info = json.loads(redis_client.hget(run_key, stage_key))

            # Compute runtime if stage is completed
            runtime = None
            if stage_info['start_time'] and stage_info['end_time']:
                start = datetime.fromisoformat(stage_info['start_time'])
                end = datetime.fromisoformat(stage_info['end_time'])
                runtime = f"{(end - start).total_seconds():.2f}s"

            stages.append({
                'name': stage_name,
                'status': stage_info['status'],
                'timestamp': stage_info['timestamp'],
                'runtime': runtime
            })
        except (TypeError, json.JSONDecodeError) as e:
            logging.error(f"Error processing stage {stage_name}: {e}")

    # Prepare response
    response_data = {
        'run_id': run_id,
        'overall_status': run_data.get('status', 'running').decode('utf-8'),
        'stages': stages
    }

    return jsonify(response_data)

@app.route("/about")
def about():
    return render_template('about.html', title='About')
