"""
DGIDB MODULE FOR DISTRIBUTED SERVICES: MAKES API CALL AND PROCESSES DRUG INTERACTIONS
Created on JANUARY 30TH 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import sys
import json
import logging
from typing import Dict, Any, List, Tuple

# Import the core DGIdb functionality
from .dgidb_legacy import run_dgidb as legacy_run_dgidb
from .filepaths import setup_service_globals, initialise_base_directories

def run_dgidb(
    monarch_input: str = 'MONDO:0007739',
    base_directory: str = None,
    **kwargs
) -> Tuple[List, List, List]:
    """
    Wrapper for legacy_run_dgidb that adapts to microservice environment

    :param monarch_input: Monarch input seed (disease identifier)
    :param base_directory: Base directory for data storage
    :param kwargs: Additional arguments to pass to legacy function
    :return: Tuple of unique nodes, unique edges, and all drugs
    """
    # Use environment variables if base_directory not provided
    if base_directory is None:
        base_directory = os.environ.get('BASE_DIRECTORY', '/app/data')

    # Extract date from environment or use current date
    date = os.environ.get('DATE', datetime.now().strftime('%Y-%m-%d'))

    # Setup global variables for the service
    today_directory, date_str, input_seed, input_file_path = setup_service_globals(
        base_directory, date, monarch_input, 'dgidb'
    )

    # Run the legacy function
    return legacy_run_dgidb(monarch_input, date, **kwargs)

def process_service_request(env_vars: Dict[str, str]) -> Dict[str, Any]:
    """
    Process DGIdb service request based on environment variables

    :param env_vars: Dictionary of environment variables
    :return: Dictionary with service results
    """
    # Extract required parameters
    monarch_input = env_vars.get('MONARCH_INPUT', 'MONDO:0007739')
    date = env_vars.get('DATE')
    base_directory = env_vars.get('BASE_DIRECTORY', '/app/data')
    layers = int(env_vars.get('LAYERS', 3))

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        # Run DGIdb discovery
        nodes, edges, all_drugs = run_dgidb(
            monarch_input=monarch_input,
            base_directory=base_directory,
            layers=layers
        )

        return {
            'status': 'success',
            'nodes': nodes,
            'edges': edges,
            'all_drugs': all_drugs
        }

    except Exception as e:
        logging.error(f"DGIdb service processing failed: {e}")
        return {
            'status': 'error',
            'message': str(e)
        }

def main():
    # Read environment variables
    env_vars = {
        key: os.environ.get(key, '')
        for key in [
            'MONARCH_INPUT', 'DATE', 'BASE_DIRECTORY', 'LAYERS'
        ]
    }

    # Process the request
    result = process_service_request(env_vars)

    # Output results (for Docker logs or future inter-service communication)
    print(json.dumps(result))

    # Exit with appropriate status code
    sys.exit(0 if result['status'] == 'success' else 1)

if __name__ == "__main__":
    main()
