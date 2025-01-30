"""
MONARCH MODULE FOR DISTRIBUTED SERVICES: MAKES API CALL AND ITERATES THROUGH DEGREES OF DISTANCE
Created on JANUARY 27TH 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import sys
import json
import logging
from typing import Dict, Any, List

# Import the core Monarch functionality
from .Monarch_legacy import run_monarch as legacy_run_monarch
from .filepaths import setup_service_globals, initialise_base_directories

def run_monarch(
    input_id: str = 'MONDO:0007739',
    base_directory: str = None,
    **kwargs
) -> Tuple[List, List, str, str, Dict]:
    """
    Wrapper for legacy_run_monarch that adapts to microservice environment

    :param input_id: MONDO disease identifier
    :param base_directory: Base directory for data storage
    :param kwargs: Additional arguments to pass to legacy function
    :return: Same return structure as legacy run_monarch
    """
    # Use environment variables if base_directory not provided
    if base_directory is None:
        base_directory = os.environ.get('BASE_DIRECTORY', '/app/data')

    # Extract date from environment or use current date
    date = os.environ.get('DATE', datetime.now().strftime('%Y-%m-%d'))

    # Setup global variables for the service
    today_directory, date_str, input_seed, input_file_path = setup_service_globals(
        base_directory, date, input_id, 'monarch'
    )

    # Run the legacy function
    return legacy_run_monarch(input_id, **kwargs)

def process_service_request(env_vars: Dict[str, str]) -> Dict[str, Any]:
    """
    Process Monarch service request based on environment variables

    :param env_vars: Dictionary of environment variables
    :return: Dictionary with service results
    """
    # Extract required parameters
    input_seed = env_vars.get('INPUT_SEED', 'MONDO:0007739')
    date = env_vars.get('DATE')
    base_directory = env_vars.get('BASE_DIRECTORY', '/app/data')
    layers = int(env_vars.get('LAYERS', 3))

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        # Run Monarch discovery
        nodes, edges, disease_id, disease_name, disease_dir = run_monarch(
            input_id=input_seed,
            base_directory=base_directory
        )

        return {
            'status': 'success',
            'nodes': nodes,
            'edges': edges,
            'disease_id': disease_id,
            'disease_name': disease_name,
            'disease_dir': str(disease_dir)  # Convert potential complex objects to string
        }

    except Exception as e:
        logging.error(f"Monarch service processing failed: {e}")
        return {
            'status': 'error',
            'message': str(e)
        }

def main():
    # Read environment variables
    env_vars = {
        key: os.environ.get(key, '')
        for key in [
            'INPUT_SEED', 'DATE', 'BASE_DIRECTORY', 'LAYERS'
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
