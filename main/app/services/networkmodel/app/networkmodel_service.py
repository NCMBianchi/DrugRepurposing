"""
NETWORK MODEL MODULE FOR DISTRIBUTED SERVICES
Created on January 30th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import sys
import json
import logging
from typing import Dict, Any, List, Tuple
import datetime

# Import the core Network Model functionality
from .networkmodel_legacy import run_networkmodel as legacy_run_networkmodel
from .filepaths import setup_service_globals, initialise_base_directories

def run_networkmodel(
    monarch_input: str = 'MONDO:0007739',
    base_directory: str = None,
    **kwargs
) -> Tuple[List, List, Dict]:
    """
    Wrapper for legacy_run_networkmodel that adapts to microservice environment

    :param monarch_input: Monarch input seed (disease identifier)
    :param base_directory: Base directory for data storage
    :param kwargs: Additional arguments to pass to legacy function
    :return: Tuple of network nodes, network edges, and network metrics
    """
    # Use environment variables if base_directory not provided
    if base_directory is None:
        base_directory = os.environ.get('BASE_DIRECTORY', '/app/data')

    # Extract date from environment or use current date
    date = os.environ.get('DATE', datetime.now().strftime('%Y-%m-%d'))

    # Setup global variables for the service
    today_directory, date_str, input_seed, input_file_path = setup_service_globals(
        base_directory, date, monarch_input, 'networkmodel'
    )

    # Run the legacy function
    return legacy_run_networkmodel(monarch_input, date, **kwargs)

def process_service_request(env_vars: Dict[str, str]) -> Dict[str, Any]:
    """
    Process Network Model service request based on environment variables

    :param env_vars: Dictionary of environment variables
    :return: Dictionary with service results
    """
    # Extract required parameters
    monarch_input = env_vars.get('MONARCH_INPUT', 'MONDO:0007739')
    date = env_vars.get('DATE')
    base_directory = env_vars.get('BASE_DIRECTORY', '/app/data')
    network_type = env_vars.get('NETWORK_TYPE', 'standard')

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        # Run Network Model generation
        network_nodes, network_edges, network_metrics = run_networkmodel(
            monarch_input=monarch_input,
            base_directory=base_directory,
            network_type=network_type
        )

        return {
            'status': 'success',
            'network_nodes': network_nodes,
            'network_edges': network_edges,
            'network_metrics': network_metrics
        }

    except Exception as e:
        logging.error(f"Network Model service processing failed: {e}")
        return {
            'status': 'error',
            'message': str(e)
        }

def main():
    # Read environment variables
    env_vars = {
        key: os.environ.get(key, '')
        for key in [
            'MONARCH_INPUT', 'DATE', 'BASE_DIRECTORY', 'NETWORK_TYPE'
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
