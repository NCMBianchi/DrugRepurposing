"""
DRUG SIMILARITY MODULE FOR DISTRIBUTED SERVICES
Created on January 30th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import sys
import json
import logging
from typing import Dict, Any, List, Tuple
import datetime

# Import the core Drug Similarity functionality
from .drugsimilarity_legacy import run_drugsimilarity as legacy_run_drugsimilarity
from .filepaths import setup_service_globals, initialise_base_directories

def run_drugsimilarity(
    monarch_input: str = 'MONDO:0007739',
    base_directory: str = None,
    **kwargs
) -> Tuple[List, List, Dict]:
    """
    Wrapper for legacy_run_drugsimilarity that adapts to microservice environment

    :param monarch_input: Monarch input seed (disease identifier)
    :param base_directory: Base directory for data storage
    :param kwargs: Additional arguments to pass to legacy function
    :return: Tuple of similarity matrix, similarity graph, and drug info
    """
    # Use environment variables if base_directory not provided
    if base_directory is None:
        base_directory = os.environ.get('BASE_DIRECTORY', '/app/data')

    # Extract date from environment or use current date
    date = os.environ.get('DATE', datetime.now().strftime('%Y-%m-%d'))

    # Setup global variables for the service
    today_directory, date_str, input_seed, input_file_path = setup_service_globals(
        base_directory, date, monarch_input, 'drugsim'
    )

    # Run the legacy function
    return legacy_run_drugsimilarity(monarch_input, date, **kwargs)

def process_service_request(env_vars: Dict[str, str]) -> Dict[str, Any]:
    """
    Process Drug Similarity service request based on environment variables

    :param env_vars: Dictionary of environment variables
    :return: Dictionary with service results
    """
    # Extract required parameters
    monarch_input = env_vars.get('MONARCH_INPUT', 'MONDO:0007739')
    date = env_vars.get('DATE')
    base_directory = env_vars.get('BASE_DIRECTORY', '/app/data')

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    try:
        # Run Drug Similarity analysis
        similarity_matrix, similarity_graph, drug_info = run_drugsimilarity(
            monarch_input=monarch_input,
            base_directory=base_directory
        )

        return {
            'status': 'success',
            'similarity_matrix': similarity_matrix,
            'similarity_graph': similarity_graph,
            'drug_info': drug_info
        }

    except Exception as e:
        logging.error(f"Drug Similarity service processing failed: {e}")
        return {
            'status': 'error',
            'message': str(e)
        }

def main():
    # Read environment variables
    env_vars = {
        key: os.environ.get(key, '')
        for key in [
            'MONARCH_INPUT', 'DATE', 'BASE_DIRECTORY'
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
