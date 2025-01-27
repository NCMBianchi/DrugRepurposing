"""
SERVICE RUN MODULE: FUNCTIONS NECESSARY TO LAUNCH DOCKER SERVICES
Created on January 22nd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import docker
import uuid
import logging
import os
import requests
import json
from typing import Dict, Any, Tuple, List

import queue_manager

import Monarch_legacy
import dgidb_legacy
import drugsimilarity_legacy
import negsamples_legacy
import networkmodel_legacy

class BaseServiceOrchestrator:
    def __init__(self, input_seed: str, date: str, base_directory: str):
        self.docker_client = docker.from_env()
        self.input_seed = input_seed
        self.date = date
        self.base_directory = base_directory
        self.service_name = f"{self.__class__.__name__.lower().replace('serviceorchestrator', '')}-{uuid.uuid4().hex[:8]}"
        self.network_name = 'drugrepurposing_default'

    def create_service(self, image: str, env: Dict[str, str] = None) -> Dict[str, Any]:
        """
        Dynamically create a Docker microservice

        :param image: Docker image name
        :param env: Environment variables dictionary
        :return: Dictionary with service details
        """
        try:
            default_env = {
                'INPUT_SEED': self.input_seed,
                'DATE': self.date,
                'BASE_DIRECTORY': self.base_directory
            }

            if env:
                default_env.update(env)

            service = self.docker_client.services.create(
                image=image,
                name=self.service_name,
                env=default_env,
                networks=[self.network_name]
            )

            logging.info(f"{self.__class__.__name__} service {self.service_name} created")

            return {
                'id': service.id,
                'name': self.service_name
            }

        except Exception as e:
            logging.error(f"Failed to create service: {e}")
            raise

    def is_service_mode(self) -> bool:
        """
        Determine if microservice mode is available
        """
        try:
            self.docker_client.ping()
            return True
        except:
            return False

class MonarchServiceOrchestrator(BaseServiceOrchestrator):
    def run(self) -> Tuple[List, List, str, str, Dict]:
        """
        Run Monarch service discovery with queue management
        """
        try:
            # Check if service can be launched immediately
            launch_status = queue_manager.service_queue_manager.request_service_instance(
                'monarch', 
                {
                    'input_seed': self.input_seed,
                    'date': self.date,
                    'base_directory': self.base_directory
                }
            )
            
            if not launch_status['can_launch']:
                raise Exception(f"Monarch service is queued. Position: {launch_status.get('queue_position', 'unknown')}")

            try:
                # Local fallback
                if not self.is_service_mode():
                    results = run_monarch(self.input_seed)
                    queue_manager.service_queue_manager.release_service_instance('monarch')
                    return results

                # Microservice mode
                service = self.create_service(
                    image='drugrepurposing/monarch:latest',
                    env={'LAYERS': '3'}
                )

                # Implement result retrieval mechanism
                results = self._wait_and_fetch_results(service)
                
                # Always release the service instance
                queue_manager.service_queue_manager.release_service_instance('monarch')
                
                return results

            except Exception as e:
                queue_manager.service_queue_manager.release_service_instance('monarch')
                logging.error(f"Monarch service run failed: {e}")
                return run_monarch(self.input_seed)

        except Exception as e:
            logging.error(f"Monarch service orchestration failed: {e}")
            return run_monarch(self.input_seed)

    def _wait_and_fetch_results(self, service):
        """
        Wait for service completion and retrieve results
        Placeholder method - needs actual implementation
        """
        # Potential implementations:
        # 1. Check service logs
        # 2. Use shared volume
        # 3. Use Redis/message queue
        # 4. Implement API endpoint in service to fetch results
        pass

class DGIdbServiceOrchestrator(BaseServiceOrchestrator):
    def run(self, nodes: List, run_layers: int) -> Tuple[List, List]:
        try:
            # Check if service can be launched immediately
            launch_status = queue_manager.service_queue_manager.request_service_instance(
                'dgidb', 
                {
                    'input_seed': self.input_seed,
                    'date': self.date,
                    'base_directory': self.base_directory,
                    'nodes': nodes,
                    'run_layers': run_layers
                }
            )
            
            if not launch_status['can_launch']:
                raise Exception(f"DGIdb service is queued. Position: {launch_status.get('queue_position', 'unknown')}")

            try:
                if not self.is_service_mode():
                    results = run_dgidb(self.input_seed, self.date, layers=run_layers)
                    queue_manager.service_queue_manager.release_service_instance('dgidb')
                    return results

                service = self.create_service(
                    image='drugrepurposing/dgidb:latest',
                    env={
                        'LAYERS': str(run_layers),
                        'NODES': json.dumps(nodes)
                    }
                )

                results = self._wait_and_fetch_results(service)
                queue_manager.service_queue_manager.release_service_instance('dgidb')
                return results

            except Exception as e:
                queue_manager.service_queue_manager.release_service_instance('dgidb')
                logging.error(f"DGIdb service run failed: {e}")
                return run_dgidb(self.input_seed, self.date, layers=run_layers)

        except Exception as e:
            logging.error(f"DGIdb service orchestration failed: {e}")
            return run_dgidb(self.input_seed, self.date, layers=run_layers)

    def _wait_and_fetch_results(self, service):
        # Placeholder implementation
        pass

class DrugSimilarityServiceOrchestrator(BaseServiceOrchestrator):
    def run(self, nodes: List, input_min_simil: float) -> Tuple[List, List]:
        try:
            # Check if service can be launched immediately
            launch_status = queue_manager.service_queue_manager.request_service_instance(
                'drugsimilarity', 
                {
                    'input_seed': self.input_seed,
                    'date': self.date,
                    'base_directory': self.base_directory,
                    'nodes': nodes,
                    'input_min_simil': input_min_simil
                }
            )
            
            if not launch_status['can_launch']:
                raise Exception(f"Drug Similarity service is queued. Position: {launch_status.get('queue_position', 'unknown')}")

            try:
                if not self.is_service_mode():
                    results = run_drugsimilarity(self.input_seed, self.date, min_simil=input_min_simil)
                    queue_manager.service_queue_manager.release_service_instance('drugsimilarity')
                    return results

                service = self.create_service(
                    image='drugrepurposing/drugsimilarity:latest',
                    env={
                        'MIN_SIMILARITY': str(input_min_simil),
                        'NODES': json.dumps(nodes)
                    }
                )

                results = self._wait_and_fetch_results(service)
                queue_manager.service_queue_manager.release_service_instance('drugsimilarity')
                return results

            except Exception as e:
                queue_manager.service_queue_manager.release_service_instance('drugsimilarity')
                logging.error(f"Drug Similarity service run failed: {e}")
                return run_drugsimilarity(self.input_seed, self.date, min_simil=input_min_simil)

        except Exception as e:
            logging.error(f"Drug Similarity service orchestration failed: {e}")
            return run_drugsimilarity(self.input_seed, self.date, min_simil=input_min_simil)

    def _wait_and_fetch_results(self, service):
        # Placeholder implementation
        pass

class NegativeSamplesServiceOrchestrator(BaseServiceOrchestrator):
    def run(self, edges: List, sim_threshold: float) -> List:
        try:
            # Check if service can be launched immediately
            launch_status = queue_manager.service_queue_manager.request_service_instance(
                'negsample', 
                {
                    'input_seed': self.input_seed,
                    'date': self.date,
                    'base_directory': self.base_directory,
                    'edges': edges,
                    'sim_threshold': sim_threshold
                }
            )
            
            if not launch_status['can_launch']:
                raise Exception(f"Negative Samples service is queued. Position: {launch_status.get('queue_position', 'unknown')}")

            try:
                if not self.is_service_mode():
                    results = generate_negative_samples(edges, similarity_threshold=sim_threshold)
                    queue_manager.service_queue_manager.release_service_instance('negsample')
                    return results

                service = self.create_service(
                    image='drugrepurposing/negsamples:latest',
                    env={
                        'SIMILARITY_THRESHOLD': str(sim_threshold),
                        'EDGES': json.dumps(edges)
                    }
                )

                results = self._wait_and_fetch_results(service)
                queue_manager.service_queue_manager.release_service_instance('negsample')
                return results

            except Exception as e:
                queue_manager.service_queue_manager.release_service_instance('negsample')
                logging.error(f"Negative Samples service run failed: {e}")
                return generate_negative_samples(edges, similarity_threshold=sim_threshold)

        except Exception as e:
            logging.error(f"Negative Samples service orchestration failed: {e}")
            return generate_negative_samples(edges, similarity_threshold=sim_threshold)

    def _wait_and_fetch_results(self, service):
        # Placeholder implementation
        pass

class NetworkModelServiceOrchestrator(BaseServiceOrchestrator):
    def run(self, nodes: List, edges: List, drug_nodes: List, drug_edges: List,
            negs_toggle: bool, run_depth: str, num_jobs: int, seed_input: str,
            input_min_simil: float, sim_threshold: float) -> Tuple[List, List, List, List]:
        try:
            # Check if service can be launched immediately
            launch_status = queue_manager.service_queue_manager.request_service_instance(
                'networkmodel', 
                {
                    'input_seed': self.input_seed,
                    'date': self.date,
                    'base_directory': self.base_directory,
                    'nodes': nodes,
                    'edges': edges,
                    'drug_nodes': drug_nodes,
                    'drug_edges': drug_edges,
                    'negs_toggle': negs_toggle,
                    'run_depth': run_depth,
                    'num_jobs': num_jobs,
                    'seed_input': seed_input,
                    'input_min_simil': input_min_simil,
                    'sim_threshold': sim_threshold
                }
            )
            
            if not launch_status['can_launch']:
                raise Exception(f"Network Model service is queued. Position: {launch_status.get('queue_position', 'unknown')}")

            try:
                if not self.is_service_mode():
                    if negs_toggle:
                        results = run_network_model_with_NS(
                            self.input_seed, self.date, run_jobs=num_jobs,
                            run_depth=run_depth, run_seed=seed_input
                        )
                    else:
                        results = run_network_model(
                            self.input_seed, self.date, run_jobs=num_jobs,
                            run_depth=run_depth, run_seed=seed_input
                        )
                    queue_manager.service_queue_manager.release_service_instance('networkmodel')
                    return results

                env = {
                    'NEGATIVE_SAMPLES': str(int(negs_toggle)),
                    'DEPTH': run_depth,
                    'JOBS': str(num_jobs),
                    'SEED': str(seed_input),
                    'MIN_SIMILARITY': str(input_min_simil),
                    'SIMILARITY_THRESHOLD': str(sim_threshold),
                    'NODES': json.dumps(nodes),
                    'EDGES': json.dumps(edges),
                    'DRUG_NODES': json.dumps(drug_nodes),
                    'DRUG_EDGES': json.dumps(drug_edges)
                }

                service = self.create_service(
                    image='drugrepurposing/networkmodel:latest',
                    env=env
                )

                results = self._wait_and_fetch_results(service)
                queue_manager.service_queue_manager.release_service_instance('networkmodel')
                return results

            except Exception as e:
                queue_manager.service_queue_manager.release_service_instance('networkmodel')
                logging.error(f"Network Model service run failed: {e}")

                if negs_toggle:
                    return run_network_model_with_NS(
                        self.input_seed, self.date, run_jobs=num_jobs,
                        run_depth=run_depth, run_seed=seed_input
                    )
                else:
                    return run_network_model(
                        self.input_seed, self.date, run_jobs=num_jobs,
                        run_depth=run_depth, run_seed=seed_input
                    )

        except Exception as e:
            logging.error(f"Network Model service orchestration failed: {e}")

            if negs_toggle:
                return run_network_model_with_NS(
                    self.input_seed, self.date, run_jobs=num_jobs,
                    run_depth=run_depth, run_seed=seed_input
                )
            else:
                return run_network_model(
                    self.input_seed, self.date, run_jobs=num_jobs,
                    run_depth=run_depth, run_seed=seed_input
                )

    def _wait_and_fetch_results(self, service):
        # Placeholder implementation
        pass
