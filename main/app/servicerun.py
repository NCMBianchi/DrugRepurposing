"""
SERVICE RUN MODULE: FUNCTIONS NECESSARY TO LAUNCH DOCKER SERVICES
Created on January 22nd 2024
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import json
import time
import logging
from typing import Dict, Any, Optional, Tuple
import uuid
import redis
import multiprocessing
import importlib
import traceback

import queue_manager

import Monarch_fallback
import dgidb_fallback
import drugsimilarity_fallback
import negsamples_fallback
import networkmodel_fallback

class ServiceOrchestrator:
    def __init__(self, redis_host='localhost', redis_port=6379):
        """
        Initialize the service orchestrator with Redis connection

        :param redis_host: Redis server host
        :param redis_port: Redis server port
        """
        self.redis_client = redis.Redis(host=redis_host, port=redis_port, decode_responses=True)
        self.logger = logging.getLogger(__name__)
        self.queue_manager = queue_manager.service_queue_manager

    def _generate_request_id(self) -> str:
        """
        Generate a unique request identifier

        :return: Unique request ID
        """
        return str(uuid.uuid4())

    def enqueue_service_request(self, service_name: str, request_params: Dict[str, Any]) -> str:
        """
        Enqueue a service request and return a request ID

        :param service_name: Name of the service (e.g., 'monarch', 'dgidb')
        :param request_params: Dictionary of request parameters
        :return: Unique request ID
        """

        # First, check with queue manager if service can be launched
        launch_status = self.queue_manager.request_service_instance(
            service_name,
            request_params
        )

        if not launch_status['can_launch']:
            # If cannot launch immediately, add additional queueing logic
            self.logger.info(f"Service {service_name} queued. Position: {launch_status.get('queue_position', 'N/A')}")

        request_id = self._generate_request_id()

        # Prepare request payload
        request_payload = {
            'request_id': request_id,
            'service_name': service_name,
            'params': request_params,
            'status': 'pending',
            'timestamp': time.time(),
            'attempts': 0
        }

        # Store request in Redis queue
        queue_key = f'service_queue:{service_name}'
        self.redis_client.rpush(queue_key, json.dumps(request_payload))

        # Store request details for tracking
        request_key = f'request:{request_id}'
        self.redis_client.hmset(request_key, {
            'service_name': service_name,
            'status': 'pending',
            'timestamp': time.time()
        })

        self.logger.info(f"Enqueued {service_name} service request: {request_id}")
        return request_id

    def _wait_and_fetch_results(
        self,
        request_id: str,
        timeout: int = 300,
        polling_interval: float = 1.0
    ) -> Dict[str, Any]:
        """
        Wait for and fetch results for a given request ID

        :param request_id: Unique request identifier
        :param timeout: Maximum wait time in seconds
        :param polling_interval: Interval between result checks
        :return: Service request results or error information
        """
        request_key = f'request:{request_id}'
        start_time = time.time()
        service_name = None

        try:
            # Retrieve service name before entering the loop
            service_name = self.redis_client.hget(request_key, 'service_name')

            while time.time() - start_time < timeout:
                # Check request status
                request_status = self.redis_client.hget(request_key, 'status')

                if request_status == 'completed':
                    # Fetch and parse results
                    result_key = f'result:{request_id}'
                    result_json = self.redis_client.get(result_key)

                    if result_json:
                        result = json.loads(result_json)

                        # Clean up Redis keys
                        self.redis_client.delete(request_key, result_key)

                        # Release service instance
                        if service_name:
                            self.queue_manager.release_service_instance(service_name)

                        return result

                elif request_status == 'failed':
                    # Fetch error information
                    error_key = f'error:{request_id}'
                    error_info = self.redis_client.get(error_key)

                    if error_info:
                        error_details = json.loads(error_info)

                        # Clean up Redis keys
                        self.redis_client.delete(request_key, error_key)

                        # Raise error after releasing service instance
                        if service_name:
                            self.queue_manager.release_service_instance(service_name)

                        raise RuntimeError(f"Service request failed: {error_details}")

                # Wait before next poll
                time.sleep(polling_interval)

            # Timeout occurred
            if service_name:
                self.queue_manager.release_service_instance(service_name)

            raise TimeoutError(f"Request {request_id} timed out after {timeout} seconds")

        except Exception as e:
            # Ensure service instance is released even if there's an error
            if service_name:
                self.queue_manager.release_service_instance(service_name)
            raise

    def process_service_request(
        self,
        service_name: str,
        request_params: Dict[str, Any],
        timeout: int = 300
    ) -> Dict[str, Any]:
        """
        Comprehensive method to enqueue, wait for, and retrieve service results

        :param service_name: Name of the service
        :param request_params: Request parameters
        :param timeout: Maximum wait time for request
        :return: Service request results
        """
        try:
            # Enqueue the request
            request_id = self.enqueue_service_request(service_name, request_params)

            # Wait for and fetch results
            return self._wait_and_fetch_results(request_id, timeout)

        except (TimeoutError, RuntimeError) as e:
            self.logger.error(f"Service request error: {e}")
            raise

class ServiceWorker:
    def __init__(self, redis_host='localhost', redis_port=6379):
        """
        Initialize the service worker for processing queued requests

        :param redis_host: Redis server host
        :param redis_port: Redis server port
        """
        self.redis_client = redis.Redis(host=redis_host, port=redis_port, decode_responses=True)
        self.logger = logging.getLogger(__name__)

    def _process_single_request(self, request_data: Dict[str, Any]) -> None:
        """
        Process a single service request

        :param request_data: Request payload dictionary
        """
        request_id = request_data['request_id']
        service_name = request_data['service_name']
        request_params = request_data['params']

        try:
            # Dynamically import service request processor
            service_module = __import__(f'{service_name}_service', fromlist=['process_service_request'])
            result = service_module.process_service_request(request_params)

            # Store successful result
            result_key = f'result:{request_id}'
            self.redis_client.set(result_key, json.dumps(result))

            # Update request status
            request_key = f'request:{request_id}'
            self.redis_client.hset(request_key, 'status', 'completed')

            self.logger.info(f"Processed {service_name} service request: {request_id}")

        except Exception as e:
            # Handle and store error
            error_key = f'error:{request_id}'
            error_info = {
                'service': service_name,
                'error_type': type(e).__name__,
                'error_message': str(e)
            }

            self.redis_client.set(error_key, json.dumps(error_info))

            # Update request status
            request_key = f'request:{request_id}'
            self.redis_client.hset(request_key, 'status', 'failed')

            self.logger.error(f"Failed to process {service_name} service request: {request_id}")

    def run(self, services: Optional[list] = None):
        """
        Continuously run workers for specified services

        :param services: List of service names to process
        """
        if services is None:
            services = ['monarch', 'dgidb', 'drugsimilarity', 'negsamples', 'networkmodel']

        workers = []

        for service in services:
            queue_key = f'service_queue:{service}'

            def worker_process(queue_key=queue_key):
                while True:
                    # Blocking pop from Redis queue
                    _, request_json = self.redis_client.blpop(queue_key)
                    request_data = json.loads(request_json)
                    self._process_single_request(request_data)

            # Start worker process
            process = multiprocessing.Process(target=worker_process)
            process.start()
            workers.append(process)

        # Wait for all workers
        for worker in workers:
            worker.join()

class ServiceRunner:
    def __init__(self, base_directory='/app/data'):
        """
        Initialize the service runner with fallback capabilities

        :param base_directory: Base directory for data storage
        """
        self.logger = logging.getLogger(__name__)
        self.base_directory = base_directory

    def _try_microservice_request(self, service_name: str, request_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Attempt to run service request via microservice

        :param service_name: Name of the service
        :param request_params: Request parameters
        :return: Service request results
        """
        try:
            # Dynamically import service request processor
            service_module = importlib.import_module(f'{service_name}_service')
            return service_module.process_service_request(request_params)
        except ImportError:
            self.logger.warning(f"Microservice for {service_name} not found. Falling back to legacy method.")
            return None
        except Exception as e:
            self.logger.error(f"Microservice request failed: {e}")
            return None

    def _run_fallback_method(self, service_name: str, request_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run legacy method for the service

        :param service_name: Name of the service
        :param request_params: Request parameters
        :return: Service request results
        """
        try:
            # Map service names to legacy function names
            legacy_map = {
                'monarch': 'run_monarch',
                'dgidb': 'run_dgidb',
                'drugsimilarity': 'run_drugsimilarity',
                'negsamples': 'run_negsamples',
                'networkmodel': 'run_networkmodel'
            }

            # Dynamically import legacy module and function
            legacy_module = importlib.import_module(f'{service_name}_fallback')
            legacy_func = getattr(legacy_module, legacy_map[service_name])

            # Prepare arguments for legacy function
            monarch_input = request_params.get('monarch_input', 'MONDO:0007739')
            date = request_params.get('date', time.strftime('%Y-%m-%d'))

            # Call legacy function with appropriate arguments
            results = legacy_func(monarch_input, date)

            # Return results in a standardized format
            return {
                'status': 'success',
                'results': results
            }
        except Exception as e:
            self.logger.error(f"Legacy method failed: {traceback.format_exc()}")
            return {
                'status': 'error',
                'message': str(e)
            }

    def process_service_request(
        self,
        service_name: str,
        request_params: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Attempt to run service request with microservice fallback to legacy method

        :param service_name: Name of the service
        :param request_params: Request parameters
        :return: Service request results
        """
        # First, try microservice approach
        microservice_result = self._try_microservice_request(service_name, request_params)

        if microservice_result:
            return microservice_result

        # If microservice fails, fall back to legacy method
        self.logger.warning(f"Falling back to legacy method for {service_name}")
        return self._run_fallback_method(service_name, request_params)

def main():
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Choose one primary purpose: either run workers or demonstrate service request
    # Typically, this would be handled by separate scripts or deployment configs

    # Option 1: Run service workers
    worker = ServiceWorker()
    worker.run()

    # Option 2: Example service request (not typically in the same script as workers)
    # runner = ServiceRunner()
    # result = runner.process_service_request(
    #     'monarch',
    #     {'monarch_input': 'MONDO:0007739', 'date': time.strftime('%Y-%m-%d')}
    # )
    # print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
