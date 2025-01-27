"""
QUEUE HELPER FOR DISTRIBUTED SERVICES
Created on January 20th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import os
import logging
from typing import Dict, Any
from collections import deque
import threading
import time
import json

class ServiceQueueManager:
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, parameters_file='Docker_parameters.txt'):
        """
        Singleton implementation to ensure only one queue manager exists
        """
        if not cls._instance:
            with cls._lock:
                if not cls._instance:
                    cls._instance = super(ServiceQueueManager, cls).__new__(cls)
                    cls._instance._initialize(parameters_file)
        return cls._instance

    def _initialize(self, parameters_file):
        """
        Initialize the service queue manager
        
        :param parameters_file: Path to the Docker parameters configuration file
        """
        self.service_queues: Dict[str, deque] = {}
        self.service_limits: Dict[str, int] = {}
        self.active_instances: Dict[str, int] = {}
        self.lock = threading.Lock()
        
        # Load service instance limits from Docker_parameters.txt
        self._load_service_limits(parameters_file)

    def _load_service_limits(self, parameters_file):
        """
        Read service instance limits from Docker_parameters.txt
        
        :param parameters_file: Path to the Docker parameters configuration file
        """
        try:
            with open(parameters_file, 'r') as f:
                for line in f:
                    if '=' in line:
                        service, limit = line.strip().split('=')
                        service_name = service.replace('_INSTANCES', '').lower()
                        self.service_limits[service_name] = int(limit)
                        self.service_queues[service_name] = deque()
                        self.active_instances[service_name] = 0
        except FileNotFoundError:
            logging.error(f"Docker parameters file {parameters_file} not found")
            # Default to a single instance if file is not found
            default_services = ['main', 'monarch', 'dgidb', 'drugsimilarity', 'negsample', 'networkmodel']
            for service in default_services:
                self.service_limits[service] = 1
                self.service_queues[service] = deque()
                self.active_instances[service] = 0

    def request_service_instance(self, service_name: str, request_details: Dict[str, Any]) -> Dict[str, Any]:
        """
        Request a service instance, queuing if maximum instances reached
        
        :param service_name: Name of the service
        :param request_details: Details of the service request
        :return: Dictionary with launch status and additional info
        """
        with self.lock:
            # Normalize service name
            service_name = service_name.lower().replace('serviceorchestrator', '')
            
            # Check if service is in known services
            if service_name not in self.service_limits:
                logging.warning(f"Unknown service: {service_name}. Using default single instance.")
                self.service_limits[service_name] = 1
                self.service_queues[service_name] = deque()
                self.active_instances[service_name] = 0

            # If under instance limit, launch immediately
            if self.active_instances[service_name] < self.service_limits[service_name]:
                self.active_instances[service_name] += 1
                return {
                    'can_launch': True,
                    'active_instances': self.active_instances[service_name],
                    'max_instances': self.service_limits[service_name]
                }
            
            # Otherwise, queue the request
            queue_position = len(self.service_queues[service_name]) + 1
            self.service_queues[service_name].append(request_details)
            
            return {
                'can_launch': False,
                'queue_position': queue_position,
                'active_instances': self.active_instances[service_name],
                'max_instances': self.service_limits[service_name]
            }
    
    def release_service_instance(self, service_name: str):
        """
        Release a service instance and check queue for next request
        
        :param service_name: Name of the service being released
        """
        with self.lock:
            # Normalize service name
            service_name = service_name.lower().replace('serviceorchestrator', '')
            
            # Decrement active instances
            self.active_instances[service_name] = max(0, self.active_instances[service_name] - 1)
            
            # Check if queued requests exist
            if self.service_queues[service_name]:
                next_request = self.service_queues[service_name].popleft()
                self.active_instances[service_name] += 1
                
                # Log the queued request being processed
                logging.info(f"Processing queued {service_name} service request")
                
                # In a more advanced implementation, you might want to:
                # 1. Publish a message to a message queue
                # 2. Call a method to actually launch the service
                # 3. Use an event-driven architecture to trigger service launch
    
    def get_service_status(self, service_name: str) -> Dict[str, Any]:
        """
        Get detailed status for a service
        
        :param service_name: Name of the service
        :return: Dictionary with service status details
        """
        service_name = service_name.lower().replace('serviceorchestrator', '')
        
        return {
            'service_name': service_name,
            'max_instances': self.service_limits.get(service_name, 1),
            'active_instances': self.active_instances.get(service_name, 0),
            'queued_requests': len(self.service_queues.get(service_name, [])),
            'can_launch_new_instance': (
                self.active_instances.get(service_name, 0) < 
                self.service_limits.get(service_name, 1)
            )
        }

# Global queue manager instance
service_queue_manager = ServiceQueueManager()