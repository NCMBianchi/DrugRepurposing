"""
QUEUE HELPER FOR DOCKER SWARM SERVICES
Created on January 20th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

from collections import deque
import docker
import threading
import time
import logging
from typing import Dict, Any
import os
import json

class ServiceQueueManager:
    def __init__(self):
        # Initialize logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Initialize Docker client
        self.client = docker.from_env()
        
        # Initialize queues and tracking
        self.services = ['monarch', 'dgidb', 'drugsim', 'negsample', 'networkmodel']
        self.job_queues = {service: deque() for service in self.services}
        self.active_jobs = {service: 0 for service in self.services}
        self.job_status = {}  # Track status of all jobs
        
        # Read configuration
        self.max_instances = self._read_max_instances()
        
        # Start queue processor
        self._start_queue_processor()
        
    def _read_max_instances(self) -> Dict[str, int]:
        """Read maximum instances from Docker_parameters.txt"""
        try:
            params = {}
            with open('Docker_parameters.txt', 'r') as f:
                for line in f:
                    if '=' in line:
                        key, value = line.strip().split('=')
                        if key.endswith('_INSTANCES'):
                            service = key.replace('_INSTANCES', '').lower()
                            params[service] = int(value)
            return params
        except Exception as e:
            self.logger.error(f"Error reading Docker_parameters.txt: {e}")
            # Default to 1 instance if file can't be read
            return {service: 1 for service in self.services}

    def submit_job(self, job_id: str, service_name: str, job_data: Dict[str, Any]) -> Dict[str, Any]:
        """Submit a job to the queue"""
        if service_name not in self.services:
            raise ValueError(f"Invalid service name: {service_name}")
        
        queue_position = len(self.job_queues[service_name]) + 1
        
        job_info = {
            'job_id': job_id,
            'service': service_name,
            'status': 'queued',
            'queue_position': queue_position,
            'data': job_data,
            'submitted_at': time.time(),
            'started_at': None,
            'completed_at': None,
            'error': None
        }
        
        self.job_queues[service_name].append(job_info)
        self.job_status[job_id] = job_info
        
        self.logger.info(f"Job {job_id} queued for {service_name} at position {queue_position}")
        return job_info

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Get the current status of a job"""
        return self.job_status.get(job_id, {'status': 'not_found'})

    async def _process_job(self, job_info: Dict[str, Any]):
        """Process a single job"""
        job_id = job_info['job_id']
        service_name = job_info['service']
        
        try:
            # Update job status
            job_info['status'] = 'running'
            job_info['started_at'] = time.time()
            
            # Scale up the service
            service = self.client.services.get(f'drugapp_{service_name}')
            current_replicas = service.attrs['Spec']['Mode']['Replicated']['Replicas']
            service.scale(current_replicas + 1)
            
            # Here you would typically wait for the job to complete
            # This might involve monitoring the service's logs or waiting for a callback
            
            # For now, we'll simulate job completion after a delay
            await asyncio.sleep(10)
            
            # Update job status
            job_info['status'] = 'completed'
            job_info['completed_at'] = time.time()
            
            # Scale down the service
            service.scale(current_replicas)
            
        except Exception as e:
            job_info['status'] = 'failed'
            job_info['error'] = str(e)
            self.logger.error(f"Error processing job {job_id}: {e}")
        
        finally:
            self.active_jobs[service_name] -= 1

    async def _process_queues(self):
        """Main queue processing loop"""
        while True:
            for service_name in self.services:
                # Check if we can process more jobs for this service
                while (self.job_queues[service_name] and 
                       self.active_jobs[service_name] < self.max_instances[service_name]):
                    
                    # Get next job
                    job_info = self.job_queues[service_name].popleft()
                    self.active_jobs[service_name] += 1
                    
                    # Process job in background
                    asyncio.create_task(self._process_job(job_info))
                    
                    # Update queue positions for remaining jobs
                    for idx, remaining_job in enumerate(self.job_queues[service_name]):
                        remaining_job['queue_position'] = idx + 1
            
            await asyncio.sleep(1)  # Check queues every second

    def _start_queue_processor(self):
        """Start the queue processor in a separate thread"""
        loop = asyncio.new_event_loop()
        thread = threading.Thread(
            target=self._run_async_loop,
            args=(loop, self._process_queues()),
            daemon=True
        )
        thread.start()

    @staticmethod
    def _run_async_loop(loop, coro):
        """Run an async coroutine in the given event loop"""
        asyncio.set_event_loop(loop)
        loop.run_until_complete(coro)

# Singleton instance
queue_manager = ServiceQueueManager()