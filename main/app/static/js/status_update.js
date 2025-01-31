document.addEventListener('DOMContentLoaded', function() {
    const runStatusContainer = document.querySelector('.run-progress');
    const runId = document.querySelector('[data-run-id]').dataset.runId;
    const startTime = new Date(document.getElementById('start-time').dataset.startTime);

    // Color mapping for stage status
    const statusColorMap = {
        'pending': '#6c757d',      // secondary
        'running': '#17a2b8',      // info
        'completed': '#28a745',    // success
        'stopped': '#dc3545',      // danger
        'failed': '#dc3545',       // danger
        'skipped': '#ffc107'       // warning
    };

    // Update progress bar and stages
    function updateRunStatus() {
        fetch(`/run_status_update/${runId}`)
            .then(response => response.json())
            .then(data => {
                // Update total run time
                updateRunTime(startTime);

                // Update progress stages
                updateProgressStages(data.stages);

                // Update overall status
                updateOverallStatus(data.overall_status);

                // Conditionally show/hide stop/delete buttons
                updateActionButtons(data.overall_status);
            })
            .catch(error => {
                console.error('Error updating run status:', error);
            });
    }

    // Update run time display
    function updateRunTime(startTime) {
        const now = new Date();
        const diffMs = now - startTime;
        const runTimeElement = document.getElementById('run-time');

        if (diffMs < 600000) { // Less than 10 minutes
            runTimeElement.textContent = `${Math.floor(diffMs / 1000)}s`;
        } else {
            runTimeElement.textContent = `${Math.floor(diffMs / 60000)} min`;
        }
    }

    // Update progress stages
    function updateProgressStages(stages) {
        const stageOrder = ['monarch', 'dgidb', 'drugsimilarity', 'negsample', 'networkmodel'];
        
        stageOrder.forEach((stageName, index) => {
            const stageMarker = document.querySelector(`.stage-marker[data-stage="${stageName}"]`);
            const stageProgressBar = stageMarker.closest('.progress-stage');
            
            // Find the corresponding stage data
            const stageData = stages.find(s => s.name === stageName);

            if (stageData) {
                // Update stage label and color
                const stageLabel = stageMarker.querySelector('.stage-label');
                
                switch (stageData.status) {
                    case 'pending':
                        stageLabel.textContent = stageName.replace('_', ' ');
                        stageProgressBar.style.backgroundColor = statusColorMap['pending'];
                        break;
                    case 'running':
                        stageLabel.textContent = `Running: ${stageName.replace('_', ' ')}`;
                        stageProgressBar.style.backgroundColor = statusColorMap['running'];
                        break;
                    case 'completed':
                        stageLabel.textContent = `${stageName.replace('_', ' ')} (${stageData.runtime})`;
                        stageProgressBar.style.backgroundColor = statusColorMap['completed'];
                        break;
                    case 'skipped':
                        if (stageName === 'negsample') {
                            stageLabel.textContent = 'SKIPPED';
                            stageProgressBar.style.backgroundColor = statusColorMap['skipped'];
                        }
                        break;
                }
            }
        });
    }

    // Update overall status badge
    function updateOverallStatus(status) {
        const statusBadge = document.getElementById('overall-status');
        statusBadge.textContent = status.toUpperCase();
        
        // Update badge class based on status
        statusBadge.className = statusBadge.className.replace(/badge-\w+/, `badge-${getStatusClass(status)}`);
    }

    // Get Bootstrap badge class for status
    function getStatusClass(status) {
        switch(status) {
            case 'running': return 'info';
            case 'completed': return 'success';
            case 'stopped': return 'warning';
            case 'failed': return 'danger';
            default: return 'secondary';
        }
    }

    // Update action buttons visibility
    function updateActionButtons(status) {
        const stopRunButton = document.getElementById('stop-run-button');
        const deleteRunButton = document.getElementById('delete-run-button');

        if (status === 'running') {
            stopRunButton.style.display = 'inline-block';
        } else {
            stopRunButton.style.display = 'none';
            deleteRunButton.textContent = 'Get Results';
            deleteRunButton.classList.remove('btn-danger');
            deleteRunButton.classList.add('btn-primary');
        }
    }

    // Initial update
    updateRunStatus();

    // Periodic updates every 5 seconds
    const updateInterval = setInterval(updateRunStatus, 5000);

    // Stop updates when run is no longer running
    function stopUpdates() {
        clearInterval(updateInterval);
    }
});