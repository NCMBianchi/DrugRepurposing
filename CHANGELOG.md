# Changelog
All notable changes to this project will be documented in this file.

## [2025.0.9] - 2025-02-14
### Released
- jupyter-notebook-2025.tar.gz: updated pipeline up to the negative samples step

### Changed
- .gitignore file to avoid also uncompressed folders

## [2025.0.8] - 2025-02-05
### Released
- jupyter-notebook-2025.tar.gz: updated pipeline up to the DGIdb step

### Changed
- Consolidated CSS styles from individual HTML files into main.css
- Fixed advancement bar height in current_runs.html
- Fixed OMIM-conversion wrapper uncollapse button positioning in home.html

### Removed
- Unused df_style.css from static/css
- Unused about.html template 
- Cached Python files from app/__pycache__

### Known Issues
- Certain Python files in the 2025.0.2-legacy release

## [2025.0.2-legacy] - 2025-01-31
### Released
- DrugRepurposing_legacy.tar.gz: Docker image with working base UI
- jupyter-notebook.tar.gz: Jupyter Notebook with a full explanation of the pipeline

## [2025.0.7] - 2025-01-31
### Added
- CHANGELOG.md to track changes (also retroactively)
- run_status.html to display specific run status
- current_runs.html to track ongoing runs
- status_update.js for realtime status updates

### Changed
- Modified routes.py for advanced UI operations (run cancellation, status checking)
- Modified routes.py to allow for and track concurrent runs

### Removed
- Jupyter Notebook implementation from 2024.1.5 commit (now available in 2025.0.2-legacy release)
- Single Docker structure from 2025.0.2 commit (now available in 2025.0.2-legacy release)

### Known Issues
- Distributed services implementation still has to be tested
- Some issues with service orchestration and queuing remain to be fixed
- UI errors in the Current Runs page (i.e. not-centered title, missing right button) 

## [2025.0.6] - 2025-01-30
### Added
- Full implementation of distributed services files (Dockerfile, requirements.txt, app directories)
- New '--parameters "(n,n,n,n,n,n)"' command for Docker instance control

### Changed
- Updated routes.py, servicerun.py and queue_manager.py for distributed services tracking and activation

## [2025.0.5] - 2025-01-27
### Added
- New servicerun.py for distributed service Images
- New queue_manager.py for service instance management based on Docker_parameters.txt

### Changed
- Preserved Flask App functionality
- Renamed platform.py to logger_utils.py to avoid package clashes
- Updated drugapp-launcher.sh for simpler commands
- Updated docker-compose.yml for shared data folder
- Split service files into '_legacy' and '_service' versions
- Modified routes.py for distributed services and queueing

## [2025.0.4] - 2025-01-24
### Changed
- Reorganized structure for proper Docker image setup

## [2025.0.3] - 2025-01-20
### Added
- Dockerfile.[service] files for Docker Swarm architecture
- docker-compose.yml for Docker Swarm architecture
- entrypoint.sh for Docker Swarm architecture
- drugapp-launcher.sh for ease-of-use
- unique.py for Docker Swarm architecture
- queue-manager.py for Docker Swarm architecture and scalability

### Changed
- Updated requirements.txt for Docker Swarm architecture
- Modified module files for Drug Repurposing in micro-services

## [2025.0.2] - 2025-01-17
### Added
- Additional JavaScript for input validation

### Changed
- Made OMIM converter collapsable
- UI improvements

## [2025.0.1] - 2025-01-16
### Added
- OMIM to MONDO converter box

### Fixed
- Minor UI issues

## [2025.0.0] - 2025-01-03
### Changed
- Updated README.md ahead of future code changes

## [2024.1.5] - 2024-09-20
### Changed
- Finalized Flask app implementation

### Known Issues
- Some input handling issues remain to be fixed

## [2024.1.4] - 2024-08-26
### Changed
- Minor updates to network_model.py
- Working pipeline in Jupyter notebook

### Added
- Draft Flask browser app in Docker image

## [2024.1.3] - 2024-08-09
### Changed
- Complete rewrite of rdf2vec.py (renamed to network_model.py)
- Overhauled ML approach

### Added
- Steps to compute negative samples

## [2024.1.2] - 2024-05-23
### Changed
- Complete rewrite of rdf2vec.py (renamed to network_model.py)
- Fully working pipeline in Jupyter notebook

## [2024.1.1] - 2024-05-10
### Changed
- Complete rewrite of drugsimilarity.py and combine_graphs.py
- Working pipeline in Jupyter notebook

## [2024.1.0] - 2024-05-04
### Changed
- Complete rewrite of Monarch.py and DGIdb.py functions
- Simplified code while maintaining equivalent results
- Added Jupyter Notebook implementation

### Removed
- Temporary metadata addition

## [2024.0.5] - 2024-04-17
### Added
- keep_node_type_mock() function in Monarch.py
- filter_edges_mock() function in Monarch.py

### Changed
- Reinstated original run_monarch() function
- Fixed orthopheno network file naming
- Temporarily disabled file removal in rdf2vec.py

## [2024.0.4] - 2024-04-15
### Added
- run_monarch_mock() function for testing
- CSV file saving for intermediate objects

### Known Issues
- First intermediate file missing 'relations' element
- Tool breaks after second run

## [2024.0.3] - 2024-04-05
### Fixed
- Coding error in build_edges() and build_nodes() functions
- Return statement now correctly returns lists instead of dataframes

### Removed
- Unused lines from original bioknowledge reviewer tool

## [2024.0.2] - 2024-03-28
### Changed
- Modified run_monarch() for MONDO URIs instead of OMIM
- Updated query URL for V3 API
- Updated JSON parsing for V3 API structure
- Renamed MIMPhenotype class to URIphenotype

### Fixed
- SSLError in get_disease_name_id()

## [2024.0.1] - 2024-03-27
### Fixed
- Dependencies in Dockerfile
- API format compatibility issues