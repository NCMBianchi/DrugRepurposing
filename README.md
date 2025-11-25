# DrugRepurposing

![BANNER](https://github.com/user-attachments/assets/e237ef0f-4b7c-42f3-b7e1-724e1b68aacf)

This ongoing project focuses on developing a <b>modular</b> and <b>automated</b> drug-repurposing pipeline that builds on the huge potential of [Knowledge Graphs](https://en.wikipedia.org/wiki/Knowledge_graph) and Machine Learning (ML) predictions.

It is based on my [MSc thesis](https://www.dropbox.com/scl/fi/6vzgfld7riqb19hm5wj6u/DRUG_REPURPOSING-thesis.pdf?rlkey=y9xm7zuxm4q035byvhe496m9k&dl=0) ([datasets from a reference full run](https://www.dropbox.com/scl/fi/prvqajjau227741z5ve91/DRUG_REPURPOSING-data.7z?rlkey=qjumdz9r93y0yv6mhc21a7bir&dl=0)), forked from [Carmen Reep's work](https://github.com/carmenreep/DrugRepurposing) which elaborated over [Núria Queralt-Rosinach's 'bioknowledge-reviewer' tool](https://github.com/NuriaQueralt/bioknowledge-reviewer). The `python_scripts` folder contains the latest versions of the functions for the pipeline. My main contrinutions to the project are:
- adding semantically-valid negative samples via the `negsamples.py` scripts, inspired by [Nicolas Hubert's 'semantic-lossfunc' project](https://github.com/nicolas-hbt/semantic-lossfunc) ([Hubert et al. 2023](https://arxiv.org/abs/2301.05601)).
- fixing minor internal and external coding issues (<i>i.e.</i> database API calls, unsupported package [pyRDF2vec](https://pypi.org/project/pyrdf2vec/), computation of drug-to-drug edges based on [RDkit](https://www.rdkit.org/docs/index.html) and 2048-bit feature similarity
- adding Torch-geomtric as an alternative ML model (<i>i.e.</i> `networkmodel_pyg.py`) for gene-to-drug edge prediction, along with XGBoost (<i>i.e.</i> `networkmodel_xgb.py`)
- adapting the node embedding via [Node2vec](https://github.com/eliorc/node2vec) step to accomodate for different ML methods (<i>i.e.</i> `embeddings_xgb.py` and `embeddings_pyg.py`)
- implementing the iteration of ML predictions with different seed inputs to improve prediction efficiency and consistency by pooling the results and ranking them based on prediction probability and occurrence across runs, which is the <b>topic of an upcoming publication</b> ([draft abstract](https://www.dropbox.com/scl/fi/em4b9mfybq5q7txjpb6ug/ABSTRACT-drug_repurposing.pdf?rlkey=5vi8ipmi8pub5rfgedmdfhvar&dl=0))


## Relevant Releases
Different aspects of the project are available in the releases sections. Below summaries of the most relevant ones.

### [COMPLETED] [Fully functioning drug-repurposing pipeline in Jupyter Notebook](https://github.com/NCMBianchi/DrugRepurposing/releases/tag/v2025.0.11)
The `jupyter-notebook-2025.ipynb` notebook takes a set of parameters (<i>e.g.</i> disease ID, similarity threshold, ML method) and then automatically performs all the steps returning a list of repurposed drugs based on predicted gene-to-drug edges.

![IMAGE-release_20250307](https://github.com/user-attachments/assets/314d7bc3-33ee-473a-8780-24faf3962253)

### [COMPLETED] [ML iterative predictions in Jupyter Notebook: testing and comparison](https://github.com/NCMBianchi/DrugRepurposing/releases/tag/v2025.0.12)
The `network-generator.ipynb` notebook builds a smaller sub-network provided of a bigger network based on centrality measure distribution, while `test-runs.ipynb` and `parse-runs.ipynb` are used to evaluate prediction efficiency and consistency across runs with different parameters.

![parse-runs-2](https://github.com/user-attachments/assets/60b54529-8b3f-458c-a3a0-489bf8920754)

### [COMPLETED] [Additional validation Notebooks for the core scripts](https://github.com/NCMBianchi/DrugRepurposing/releases/tag/v2025.1.0)
We implemented new notebooks to validate the subset network (<i>i.e.</i> `network-generator-validation.ipynb`) used in iteration-consensus runs, and  our approach to generate negative samples compared to other methods (<i>i.e.</i> `GCN_negsamples-ablation.ipynb`). We also implemented new notebooks to assess the prediction metrics (<i>e.g.</i> `XGB_raw_10it-evalmetrics.ipynb`) and the overlap of predictions (<i>e.g.</i> `networkmodel-convergence-full.ipynb`) between different models and methods.

![release_notebooks-1](https://github.com/user-attachments/assets/937444a9-6cd3-4a29-8b15-17067d0bfb31)

### [COMPLETED] [No-UI Docker Image to run all the notebooks]()
This Docker Image (~110MB) allows to create a Container via Docker Compose that installs the correct version of Python and all dependencies to run all the notebooks and core scripts in a fully isolated environment.
Datasets from our own runs of the notebooks and scripts (~45GB) are available on Zenodo: []().

![]()

### [ONGOING] [Scalable and modular iterative pipeline in Docker-compose](https://github.com/NCMBianchi/DrugRepurposing/releases/tag/v2025.0.12-docker-compose)
Pipeline functionalities are split in several modules that manage 1-2 `.py` scripts, and such modules are then implemented in multiple instances of Docker containers built by Docker-compose and managed by a `Celery_app` container for cueing and data transfer. This also allows for further parallelisation and resource optimisation. The end result will be merged with the FlaskApp-based web UI docker.

![GITHUB_update-20250407-1](https://github.com/user-attachments/assets/2444141e-0b9c-4b5b-afff-f819e9b6444e)

### [PAUSED] [FlaskApp-based Web UI in Docker](https://github.com/NCMBianchi/DrugRepurposing/releases/tag/v2025.0.12-docker-ui)
All the files within `app/services/`, `app/static/` and `app/templates/` are fully working, but the pipeline functionalities are outdated and will be updated once the Docker-compose modular version is finished.

![DRUGAPP_home-4-EDIT](https://github.com/user-attachments/assets/74c5c9e7-d4f9-4fc2-8c9c-483c72c075d7)

## Aknowledgements
<b>LUMC and LIACS (NL)</b> [Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198); <b>Politecnico di Milano, DEIB (IT)</b> [Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ); <b>non-acaedemic</b> [Nicolas Hubert](https://scholar.google.com/citations?user=nHtB06wAAAAJ).

## License
This project is licensed under the MIT License. Feel free to use and modify the code as per your needs.

