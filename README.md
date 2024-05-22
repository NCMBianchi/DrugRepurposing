# DrugRepurposing

## Notes on updating the pipeline (Niccolò Bianchi)
Due to many issues with the tool –further notes available in the README.md files in previous commits– I developed new code that is aimed at performing the same tasks as in the python files (<i>i.e.</i> ```Monarch.py```, ```DGIdb.py```, ```drugsimilarity.py```, ```combine_graphs.py``` and ```rdf2vec.py```) in a Jupyter Notebook.

After writing the whole code, functions will once again be divided in python scripts that can be called by the code in <b>```routes.py```</b> which manages the code after input in the flaskApp: I indeed still want this tool to run as a self-contained docker that can be:
	- either installed locally by tech-savy users that want to customise operations and tweak the code,
	- or accessed remotely online through a client-server connection while the flaskApp runs in a server.

![README_image-1_bis](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/7db18469-0998-42f7-8bb2-23ee2af35f3c)

### ```Monarch.py``` FIXED
Issues with the new V3 API by Monarch Initiative are solved. Runtime is about 15min.

Some of the metadata information (<i>e.g.</i> ```RO:``` for relational ontology) are no longer available, and the new API doesn't accept ```OMIM:``` IDs as inputs. Nonetheless these information were just metadata, and were not necessary.

The new version uses <b>```MONDO:```</b> IDs as input for the diseases –mostly focusing on Huntington's Disease, ```MONDO:0007743```. Most of the redundant steps were removed, and current output ```.csv``` files have these structures:

| node_id       | node_name     |
| ------------- | ------------- |
|               |               |

The <b>```nodes```</b> object is therefore a list of dictionaries (<i>i.e.</i> ```{id: label}```) for each node.

| subj_id       | subj_name     | relation      | obj_name      | obj_name      | notes         |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|               |               |               |               |               |               |

The <b>```edges```</b> object instead is a list of lists for each edge, some elements of which are dictionaries (<i>i.e.</i> ```[{subj}{rel}{obj}{notes}]```).

These data-structures are retained also in the next steps.

### ```DGIdb.py``` FIXED
This code section manages calls to the Entrez API to convert gene IDs from the ones used by Monarch Initiative (<i>i.e.</i> ```HGNC:```, ```ZFIN:```), to the Entrez ones accepted by DGIdb's API to obtain interactions between genes and drugs.

Running this section can take up to 4h if fetching gene-to-drug associations for all 3 layers in the Monarch API call, so the script in the Jupyter notebook only runs layers 1 and 2 –thanks to an optional variable in main function ```run_dgidb()``` that allows to select up to which layer to consider. Runtime gets to about 30min.

Furthermore, I also already made functions that would use DGIdb's new V5 GraphQL API –as the old V2 API will be no longer supported following June 1st 2024. This API does not require Entrez API anymore, as it uses gene IDs found in Monarch. The result is that this step is now much faster: a full 3-layer run now takes about 35-40min instead of almost 4h.

### ```drugsimilarity.py``` FIXED
Similarity between drugs found in the previous step is calculated by first converting drug IDs into SMILES notations (<i>i.e.</i> strings corresponding to molecular structures), then RDKit objects (<i>i.e.</i> 2D/3D molecular structure objects) and finally in ECFP bitvector (<i>i.e.</i> Extended Connectivity FingerPrints) which store information on whether certain features are present or not.

Based on that, similarity is computed based on Tanimoto's coefficient.

### ```combine_graphs.py``` SKIPPED
This section is no longer necessary, as data objects have been structured coherently throughout the entire pipeline.

### ```rdf2vec.py``` FIXED (renamed as ```network_model.py```)
Due to issues with PyPi package ```pyrdf2vec``` requireing a version of ```torch``` that is no longer supported, I switched to a non-RDF approach via ```node2vec```: data is still structured as a knowledge graph. I changed the logic of which data the model is trained (<i>i.e.</i> no longer on a subset of the drug-to-gene associations, but their entirety) and of which the prediction is made on (<i>i.e.</i> no longer the other subset of drug-to-gene associations, but the entire biological and drug associations network, plus all other drugs on DGIdb) in order to effectively 'discover' new repurposed drugs outside of those already related to the disease of interest.

![Huntington disease_2024-05-22_full_network_plot](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/04788652-457c-4a7a-977a-7b2debe3b391)

## Notes on the Docker Image / FlaskApp (Niccolò Bianchi)

Later on I will also work on more changes in the <b>.html</b> and <b>.css</b> files used by the flask App, in order to enhance the UI of the tool.

![README_image-4b](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/87091903-4416-40b4-a6a1-9c1f6b7334a3)

## Original description (Carmen Reep)
Automated drug repurposing workflow for rare diseases using Flask.

This drug repurposing app finds candidate compounds for a symptom of a specific rare disease.
For this, it uses a workflow that has four main steps:
1. Creation of an RDF knowledge graph by selecting a subnetwork from Monarch Initiative, adding drugs using the DGIdb database, creating drug-drug edges based on drug compound structure similarity.
2. Graph embedding using RDF2Vec.
3. Creation of edge representations (=training/test data) by fusing drug and gene feature vectors.
4. Training an XGBoost machine learning model which is evaluated using repeated stratified 10-fold cross validation, where the best model is used to predict unknown drug-gene interactions of interest.

(...)

## INSTRUCTIONS: how to use the tool
### Docker
Build docker image: <code>docker build --tag drugapp .</code>

This might take a while (5-20min depending on your internet connection) due to the entire image size being almost 8GB due to some very big python packages.

Run the image: <code>docker run -d -p 5000:5000 --name drugapp drugapp</code>

## Aknowledgements
### LUMC and LIACS (NL)
[Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198)

### Politecnico di Milano (IT)
[Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ)

### SWAT4HCLS
Presentation by Carmen Reep on the original work this fork repository is based on: [Automated drug repurposing workflow for rare diseases](https://youtu.be/RsfUrRhZAso?si=Og1z1RdPaukpPIbP)
