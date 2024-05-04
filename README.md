# DrugRepurposing

## Notes on Pipeline Update, by Niccolò Bianchi
Due to many issues with the tool –further notes available in the README.md files in previous commits– I developed new code that is aimed at performing the same tasks as in the python files (<i>i.e.<i/> ```Monarch.py```, ```DGIdb.py```, ```drugsimilarity.py```, ```combine_graphs.py``` and ```rdf2vec.py```) in a Jupyter Notebook.

After writing the whole code, functions will once again be divided in python scripts that can be called by the code in <b>```routes.py```<b/> which manages the code after input in the flaskApp: I indeed still want this tool to run as a self-contained docker that can be:
	- either installed locally by tech-savy users that want to customise operations and tweak the code,
	- or accessed remotely online through a client-server connection while the flaskApp runs in a server.

![README_image-1_bis](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/7db18469-0998-42f7-8bb2-23ee2af35f3c)

So far I solved issues with the new V3 API by Monarch Initiative –replacing the old Biolink API: some of the metadata information (<i>e.g.<i/> ```RO:``` for relational ontology) are no longer available, nor does the new API accept ```OMIM:``` IDs as inputs. The new version uses <b>```MONDO:```<b/> IDs as input for the diseases –mostly focusing on Huntington's Disease, ```MONDO:0007743```.

This section originally was in file <b>```Monarch.py```<b/> that managed API calls to build the network of associations between disease, phenotypes and genes.

![README_image-2](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/86c5166f-fac8-42b5-8e62-5a00be65c074)

I also skipped most the redundant steps that were in the original code: the results is a more linear pipeline. Yet by removing some of the metadata-adding steps –which broke the original code, though– output information is slightly less wide. Now the output files have these structures:

| col1 NODEs    | col2 NODEs    |
| node_id       | node_name     |
| ------------- | ------------- |
|               |               |

The <b>```nodes```<b/> object is therefore a list of dictionaries (<i>i.e.<i/> {id: label}) for each node.

| col1 EDGEs    | col2 EDGEs    | col3 EDGEs    | col4 EDGEs    | col5 EDGEs    |
| subj_id       | subj_name     | relation      | obj_name      | obj_name      |
| ------------- | ------------- | ------------- | ------------- | ------------- |
|               |               |               |               |               |

The <b>```edges```<b/> object instead is a list of lists for each edge, some elements of which are dictionaries (<i>i.e.<i/> {subj}{rel}{obj}). 

![IMAGE_placeholder]

I also implemented the functions originally in <b>```DGIdb.py```<b/>, which manages calls to the Entrez API to convert gene IDs from the ones used by Monarch Initiative (<i>i.e.<i/> ```HGNC:```, ```ZFIN:```), to the Entrez ones accepted by DGIdb's API to obtain interactions between genes and drugs.

Those calls took quite a lot: ~15-20min for Monarch Initiative (65MB) and ~200min for ENTREZ and DGIdb (7MB). Initial seed was <b>```MONDO:0007743```<b/>, three layers of neighbours for the association network:

| Nodes         | Edges         |
| ------------- | ------------- |
| 7.375         | 12.726        |

Following association also to drugs starting from about 3K genes:

| Nodes         | Edges         |
| ------------- | ------------- |
| 9.975         | 16.576        |

![README_image-4b](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/87091903-4416-40b4-a6a1-9c1f6b7334a3)

Later on I will also work on more changes in the <b>.html</b> and <b>.css</b> files used by the flask App, in order to enhance the UI of the tool.

## Original Description, by Carmen Reep
Automated drug repurposing workflow for rare diseases using Flask.

This drug repurposing app finds candidate compounds for a symptom of a specific rare disease.
For this, it uses a workflow that has four main steps:
1. Creation of an RDF knowledge graph by selecting a subnetwork from Monarch Initiative, adding drugs using the DGIdb database, creating drug-drug edges based on drug compound structure similarity.
2. Graph embedding using RDF2Vec.
3. Creation of edge representations (=training/test data) by fusing drug and gene feature vectors.
4. Training an XGBoost machine learning model which is evaluated using repeated stratified 10-fold cross validation, where the best model is used to predict unknown drug-gene interactions of interest.

(...)

## Getting started
### Docker
Build docker image: <code>docker build --tag drugapp .</code>

This might take a while (5-20min depending on your internet connection) due to the entire image size being almost 8GB due to some very big python packages.

Run the image: <code>docker run -d -p 5000:5000 --name drugapp drugapp</code>

### <s>Without docker</s> (no longer suggested: requires out-of-date python packages)
<s>Download the `app` folder, then run the Flask app with <code>python run.py</code>. </s>

## Aknowledgements
### LUMC and LIACS (NL)
[Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198)

### Politecnico di Milano (IT)
[Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ)

### SWAT4HCLS
Presentation by Carmen Reep on the original work this fork repository is based on: [Automated drug repurposing workflow for rare diseases](https://youtu.be/RsfUrRhZAso?si=Og1z1RdPaukpPIbP)
