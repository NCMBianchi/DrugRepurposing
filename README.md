# DrugRepurposing

## Notes on Pipeline Update, by Niccolò Bianchi
I am currently making modifications to the code in this tool due to changes to Monarch API, which is an integral part of the tool's pipeline.

The tool is <b>currently not working</b> neither in the old version (support to the old Biolink API was discontinued on March 27th 2024, although it was even planned to be discontinued since March 20th 2023 –so it had one extra year of life), nor the new one as some minor changes are still required.

![README_image-1_bis](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/7db18469-0998-42f7-8bb2-23ee2af35f3c)

Moreover, the new API <b>no longer supports OMIM</b> as an input seed for queries, only their own URIs (<i>e.g.</i> MONDO). Changes were made, and the tool now is able to make an API call to the new V3 version.

![README_image-2](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/86c5166f-fac8-42b5-8e62-5a00be65c074)

The tool is currently still not working: some functions are not correctly handling some of the information.

![README_image-5b py](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/cdb9be49-791b-4e74-a5f6-875e6547e53f)

At first I wanted to bypass the <b>`run_monarch()`</b> macro-function by fetching older outputs from Carmen Reep's own runs in 2022 –by making a <b>`run_monarch_mock()`</b> function that only copies and renames reference files in a dedicated subdirectory. Yet this doesn't work, as the same issue stopping `run_monarch()` also halts <b>`run_monarch_symptom()`</b>. So I went back at fixing the steps that lead to an empty dataframe.

![README_image-3a](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/5a3af0d8-8ca1-4e94-8298-76a556d36c80)

Parts of the <b>`keep_node_type()`</b> function within `get_orthopheno_list()` are not working as the <i>RO:</i> URIs listed no longer are supported and have to be updated. In order to temporarily bypass this, I made a <b>`keep_node_type_mock()`</b> function that returns the entire set.

![README_image-6a](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/f1e4f9a2-31ab-421e-9808-db14ef70d192)

Also parts of the <b>`filter_edges()`</b> function are not working, as matches cannot be found and the resulting `keep` set was empty. In order to temporarily bypass this, I made a <b>`filter_edges_mock()`</b> function that returns the entire set.

![README_image-7a](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/bb3a8810-4e4a-4653-8caa-a2ec16a1491c)

Finally when running the <b>`add_attributes()`</b> function within the `get_orthopheno_list()` the resulting `metaedges` set would also be empty. Apparently no match between IDs can be found between those in the `edges` list and those in lists of subjects and objects. I could not bypass this, so this is where I am currently stuck.

![README_image-9d](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/cd4f7078-2606-4a74-8fbb-5e2fefe0956c)

I added a few lines that create <b><i>.csv</b></i> output files and based on that matching IDs seem to be present, although not in the same order.

![README_image-8a](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/d8c6e7c6-90c0-47c0-9fb9-dea8044f910c)

Since the original code deleted all temporary <i>.csv</i> files when durring the final <b>`rdf2vec.py`</b> script, and those are now not available for checks and workarounds, I commented out the lines that delete such files in order to then check them and store them at least for a few test runs when I will be able to get to that point.

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

The whole drug repurposing workflow is saved in five separate Python files: <code>Monarch.py</code>, <code>DGIdb.py</code>, <code>drugsimilarity.py</code>, <code>combine_graphs.py</code>, and <code>rdf2vec.py</code>. To run this whole process with Flask, three extra Python files were created: <code>\_\_init\_\_.py</code>, `routes.py`, and `forms.py`. 

The `__init__.py` file imports the `Flask` class from the `flask` package version 1.1.2 and instantiates the Flask app. The `routes.py` file specifies the routes of the app, which are acts that bind a URL to a view function (a function that responds to a request). The `forms.py` file includes form classes for user input and uses WTForms to render and validate the forms.

Each time the user makes a request (by clicking on a link or giving input via a form), the `render_template()` function is run, which generates output from a template file that is found in the \`templates\' folder of the app. For a nice general style of the app, each template file extends the `layout.html` template, which provides a site header and title and a section with links of all previously run predictions.

To start the server, the `run()` method of the `Flask` object is called. It returns the URL where the server is available. When opening the URL, the home page of the app is shown, which is rendered from the `home.html` template file. On this home page, there are several options for the user. The user can choose an existing disease graph, which has been previously created using a disease seed as input for Monarch, where the date of creation of the graph is specified. If the user is interested in a disease that is not in this list of previously created graphs, or the user wants a newer version of a graph, the user can create a new disease graph by specifying the phenotype MIM number of the disease of interest (e.g. \`143100\' for Huntington's disease). This enables the Python `Monarch.py` file to run the function that creates a new network from Monarch using this input number as seed. Creating a new network form Monarch takes some time, so the user is asked to be patient. After running a new disease graph, the graph is saved and added to the existing disease graphs on the homepage.

After selecting an existing disease graph, or creating a new disease graph, the `symptoms.html` template is rendered. Here, all symptoms of the disease of interest are listed, and the user is asked to select one symptom of interest. When a user selects a symptom, again, the Python `Monarch.py` file is called to run the function that creates a new network from Monarch, with this time only the chosen symptom as seed. After creating this symptom Monarch graph, this graph is merged with the disease Monarch graph to create the final Monarch graph. Then `DGIdb.py`, `drugsimilarity.py`, `combine_graphs.py`, and `rdf2vec.py` are run in this order. The last Python file saves the ranked predicted drug-gene interactions as an HTML file in the `templates` folder of the app. When everything is run, this predictions HTML file is rendered and the list of ranked predictions is shown on screen. This list includes the drugs, the genes each drug interacts with, together with the interaction type and confidence. Each drug and gene is presented with its label. Clicking on a label will direct the user to the URI link of the drug or gene, which shows all up-to-date information about that drug or gene.

On every page (home, symptom, prediction), we added the option to select a previously run prediction, to avoid long running times.

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
