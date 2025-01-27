# DrugRepurposing

![README_image-1_bis](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/7db18469-0998-42f7-8bb2-23ee2af35f3c)

## 2025 UPDATES: scalability (iterations of dockerised modules), additional neural network ML option, batch runs

Here's how the tool works in its 2024.1.5 release:

![GITHUB_update-20250103-1](https://github.com/user-attachments/assets/88773663-84a1-422f-bf92-1e0a926bcf2d)

I am currently developing a new version of the tool that allows for:
- **scalability**, so that if this tool would ever be implemented in a server and multiple users would access it, they could all run their queries cuncurrently and avoid long cues as the modules that require more time would exist in multiple copies with their own dedicated core
- a new ML approach based on **graphic neural networks** (GNNs) so that the user could choose from a dropdown menu whether to run a gradient boosting method (*i.e.* with Python package `XGboost`) or a neural network method (*i.e.* with Python package `PyG`, a.k.a. PyTorch Geometric)
- **batch runs**, so that the user can launch a series of *n* runs with different seeds and obtain a distribution of results from the gene-to-drug prediction, and determine which of the repurposed drugs were returned by the most runs: this is meant to solve the issue of having a certain range of variability based on the seed for the ML step, altough running times will then be much longer

Here's how the new tool will be in its future 2025.1.0 release:

![GITHUB_update-20250103-2](https://github.com/user-attachments/assets/cd844433-c210-4791-94b4-ec6acc236cf3)

These changes will be tested on a smaller artificial network –based on Huntington's Disease, but with limited number of nodes compared to the original 61K nodes in the network in my [2024 MSc thesis](https://www.dropbox.com/scl/fi/6vzgfld7riqb19hm5wj6u/DRUG_REPURPOSING_thesis_Premium.pdf?rlkey=y9xm7zuxm4q035byvhe496m9k&dl=0)– to determine if different methods lead to a higher prediction rate. Other than the two ML appraoches described above, I will also test whether other methods of generating semantically valid negative samples lead to better results.

### Intermediate Updates
So far, **versions 2025.0.1** and **2025.0.2** add the [OMIM conversion tool](https://github.com/NCMBianchi/OMIM-converter) I developed to turn Monarch Initiative's IDs for diseases, genes and phenotypes (*i.e* respectively `MONDO:`, `HGNC:` and `PH:`) into OMIM IDs and viceversa. For the purpose of this tool, I only implemented the OMIM-to-MONDO conversion: if a user inputs the OMIM for their disease of interest, the conversion tool will fill the input box with the desired Monarch Initiative ID that is accepted by Monarch V3 API as queries. Such conversion box is collapsable. Moreover, a custom .js now checks that manual inputs (*i.e.* `Disease URI` and `ML seed`) are in the correct format before launching the run.

A backup version of the single Docker app is available in the **back-2025.0.2 branch**.

**Versions 2025.0.3**, **2025.0.4** and **2025.0.5** introduced the Docker Swarm directory architecture for distributed services, that will allow for scalability and queues in runs. A `drugapp-launcher.sh` script was be added. It allows to either build the docker images and launch the containers (*i.e.* ```druapp-launcher.sh --build```) or remove them (*i.e.* ```drugapp-launcher.sh --remove```). Such command has to be executable beforehand with ```chmod +x drugapp-launcher.sh```.

The next step to complete the transition to Docker Swarm distributed services require the update of specific services' Python files, as well as adding *\_wait_and_fetch_results()* function to fully pass data between services within the ```servicerun.py``` file.

## PIPELINE: Jupyter Notebook (Niccolò Bianchi)
The entire pipeline can be launched within the Jupyter Notebook. A full run dataset is available on [dropbox](https://www.dropbox.com/scl/fi/prvqajjau227741z5ve91/data.7z?rlkey=qjumdz9r93y0yv6mhc21a7bir&st=jqsfgypj&dl=0).

My major contributions to the project are:
1. updating scripts that interacted with database APIs (<i>i.e.</i> Monarch Initiative, DGIdb)
2. fixing issue with defunct Python packages (<i>e.g.</i> pyRDF2vec)
3. changing the code for the training, testing ang prediction steps in order to make the prediction on drugs that are not only those already part of the network of nodes associated to the disease of interest, in ```networkmodel.py```
4. adding a step that computes negative triples, in ```negsamples.py```
5. updating the training, testing and prediction steps accordingly, in ```networkmodel.py```

![Huntington disease_2024-07-31_full_network_plot_with_NS](https://github.com/user-attachments/assets/eed9bbfc-5168-4a6b-a19c-ba7b2e91c25e)

A link to my dissertation will be posted here in the short future, where I discussed the results and future research.

## BROWSER APP: Flask and Docker (Niccolò Bianchi)

All the scripts for the app to run are complete (<i>i.e.</i> ```routes.py```, ```monarch.py```, ```DGIdb.py```, ```drugsimilarity.py```, ```networkmodel.py```, ```negsamples.py``` and ```networkmodel.py```), as well as the .html, .css and .js files for the web app to render.

I have to solve a few minor issues with the input formatting and the exceptions –which now cause the app to halt at start when running 'random' as seed. Afterwards, I want to implement an advancement bar that would appear when launching the app, as shown below:

![DRUGAPP_home-4-EDIT](https://github.com/user-attachments/assets/74c5c9e7-d4f9-4fc2-8c9c-483c72c075d7)

## How to use

Launch the Docker Client on your device, then ```docker build [--no-cache] --tag drugapp .``` to build the Docker image based on a provided `Dockerfile`, then launch the Docker container with ```docker run -d -p 5001:5000 --name drugapp drugapp```. You can choose any door other than `5001:5000`, based on your requirements.

## Aknowledgements
### LUMC and LIACS (NL)
[Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198)

### Politecnico di Milano (IT)
[Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ)

### SWAT4HCLS
Presentation by Carmen Reep on the original work this fork repository is based on: [Automated drug repurposing workflow for rare diseases](https://youtu.be/RsfUrRhZAso?si=Og1z1RdPaukpPIbP)
