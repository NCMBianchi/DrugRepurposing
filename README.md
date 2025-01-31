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

Both the Jupyter Notebook from version 2024.1.5 and the single Docker app are available in the **v2025.0.2-legacy** release.

**Versions 2025.0.3** through **2025.0.7** introduce the Docker Swarm directory architecture for distributed services, that will allow for scalability and queues in runs. A `drugapp-launcher.sh` script was be added. It allows to either build the docker images and launch the containers (*i.e.* ```druapp-launcher.sh --build```) or remove them (*i.e.* ```drugapp-launcher.sh --remove```). Such command has to be executable beforehand with ```chmod +x drugapp-launcher.sh```. Additional static pages `run_status.html` and `current_runs.html` have been added to allow for UI navigation of concurrent runs in the app. The `status_update.js` script handles runtime display.

The next step is to fully test the functionality of distributed services with a small (real or artificial) network ahead of adding the new alternatives for the algorithm for negative samples generation based on centrality measures, and for the ML prediction based on graph neural networks (_i.e._ GNN, via the `PyG` package).

## PIPELINE: Jupyter Notebook (Niccolò Bianchi)
The entire pipeline can be launched within the Jupyter Notebook. A full run dataset is available on [dropbox](https://www.dropbox.com/scl/fi/prvqajjau227741z5ve91/data.7z?rlkey=qjumdz9r93y0yv6mhc21a7bir&st=jqsfgypj&dl=0).

My major contributions to the project are:
1. updating scripts that interacted with database APIs (<i>i.e.</i> Monarch Initiative, DGIdb)
2. fixing issue with defunct Python packages (<i>e.g.</i> pyRDF2vec)
3. changing the code for the training, testing ang prediction steps in order to make the prediction on drugs that are not only those already part of the network of nodes associated to the disease of interest, in ```networkmodel.py```
4. adding a step that computes negative triples, in ```negsamples.py```
5. updating the training, testing and prediction steps accordingly, in ```networkmodel.py```

![Huntington disease_2024-07-31_full_network_plot_with_NS](https://github.com/user-attachments/assets/eed9bbfc-5168-4a6b-a19c-ba7b2e91c25e)

## Aknowledgements
### LUMC and LIACS (NL)
[Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198)

### Politecnico di Milano (IT)
[Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ)

### SWAT4HCLS
Presentation by Carmen Reep on the original work this fork repository is based on: [Automated drug repurposing workflow for rare diseases](https://youtu.be/RsfUrRhZAso?si=Og1z1RdPaukpPIbP)
