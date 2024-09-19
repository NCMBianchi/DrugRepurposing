# DrugRepurposing

![README_image-1_bis](https://github.com/NCMBianchi/DrugRepurposing/assets/111352723/7db18469-0998-42f7-8bb2-23ee2af35f3c)

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

## Aknowledgements
### LUMC and LIACS (NL)
[Carmen Reep](https://www.researchgate.net/profile/Carmen-Reep), [Núria Queralt-Rosinach](https://www.researchgate.net/scientific-contributions/Nuria-Queralt-Rosinach-2198951627), [Katy J. Wolstencroft](https://www.researchgate.net/profile/Katy-Wolstencroft), [Armel Lefebvre](https://0-scholar-google-com.brum.beds.ac.uk/citations?user=O363fEMAAAAJ&hl=en), [Marco Spruit](https://scholar.google.com/citations?user=GFvyyeAAAAAJ), [Mireia Palou Tort](https://nl.linkedin.com/in/mireia-palou-tort-295909198)

### Politecnico di Milano (IT)
[Rosario M. Piro](https://scholar.google.com/citations?user=HuNyLrcAAAAJ)

### SWAT4HCLS
Presentation by Carmen Reep on the original work this fork repository is based on: [Automated drug repurposing workflow for rare diseases](https://youtu.be/RsfUrRhZAso?si=Og1z1RdPaukpPIbP)