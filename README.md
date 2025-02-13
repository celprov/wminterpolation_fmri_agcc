# Project Description
Agenesis of the corpus callosum (AgCC) is characterized by a developmental absence of the corpus callosum, the largest white matter bundle connecting the two cerebral hemispheres. 
A big question regarding this condition is how these patients can show functional connectivity comparable to a healthy brain, despite the lack of a key white matter structure. 
In this project, we aimed to identify which white matter structures compensate for the absence of the corpus callosum in order to sustain healthy brain function. 
We used a recently introduced graph signal recovery framework to interpolate signals from the gray matter into the white matter using the anatomical information as a guide. 
This allowed us to capture how temporally correlated gray matter areas are mediated by the underlying structure. 
The co-activation patterns related to the default mode network obtained presented not only different patterns of activation and deactivation but also different white matter structure visibility and different specificity to control or AgCC subjects.
Furthermore, we identified the most important white matter bundles and gray matter regions that discriminate between AgCC and a healthy brain using support vector machine recursive feature elimination (SVM-RFE).

# Repository Structure 
- All the code necessary for the preprocessing is situated in the folder `Preprocessing-Hub`

- `connectivity-pipeline` contains the code necessary to perform the inpainting of the grey matter signal
	into the white matter

- `kmeans` contains the code in order to perform the clustering and obtain the PCC-seed CAPs
- `WholeBrainCaps` contains the code to obtain whole-brain CAPs

- All the analysis performed on the CAPs are contained in `analysis_after_clustering`
	including the feature selection using SVM-RFE and the classification using SVM

- The folder `data` contains the preprocessed data, the inpainted volumes and the CAPs

- Each subfolder contains a text file that explain the pipeline to follow and in which order to run the different code files.