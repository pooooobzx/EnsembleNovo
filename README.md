
Ensemble by running:

python3 Ensemble.py

file structures:

scores_*: confidence scores output by each model, with * replaced for each model (forward is clip, backward is reverse and non-auto is prime)

denovo_*_new: denovo pepetide resutls by each model

clipmodel.py: model to encode peptide and spectrum to calculate the pairwise cosine similarity scores

Model_NAT: folder containing model to run non-autoregressive sequencing. 


venn.py: get statistics of how much peptide are correctly predicted for each model for drawing the venn diagram.

