# Modeling protein Kunitz domains with a Hidden Markov model
## Overview 
Kunitz domains are compact cysteine-rich protein domains involved in protease inhibition, toxin activity, coagulation, and disease-related pathways. Because sequence divergence can obscure their detection, accurate computational annotation remains important for comparative protein analysis.

This repository details a comprehensive bioinformatics workflow for the characterization of Kunitz domains. Data was retrieved from the Protein Data Bank (`PDB`) and cleaned. The obtain sequence were cluster avoid data redundancy and biases, and structural alignment with `mTMalign` was carried out to identify conserved motifs. A Hidden Markov Model (HMM) is then built to capture the unique sequence patterns of Kunitz domains. Finally, the model's predictive performance is thoroughly evaluated using both positive and negative datasets derived from UniProt, culminating in the determination of an optimal E-value threshold and assessment of its classification accuracy.

## Objectives

- Construct profile HMMs from structurally aligned Kunitz domain sequences 
    - Retrieve sequence and structure form the PDB
    - Filter and cluster data to avoid biases and redundancy
    - Train an HMM with the multiple sequence alignment for the Kunitz domain
    - Build positive and negative set from UniProt Search
    - Evaluate model performance on positive and negative sets

## Methodology
### [1. Data collection](./01_DataCollection/DataCollection.ipynb)
Protein structures and features were selected from the PDB. The obtain results were filtered to retain only Pfam-annotated entries, and duplicated were remove base on PDB ID.

### [2. Data Processing](./02_DataProcessing/DataProcessing.ipynb)
Representative sequences were identified with `MMseqs2`. CIF files for the representative PDB entries were downloaded, and the corresponding chains were extracted using the identified chain IDs. Structural alignments were generated with `mTM-align` and `PDBeFold` and exported as FASTA files. Sequences that did not align consistently with the majority of sequences were excluded. Sequence logos were compared with `Skylign` to select the alignment that best captured the conserved cysteine pattern and was therefore most suitable for HMM training.

### [3. Model Training](./03_ModelTraning/ModelTraning.ipynb)
To standardize alignment lengths, sequences were truncated to 139 residues for mTM-align. The HMM profile was then constructed from the aligned sequences using `HMMER’s hmmbuild` function. 

### [4. Model Testing](./04_ModelEvaluation/ModelEvaluation.ipynb)
Positive and negative datasets were obtained from `UniProt`: reviewed entries annotated with Pfam PF00014 were labeled as the positive set, whereas reviewed entries without Kunitz domain annotation were used as the negative set. `hmmsearch` was run on both sets using the train model, and the outputs were parsed into DataFrames containing target names, E-values, scores, biases, and labels (1 for positive and 0 for negative). Non-matching negative entries were assigned a nominal E-value of 100. The datasets were shuffled and split into two halves for cross-validation. Performance was assessed using accuracy and the Matthews correlation coefficient (MCC) across E-value thresholds from $10^{−1}$ to $10^{−14}$. The optimal threshold was then applied to the second set to compute the confusion matrix, ROC curve, and AUC.

## Results
A total of 135 unique entries were obtain from PDB annotated as Pfam (PF00014), with resolution lower than 3.5 an of length between 40-80 residue. Then the sequence were cluster in 20 groups with 95% of identity and 85% of coverage and a representative sequence for each was selected. This sequence were structurally align to derive the multiple sequence aligment. Problematic sequence were excluded and sequence logo was revise to ensure the caption of the domain representatives cysteines

The HMM model was trained with the multiple sequence aligment. For the sequence evaluation a set of positive entries ( Pfam and review) and a set of negative entries (Not Pfam and review) were downloaded from UniProt and run with the model, so a E-value, fix to a dimension of 1000 sequences, was computed for each sequence in both set. For the positive set 3 sequence were excluded to avoid overfitting. Both sets were shuffle and merge into one data frame to create the cross validation sets. 

| |Positive| Negative| Total|
|--|-------|--------|------|
|sequences| 395|574229|574624|

The optimal threshold obtain with 2 cross validation set was $10^{-5}$ giving an efficient model for identifying the Kunitz domain in protein sequences. 

||MCC|Accuracy|F1|
|----|---|----|---|
|Score|  0.991|0.999|0.991|


As evidence by the metrics and the AUC, the model was capable of identifying the presence on the Kunitz domain in the protein sequence with high accuracy.

Confusion matrix was also computed 
|                | predicted negative | predicted positive|
|---|----|---|
|actual negative     |         287115       |            0|
|actual positive    |               6        |         191|

6 false negative sequence were identify, caused mainly by by the loop structure of the protein or the atypical cysteine localization in the sequence.


## Environment setup and tools  

To run this pipeline, the following command-line tools, packages and web resources must be available

#### Command-line tools
- [hmmer](https://hmmer.org/documentation.html): for building from a multiple sequence aligment HMM and querying protein sequences against a trained model.
- [mmseqs2](https://github.com/soedinglab/mmseqs2?tab=readme-ov-file): use for sequence clustering to reduce redundancy
- [mtm-align](https://yanglab.qd.sdu.edu.cn/mTM-align/): for multiple sequence aligment based on structural aligment 

#### Web tools
- [PDB](https://www.rcsb.org/): database to query for obtaining the protein sequence, the protein structure and other features for the model training
- [UniProt](https://www.uniprot.org/): database tp query for obtaining the positive and negative sets
- [Skylign](https://skylign.org/): for generate sequence logo visualizations from the multiple sequence alignment
- [NCBI MSA Viewer](https://www.ncbi.nlm.nih.gov/projects/msaviewer/?appname=ncbi_msav&openuploaddialog): foe visualization of multiple sequence aligment
- [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver): for multiple sequence aligment based on structural aligment
- [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/#table): for reference of biological information on the Kunitz-type domain
- [InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/PF00014/): for confirming domain annotations

#### Python dependencies 

|Library|	Purpose|
|-------|------|
|`pandas`|	Manage tabular data and data frames or read/write CSVs|
|`requests`| Queries the PDB throw the API|
|`pathlib`| Manages the files paths|
|`numpy`|	Perform numerical and array operations|
| `sklearn`| Model evaluation and metrics|
| matplotlib| Plots and graphs|

### Installation via Conda
For installing the command-line tools a conda enviorment was created and all the programs install 

```bash
# Create and activate a dedicated conda environment
conda create -n hmm_kunitz python=3.10
conda activate hmm_kunitz

# Install core bioinformatics tools
conda install -c bioconda hmmer
conda install -c conda-forge -c bioconda mmseqs2
conda install bioconda::mtm-align
```
