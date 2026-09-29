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
### 1. Data collection 
Protein structures and features were selected from the PDB. The obtain results were filtered to retain only Pfam-annotated entries, and duplicated were remove base on PDB ID.

### 2. Data Processing
Representative sequences were identified with `MMseqs2`. CIF files for the representative PDB entries were downloaded, and the corresponding chains were extracted using the identified chain IDs. Structural alignments were generated with `mTM-align` and `PDBeFold` and exported as FASTA files. Sequences that did not align consistently with the majority of sequences were excluded. Sequence logos were compared with `Skylign` to select the alignment that best captured the conserved cysteine pattern and was therefore most suitable for HMM training.

### 3. Model Training
To standardize alignment lengths, sequences were truncated to 139 residues for mTM-align. The HMM profile was then constructed from the aligned sequences using `HMMER’s hmmbuild` function. 

### 4. Model Testing
Positive and negative datasets were obtained from `UniProt`: reviewed entries annotated with Pfam PF00014 were labeled as the positive set, whereas reviewed entries without Kunitz domain annotation were used as the negative set. `hmmsearch` was run on both sets using the train model, and the outputs were parsed into DataFrames containing target names, E-values, scores, biases, and labels (1 for positive and 0 for negative). Non-matching negative entries were assigned a nominal E-value of 100. The datasets were shuffled and split into two halves for cross-validation. Performance was assessed using accuracy and the Matthews correlation coefficient (MCC) across E-value thresholds from $10^{−1}$ to $10^{−14}$. The optimal threshold was then applied to the second set to compute the confusion matrix, ROC curve, and AUC.

## Results

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
- [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver): for multiple sequence aligment based on structural aligment
- [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/#table): for reference of biological information on the Kunitz-type domain
- [InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/PF00014/): for confirming domain annotations

#### Python dependencies 

|Library|	Purpose|
|-------|------|
|`pandas`|	Manage tabular data and data frames or read/write CSVs|
|`requests`| for querying the PDB throw the API|
|`numpy`|	Perform numerical and array operations|

### Installation via Conda
For installing the command-line tools a conda enviorment was created and a;; the programs instal 

```bash
# Create and activate a dedicated conda environment
conda create -n hmm_kunitz python=3.10
conda activate hmm_kunitz

# Install core bioinformatics tools
conda install -c bioconda hmmer
conda install -c conda-forge -c bioconda mmseqs2
conda install bioconda::mtm-align
```