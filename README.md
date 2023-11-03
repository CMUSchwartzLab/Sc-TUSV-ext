# Sc-TUSV-ext
Sc-TUSV-ext: Tumor evolutionary tree reconstruction from single-cell DNA-seq data of single nucleotide variants (SNV), copy number alterations (CNA) and structural variants (SV).

## Installation
To run Sc-TUSV-ext, two separate environments are required. 
- The first environment is for MEDICC2. You can create the medicc2 environment by running the following commands: <br>
    ```
    conda create -n medicc_env
    
    conda activate medicc_env
    
    conda install -c bioconda -c conda-forge medicc2
    ```
  Or, please follow the instructions as described in [https://bitbucket.org/schwarzlab/medicc2/src/master/](https://bitbucket.org/schwarzlab/medicc2/src/master/).

- The second environment is for running the program with python 2.7. This will need the following commands: <br>
    ```
    conda create -n sctusvext python=2.7
    conda activate sctusvext
    conda config --add channels conda-forge
    conda config --add channels bioconda
    ```
  Then, you will need the following packages in the  `sctusvext` environment.
      - `numpy` <br>
      - `pandas` <br>
      - `ete2` <br>
      - `gurobipy` <br>
      - `graphviz` <br>
      - `biopython` <br>
      - `PyVCF`
- We use the Gurobi optimzer for Sc-TUSV-ext. To acquire Gurobi license, you can sign up as an academic user in the Gurobi website - [https://www.gurobi.com/downloads/end-user-license-agreement-academic/](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). 
  
## Inputs and Outputs
### Input
The input folder should contain the processed variant called scDNAseq files in VCF format. An example can be found in the `example/sample/` folder. 

### Outputs
- Medicc2 output folder: Output of the MEDICC2 method.
- clusters.tsv: The clone assignments for the single-cells according to the MEDICC2 distances. This file is saved inside the MEDICC2 output folder.
- T.dot: Output tree with the `clone assignments` in the nodes and  `phylogenetic cost/number of SNV and SV mapped` in the branches.
- Z.tsv: Clone assignment matrix.
- C.tsv: The estimated copy numbers of the clones.


## Instructions for running
We suggest to run and change the `sctusvext.sh` file as per user's need. The command to run the file:
```
./sctusvext.sh input_folder output_folder number_of_leaves
```
For example, if you wish to have 3 leaves in the tree, i.e. 5 clones, the following command should be run - 
```
./sctusvext.sh example/sample/ example/sample_output 3
```

