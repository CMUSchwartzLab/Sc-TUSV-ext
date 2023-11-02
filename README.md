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
  
