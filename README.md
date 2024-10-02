# Sc-TUSV-ext

Sc-TUSV-ext is an integer linear programming (ILP) based tumor clonal lineage inference method using single nucleotide variants (SNV), copy number alterations (CNA) and structural variants (SV) from single-cell DNA sequencing data.


![image](https://github.com/CMUSchwartzLab/Sc-TUSV-ext/assets/38611728/722ef024-b9d9-44a8-9215-542245c3baac)



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
  Then, you will need the following packages in the  `sctusvext` environment. <br>
      - `numpy` <br>
      - `pandas` <br>
      - `ete2` <br>
      - `gurobipy` <br>
      - `graphviz` <br>
      - `biopython=1.76` <br>
      - `scipy` <br>
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

---
In addition to this, we support the framework to run with a set of different settings accessible through the line 13 of `sctusvext.sh` file. The settings are the following:

```
python sc-tusv-ext.py -i example/sample -o "example/sample_output/output_sctusvext" -c2cl "example/sample_output/medicc2_output/clusters.tsv" -n 3 -c 10 -t 2 -r 2 -p 8 -m 1000 -col -b -sv_ub 80 -C 120 
```
Following inputs are mandatory:
- `-i` : input folder
- `-o` : output folder
- `-c2cl` : clone assignment file.
- `-n` : number of leaves.
- `-c` : maximum copy number allowed for any breakpoint or segment on any node
- `-t` : maximum number of coordinate-descent iterations
- `-r` : number of random initializations of the coordinate-descent algorithm
- `-col` : binary flag whether to collapse the redundant nodes
- `-sv_ub` : the number of subsampled SV breakpoints 
- `-const` : number of total subsampled breakpoints and SNVs
- `-m` : maximum time (seconds) in each coordinate descent iteration

Optional parameters:
- `-b` : binary flag for the regularization parameters to be set automatically
- `-l` : lambda regularization parameter for weighting the phylogenetic cost
- `-p` : number of processors to use (uses all the available cores by default)
- `-s` : number of segments (in addition to those containing breakpoints) that are randomly kept (default keeps all the segments)
  
