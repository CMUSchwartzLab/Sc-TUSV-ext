conda init bash
#------------------
# Command to run Sc-TUSV-ext
# ./sctusvext.sh ../SimGenomeSc/sims/test_code/patient1/sample ../SimGenomeSc/sims/test_code/patient1/outputs 
#------------------
echo "Command: ./sctusvext.sh input_folder output_folder number_of_clones"
input_folder=$1
output_folder=$2
n_clone=$3
#tusv_time=$3 # 1000/1500 is a good number for a small dataset with 100s of mutations.


rm -r $output_folder
mkdir $output_folder
# create medicc2 input
mkdir "$output_folder/medicc2_input"
python vcf_to_medicc2.py "$input_folder/" "$output_folder/medicc2_input" # creates a medicc2 format input from the vcf files.


# run medicc2
source /home/nbristy/miniconda3/bin/activate medicc_env # you might need to change 'medicc_env' with your medicc2 environment's name.
echo $CONDA_DEFAULT_ENV
medicc2 "$output_folder/medicc2_input/input.tsv" "$output_folder/medicc2_output" --plot 'bars' # runs medicc2 with default settings.
conda /home/nbristy/miniconda3/bin/deactivate

# run clustering
rm "$output_folder/medicc2_output/clusters.tsv" "$output_folder/medicc2_output/cluster.png"
python cluster.py "$output_folder/medicc2_output" "$input_folder/" $n_clone # clusters single cells using hierarchical clustering.

# runs sctusvext 
source /home/nbristy/miniconda3/bin/activate tusv2 # you might need to change 'tusv2' with your Sc-TUSV-ext environment's name.
echo $CONDA_DEFAULT_ENV
module load gurobi902  # Update with your gurobi module. If it is auto uploaded, you do not need this line.
rm -r "$output_folder/output_sctusvext"
mkdir "$output_folder/output_sctusvext"
# you can change the parameters according to your need.
python sc-tusv-ext.py -i $input_folder -o "$output_folder/output_sctusvext" -c2cl "$output_folder/medicc2_output/clusters.tsv" -n $n_clone -c 10 -t 2 -r 2 -p 8 -m 500 -col -b -sv_ub 80 -C 120 
conda /home/nbristy/miniconda3/bin/deactivate 


