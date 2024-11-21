conda init bash
#------------------
# Command to run Sc-TUSV-ext
# ./sctusvext.sh ../SimGenomeSc/sims/test_code/patient1/sample ../SimGenomeSc/sims/test_code/patient1/outputs n_leaves percentage
#------------------
echo "Command: ./sctusvext.sh input_folder output_folder number_of_leaves"
input_folder=$1
output_folder=$2
n_leaf=$3
percentage=$4
clusters_file="$output_folder/pred_kmeans_clusters.tsv"

rm -r $output_folder
mkdir $output_folder

# runs sctusvext 
source /home/nbristy/miniconda3/bin/activate tusv2 # you might need to change 'tusv2' with your Sc-TUSV-ext environment's name.
echo $CONDA_DEFAULT_ENV
module load gurobi902  # Update with your gurobi module. If it is auto uploaded, you do not need this line.
# you can change the parameters according to your need.
python sc-tusv-ext.py -i $input_folder -o "$output_folder" -n $n_leaf -c 10 -t 2 -r 2 -p 8 -m 1000 -col -b -sv_ub 80 -C 120 -x $percentage
conda /home/nbristy/miniconda3/bin/deactivate


