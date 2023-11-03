import pandas as pd
import os
import sys

folder = sys.argv[1] #"../DG1134/"
output_folder = sys.argv[2]
cells = sorted(os.listdir(folder))#[1:]
print('Number of cells', len(cells))
fw = open(output_folder+'/input.tsv','w')
print('created file ',output_folder,'/input.tsv')
fw.write("sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n")
print(cells)
for cell in cells:
    cellname = cell.split('.')[0]
    # Split the data by lines
    vcf_data = open(folder+'/'+cell,'r')
    lines = [line.rstrip() for line in vcf_data.readlines()]

    # Parse metadata, header, and variant data
    metadata = [line for line in lines if line.startswith("##")]
    header = [line for line in lines if line.startswith("#") and not line.startswith("##")]
    variant_data = [line.split('\t') for line in lines if not line.startswith("#")]
    
    for variant in variant_data:
        if 'cnv' in variant[2]:
            fw.write(cellname+'\t'+'chrom'+variant[0]+'\t'+variant[1]+'\t'+variant[7].split(';')[0].split('=')[1]+'\t'+str(round(float(variant[9].split(':')[1].split(',')[0])))+'\t'+str(round(float(variant[9].split(':')[1].split(',')[1])))+'\n')

fw.close()