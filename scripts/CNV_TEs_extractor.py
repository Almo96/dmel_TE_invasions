# This script will create a .csv for every sample in the input directory.
# Every .csv will contain the copy number of every trapsoson provided for the sample.
# The input files are produced with deviaTE
#Creator Almo

import os

def extract_tes(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
        num_strs = lines[1].split(':')[1].strip().split(' or ')
        num1 = float(num_strs[0].split()[0])
        num2 = float(num_strs[1].split()[0])
        
    return num1, num2

input_folder = '/Volumes/INTENSO/merged/deviaTE_analysis'
output_folder = '/Volumes/INTENSO/merged/CSV'

# create a dictionary to store the data
data = {}

# iterate over the files in the input folder
for filename in os.listdir(input_folder):
    if '.bam' in filename:
        sample = filename.split('.')[0]
        te = filename.split('.')[-1]
        input_file = os.path.join(input_folder, filename)
        
        # extract the data from the file
        all_reads, hq_reads = extract_tes(input_file)

        # add the data to the dictionary
        if sample not in data:
            data[sample] = {}
        data[sample][te] = (all_reads, hq_reads)

# write the data to separate CSV files
for sample in data:
    output_file = os.path.join(output_folder, f'{sample}.csv')
    with open(output_file, 'w') as f:
        f.write('Sample,TE,All_reads,HQ_reads\n')
        for te in data[sample]:
            all_reads, hq_reads = data[sample][te]
            f.write(f'{sample},{te},{all_reads:.2f},{hq_reads:.2f}\n')
