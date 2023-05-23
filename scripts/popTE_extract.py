import os
import csv

input_directory = '/Volumes/HD_Almo/stat_POP2'
output_file = '/Volumes/HD_Almo/Pop2_reads.csv'

# Define the CSV header
header = ['run_accession', 'te', 'reads', 'reads_in_file', 'reads_mapped']


# Open the output CSV file
with open(output_file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(header)

    # Iterate over the text files in the input directory
    for file_name in os.listdir(input_directory):
        if file_name.endswith('.txt'):
            run_accession = file_name.split(".")[0]
            file_path = os.path.join(input_directory, file_name)
            with open(file_path, 'r') as input_file:
                lines = input_file.readlines()
                for i in range(len(lines)):
                    line = lines[i].strip()
                    if i == 0:
                        data = line.split('\t')
                        reads_in_file = data[-1]
                    elif i == 1:
                        data = line.split('\t')
                        reads_mapped = data[-1]
                    elif line.startswith('te') or line.startswith('refchr'):
                        data = line.split('\t')
                        te = data[1]
                        reads = data[2]
                        row = [run_accession, te, reads, reads_in_file, reads_mapped]
                        writer.writerow(row)