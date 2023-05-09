'''
Author: RICCARDO PIANEZZA

mpileup2PCA.py

This script takes as input the path to a folder containing multiple mpileup files (es. file.RAW output from DeviaTE).
Each file has the bases count for each position in a sequence to which the FASTQ file has been previously mapped to (es. a TE).
It merges all the files, converts the raw counts into allele frequencies and select only the major allele frequency for each position.

The output (please also provide the path to the output folder!) can be directly used to perfmorm a PCA: one row per sample, one column per SNP.
'''

import argparse
import os
import statistics

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Combine VCF files into a single TSV file.")
parser.add_argument("folder", help="Path to the folder containing the VCF files")
parser.add_argument("output", help="Path to the output folder")
parser.add_argument("--min-freq", type=float, help="Minimum frequency to be considered an allele")
parser.add_argument("--min-count", type=int, help="Minimum number of samples that overcome the min freq to be considered a SNP")
args = parser.parse_args()

file_dir = args.folder
files = os.listdir(file_dir)

'''
Takes as input the path to a mpileup file (.RAW files from DeviaTE), remove all the unnecessary
information for this script. Can be runned in a loop for all the files in the folder.
'''
def clean_mpileup(path):
    filename = os.path.splitext(os.path.basename(path))[0].split('.')[0]
    with open(file_dir+path, 'r') as f_in:
        with open(args.output+"/"+filename+".cleaned", 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                columns = line.strip().split(' ')
                new_columns = columns[0:1]+[filename]+columns[2:3]+columns[4:8]
                new_line = '\t'.join(new_columns) + '\n'
                f_out.write(new_line)


'''
Takes as input the path to the folder where a mpileup file is and the path to the file itself.
Return the header of the merged file which will be then filled in the next function.
'''
def sync_header(folder, path):
    columns = []
    with open(folder+"/"+path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            colname = line[0] + '_' + line[2]
            columns.append(colname)

    column_names = '\t'.join(columns)
    return "ID\t" + column_names + '\n'


'''
Takes as input the path to the folder where a mpileup file is and the path to the file itself.
Can be looped for all the .RAW files in the folder. Appends file's genetic info to the merged file.
'''
def cleaned2sync(folder, path):
    rowname = ""
    columns = []
    data = {}
    with open(folder+"/"+path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if rowname == "":
                rowname = line[1]
            colname = line[0] + '_' + line[2]
            columns.append(colname)
            values = ':'.join([line[3], line[4], line[5], line[6]])
            if colname not in data:
                data[colname] = []
            data[colname].append(values)

    output = []
    output.append(rowname)
    row_values = '\t'.join(['\t'.join(data[col]) for col in columns])

    return path + "\t" + row_values + '\n'


'''
Takes as input the path to the merged file.
Converts the raw read counts into allele frequencies and writes a new file (.af). If there are 5 or less reads
supporting a particular position, the frequencies of that position are assigned to -1, which represents
a missing value.
'''
def allele_freqs(merged_sync):
    with open(merged_sync, 'r') as infile, open(args.output+".af", 'w') as outfile:
        header = infile.readline()
        outfile.write(header)
        for line in infile:
            cells = line.strip().split('\t')
            outfile.write(cells[0]+ '\t')
            for cell in cells[1:]:
                values = cell.strip().split(':')
                total = sum([int(v) for v in values])
                if total <= 10:
                    frequencies = [-1,-1,-1,-1]
                else:
                    frequencies = [round(int(v)/total, 2) for v in values]
                outfile.write(':'.join([str(f) for f in frequencies]) + '\t')
            outfile.write('\n')

'''
Parse the .af file by column (position) and investigate if it's a SNP or not. Return a list of SNP indexes.
SNP calling strategy: if a position has at least 2 alleles and each allele is present in at least n individuals,
then it's a SNP. The minimum frequency to be considered an allele is defined by the args.min-freq,
the minimum number of samples that has to carry the alleles is defined by args.min-count.
'''
#def detect_SNP(af_sync): OLD VERSION based on frequencies
'''
    with open(af_sync, 'r') as infile:
        header = infile.readline()
        num = header.strip().split('\t')
        snps = [0]
        for i in range(1, len(num)):
            infile.seek(0)
            infile.readline()
            maxfreqs = [0, 0, 0, 0]
            alleles = 0
            for line in infile:
                cell = line.strip().split('\t')[i]
                values = cell.strip().split(':')
                for x,value in enumerate(values):
                    if float(value) > float(maxfreqs[x]):
                        maxfreqs[x] = value
            maxfreqs = [float(x) for x in maxfreqs]
            for freq in maxfreqs:
                if freq > 0.4:
                    alleles += 1
            if alleles >= 2:
                snps.append(i)
            i+=1
    return(snps)

'''
#def detect_SNP(af_sync): OLD VERSION based on variance
'''
    with open(af_sync, 'r') as infile:
        header = infile.readline()
        num = header.strip().split('\t')
        snps = [0]
        for i in range(1, len(num)):
            infile.seek(0)
            infile.readline()
            freqs = [[],[],[],[]]
            for line in infile:
                cell = line.strip().split('\t')[i]
                values = cell.strip().split(':')
                for x,value in enumerate(values):
                    freqs[x].append(float(value))
            variances = [statistics.variance(allele) for allele in freqs]
            maxi = max(variances)
            if maxi > 0.1388:
                print(str(maxi)+": "+str(num[i]))
                snps.append(i)
            i+=1
            
    return(snps)
'''
#def detect_SNP(af_sync): OLD VERSION based on counts
'''
    with open(af_sync, 'r') as infile:
        header = infile.readline()
        num = header.strip().split('\t')
        snps = [0]
        for i in range(1, len(num)):
            infile.seek(0)
            infile.readline()
            alleles = 0
            counts = [0,0,0,0]
            freqs = [[],[],[],[]]
            for line in infile:
                cell = line.strip().split('\t')[i]
                values = cell.strip().split(':')
                for x,value in enumerate(values):
                    freqs[x].append(float(value))
            for a, base in enumerate(freqs):
                for sample in base:
                    if sample > 0.5:
                        counts[a]+=1
            for count in counts:
                if count > 10 and count < 80:
                    alleles += 1
            if alleles >= 1:
                snps.append(i)
            i+=1
            
    return(snps)
'''
def detect_SNP(af_sync, min_freq, min_count):

    with open(af_sync, 'r') as infile:
        header = infile.readline()
        num = header.strip().split('\t')
        snps = [0]
        for i in range(1, len(num)):
            infile.seek(0)
            infile.readline()
            maxfreqs = [0, 0, 0, 0]
            freqs = [[],[],[],[]]
            counts = [0, 0, 0, 0]
            alleles = 0
            ids = 0
            for line in infile:
                cell = line.strip().split('\t')[i]
                values = cell.strip().split(':')
                for x,value in enumerate(values):
                    freqs[x].append(float(value))
                    if float(value) > float(maxfreqs[x]):
                        maxfreqs[x] = value
            maxfreqs = [float(x) for x in maxfreqs]
            for a, base in enumerate(freqs):
                for sample in base:
                    if sample > min_freq:
                        counts[a]+=1
            for count in counts:
                if count >= min_count:
                    ids += 1
            for freq in maxfreqs:
                if freq > min_freq:
                    alleles += 1
            if (alleles >= 2) and (ids >= 2):
                snps.append(i)
            i+=1
    return(snps)




'''
Takes as input the path to the merged.af file and a vector of indexes representing the allele to be selected
for each SNP. It filters the .af file keeping only the values of the indexes in the vector (es. only major alleles).
'''
def selectSNP(af_sync, SNPvector):
    with open(af_sync, 'r') as infile, open(args.output+".SNPs", 'w') as outfile:
        for line in infile:
            new_line = ""
            cells = line.strip().split('\t')
            for i, cell in enumerate(cells):
                if i in SNPvector:
                    new_line = new_line + cell + "\t"
            new_line = new_line.rstrip('\t')
            new_line = new_line + "\n"
            outfile.write(new_line)
                
                
        

'''
Takes as input the path to the merged.af file.
For every position, it detects the major allele index (0,1,2,3 for A,T,C,G) and store the info in a vector.
'''
#def minor_allele_vector(af_sync):
'''
    with open(af_sync, 'r') as infile:
        header = infile.readline()
        length = len(header.strip().split('\t'))
        maf_indexes = []
        for i in range(length-1):
            sums = [0, 0, 0, 0]
            next(infile)
            for line in infile:
                cells = line.strip().split('\t')[1:]
                values = cells[i].split(':')
                sums[0] += float(values[0])
                sums[1] += float(values[1])
                sums[2] += float(values[2])
                sums[3] += float(values[3])
            #maf_index = sums.index(max(sums))
            max_num = max(sums)
            temp = sums
            temp.remove(max_num)
            minor_allele = max(temp)
            for f in sums:
                if f==minor_allele:
                    maf_index = sums.index(f)
            maf_indexes.append(maf_index)
            infile.seek(0)
    return maf_indexes
'''
def major_allele_vector(af_sync):
    with open(af_sync, 'r') as infile:
        header = infile.readline()
        length = len(header.strip().split('\t'))
        maf_indexes = []
        for i in range(length-1):
            sums = [0, 0, 0, 0]
            next(infile)
            for line in infile:
                cells = line.strip().split('\t')[1:]
                values = cells[i].split(':')
                sums[0] += float(values[0])
                sums[1] += float(values[1])
                sums[2] += float(values[2])
                sums[3] += float(values[3])
            maf_index = sums.index(max(sums))
            maf_indexes.append(maf_index)
            infile.seek(0)
    return maf_indexes


'''
Takes as input the path to the merged.af file and the index vector created by the function "minor_allele_vector".
For every position, it select the major allele in all the samples. Writes a new file (.PCAable).
'''
def major_allele_selection(af_sync, vector):
    with open(af_sync, 'r') as infile, open(args.output+".PCAable", 'w') as outfile:
        header = infile.readline()
        outfile.write(header)
        for line in infile:
            i = 0
            cells = line.strip().split('\t')
            outfile.write(cells[0])
            for cell in cells[1:]:
                maf_index = vector[i]
                values = cell.split(':')
                ma = values[maf_index]
                outfile.write('\t' + str(ma))
                i += 1
            outfile.write('\n')


for file in files:                                      # Iterates over all the .RAW files
    clean_mpileup(file)                                 # Cleans all the .RAW files from unnecessary infos

cleaned_dir = args.output
files = os.listdir(cleaned_dir)                         # Stores the cleaned file names in a list

with open(args.output+"merged.sync", 'w') as o:         
        o.write(sync_header(cleaned_dir, files[0]))     # Writes the header of the merged file

for file in files:                                      # Iterates over all the "cleaned" files
    with open(args.output+"merged.sync", 'a') as o:
        o.write(cleaned2sync(cleaned_dir, file))        # Append the genetic info from each file to the merged file

allele_freqs(args.output+"merged.sync")                 # Convert the raw read counts in the merged files into allele freqs


positions = detect_SNP(args.output+".af", args.min_freq, args.min_count)
selectSNP(args.output+".af", positions)


indexes = major_allele_vector(args.output+".SNPs")        # Create the vector of major allele indexes for all the positions
major_allele_selection(args.output+".SNPs", indexes)      # Produce the final output with only the major allele for each position

