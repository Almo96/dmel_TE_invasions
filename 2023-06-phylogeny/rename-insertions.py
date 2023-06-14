import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Rename sequences in phylip file.")
parser.add_argument("phylip", help="Path to phylip file")
parser.add_argument("msafasta", help="Path to MSA.fasta file")
parser.add_argument("output", help="Path to the output file")
parser.add_argument("te", help="TE name")
args = parser.parse_args()

'''
Takes as input a phylip file and renames the TE insertions in a subsequential order by taxa.

es. Dmel_1
Dmel_2
Dsim_1
Dmel_3
Dsim_2...

Also write a file with the new names to keep track of the changes (.insertion.names.txt)
'''
def rename(path, msa, out, family):
    taxa = {}
    with open(path, 'r') as input_file, open(msa, 'r') as fasta:
        with open(out, 'w') as output_file:
            names = []
            fasta_names = []
            for line in input_file:
                if line.startswith(str(family)):
                    name = line.strip().split(' ')[0]
                    te_name = name.strip().split('_')[0]
                    species = name.strip().split('_')[1]
                    if species not in taxa:
                        taxa[species] = 1
                    else:
                        taxa[species] += 1
                    new_name = str(te_name) + "_" + str(species) + "_" + str(taxa[species])
                    new_line = new_name + " " + " ".join(line.strip().split(' ')[1:]) + "\n"
                    names.append(new_name)
                    output_file.write(new_line)
                else:
                    output_file.write(line)
            for line in fasta:
                if line.startswith(">"):
                    fasta_name = line[1:]
                    fasta_names.append(fasta_name)
                else:
                    continue
        with open(out+"_insertion-names.txt", 'w') as insertions:
            for item1, item2 in zip(names, fasta_names):
                insertions.write(f"{item1}\t{item2}")

rename(args.phylip, args.msafasta, args.output, args.te)
