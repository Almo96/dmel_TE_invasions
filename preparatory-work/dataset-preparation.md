Dataset preparation
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.4.0     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

The dataset we are going to use for our analysis is composed by:

-   Global diversity lines (**gdl**) - 85 samples
-   Old lab strains from the Tirant paper (**tirant**, ols) - 12 samples
-   Museum specimens from 1800 (**museum**) - 25 samples

In total, we have **122 samples**.

The **metadata** file preparation is described in **metadata.md**. Here,
I want to describe the bioinformatic steps that were followed to prepare
the dataset.

## Data download

The files were downloaded in *FASTQ.gz* format from the ENA website
using **wget**. The accession numbers are:

-   PRJNA268111 (**gdl**)
-   PRJNA634847 (**tirant**)
-   PRJNA945389 (**museum**)

Example of the wget command, with the file *ftp.txt* containing all the
ftp addresses for the files, one in each line.

``` bash
wget -i ftp.txt -P /Volumes/INTENSO/data-meseum_invasions/old-lab-strains/
```

After downloading the files, we checked the **md5 sum** of all the files
to assess their completeness. Example md5 for all the fastq.gz files in
the folder:

``` bash
md5 *fastq.gz > md5-downloaded.txt
```

Both the ftp and the expected md5 can be found in the ENA website report
for the project (download in tsv format after selecting the columns of
interest).

## Trimming

To be consistent with further analysis across all the samples, we
decided to **trim** the reads of all the files to a length of **100
bp**. This threshold was decided after careful analysis of the read
length distribution in all the samples.

Before starting the trimming we measured the length distribution of
every sample with the following command:

``` bash
for FILE in *; do gzip -cd $FILE|awk '{print $1}'|paste - - - -|awk '{print length($2)}'|sort|uniq -c > $FILE.txt; done
```

Then we trimmed all the reads to 100nt and we discarded all the reads
under 100nt:

``` bash
ls *gz | parallel -j 4 'gzip -cd {} | cut -c-100 | awk "NR%4==2 && length(\$0)>=100{print prev; print \$0; getline; print \$0; getline; print \$0} {prev=\$0}" | gzip -c > trimmed/{}'
```

[GNU Parallel](https://doi.org/10.5281/zenodo.7761866) was used to
parallelize the process, making it faster.

After trimming we controlled the read length distribution of the trimmed
files using the previous command, this time with parallelization to
speed it up, the expected output of this control is to have all the
reads at 100nt, with the total number of reads in the trimmed files
equal to the sum of all the reads with more than 100nt in the original
files.

``` bash
ls * | parallel --jobs 10 'gzip -cd {} | awk "{print \$1}" | paste - - - - | awk "{print length(\$2)}" | sort | uniq -c > length_post_trimming/museum_post_trimming_length/{}.txt'
```

The two trimmed fastq files for sample are then merge because deviaTE
analysis only takes single reads as input, to merge the two files we
used the following command:

``` bash
find . -name '*_1.fastq.gz' | parallel -j4 'n={}; n=${n%_1.fastq.gz}; gzip -cd {} ${n}_2.fastq.gz | gzip -c > merged/${n}.fastq.gz'
```

## Mapping to reference library (bam files)

After trimming, the next step was to map the fastq files to the TEs
reference library (*/ref/teseqs-3scg-dmel.fasta*).

``` bash
find . -name "*.fq.gz" | parallel -j 4 'n={/.}; bwa bwasw -t 10 teseqs-3scg-dmel.fasta {} | samtools sort -@ 4 -m 3G - > map/${n}.sort.bam'
```

Then indexing with samtools:

``` bash
for i in *bam;do samtools index $i;done
```

## deviaTE_analysis

We used the deviaTE_analysis in parallel, this generated a file for each
TE in the reference for every sample, we have 180 TEs, 3 single-copy
genes for the normalization and 122 samples.

``` bash
ls -1 *bam | parallel -j 20 'cat /Volumes/INTENSO/data-meseum_invasions/ref/TEnames.txt | while read TE; do deviaTE_analyse --input {} --single_copy_genes Dmel_tj,Dmel_rpl32,Dmel_rhi --library /Volumes/INTENSO/data-meseum_invasions/ref/teseqs-3scg-dmel.fasta --family $TE; done'
```

## Plotting

Finally we used deviTE_plot for generating the plots of the TEs of
interest:

``` bash
for file in *.OPUS; do deviaTE_plot --input "$file" & done
```

We also generated all the TEs for Harwich to screen for the reads
distribution and remove outliers

``` bash
for file in SRR11846560.*; do deviaTE_plot --input "$file" & done
```
