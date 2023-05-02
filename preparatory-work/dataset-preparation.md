Dataset preparation
================

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

The dataset we are going to use for our analysis is composed by:

- Global diversity lines (**gdl**) - 86 samples
- Old lab strains from the Tirant paper (**tirant**, ols) - 14 samples
- Museum specimens from 1800 (**museum**) - 25 samples

In total, we have **125 samples**.

The **metadata** file preparation is described in **metadata.md**. Here,
I want to describe the bioinformatic steps that were followed to prepare
the dataset.

## Data download

The files were downloaded in *FASTQ.gz* format from the ENA website
using **wget**. The accession numbers are:

- PRJNA268111 (**gdl**)
- PRJNA634847 (**tirant**)
- PRJNA945389 (**museum**)

Example of the wget command, with the file *ftp.txt* containing all the
ftp addresses for the files, one in each line.

    wget -i ftp.txt -P /Volumes/INTENSO/data-meseum_invasions/old-lab-strains/

After downloading the files, we checked the **md5 sum** of all the files
to assess their completeness. Example md5 for all the fastq.gz files in
the folder:

    md5 *fastq.gz > md5-downloaded.txt

Both the ftp and the expected md5 can be found in the ENA website report
for the project (download in tsv format after selecting the columns of
interest).

## Trimming

To be consistent with further analysis across all the samples, we
decided to **trim** the reads of all the files to a length of **100
bp**. This threshold was decided after careful analysis of the read
length distribution in all the samples.

Command for reads trimming:

## Mapping to reference library (bam files)

After trimming, the next step was to map the fastq files to the TEs
reference library (*/ref/teseqs-3scg-dmel.fasta*).

Command for mapping:

    bwa bwasw