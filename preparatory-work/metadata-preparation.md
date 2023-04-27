Metadata preparation
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

I created a **unique metadata file** for the three merged datasets we
are going to use:

- Global diversity lines (**gdl**) - 86 samples
- Old lab strains from the Tirant paper (**tirant**, ols) - 12 samples
  (3 files are RNA seq)
- Museum specimens from 1800 (**museum**) - 25 samples

The file contains the following information:

- **Sample identifier** (should match the name of the downloaded files)
- **Sample name**
- **Year** of collection
- **Location** of collection
- **Latitude** of collection
- **Longitude** of collection
- **Source study**

In total, we have **123 samples**. Note that for some samples we have
some missing data (location, year).

``` r
gdl <- read_tsv("/Volumes/Temp1/Dmel-stealthTEs/gdl/metadata-gdl.txt") %>% select(-sample_accession)
```

    ## Rows: 86 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample_accession, run_accession, sample
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ols <- read_tsv("/Volumes/Temp1/Dmel-stealthTEs/old-lab-strains/metadata.txt") %>% filter(!(run_accession %in% c("SRR11846566","SRR11846567","SRR12831808")))
```

    ## Rows: 14 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): run_accession, sample, year, study, location
    ## dbl (2): lat, long
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
museum <- read_tsv("/Volumes/Temp1/Dmel-stealthTEs/museum/metadata-museum.txt") %>% select(-sample_accession)
```

    ## Rows: 25 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample_accession, run_accession, sample
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
museum2 <- read_tsv("/Volumes/Temp1/Dmel-stealthTEs/museum/metadata-museum2.txt") %>% select(-attribution, -sra_accession)
```

    ## Rows: 25 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): year, location, attribution, sra_accession, sample
    ## dbl (2): lat, long
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
museum_final <- museum %>% inner_join(museum2, by="sample") %>% mutate(study="museum")
```

``` r
gdl_final <- gdl %>% 
  mutate(location = case_when(startsWith(sample, "B") ~ "Beijing, China", startsWith(sample, "I") ~ "Ithaca, USA", startsWith(sample, "N") ~ "Netherlands", startsWith(sample, "T") ~ "Tasmania, Australia", startsWith(sample, "Z") ~ "Zimbabwe"), lat = case_when(startsWith(sample, "B") ~ 40, startsWith(sample, "I") ~ 42, startsWith(sample, "N") ~ 52, startsWith(sample, "T") ~ 43, startsWith(sample, "Z") ~ 19), , long = case_when(startsWith(sample, "B") ~ 116, startsWith(sample, "I") ~ 76, startsWith(sample, "N") ~ 5, startsWith(sample, "T") ~ 147, startsWith(sample, "Z") ~ 29), study="gdl")
```

``` r
metadata <- bind_rows(museum_final, ols, gdl_final)
metadata
```

    ## # A tibble: 123 × 7
    ##    run_accession sample year        location                     lat  long study
    ##    <chr>         <chr>  <chr>       <chr>                      <dbl> <dbl> <chr>
    ##  1 SRR23876562   H9     Mid 1800s   Passau, Germany               46    13 muse…
    ##  2 SRR23876563   H10    Early 1800s Lund, Sweden                  56    13 muse…
    ##  3 SRR23876564   H4     Early 1800s Lund, Sweden                  56    13 muse…
    ##  4 SRR23876565   H5     Late 1800s  Zealand, Denmark              55    12 muse…
    ##  5 SRR23876566   H6     Early 1800s Småland, Sweden (but infe…    57    15 muse…
    ##  6 SRR23876567   H7     1933        Lund, Sweden                  56    13 muse…
    ##  7 SRR23876568   H8     1933        Lund, Sweden                  56    13 muse…
    ##  8 SRR23876569   H25    Mid 1800s   Passau, Germany               46    13 muse…
    ##  9 SRR23876570   H24    1933        Lund, Sweden                  56    13 muse…
    ## 10 SRR23876571   H23    1933        Lund, Sweden                  56    13 muse…
    ## # ℹ 113 more rows

``` r
write_tsv(metadata, "/Volumes/Temp1/Dmel-stealthTEs/dataset-metadata")
```
