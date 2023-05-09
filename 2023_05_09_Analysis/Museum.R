library(tidyverse)
theme_set(theme_bw())

#tail -q -n +2 *.csv > GDL_ols_museum.csv

#Early
#
#SRR23876564.fastq.gz.sort.bam.OPUS.pdf
#SRR23876563.fastq.gz.sort.bam.OPUS.pdf
#SRR23876566.fastq.gz.sort.bam.OPUS.pdf
#SRR23876582.fastq.gz.sort.bam.OPUS.pdf
#SRR23876583.fastq.gz.sort.bam.OPUS.pdf
#SRR23876584.fastq.gz.sort.bam.OPUS.pdf

#
#MID
#SRR23876562.fastq.gz.sort.bam.OPUS.pdf
#SRR23876569.fastq.gz.sort.bam.OPUS.pdf

#
#Late
#SRR23876565.fastq.gz.sort.bam.OPUS.pdf



#The dating of samples is not fully documented, the Early, Mid and Late 1800 sample have been dated 1800, 1850 and 1875 for the analysis

df_0 <- read.csv("/Volumes/INTENSO/merged/CSV/GDL_ols_museum.csv", header = FALSE)
names(df_0) <- c("run_accession","TE", "All_reads", "HQ_reads")

#Removing single copy genes used for the normalization
df_1 <- df_0 %>%
  filter(!(TE %in% c("Dmel_tj", "Dmel_rpl32", "Dmel_rhi")))


df_metadata <- read.table("/Users/ascarpa/Downloads/dataset-metadata.txt", sep = "\t", header = TRUE)

df_2 <- inner_join(df_1, df_metadata, by = "run_accession") 

#Removing non Dmel TEs, they only map low quality reads
df <- df_2 %>%
  filter(!grepl("^DM|^DV|^DNTOMRETA", TE))


df_museum <- subset(df, study == "museum")
df_museum_1800_vs_1933 <-subset(df_museum, estimated_year == 1800 | estimated_year == 1933)

museum_1800 <- df_museum %>%
  filter(estimated_year == 1800) %>%
  group_by(TE) %>%
  summarize(avg_TE_all = mean(All_reads),
            avg_TE_HQ = mean(HQ_reads))
names(museum_1800) <- c("TE","avg_1800_all", "avg_1800_HQ")


museum_1933 <- df_museum %>%
  filter(estimated_year == 1933) %>%
  group_by(TE) %>%
  summarize(avg_TE = mean(All_reads),
            avg_TE_HQ = mean(HQ_reads))
names(museum_1933) <- c("TE","avg_1933_all", "avg_1933_HQ")

museum_difference <- inner_join(museum_1800, museum_1933, by ="TE")
museum_difference$diff_all <- (museum_difference$avg_1933_all-museum_difference$avg_1800_all)
museum_difference$diff_norm_all <- (museum_difference$avg_1933_all-museum_difference$avg_1800_all)/museum_difference$avg_1933_all
museum_difference$diff_HQ <- (museum_difference$avg_1933_HQ-museum_difference$avg_1800_HQ)
museum_difference$diff_norm_HQ <- (museum_difference$avg_1933_HQ-museum_difference$avg_1800_HQ)/museum_difference$avg_1933_HQ


sorted_museum_difference_all <- museum_difference %>%
  arrange(desc(diff_norm_all)) %>%
  head(10)
print(sorted_museum_difference_all)

sorted_museum_difference_HQ <- museum_difference %>%
  arrange(desc(diff_norm_HQ)) %>%
  head(10)
print(sorted_museum_difference_HQ)


gg_different_TEs_all <- ggplot(sorted_museum_difference_all, aes(x = TE, y = diff_all))+ 
  geom_bar(stat = "identity", color='skyblue',fill='steelblue')+
  labs(title = "All reads - Difference between 1800 and 1933 in the museum samples")+
  xlab("TE") +
  ylab("Difference in copy number")+
  coord_cartesian(xlim = c(1, 10))
plot(gg_different_TEs_all)


gg_different_TEs_HQ <- ggplot(sorted_museum_difference_HQ, aes(x=TE, y=diff_HQ))+ 
  geom_bar(stat = "identity", color='skyblue',fill='steelblue')+
  labs(title = "HQ reads - Difference between 1800 and 1933 in the museum samples")+
  xlab("TE") +
  ylab("Difference in copy number")+
  coord_cartesian(xlim = c(1, 10))
plot(gg_different_TEs_HQ)




df_museum_412 <- subset(df, study == "museum" & TE == "412")
gg_412_all <- ggplot(df_museum_412, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "412 all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_412_all)
gg_412_HQ <- ggplot(df_museum, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "412 HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_412_HQ)


df_museum_BLOOD <- subset(df, study == "museum" & TE == "BLOOD")
gg_BLOOD_all <- ggplot(df_museum_BLOOD, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "BLOOD all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_BLOOD_all)
gg_BLOOD_HQ <- ggplot(df_museum_BLOOD, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "Blood HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_BLOOD_HQ)


df_museum_OPUS <- subset(df, study == "museum" & TE == "OPUS")
gg_OPUS_all <- ggplot(df_museum_OPUS, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "OPUS all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_all)
gg_OPUS_HQ <- ggplot(df_museum_OPUS, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "OPUS HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_HQ)


df_museum_TIRANT <- subset(df, study == "museum" & TE == "TIRANT")
gg_OPUS_all <- ggplot(df_museum_TIRANT, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "TIRANT all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_all)
gg_OPUS_HQ <- ggplot(df_museum_TIRANT, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "TIRANT HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_HQ)


df_museum_GYPSY2 <- subset(df, study == "museum" & TE == "GYPSY2")
gg_OPUS_all <- ggplot(df_museum_GYPSY2, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "GYPSY2 all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_all)
gg_OPUS_HQ <- ggplot(df_museum_GYPSY2, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "TIRANT HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_OPUS_HQ)



df_museum_KEPLER <- subset(df, study == "museum" & TE == "KEPLER")
gg_KEPLER_all <- ggplot(df_museum_KEPLER, aes(x=as.factor(estimated_year), y=All_reads)) + 
  geom_boxplot() +
  labs(title = "KEPLER all reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_KEPLER_all)
gg_KEPLER_HQ <- ggplot(df_museum_KEPLER, aes(x=as.factor(estimated_year), y=HQ_reads)) + 
  geom_boxplot() +
  labs(title = "KEPLER HQ reads")+
  xlab("Estimated year") +
  ylab("Copy number")
plot(gg_KEPLER_HQ)


