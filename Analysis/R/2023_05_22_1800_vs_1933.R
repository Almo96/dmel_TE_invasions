#The same analysis performed on reads trimmed at 50nt instead of 100nt

library(tidyverse)
library(writexl)
library(ggpubr)
theme_set(theme_bw())

df_0 <- read.csv("/Volumes/HD_Almo/dmel_museum_50/CSV/GDL_ols_50_museum.csv", header = FALSE)
names(df_0) <- c("run_accession","TE", "All_reads", "HQ_reads")


df_metadata <- read.table("/Users/ascarpa/Downloads/dataset-metadata.txt", sep = "\t", header = TRUE)

df_2 <- inner_join(df_0, df_metadata, by = "run_accession") 


df_museum_1800 <- subset(df_2, run_accession == "SRR23876582")
df_Harwich <- subset(df_2, run_accession == "SRR11846560")
df_Lausanne_S <- subset(df_2, run_accession == "SRR11846564")


#I only keep TEs present in D. Harwich, this would not show lost TEs but I checked and there are no TEs lost
df_filtering <- df_Harwich %>%
  filter(All_reads>1) %>%
  select("TE")

#The following TEs are removed because we are using all the reads, but this can cause misleading reasults if the reads
#map on a single part of the TE, like streteches of As in the LINEs, the removed TEs are picked after a screening of
#the deviaTE_plots.

df_filtering <- df_filtering %>%
  filter(!(TE %in% c("DV26847", "KEPLER", "Q", "TARTVIR", "TARTYAK")))

df_museum_1800 <- inner_join(df_museum_1800, df_filtering, by = "TE")
df_Harwich <- inner_join(df_Harwich, df_filtering, by = "TE")
df_Lausanne_S <- inner_join(df_Lausanne_S, df_filtering, by = "TE")

df_museum_1800 <- df_museum_1800 %>%
  select(TE, All_reads)
df_Harwich <- df_Harwich %>%
  select(TE, All_reads)
df_Lausanne_S <- df_Lausanne_S %>%
  select(TE, All_reads)


df_m_Har <- inner_join(df_museum_1800, df_Harwich, by = "TE")
names(df_m_Har) <- c("TE", "copies_museum", "copies_Har")
df_m_Lau <- inner_join(df_museum_1800, df_Lausanne_S, by = "TE")
names(df_m_Lau) <- c("TE", "copies_museum", "copies_Lau")


df_m_Har$fold_enrichment <- df_m_Har$copies_Har/df_m_Har$copies_museum
df_m_Lau$fold_enrichment <- df_m_Lau$copies_Lau/df_m_Lau$copies_museum

df_m_Har$log_fold_enrichment <- log2(df_m_Har$fold_enrichment)
df_m_Lau$log_fold_enrichment <- log2(df_m_Lau$fold_enrichment)


red_c <- "#c1272d"
blue_c <- "#155BE8"
TEs_c <- c("412" = red_c, "BLOOD" = red_c, "OPUS" = red_c,"TIRANT" = blue_c, "DMHFL1" = blue_c, "DMIFACA" = blue_c, "PPI251" = blue_c)


#Figure 1A
ggplot(data = df_m_Har, aes(x = TE, y = fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(data = df_m_Lau, aes(x = TE, y = fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


#Figure 1A log
gg_1A_log_Har <- ggplot(data = df_m_Har, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  coord_cartesian(ylim = c(-2, 7))+
  scale_y_continuous(breaks = c(-2,0,2,4,6))+
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))

plot(gg_1A_log_Har)


gg_1A_log_Har_s <- gg_1A_log_Har + ylab("")+xlab("")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


gg_1A_log_Lau <- ggplot(data = df_m_Lau, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  coord_cartesian(ylim = c(-2, 7))+
  scale_y_continuous(breaks = c(-2,0,2,4,6))+
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))

plot(gg_1A_log_Lau)

#Combined
gg_1A <- ggarrange(gg_1A_log_Har_s, gg_1A_log_Lau, ncol = 1, nrow = 2)
plot(gg_1A)


#Figure check Harwich and Lousanne_S

df_Har_Lau <- inner_join(df_Harwich, df_Lausanne_S, by = "TE")
names(df_Har_Lau) <- c("TE", "copies_Har", "copies_Lau")
df_Har_Lau$fold_enrichment <- df_Har_Lau$copies_Har/df_Har_Lau$copies_Lau
df_Har_Lau$log_fold_enrichment <- log2(df_Har_Lau$fold_enrichment)

ggplot(data = df_Har_Lau, aes(x = TE, y = fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
Har_Lau <- ggplot(data = df_Har_Lau, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("fold Enrichment") +
  coord_cartesian(ylim = c(-1.5,3))+
  scale_y_continuous(breaks = c(-2,0,2))+
  annotate("text", x=9, y=3, label= "Harwich (1967)", size=8) +
  annotate("text", x=10.5, y=-1.5, label= "Lausanne-S (1938)", size=8) +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

plot(Har_Lau)


#Figure 1C GDL test
df_museum <- df_2 %>%
  filter(study == "museum", estimated_year == 1800 | estimated_year == 1850)

#Check, expected output 8, 6 early and 2 mid samples
count <- df_museum %>% 
  filter(TE == "OPUS") %>% 
  summarize(count = n())
# Output the count
print(count)
df_GDL <- df_2 %>%
  filter(study == "gdl")

#Figure 1C
df_museum_a <- df_2 %>%
  filter(study == "museum", estimated_year == 1800 | estimated_year == 1933)

df_GDL_museum_8a_TEs <- subset(df_museum_a, TE == "412"|TE== "OPUS"|TE == "BLOOD"|TE =="CIRC"|TE =="INVADER4")

df_GDL_museum_8a_TEs$study <- factor(df_GDL_museum_8a_TEs$estimated_year, levels = c(1800, 1933))
df_GDL_museum_8a_TEs$TE <- factor(df_GDL_museum_8a_TEs$TE, levels = c("412", "OPUS","BLOOD","CIRC","INVADER4"))

fig_1C <- ggplot(df_GDL_museum_8a_TEs, aes(x = as.factor(estimated_year), y = All_reads)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("1800", "1933")), map_signif_level = TRUE, textsize = 2) +
  facet_wrap(~ TE, nrow = 1) +
  labs(x = "year", y = "copy number")

plot(fig_1C)

#Sanity check on single copy genes
df_GDL_museum_8b_TEs <- subset(df_museum_a, TE=="Dmel_rhi"|TE=="Dmel_rpl32"|TE=="Dmel_tj")

df_GDL_museum_8b_TEs$study <- factor(df_GDL_museum_8b_TEs$estimated_year, levels = c(1800, 1933))

ggplot(df_GDL_museum_8b_TEs, aes(x = as.factor(estimated_year), y = All_reads)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("1800", "1933")), map_signif_level = TRUE) +
  facet_wrap(~ TE, nrow = 1) +
  labs(x = "Study", y = "Copy number")


#Figure 2

#1800 early: H10, H13
#1800 mid: H9, H25
#1800 late: H5
#1933: H21, H24
#Crimea
#Lausanne-S
#Urbana-S
#Berlin-K
#Dmel68
#Harwich


TE_2 <- c("412", "OPUS", "BLOOD", "TIRANT", "DMIFACA", "DMHFL1", "PPI251")
samples_2 <- c("H10", "H13", "H9", "H25", "H5", "H21", "H24","Crimea", "Lausanne-S", "Urbana-S", "Berlin-K", "Dmel68", "Harwich")


df_fig2 <- df_2 %>%
  select(TE, sample, All_reads)%>%
  filter(TE %in% TE_2)%>%
  filter(sample %in% samples_2)

#Convert to wide format
table_fig2 <- df_fig2 %>%
  pivot_wider(names_from = sample, values_from = All_reads)

table_fig2_ordered <- table_fig2 %>%
  arrange(match(TE, TE_2)) %>%
  select(TE, all_of(samples_2))

table_fig2_ordered

# Save as XLSX
write_xlsx(table_fig2_ordered, "/Volumes/INTENSO/deviaTE_plots/Paper/Figure_2/table_fig2.xlsx")
