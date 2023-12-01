library(tidyverse)
library(writexl)
library(ggpubr)
theme_set(theme_bw())

#100nt reads
df_0 <- read.csv("/Volumes/INTENSO/merged/CSV/GDL_ols_museum.csv", header = FALSE)

#50nt reads
#df_0 <- read.csv("/Volumes/HD_Almo/dmel_museum_50/CSV/GDL_ols_50_museum.csv", header = FALSE)
names(df_0) <- c("run_accession","TE", "All_reads", "HQ_reads")


df_metadata <- read.table("/Users/ascarpa/dmel_TE_invasions/dataset-metadata", sep = "\t", header = TRUE)

df_2 <- inner_join(df_0, df_metadata, by = "run_accession") 


df_museum_1800 <- subset(df_2, run_accession == "SRR23876566") #H6
df_Harwich <- subset(df_2, run_accession == "SRR11846560") #Har
df_1933 <- subset(df_2, run_accession == "SRR23876576") #H19

#Investigation asked by the reviewer about the other TEs that show fold enrichmnent,
#TARTC is for telomeric maintainance
#recent with a line from 2004
recent <- subset(df_2, run_accession == "SRR1663548")

#I only keep TEs present in D. Harwich, this would not show lost TEs but I checked and there are no TEs lost
df_filtering <- df_Harwich %>%
  dplyr::filter(All_reads>1) %>%
  dplyr::select("TE")

# fragment library
# df_filtering <- df_filtering %>% filter(!str_detect(TE, "CM0347"))

#The following TEs are removed because we are using all the reads, but this can cause misleading reasults if the reads
#map on a single part of the TE, like streteches of As in the LINEs, the removed TEs are picked after a screening of
#the deviaTE_plots.

df_filtering <- df_filtering %>%
  dplyr::filter(!(TE %in% c("Dmel_rhi", "Dmel_rpl32", "Dmel_tj", "DV26847", "KEPLER", "Q", "TARTVIR", "TARTYAK")))

df_museum_1800 <- inner_join(df_museum_1800, df_filtering, by = "TE")
df_Harwich <- inner_join(df_Harwich, df_filtering, by = "TE")
df_1933 <- inner_join(df_1933, df_filtering, by = "TE")
recent <- inner_join(recent, df_filtering, by = "TE")

df_museum_1800 <- df_museum_1800 %>%
  dplyr::select(TE, All_reads)
df_Harwich <- df_Harwich %>%
  dplyr::select(TE, All_reads)
df_1933 <- df_1933 %>%
  dplyr::select(TE, All_reads)
recent <- recent %>%
  dplyr::select(TE, All_reads)

df_m_Har <- inner_join(df_museum_1800, df_Harwich, by = "TE")
names(df_m_Har) <- c("TE", "copies_museum", "copies_Har")
df_m1800_1933 <- inner_join(df_museum_1800, df_1933, by = "TE")
names(df_m1800_1933) <- c("TE", "copies_museum", "copies_1933")

df_m_recent <- inner_join(df_museum_1800, recent, by = "TE")
names(df_m_recent) <- c("TE", "copies_museum", "copies_recent")
df_Lau_recent <- inner_join(df_1933, recent, by = "TE")
names(df_Lau_recent) <- c("TE", "copies_1933", "copies_recent")
df_Har_recent <- inner_join(df_Harwich, recent, by = "TE")
names(df_Har_recent) <- c("TE", "copies_Har", "copies_recent")



df_m_Har$fold_enrichment <- as.numeric(df_m_Har$copies_Har)/as.numeric(df_m_Har$copies_museum)
df_m1800_1933$fold_enrichment <- as.numeric(df_m1800_1933$copies_1933)/as.numeric(df_m1800_1933$copies_museum)

df_m_Har$log_fold_enrichment <- log2(df_m_Har$fold_enrichment)
df_m1800_1933$log_fold_enrichment <- log2(df_m1800_1933$fold_enrichment)


red_c <- "#c1272d"
blue_c <- "#155BE8"
purple_c <- "#AE37FF"
TEs_c <- c("412" = red_c, "BLOOD" = red_c, "OPUS" = red_c,"TIRANT" = blue_c, "DMHFL1" = blue_c, "DMIFACA" = blue_c, "PPI251" = blue_c)
#TEs_c <- c("TARTC" = purple_c,"412", "DM23420" = purple_c,"412" = red_c, "BLOOD" = red_c, "OPUS" = red_c,"TIRANT" = blue_c, "DMHFL1" = blue_c, "DMIFACA" = blue_c, "PPI251" = blue_c)


gg_1A_log_Har <- ggplot(data = df_m_Har, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  coord_cartesian(ylim = c(-2, 7))+
  scale_y_continuous(breaks = c(-2,0,2,4,6))+
  scale_fill_manual(values = TEs_c) +
  annotate("text", x=12, y=4, label= "Harwich (1967)", size=12) +
  annotate("text", x=12, y=-2, label= "18SL6 (~1800)", size=12) +
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))

plot(gg_1A_log_Har)

# Prevent the removal of the P-element, this would remove one column and the combined figure wouldn't be alligned
df_m1800_1933 <- df_m1800_1933 %>% mutate_all(~replace(., is.nan(.), 0))

gg_2A_log_1933_1800 <- ggplot(data = df_m1800_1933, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TE") +
  ylab("Fold Enrichment") +
  coord_cartesian(ylim = c(-2, 7))+
  scale_y_continuous(breaks = c(-2,0,2,4,6))+
  scale_fill_manual(values = TEs_c) +
  annotate("text", x=12, y=4, label= "19SL19 (1933)", size=12) +
  annotate("text", x=12, y=-2, label= "18SL6 (~1800)", size=12) +
  theme(legend.position = "none", axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

plot(gg_2A_log_1933_1800)

#Figure 1C
df_museum_a <- df_2 %>%
  filter(study == "museum", estimated_year == 1800 | estimated_year == 1933)

df_GDL_museum_8a_TEs <- subset(df_museum_a, TE == "412"|TE== "OPUS"|TE == "BLOOD"|TE =="CIRC"|TE =="INVADER4")

df_GDL_museum_8a_TEs$study <- factor(df_GDL_museum_8a_TEs$estimated_year, levels = c(1800, 1933))
df_GDL_museum_8a_TEs$TE <- factor(df_GDL_museum_8a_TEs$TE, levels = c("412", "BLOOD","OPUS","CIRC","INVADER4"))

fig_1C <- ggplot(df_GDL_museum_8a_TEs, aes(x = as.factor(estimated_year), y = All_reads)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("1800", "1933")), map_signif_level = TRUE, textsize = 6) +
  facet_wrap(~ TE, nrow = 1, labeller = labeller(TE = 
                                                   c("412" = "412",
                                                     "BLOOD" = "Blood",
                                                     "OPUS" = "Opus",
                                                     "CIRC" = "Circe",
                                                     "INVADER4" = "Invader-4")))+
  labs(x = "year", y = "copy number")+
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) 

plot(fig_1C)
ggsave("/Volumes/INTENSO/deviaTE_plots/Paper/Figure_1/Fig_1C.png", fig_1C, width = 10, height = 3, dpi = 300)



# p-values Bonferroni correction for multiple testing
df_museum_b <- df_2 %>%
  filter(study == "museum", estimated_year == 1800 | estimated_year == 1850 | estimated_year == 1933)
df_museum_b <- df_museum_b %>%
  mutate(year = ifelse(estimated_year %in% c(1880, 1850), 1800, estimated_year))

df_museum_b <- df_museum_b %>%
  dplyr::select(run_accession, TE, All_reads, year)

df_t <- df_museum_b %>%
  pivot_wider(names_from = TE, values_from = All_reads)

t.test(BLOOD ~ as.factor(year), data = df_t)
t.test(INVADER4 ~ as.factor(year), data = df_t)

TE_vector <- unique(df_museum_a$TE)
results_df <- data.frame(TE = character(length(TE_vector)), p_value = numeric(length(TE_vector)))

# Loop through TE_vector and perform t-test for each variable
for (i in seq_along(TE_vector)) {
  # Extract the variable name from TE_vector
  var <- TE_vector[i]
  
  # Convert the variable to numeric (assuming it's not already numeric)
  df_t[[var]] <- as.numeric(df_t[[var]])
  
  # Perform the t-test
  t_test_result <- t.test(df_t[[var]] ~ df_t$year)
  
  # Store results in the data frame
  results_df[i, ] <- list(TE = var, p_value = t_test_result$p.value)
}

# Display the results data frame
print(results_df)


results_df <- na.omit(results_df)
results_df$p_value_Bonferroni <- p.adjust(results_df$p_value, method = "bonferroni")

ggplot(data = results_df, aes(x = TE, y = -log10(p_value_Bonferroni), color = TE)) +
  geom_point() +
  labs(x = "TE",
       y = "-log10(p_value_Bonferroni)") +
  geom_hline(yintercept = -log10(0.0001), linetype = "dashed", color = "red")+
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))+
  scale_color_manual(values = c("BLOOD" = "red", "OPUS" = "red", "412" = "red"))
  

# divide males and females
unique(df_t$run_accession)
females <- c("SRR23876564", "SRR23876566", "SRR23876567", "SRR23876570", "SRR23876571", "SRR23876572",
            "SRR23876574", "SRR23876575", "SRR23876576", "SRR23876577", "SRR23876583", "SRR23876585",
            "SRR23876562", "SRR23876569")



males <- c("SRR23876563", "SRR23876568", "SRR23876573", "SRR23876578", "SRR23876579", "SRR23876580",
          "SRR23876581", "SRR23876582", "SRR23876584", "SRR23876586")


df_t_males <- df_t %>%
  filter(run_accession %in% males)
df_t_females <- df_t %>%
  filter(run_accession %in% females)

results_df_males <- data.frame(TE = character(length(TE_vector)), p_value = numeric(length(TE_vector)))
results_df_females <- data.frame(TE = character(length(TE_vector)), p_value = numeric(length(TE_vector)))

for (i in seq_along(TE_vector)) {
  var <- TE_vector[i]
  df_t_males[[var]] <- as.numeric(df_t_males[[var]])
  t_test_result_m <- t.test(df_t_males[[var]] ~ df_t_males$year)
  results_df_males[i, ] <- list(TE = var, p_value = t_test_result_m$p.value)
}

print(results_df_males)


for (i in seq_along(TE_vector)) {
  var <- TE_vector[i]
  df_t_females[[var]] <- as.numeric(df_t_females[[var]])
  t_test_result_f <- t.test(df_t_females[[var]] ~ df_t_females$year)
  results_df_females[i, ] <- list(TE = var, p_value = t_test_result_f$p.value)
}

print(results_df_females)

results_df_males <- na.omit(results_df_males)
results_df_females <- na.omit(results_df_females)

results_df_males$p_value_Bonferroni <- p.adjust(results_df_males$p_value, method = "bonferroni")
results_df_females$p_value_Bonferroni <- p.adjust(results_df_females$p_value, method = "bonferroni")


males <- ggplot(data = results_df_males, aes(x = TE, y = -log10(p_value_Bonferroni), color = TE)) +
  geom_point() +
  ggtitle("Males")+
  labs(x = "TE",
       y = "-log10(p_value_Bonferroni)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))+
  scale_color_manual(values = c("BLOOD" = "red", "OPUS" = "red", "412" = "red"))

plot(males)

females <- ggplot(data = results_df_females, aes(x = TE, y = -log10(p_value_Bonferroni), color = TE)) +
  geom_point() +
  ggtitle("Females")+
  labs(x = "TE",
       y = "-log10(p_value_Bonferroni)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  theme(legend.position = "none", axis.text.x = element_text(size = 4, angle = 45, hjust = 1))+
  scale_color_manual(values = c("BLOOD" = "red", "OPUS" = "red", "412" = "red"))

plot(females)


#Figure check recent, museum and Lousanne_S


df_m_recent$fold_enrichment <- df_m_recent$copies_recent/df_m_recent$copies_museum
df_Lau_recent$fold_enrichment <- df_Lau_recent$copies_recent/df_Lau_recent$copies_1933
df_Har_recent$fold_enrichment <- df_Har_recent$copies_recent/df_Har_recent$copies_Har

df_m_recent$log_fold_enrichment <- log2(df_m_recent$fold_enrichment)
df_Lau_recent$log_fold_enrichment <- log2(df_Lau_recent$fold_enrichment)
df_Har_recent$log_fold_enrichment <- log2(df_Har_recent$fold_enrichment)


m_recent <- ggplot(data = df_m_recent, aes(x = TE, y = as.numeric(log_fold_enrichment), fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("fold Enrichment") +
  coord_cartesian(ylim = c(-1.5,3))+
  scale_y_continuous(breaks = c(-2,0,2))+
  annotate("text", x=9, y=3, label= "recent (2004)", size=8) +
  annotate("text", x=10.5, y=-1.5, label= "museum (~1800)", size=8) +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

plot(m_recent)

Lau_recent <- ggplot(data = df_Lau_recent, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("fold Enrichment") +
  coord_cartesian(ylim = c(-1.5,3))+
  scale_y_continuous(breaks = c(-2,0,2))+
  annotate("text", x=9, y=3, label= "recent (2004)", size=8) +
  annotate("text", x=10.5, y=-1.5, label= "Lausanne-S (1938)", size=8) +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

plot(Lau_recent)

Har_recent <- ggplot(data = df_Har_recent, aes(x = TE, y = log_fold_enrichment, fill = TE)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("fold Enrichment") +
  coord_cartesian(ylim = c(-1.5,3))+
  scale_y_continuous(breaks = c(-2,0,2))+
  annotate("text", x=9, y=3, label= "recent (2004)", size=8) +
  annotate("text", x=10.5, y=-1.5, label= "Harwich (1967)", size=8) +
  scale_fill_manual(values = TEs_c) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

plot(Har_recent)

#Figure check Harwich and Lousanne_S

df_Har_Lau <- inner_join(df_Harwich, df_1933, by = "TE")
names(df_Har_Lau) <- c("TE", "copies_Har", "copies_1933")
df_Har_Lau$fold_enrichment <- df_Har_Lau$copies_Har/df_Har_Lau$copies_1933
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
