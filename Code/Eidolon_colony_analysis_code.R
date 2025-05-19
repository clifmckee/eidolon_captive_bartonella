### Required code for Eidolon captive colony analyses
### R version 4.5.0 "How About a Twenty-Six"
### Works 19 May 2025

# Clear environment
rm(list=ls())

# Install ggtree package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree", force = TRUE)

# Load necessary packages
library(reshape2)
library(viridis)
library(permute)
library(lattice)
library(vegan)
library(cowplot)
library(segmented)
library(nlme)
library(mgcv)
library(binom)
library(stats)
library(MuMIn)
library(ggthemes)
library(tidyverse)
library(ggtree)
library(phytools)
library(magrittr)

# Print out session info
sessionInfo()

# Load stuff for analysis of slopes from gam
tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)

# PHYLOGENETIC ANALYSIS ------------

#### Sup Fig: Bartonella ITS tree ####

# Read in tree data
Bart_ITS_tree_data <- read.csv("Sequencing/Bart_ITS_tree_metadata.csv")

# Read in tree
Bart_ITS_tree <- read.tree("Sequencing/Bartonella ITS model selection/AICc/Bartonella_ITS_sequences_MAFFT-L-INS-i.fasta.contree")

# Root tree at midpoint
root_Bart_ITS_tree <- midpoint.root(Bart_ITS_tree)

# Order tree data by tip label
o_Bart_ITS_tree_data <- Bart_ITS_tree_data[match(root_Bart_ITS_tree$tip.label, Bart_ITS_tree_data$Nickname),]

# Combine tip labels and tree data
df_Bart_ITS_tree_data <- data.frame(label = root_Bart_ITS_tree$tip.label,
                                     species = o_Bart_ITS_tree_data$Bartonella.species,
                                     isolate = o_Bart_ITS_tree_data$Isolate.Strain.Clone,
                                     accession = o_Bart_ITS_tree_data$GenBank.accession.number,
                                     group = factor(o_Bart_ITS_tree_data$Group))

# Plot tree
ggtree(root_Bart_ITS_tree) %<+% df_Bart_ITS_tree_data +
    geom_point(aes(size = as.numeric(label))) +
    geom_tiplab(aes(label = paste(species, isolate, accession), color = group), size = 3) +
    geom_text2(aes(subset = !isTip & as.numeric(label) > 70, label = label),
               size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
    geom_treescale(x = 0, y = -0.3, width = 0.5, fontsize = 3) +
    xlim(0, 1.8) +
    scale_size_binned(range = c(0, 1), breaks = c(50, 80)) +
    scale_color_manual(values = colorblind_pal()(8)[c(1, 7, 2)]) +
    theme_tree() +
    theme(legend.position = "none")
ggsave("Sequencing/Bart_ITS_tree.png", height = 4, width = 5, units = "in", dpi = 300, bg = "white")

#### Sup Fig: Bartonella ftsZ tree ####

# Read in tree data
Bart_ftsZ_tree_data <- read.csv("Sequencing/Bart_ftsZ_tree_metadata.csv")

# Read in tree
Bart_ftsZ_tree <- read.tree("Sequencing/Bartonella ftsZ model selection/AICc/Bartonella_ftsZ_sequences_MAFFT-L-INS-i_GB.fasta.contree")

# Root tree at midpoint
root_Bart_ftsZ_tree <- midpoint.root(Bart_ftsZ_tree)

# Order tree data by tip label
o_Bart_ftsZ_tree_data <- Bart_ftsZ_tree_data[match(root_Bart_ftsZ_tree$tip.label, Bart_ftsZ_tree_data$Nickname),]

# Combine tip labels and tree data
df_Bart_ftsZ_tree_data <- data.frame(label = root_Bart_ftsZ_tree$tip.label,
                                    species = o_Bart_ftsZ_tree_data$Bartonella.species,
                                    isolate = o_Bart_ftsZ_tree_data$Isolate.Strain.Clone,
                                    accession = o_Bart_ftsZ_tree_data$GenBank.accession.number,
                                    group = factor(o_Bart_ftsZ_tree_data$Group))

# Plot tree
ggtree(root_Bart_ftsZ_tree) %<+% df_Bart_ftsZ_tree_data +
  geom_point(aes(size = as.numeric(label))) +
  geom_tiplab(aes(label = paste(species, isolate, accession), color = group), size = 3) +
  geom_text2(aes(subset = !isTip & as.numeric(label) > 80, label = label),
             size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
  geom_treescale(x = 0, y = -0.5, width = 0.3, fontsize = 3) +
  xlim(0, 2.2) +
  scale_size_binned(range = c(0, 1), breaks = c(50, 80)) +
  scale_color_manual(values = colorblind_pal()(8)[c(1, 7, 2)]) +
  theme_tree() +
  theme(legend.position = "none")
ggsave("Sequencing/Bart_ftsZ_tree.png", height = 5, width = 5, units = "in", dpi = 300, bg = "white")

#### Sup Fig: Bartonella gltA tree ####

# Read in tree data
Bart_gltA_tree_data <- read.csv("Sequencing/Bart_gltA_tree_metadata.csv")

# Read in tree
Bart_gltA_tree <- read.tree("Sequencing/Bartonella gltA model selection/AICc/Bartonella_gltA_sequences_MAFFT-L-INS-i_GB.fasta.contree")

# Root tree at midpoint
root_Bart_gltA_tree <- midpoint.root(Bart_gltA_tree)

# Order tree data by tip label
o_Bart_gltA_tree_data <- Bart_gltA_tree_data[match(root_Bart_gltA_tree$tip.label, Bart_gltA_tree_data$Nickname),]

# Combine tip labels and tree data
df_Bart_gltA_tree_data <- data.frame(label = root_Bart_gltA_tree$tip.label,
                                     species = o_Bart_gltA_tree_data$Bartonella.species,
                                     isolate = o_Bart_gltA_tree_data$Isolate.Strain.Clone,
                                     accession = o_Bart_gltA_tree_data$GenBank.accession.number,
                                     group = factor(o_Bart_gltA_tree_data$Group))

# Plot tree
ggtree(root_Bart_gltA_tree) %<+% df_Bart_gltA_tree_data +
  geom_point(aes(size = as.numeric(label))) +
  geom_tiplab(aes(label = paste(species, isolate, accession), color = group), size = 2.5) +
  geom_text2(aes(subset = !isTip & as.numeric(label) > 80, label = label),
             size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
  geom_treescale(x = 0, y = -1, width = 0.1, fontsize = 3) +
  xlim(0, 0.8) +
  scale_size_binned(range = c(0, 1), breaks = c(50, 80)) +
  scale_color_manual(values = colorblind_pal()(8)[c(1, 7, 2)]) +
  theme_tree() +
  theme(legend.position = "none")
ggsave("Sequencing/Bart_gltA_tree.png", height = 9, width = 9, units = "in", dpi = 300, bg = "white")

#### Sup Fig: Bartonella conc ftsZ + gltA tree ####

# Read in tree data
Bart_conc_tree_data <- read.csv("Sequencing/Bart_conc_tree_metadata.csv")

# Read in tree
Bart_conc_tree <- read.tree("Sequencing/Bartonella conc model selection/BIC/Bartonella_conc_sequences.fasta.contree")

# Root tree at midpoint
root_Bart_conc_tree <- midpoint.root(Bart_conc_tree)

# Order tree data by tip label
o_Bart_conc_tree_data <- Bart_conc_tree_data[match(root_Bart_conc_tree$tip.label, Bart_conc_tree_data$Nickname),]

# Combine tip labels and tree data
df_Bart_conc_tree_data <- data.frame(label = root_Bart_conc_tree$tip.label,
                                     species = o_Bart_conc_tree_data$Bartonella.species,
                                     isolate = o_Bart_conc_tree_data$Isolate.Strain.Clone,
                                     group = factor(o_Bart_conc_tree_data$Group))

# Plot tree
ggtree(root_Bart_conc_tree) %<+% df_Bart_conc_tree_data +
  geom_point(aes(size = as.numeric(label))) +
  geom_tiplab(aes(label = paste(species, isolate), color = group), size = 2.5) +
  geom_text2(aes(subset = !isTip & as.numeric(label) > 80, label = label),
             size = 2, hjust = 1.1, vjust = -0.3, color = "grey50") +
  geom_treescale(x = 0, y = -2, width = 0.1, fontsize = 3) +
  xlim(0, 0.9) +
  scale_size_binned(range = c(0, 1), breaks = c(50, 80)) +
  scale_color_manual(values = colorblind_pal()(8)[c(1, 6, 7, 2)]) +
  theme_tree() +
  theme(legend.position = "none")
ggsave("Sequencing/Bart_conc_tree.png", height = 10, width = 10, units = "in", dpi = 300, bg = "white")

# PREVALENCE AND INFECTION LOAD -------------------------------------------

#### Calculating confidence intervals ####

# Captive bats positive markers
one <- c(9, 13, 62, 57, 37, 44, 44, 40, 10, 18, 27, 33, 35, 39)
two <- c(6, 11, 55, 46, 27, 37, 36, 27, 3, 12, 14, 24, 22, 36)
three <- c(4, 8, 45, 38, 23, 27, 27, 23, 2, 5, 9, 17, 15, 30)
four <- c(3, 6, 18, 27, 17, 10, 12, 7, 0, 0, 2, 8, 9, 21)
total <- c(12, 22, 73, 61, 72, 63, 69, 67, 58, 81, 87, 80, 85, 81)
# Wilson score confidence intervals
prev.CI1 <- binom.confint(x=c(one, two, three, four), n=rep(total, 4), method="wilson")
write.csv(prev.CI1, "Results/confint.csv")

# Flies and wild bats positive markers
flies1 <- c(24, 14, 47)
flies2 <- c(23, 12, 42)
flies3 <- c(18, 11, 33)
flies4 <- c(17, 11, 21)
flies.mixed <- c(2, 5, 28)
flies.total <- c(28, 18, 50)
# Wilson score confidence intervals
prev.CI2 <- binom.confint(x=c(flies1, flies2, flies3, flies4, flies.mixed), n=rep(flies.total, 5), method="wilson")
write.csv(prev.CI2, "Results/confint_flies.csv")

# Co-infections
mixed.pos <- c(4, 3, 24, 17, 7, 7, 5, 7, 0, 2, 7, 11, 8, 12)
mixed.total <- c(12, 22, 73, 61, 72, 63, 69, 67, 58, 81, 87, 80, 85, 81)
mixed.CI <- binom.confint(mixed.pos, mixed.total, method="wilson")
write.csv(mixed.CI, "Results/confint_mixed.csv")

#### Figure: Prevalence ####

# Read in captive prevalence data
prev.CI <- read.csv("Data/captive_colony_prevalence.csv", head=T)
# Read in fly and wild bat prevalence data
flies.CI <- read.csv("Data/fly_comparisons_prevalence.csv", head=T)
# Join data
CI_join <- full_join(prev.CI, flies.CI) %>%
  replace_na(list(Group = "captive bats"))

# Line plot of captive prevalence data with points for fly and wild bat prevalence data
ggplot(data = CI_join, aes(x = Day, y = prev.mean, color = Group)) +
  geom_pointrange(aes(x = Day, y = prev.mean, ymin = prev.lower, ymax = prev.upper, color = Group), alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = -0.05,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                   "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Prevalence", breaks = seq(0, 1, 0.25), limits = c(-0.05, 1)) +
  theme_cowplot()
ggsave("Results/Fig_prevalence.png", height =  3.5, width = 8, units = "in", dpi = 300, bg = "white")

# Read in captive RT-PCR data
inf.dat <- read.csv("Data/captive_colony_infection.csv", head=T)
# Read in fly and wild bats RT-PCR data
inf.fly <- read.csv("Data/fly_comparisons_infection.csv", head=T)
# Join data and filter RT-PCR data to remove negative results
inf_join_filter <- full_join(inf.dat, inf.fly) %>%
  replace_na(list(Group = "captive bats")) %>%
  filter(!is.na(Ct), Ct<40)
# Statistics
inf_join_filter_stats <- inf_join_filter %>%
  group_by(Group, Day) %>%
  summarize(n = n(),
            mean = mean(Ct, na.rm = TRUE),
            min = fivenum(Ct, na.rm = TRUE)[1],
            lower = fivenum(Ct, na.rm = TRUE)[2],
            median = fivenum(Ct, na.rm = TRUE)[3],
            upper = fivenum(Ct, na.rm = TRUE)[4],
            max = fivenum(Ct, na.rm = TRUE)[5])

# Scatter plot of captive RT-PCR data with boxes for fly and wild bat prevalence data
ggplot(data = inf_join_filter, aes(x = Day, y = Ct, color = Group)) +
  geom_jitter(size = 0.5, shape = 21, width = 5, height = 0) +
  geom_point(data = inf_join_filter_stats, aes(x = Day, y = median, color = Group), size = 2.5, alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = 45,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                   "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "Day", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_reverse(name = "Ct value", breaks = seq(45, 20, -5), limits = c(45, 20)) +
  theme_cowplot()
ggsave("Results/SuppFig_infection_load.png", height =  3.5, width = 8, units = "in", dpi = 300, bg = "white")

#### Test for effect of flies on prevalence ####

# Effect on prevalence for all age groups after fly reintroduction
prop.test(x=c(32, 0), n=c(53, 53), alternative="greater")

# Difference in effect of flies on prevalence for late vs. early cohort bats
prop.test(x=c(16, 16), n=c(17, 36), alternative="greater")

# Effect on species change for all age groups after fly reintroduction
prop.test(x=c(48, 0), n=c(84, 84), alternative="greater")

# Difference in effect of flies on species change for late vs. early cohort bats
prop.test(x=c(22, 26), n=c(28, 56), alternative="greater")

# All age groups together after fly reintroduction
# Pos/neg test for Jul11 and 17Jan12 neg before
prop.test(x=c(17, 14), n=c(24, 22), alternative="less")

# Pos/neg test for 17Jan12 neg before
prop.test(x=c(17, 15), n=c(27, 26), alternative="less")

# Became pos or changed sp. for Jul11 and 17Jan12
prop.test(x=c(19, 28), n=c(36, 41), alternative="less")

# Became pos or changed sp. for 17Jan12
prop.test(x=c(19, 29), n=c(39, 45), alternative="less")

# Sp in wild bat or fly for Jul11 and 17Jan12
prop.test(c(13, 13), n=c(26, 26), alternative="greater")

# Sp in wild bat or fly for 17Jan12
prop.test(c(13, 13.5), n=c(27, 27), alternative="greater")

# Adult and sexually immature bats
# Pos/neg test for Jul11 and 17Jan12 neg before
prop.test(x=c(5, 10), n=c(11, 18), alternative="less")

# Pos/neg test for 17Jan12 neg before
prop.test(x=c(5, 11), n=c(14, 22), alternative="less")

# Became pos or changed sp. for Jul11 and 17Jan12
prop.test(x=c(6, 19), n=c(18, 31), alternative="less")

# Became pos or changed sp. for 17Jan12
prop.test(x=c(6, 20), n=c(21, 35), alternative="less")

# Sp in wild bat or fly for Jul11 and 17Jan12
prop.test(x=c(5, 6), n=c(12, 12), alternative="greater")

# Sp in wild bat or fly for 17Jan12
prop.test(x=c(5, 6.5), n=c(13, 13), alternative="greater")

# Neonate and juvenile bats
# Pos/neg test for Jul11 and 17Jan12 neg before
prop.test(x=c(12, 4), n=c(13, 4), alternative="less")
fisher.test(x=matrix(c(12, 1, 4, 0), nrow=2), alternative="less")

# Pos/neg test for 17Jan12 neg before
prop.test(x=c(12, 4), n=c(13, 4), alternative="less")
fisher.test(x=matrix(c(12, 1, 4, 0), nrow=2), alternative="less")

# Became pos or changed sp. for Jul11 and 17Jan12
prop.test(x=c(13, 9), n=c(18, 10), alternative="less")
fisher.test(x=matrix(c(13, 5, 9, 1), nrow=2), alternative="less")

# Became pos or changed sp. for 17Jan12
prop.test(x=c(13, 9), n=c(18, 10), alternative="less")
fisher.test(x=matrix(c(13, 5, 9, 1), nrow=2), alternative="less")

# Sp in wild bat or fly for Jul11 and 17Jan12
prop.test(x=c(8, 7), n=c(14, 14), alternative="greater")

# Sp in wild bat or fly for 17Jan12
prop.test(x=c(8, 7), n=c(14, 14), alternative="greater")

#### Regression for March ####

# Sp in colony bat and sampled fly for Mar10
prop.test(x=c(9, 13), n=c(26, 26), alternative="greater")

# March, both sexes
mar.pos <- data.frame(read.csv("Data/Mar10_all.csv", header=T))
cor.test(mar.pos$flies, mar.pos$pos)
mar.mod1 <- glm(formula=pos~sex+age+flies, data=mar.pos, family=binomial(link="logit"), na.action="na.fail")
summary(mar.mod1)
summary(aov(mar.mod1))
mar.mod1.sel <- dredge(mar.mod1); mar.mod1.sel

# March, females only
marF.pos <- data.frame(read.csv("Data/Mar10_F.csv", header=T))
cor.test(marF.pos$flies, marF.pos$pos)
cor.test(as.numeric(marF.pos$preg=="P"), marF.pos$pos)
mar.mod3 <- glm(formula=pos~age+preg+flies, data=marF.pos, family=binomial, na.action="na.fail")
summary(mar.mod3)
summary(aov(mar.mod3))
mar.mod3.sel <- dredge(mar.mod3); mar.mod3.sel

# CHANGES IN BARTONELLA DIVERSITY -----------------------------------------

#### Calculate diversity indices ####

# Read in captive relative count data
rel.count <- as.data.frame(read.csv("Data/captive_colony_relative_counts.csv", head=T))
# Add Day as column
rel.count$Day <- prev.CI$Day
# Reshape relative count data into long format
m.rel.count <- melt(rel.count, id.vars=c("Date", "Day"))
# Read in raw counts
bat.raw.count.data <- as.data.frame(read.csv("Data/captive_colony_counts.csv", head=T))

# Read in fly and wild bat relative count data
fly.counts <- as.data.frame(read.csv("Data/fly_comparisons_relative_counts.csv", head=T))
# Reshape relative count data into long format
m.fly.counts <- melt(fly.counts)
# Read in fly abundance data
fly.raw.count.data <- as.data.frame(read.csv("Data/fly_comparisons_counts.csv", head=T))

# Read in captive relative count data
count.comp <- as.data.frame(read.csv("Data/captive_colony_counts_comparison.csv", head=T))
# Reshape relative count data into long format
m.count.comp <- melt(count.comp)

# Recalculate infection diversity data
inf.div <- rel.count[, c(1, 10)]
inf.div$Species.richness <- rowSums(bat.raw.count.data[, 3:10] != 0)
inf.div$Shannon.number <- exp(diversity(x=bat.raw.count.data[, 3:10], index="shannon"))
inf.div$InvSimpson <- diversity(x=bat.raw.count.data[, 3:10], index="invsimpson")
inf.div$Total.count <- apply(bat.raw.count.data[, 3:10], 1, sum)

# Read in fly and wild bat diversity data
fly.div <- fly.counts[c(1, 3, 4), c(1, 2)]
fly.div$Species.richness <- rowSums(fly.raw.count.data[, 4:11] != 0)
fly.div$Shannon.number <- exp(diversity(x=fly.raw.count.data[, 4:11], index="shannon"))
fly.div$InvSimpson <- diversity(x=fly.raw.count.data[, 4:11], index="invsimpson")
fly.div$Total.count <- apply(fly.raw.count.data[, 4:11], 1, sum)

# Calculate error around species richness and Shannon numbers
# Define function
div.samp <- function(diversity.table, relative.counts, iter, lower, upper){
  div.sum <- array(0, c(dim(diversity.table)[1], 6))
  for(i in 1:dim(diversity.table)[1]){
    dist <- array(0, c(iter, 3))
    for(j in 1:iter){
      samp <- rmultinom(n=1, size=diversity.table$Total.count[i],
                       prob=relative.counts[i, (ncol(relative.counts)-7):ncol(relative.counts)])
      samp.sp <- colSums(samp != 0)
      samp.rel <- samp/sum(samp)
      samp.shan <- exp(diversity(x=samp.rel, index="shannon"))
      samp.simp <- diversity(x=samp.rel, index="invsimpson")
      dist[j, 1:3] <- c(samp.sp, samp.shan, samp.simp)
    }
    div.sum[i, ] <- c(quantile(x=dist[, 1], probs=c(lower, upper)),
                     quantile(x=dist[, 2], probs=c(lower, upper)),
                     quantile(x=dist[, 3], probs=c(lower, upper)))
  }
  colnames(div.sum) <- c("lower.richness", "upper.richness",
                        "lower.Shannon", "upper.Shannon",
                        "lower.Simpson", "upper.Simpson")
  return(div.sum)
}

# Calculate diversity for captive bats
inf.div.samp <- div.samp(diversity.table=inf.div, relative.counts=bat.raw.count.data[, 3:10],
                        iter=1000, lower=.025, upper=.975)
inf.div3 <- cbind(inf.div, inf.div.samp)

# Calculate diversity for flies and wild bats
fly.div.samp <- div.samp(diversity.table=fly.div, relative.counts=fly.raw.count.data[, 4:11],
                        iter=1000, lower=.025, upper=.975)
fly.div3 <- cbind(fly.div, fly.div.samp)

# Join diversity data
div3_join <- full_join(inf.div3, fly.div3) %>%
  replace_na(list(Group = "captive bats")) %>%
  mutate(Day = case_when(Group == "captive flies" ~ 221,
                         Group == "wild flies" ~ 908,
                         Group == "wild bats" ~ 898,
                         TRUE ~ Day))

#### Figure: Bartonella species frequencies ####

# Line plot of captive relative count data by Bartonella species
SuppFigRelCounta <- ggplot(data=m.rel.count, aes(x=Day, y=value, fill=variable)) +
  geom_col(position = "stack") +
  annotate(geom = "text", x = prev.CI$Day, y = -0.05,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_fill_colorblind(name = "Species") +
  scale_x_continuous(name = "Date", breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Relative counts", breaks = seq(0, 1, 0.25)) +
  theme_cowplot()

# Stacked bar plot of captive relative count data by Bartonella species
SuppFigRelCountb1 <- ggplot(data=m.count.comp, aes(x=Group, y=value, fill=variable)) +
  geom_col(position = "stack", width = 0.75) +
  xlab("") +
  scale_y_continuous(name = "Relative counts", breaks = seq(0, 1, 0.25)) +
  scale_fill_colorblind(name="Species") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none") +
  ggtitle("Before J12\nvs. after J12")

# Stacked bar plot of Mar10 bat and fly relative count data
SuppFigRelCountb2 <- ggplot(data=m.fly.counts[which(m.fly.counts$Split=="10-Mar"),], aes(x=Group, y=value, fill=variable)) +
  geom_col(position="stack", width = 0.75) +
  xlab("Group") +
  scale_y_continuous(name = "", breaks = seq(0, 1, 0.25)) +
  scale_fill_colorblind(name="Species") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none") +
  ggtitle("M10 bats\nvs. flies")

# Stacked bar plot of Mar10 bat and fly relative count data
SuppFigRelCountb3 <- ggplot(data=m.fly.counts[which(m.fly.counts$Split=="17-Jan-12"),], aes(x=Group, y=value, fill=variable)) +
  geom_col(position="stack", width = 0.75) +
  scale_fill_colorblind(name="Species") +
  theme_cowplot() +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none") +
  ggtitle("After J12 bats\nvs. flies")

# Combine plots
SuppFigRelCountb <- plot_grid(SuppFigRelCountb1, SuppFigRelCountb2, SuppFigRelCountb3, labels=c("B", "C", "D"), ncol=3, align="hv")
plot_grid(SuppFigRelCounta, SuppFigRelCountb, labels=c("A", ""), nrow=2, rel_heights = c(0.625, 0.375))
ggsave("Results/Fig_rel_count.png", height = 7, width = 7, units = "in", dpi = 300, bg = "white")

#### Sup Fig: Measures of Bartonella diversity ####

# Line plot of captive species richness data with points for fly and wild bat data
FigDiva <- ggplot(data = div3_join, aes(x = Day, y = Species.richness, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = 0,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                   "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Species richness", breaks = seq(0, 8, 2), limits = c(0, 8)) +
  theme_cowplot()

# Line plot of captive Shannon diversity data with points for fly and wild bat data
FigDivb <- ggplot(data = div3_join, aes(x = Day, y = Shannon.number, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_linerange(aes(ymin = lower.Shannon, ymax = upper.Shannon)) +
  annotate(geom = "text", x = prev.CI$Day, y = 0,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                   "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Shannon number", breaks = seq(0, 7, 1), limits = c(0, 7)) +
  theme_cowplot()

# Line plot of captive inverse Simpson diversity data with points for fly and wild bat data
FigDivc <- ggplot(data = div3_join, aes(x = Day, y = InvSimpson, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_linerange(aes(ymin = lower.Simpson, ymax = upper.Simpson)) +
  annotate(geom = "text", x = prev.CI$Day, y = 0,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "Day", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Inv. Simpson index", breaks = seq(0, 6, 1), limits = c(0, 6)) +
  theme_cowplot()

# Read in captive RT-PCR data
inf.dat <- read.csv("Data/captive_colony_infection.csv", head=T)
# Read in fly and wild bats RT-PCR data
inf.fly <- read.csv("Data/fly_comparisons_infection.csv", head=T)
# Join data
inf_join <- full_join(inf.dat, inf.fly) %>%
  replace_na(list(Group = "captive bats"))
# Statistics
inf_join_stats_Numspecies <- inf_join %>%
  group_by(Group, Day) %>%
  summarize(n = n(),
            mean = mean(NumSpecies, na.rm = TRUE),
            min = fivenum(NumSpecies, na.rm = TRUE)[1],
            lower = fivenum(NumSpecies, na.rm = TRUE)[2],
            median = fivenum(NumSpecies, na.rm = TRUE)[3],
            upper = fivenum(NumSpecies, na.rm = TRUE)[4],
            max = fivenum(NumSpecies, na.rm = TRUE)[5])

# Line plot of number of Bartonella species in a sample with points for fly and wild bat data
FigDivd <- ggplot(data = inf_join, aes(x = Day, y = NumSpecies, color = Group)) +
  geom_count(shape = 21) +
  geom_point(data = inf_join_stats_Numspecies, aes(x = Day, y = mean, color = Group), size = 3, alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = 0,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_size_area(name="Counts", breaks = c(1, 20, 40)) +
  scale_x_continuous(name = "Day", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Species in sample", breaks = seq(0, 4, 1), limits = c(0, 4)) +
  theme_cowplot()

# Combine plots
plot_grid(FigDiva, FigDivb, FigDivc, FigDivd, labels=c("A", "B", "C", "D"), nrow = 2, ncol=2, align="hv")
ggsave("Results/SuppFig_diversity.png", height = 7, width = 15, units = "in", dpi = 300, bg = "white")

#### Fig: Bartonella beta diversity ####

## Read in data for each date and calculate beta diversity
# July 2009
Jul09 <- read.csv("Data/Jul09_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- vegdist(Jul09, "jaccard")
dates <- rep("2009-07-28", length(index1))
days <- rep(0, length(vegdist(Jul09, "jaccard")))

# November 2009
Nov09 <- read.csv("Data/Nov09_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Nov09, "jaccard"))
dates <- c(dates, rep("2009-11-05", length(as.vector(vegdist(Nov09, "jaccard")))))
days <- c(days, rep(100, length(vegdist(Nov09, "jaccard"))))

# January 2010
Jan10 <- read.csv("Data/Jan10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Jan10, "jaccard"))
dates <- c(dates, rep("2010-01-28", length(as.vector(vegdist(Jan10, "jaccard")))))
days <- c(days, rep(184, length(vegdist(Jan10, "jaccard"))))

# March 2010
Mar10 <- read.csv("Data/Mar10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Mar10, "jaccard"))
dates <- c(dates, rep("2010-03-01", length(as.vector(vegdist(Mar10, "jaccard")))))
days <- c(days, rep(216, length(vegdist(Mar10, "jaccard"))))

# May 2010
May10 <- read.csv("Data/May10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(May10, "jaccard"))
dates <- c(dates, rep("2010-05-21", length(as.vector(vegdist(May10, "jaccard")))))
days <- c(days, rep(297, length(vegdist(May10, "jaccard"))))

# July 2010
Jul10 <- read.csv("Data/Jul10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Jul10, "jaccard"))
dates <- c(dates, rep("2010-07-14", length(as.vector(vegdist(Jul10, "jaccard")))))
days <- c(days, rep(351, length(vegdist(Jul10, "jaccard"))))

# September 2010
Sep10 <- read.csv("Data/Sep10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Sep10, "jaccard"))
dates <- c(dates, rep("2010-09-23", length(as.vector(vegdist(Sep10, "jaccard")))))
days <- c(days, rep(422, length(vegdist(Sep10, "jaccard"))))

# November 2010
Nov10 <- read.csv("Data/Nov10_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Nov10, "jaccard"))
dates <- c(dates, rep("2010-11-05", length(as.vector(vegdist(Nov10, "jaccard")))))
days <- c(days, rep(465, length(vegdist(Nov10, "jaccard"))))

# March 2011
Mar11 <- read.csv("Data/Mar11_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Mar11, "jaccard"))
dates <- c(dates, rep("2011-03-04", length(as.vector(vegdist(Mar11, "jaccard")))))
days <- c(days, rep(584, length(vegdist(Mar11, "jaccard"))))

# July 2011
Jul11 <- read.csv("Data/Jul11_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Jul11, "jaccard"))
dates <- c(dates, rep("2011-07-13", length(as.vector(vegdist(Jul11, "jaccard")))))
days <- c(days, rep(715, length(vegdist(Jul11, "jaccard"))))

# 17 January 2012
Jan1712 <- read.csv("Data/Jan1712_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Jan1712, "jaccard"))
dates <- c(dates, rep("2012-01-17", length(as.vector(vegdist(Jan1712, "jaccard")))))
days <- c(days, rep(903, length(vegdist(Jan1712, "jaccard"))))

# 31 January 2012
Jan3112 <- read.csv("Data/Jan3112_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Jan3112, "jaccard"))
dates <- c(dates, rep("2012-01-31", length(as.vector(vegdist(Jan3112, "jaccard")))))
days <- c(days, rep(917, length(vegdist(Jan3112, "jaccard"))))

# February 2012
Feb12 <- read.csv("Data/Feb12_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Feb12, "jaccard"))
dates <- c(dates, rep("2012-02-14", length(as.vector(vegdist(Feb12, "jaccard")))))
days <- c(days, rep(931, length(vegdist(Feb12, "jaccard"))))

# March 2012
Mar12 <- read.csv("Data/Mar12_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index1 <- c(index1, vegdist(Mar12, "jaccard"))
dates <- c(dates, rep("2012-03-15", length(as.vector(vegdist(Mar12, "jaccard")))))
days <- c(days, rep(961, length(vegdist(Mar12, "jaccard"))))

# March 2010 flies
flies.Mar10 <- read.csv("Data/Mar10_flies_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index2 <- as.vector(vegdist(flies.Mar10, "jaccard"))
date <- rep("2010-03-11", length(as.vector(vegdist(flies.Mar10, "jaccard"))))
day <- rep(226, length(as.vector(vegdist(flies.Mar10, "jaccard"))))
group <- rep("captive flies", length(as.vector(vegdist(flies.Mar10, "jaccard"))))

# 17 January 2012 wild bats
wild.Jan1712 <- read.csv("Data/Jan1712_wild_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index2 <- c(index2, vegdist(wild.Jan1712, "jaccard"))
date <- c(date, rep("2012-01-12", length(as.vector(vegdist(wild.Jan1712, "jaccard")))))
day <- c(day, rep(898, length(as.vector(vegdist(wild.Jan1712, "jaccard")))))
group <- c(group, rep("wild bats", length(as.vector(vegdist(wild.Jan1712, "jaccard")))))

# 17 January 2012 wild flies
flies.Jan1712 <- read.csv("Data/Jan1712_flies_beta_diversity.csv", head=T) %>% 
  mutate(across(everything(), \(x) as.numeric(x > 0, na.rm = TRUE)))
index2 <- c(index2, vegdist(flies.Jan1712, "jaccard"))
date <- c(date, rep("2012-01-22", length(as.vector(vegdist(flies.Jan1712, "jaccard")))))
day <- c(day, rep(908, length(as.vector(vegdist(flies.Jan1712, "jaccard")))))
group <- c(group, rep("wild flies", length(as.vector(vegdist(flies.Jan1712, "jaccard")))))

# Output captive calculations into data frame
dissim <- as.data.frame(dates)
dissim$days <- days
dissim$index1 <- index1
dissim$groups <- "captive bats"

# Output fly and wild bat calculations into data frame
flysim <- as.data.frame(date)
flysim$day <- day
flysim$index2 <- index2
flysim$group <- group

# Join data
sim_join <- full_join(dissim, flysim,
                      by = c("dates" = "date",
                             "days" = "day",
                             "index1" = "index2",
                             "groups" = "group"))
# Statistics
sim_stats <- sim_join %>%
  group_by(groups, days) %>%
  summarize(n = n(),
            mean = mean(index1, na.rm = TRUE),
            min = fivenum(index1, na.rm = TRUE)[1],
            lower = fivenum(index1, na.rm = TRUE)[2],
            median = fivenum(index1, na.rm = TRUE)[3],
            upper = fivenum(index1, na.rm = TRUE)[4],
            max = fivenum(index1, na.rm = TRUE)[5])

# Line plot for captive beta diversity data with points for fly and wild bat data
BetaDiv <- ggplot(data = sim_join, aes(x = days, y = index1, color = groups)) +
  geom_count(shape = 21) +
  geom_point(data = sim_stats, aes(x = days, y = mean, color = groups), size = 3, alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = -0.15,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_size_area(name="Counts", breaks = c(1, 100, 500)) +
  scale_x_continuous(name = "Day", breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Beta diversity", breaks = seq(0, 1, 0.2)) +
  theme_cowplot()
ggsave("Results/SuppFig_beta_div.png", height = 3.5, width = 7, units = "in", dpi = 300, bg = "white")

#### Fig: Infection durations ####

# Read in table
inf.dur <- read.csv("Data/infection_durations.csv", head=T)

# Test differences between species
kruskal.test(Duration~Species, inf.dur)
pairwise.wilcox.test(inf.dur$Duration, inf.dur$Species, p.adjust.method="bonferroni")

# Summary statistics
dur_stats <- inf.dur %>%
  group_by(Species) %>%
  summarize(n = n(),
            mean = mean(Duration, na.rm = TRUE),
            min = fivenum(Duration, na.rm = TRUE)[1],
            lower = fivenum(Duration, na.rm = TRUE)[2],
            median = fivenum(Duration, na.rm = TRUE)[3],
            upper = fivenum(Duration, na.rm = TRUE)[4],
            max = fivenum(Duration, na.rm = TRUE)[5])
  
# Plot
ggplot(data = inf.dur, aes(x = Species, y = Duration)) +
  geom_count(shape = 21) +
  geom_point(data = dur_stats, aes(x = Species, y = mean), size = 3, alpha = 0.8) +
  annotate("text", x = seq(1, 8, 1), y = -25, label=dur_stats$n[c(1, 2, 3, 4, 5, 8, 6, 7)], size = 4) +
  scale_x_discrete(limits=c("E1", "E2", "E3", "E4", "E5", "Ew", "Eh6", "Eh7")) +
  scale_y_continuous(name="Infection duration", breaks = seq(0, 600, 100)) +
  scale_size_area(name="Counts", breaks = c(1, 5)) +
  theme_cowplot()
ggsave("Results/Fig_durations.png", height = 3.5, width = 7, units = "in", dpi = 300, bg = "white")

# COINFECTION FREQUENCY ANALYSIS ------------------------------------------

#### Functions ####

# Create function to perform multinomial test of co-infection frequency for a focal infection species
# Citation for analysis is Pepin KM et al. (2012) Influenza
# "Multiannual patterns of influenza A transmission in Chinese live bird market systems"
multi.test <- function(observed, # Observed counts for each co-infection with a focal species
                      coinfections, # Total count of co-infections with a focal species
                      prob.obs, # Vector for the observed frequency of each co-infection out of the total
                      prob.exp) # Vector for the expected frequency of each co-infection
{
  # Calculate multinomial likelihood of co-infection counts based on observed co-infection frequencies
  test.obs <- log(dmultinom(x=observed, size=coinfections, prob=prob.obs))
  # Calculate multinomial likelihood of co-infection counts based on expected co-infection frequencies
  test.exp <- log(dmultinom(x=observed, size=coinfections, prob=prob.exp))
  # Calculate likelihood ratio of expected over observed
  test.LR <- -2*(test.exp-test.obs)
  # Make a correction for the difference between the moments of the
  # likelihood ratio statistic and the chi-square distribution
  test.corr <- 1+((sum((1/prob.exp)-1))/(6*coinfections*(length(prob.obs)-1)))
  # Adjust likelihood ratio with correction
  test.LR.corr <- test.LR/test.corr
  # Calculate right-tailed p-value for likelihood ratio from chi-square distribution
  test.p <- 1-pchisq(q=test.LR.corr, df=length(prob.obs)-1)
  return(list("Corrected likelihood ratio"=test.LR.corr, "Test p-value"=test.p))
}

# Create function to perform binomial test of co-infection frequency for a focal infection species
bin.test <- function(observed, # Observed counts for each co-infection with a focal species
                    coinfections, # Total count of co-infections with a focal species
                    prob.obs, # Vector for the observed frequency of each co-infection out of the total
                    prob.exp) # Vector for the expected frequency of each co-infection
{
  # Calculate binomial likelihoods of co-infection counts based on observed co-infection frequencies
  test.obs <- log(dbinom(x=observed, size=coinfections, prob=prob.obs))
  # Calculate binomial likelihoods of co-infection counts based on expected co-infection frequencies
  test.exp <- log(dbinom(x=observed, size=coinfections, prob=prob.exp))
  # Calculate likelihood ratios of expected over observed for all co-infections
  test.LR <- -2*(test.exp-test.obs)
  # Make a correction for the difference between the moments of the
  # likelihood ratio statistic and the chi-square distribution
  test.corr <- 1+((sum((1/prob.exp)-1))/(6*coinfections*(length(prob.obs)-1)))
  # Adjust likelihood ratios with correction
  test.LR.corr <- test.LR/test.corr
  # Calculate right-tailed p-values for likelihood ratios from chi-square distribution
  test.p <- 1-pchisq(q=test.LR.corr, df=1)
  return(list("Corrected likelihood ratio"=test.LR.corr, "Test p-values"=test.p))
}

#### Before/after infection frequency ####

# Before/after species counts
multi.test(observed=c(23, 22, 37, 6, 52, 37, 1, 0),
           coinfections=178,
           prob.obs=c(0.1292134831, 0.1235955056, 0.2078651685, 0.03370786517, 0.2921348315, 0.2078651685, 0.005617977528, 0),
           prob.exp=c(0.02600472813, 0.02600472813, 0.1560283688, 0.1371158392, 0.1111111111, 0.4326241135, 0.09456264775, 0.01654846336))
bin.test(observed=c(23, 22, 37, 6, 52, 37, 1, 0),
         coinfections=178,
         prob.obs=c(0.1292134831, 0.1235955056, 0.2078651685, 0.03370786517, 0.2921348315, 0.2078651685, 0.005617977528, 0),
         prob.exp=c(0.02600472813, 0.02600472813, 0.1560283688, 0.1371158392, 0.1111111111, 0.4326241135, 0.09456264775, 0.01654846336))

# Before/after species abundance
multi.test(observed=c(34, 34, 53, 8, 99, 57, 1, 0),
           coinfections=286,
           prob.obs=c(0.1188811189, 0.1188811189, 0.1853146853, 0.02797202797, 0.3461538462, 0.1993006993, 0.003496503497, 0),
           prob.exp=c(0.01610738255, 0.01879194631, 0.1395973154, 0.1181208054, 0.1248322148, 0.511409396, 0.06040268456, 0.01073825503))
bin.test(observed=c(34, 34, 53, 8, 99, 57, 1, 0),
         coinfections=286,
         prob.obs=c(0.1188811189, 0.1188811189, 0.1853146853, 0.02797202797, 0.3461538462, 0.1993006993, 0.003496503497, 0),
         prob.exp=c(0.01610738255, 0.01879194631, 0.1395973154, 0.1181208054, 0.1248322148, 0.511409396, 0.06040268456, 0.01073825503))

#### Bats/flies infection frequency ####

# Mar10 bats and flies
multi.test(observed=c(0, 1, 2, 0, 7, 15, 0, 1),
           coinfections=26,
           prob.obs=c(0, 1, 2, 0, 7, 15, 0, 1)/26,
           prob.exp=c(0.023, 0.045, 0.14, 0.13, 0.083, 0.55, 0.023, 0.015))
bin.test(observed=c(0, 1, 2, 0, 7, 15, 0, 1),
         coinfections=26,
         prob.obs=c(0, 1, 2, 0, 7, 15, 0, 1)/26,
         prob.exp=c(0.023, 0.045, 0.14, 0.13, 0.083, 0.55, 0.023, 0.015))

# 17Jan12 bats and flies
multi.test(observed=c(0, 1, 5, 2, 6, 6, 0, 0),
           coinfections=20,
           prob.obs=c(0, 1, 5, 2, 6, 6, 0, 0)/20,
           prob.exp=c(1, 5, 21, 6, 17, 23, 5, 6)/sum(c(1, 5, 21, 6, 17, 23, 5, 6)))
bin.test(observed=c(0, 1, 5, 2, 6, 6, 0, 0),
         coinfections=20,
         prob.obs=c(0, 1, 5, 2, 6, 6, 0, 0)/20,
         prob.exp=c(1, 5, 21, 6, 17, 23, 5, 6)/sum(c(1, 5, 21, 6, 17, 23, 5, 6)))

# 17Jan12 captive and wild bats
multi.test(observed=c(1, 5, 21, 6, 17, 23, 5),
           coinfections=78,
           prob.obs=c(1, 5, 21, 6, 17, 23, 5)/78,
           prob.exp=c(34, 34, 53, 8, 99, 57, 1)/sum(c(34, 34, 53, 8, 99, 57, 1)))
bin.test(observed=c(1, 5, 21, 6, 17, 23, 5),
         coinfections=78,
         prob.obs=c(1, 5, 21, 6, 17, 23, 5)/78,
         prob.exp=c(34, 34, 53, 8, 99, 57, 1)/sum(c(34, 34, 53, 8, 99, 57, 1)))

#### Before after flies coinfection frequency ####

# For these tests, species new2 was not included because it was detected in any of these timepoints
# Adjusted expected frequencies before flies to exclude new2

# Species co-infections
multi.RLR <- abs(c(-0.7869727436, -2.299199917, -1.532922535, -1.905969218, -0.5968339487, -0.2217247184, -1.584272974))
2*pnorm(q=multi.RLR, sd=sqrt(12+10), lower.tail=F)
bin.RLR <- abs(c(-1.0528832, -0.9882279369, -1.252403006, 0.2869797434, 1.295475583, -1.988336027, -1.586740833, -4.543580793, 7.080997793, -0.6664980442, -4.580046852, -1.114218595, -2.857970767, -2.249298051, 2.74421875, -1.474793372, -10.74048803, 3.113515282, 9.475720802, -1.732867764, -1.887819224, -2.91757545, 2.319958218, -6.153212065, 1.097931146, 2.113546777, -3.416579303, -1.904410177, 0.4360878794, 0.658828829, -0.6091411429, 0.4435645137, 2.043782012, -0.2575337067, 0.07330337225, -1.108898354, -2.071580539, -3.143537266, -0.08350079776, -3.829263839, 4.330446999, -0.35267413))
2*pnorm(q=bin.RLR, sd=sqrt(2+2), lower.tail=F)

# EXTRA FIGURES -----------------------------------------------------------

#### Sup Fig: Change in infection between groups (CIG) ####

# Calculate proportion of age groups and sexes positive at start of study and by end of study
# Neonate (NEO)/juvenile (JUV), sexually immature adult (SI), sexually mature adult (A); female (F) and male (M)
change <- data.frame(age.sex=c("NEO/JUV", "SI", "A", "F", "M"),
                    group=c(rep("begin", 5), rep("end", 5)),
                    groupnum=c(1, 2, 3, 1, 2, 1, 2, 3, 1, 2),
                    start.end=c(0, 17, 42, 31, 29, 29, 17, 49, 53, 42),
                    total=c(33, 17, 52, 58, 44, 33, 17, 52, 58, 44))
change.CI <- binom.confint(change$start.end, change$total, methods="wilson")
change <- cbind(change, change.CI[,4:6])

# Test to see if there are differences between groups at start and end
# NEO/JUV start vs. end
prop.test(x=c(0, 29), n=c(33, 33), alternative="less")

# SI start vs. end
prop.test(x=c(17, 17), n=c(17, 17), alternative="less")
fisher.test(x=matrix(c(17, 0, 17, 0), nrow=2), alternative="less")

# A start vs. end
prop.test(x=c(42, 49), n=c(52, 52), alternative="less")

# NEO/JUV, SI, and A at start
prop.test(x=c(0, 17, 42), n=c(33, 17, 52), alternative="two.sided")

# NEO/JUV, SI, and A at end
prop.test(x=c(29, 17, 49), n=c(33, 17, 52), alternative="two.sided")

# F start vs. end
prop.test(x=c(31, 53), n=c(58, 58), alternative="less")

# M start vs. end
prop.test(x=c(29, 42), n=c(44, 44), alternative="less")

# F vs. M at start
prop.test(x=c(31, 29), n=c(58, 44), alternative="two.sided")

# F vs. M at end
prop.test(x=c(53, 42), n=c(58, 44), alternative="two.sided")
fisher.test(x=matrix(c(53, 5, 42, 2), nrow=2), alternative="two.sided")

# Plot proportion positive by age group
SuppFigCIGa <- ggplot(data=change[c(1:3, 6:8),], aes(x=as.factor(groupnum), y=mean, color=group)) +
  geom_pointrange(aes(y=mean, ymin=lower, ymax=upper), position=position_dodge(width=0.75)) +
  theme_cowplot() +
  scale_x_discrete(name="Age class", labels=c("1"="NEO/JUV", "2"="SI", "3"="A")) +
  scale_y_continuous(name="Proportion positive", limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  scale_color_manual(name="Group", labels=c("Start", "End"), values=c("#a6611a", "#dfc27d"))

# Plot proportion positive by sex
SuppFigCIGb <- ggplot(data=change[c(4:5, 9:10),], aes(x=age.sex, y=mean, color=group)) +
  geom_pointrange(aes(y=mean, ymin=lower, ymax=upper), position=position_dodge(width=0.75)) +
  theme_cowplot() +
  scale_x_discrete(name="Group") +
  scale_y_continuous(name="Proportion positive", limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  scale_color_manual(name="Sex", labels=c("Start", "End"), values=c("#80cdc1", "#018571"))

# Combine plots
plot_grid(SuppFigCIGa, SuppFigCIGb, labels=c("A", "B"), nrow=2, align="hv")
ggsave("Results/SuppFig_infection_change.png", height=6, width=6, units="in", dpi=300, bg="white")

#### Sup Fig: Prevalence by number of positive markers ####

# Melt data
pos_markers_melt <- melt(CI_join,
                     measure.vars = c("one", "two", "three", "four"),
                     id.vars = c("Date", "Day", "Group"))

# Line plot of captive prevalence data separated by number of positive markers
ggplot(data=pos_markers_melt, aes(x=Day, y=value, color=variable, shape=Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values=c("grey60", "grey40", "grey20", "grey0"),
                      name="Positive\nmarkers", labels=c("1", "2", "3", "4")) +
  scale_shape_manual(values=c(20, 15, 17, 18)) +
  scale_x_continuous(name = "Day", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Prevalence", breaks = seq(0, 1, 0.25)) +
  theme_cowplot()
ggsave("Results/SuppFig_pos_markers.png", height=3.5, width=7, units="in", dpi=300, bg="white")

#### Sup Fig: Mixed infection prevalence and number of positive markers ####

# Statistics
inf_join_stats_Markers <- inf_join %>%
  group_by(Group, Day) %>%
  summarize(n = n(),
            mean = mean(Markers, na.rm = TRUE),
            min = fivenum(Markers, na.rm = TRUE)[1],
            lower = fivenum(Markers, na.rm = TRUE)[2],
            median = fivenum(Markers, na.rm = TRUE)[3],
            upper = fivenum(Markers, na.rm = TRUE)[4],
            max = fivenum(Markers, na.rm = TRUE)[5])

# Line plot of number of positive markers with points for fly and wild bat data
ggplot(data = inf_join, aes(x = Day, y = Markers, color = Group)) +
  geom_count(shape = 21) +
  geom_point(data = inf_join_stats_Markers, aes(x = Day, y = mean, color = Group), size = 3, alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = 0,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_size_area(name="Counts", breaks = c(1, 10, 20)) +
  scale_x_continuous(name = "", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Positive markers", breaks = seq(0, 4, 1)) +
  theme_cowplot()
ggsave("Results/SuppFig_markers.png", height=3.5, width=7, units="in", dpi=300, bg = "white")

# Line plot for captive mixed infection data with points for fly and wild bat data
ggplot(data = CI_join, aes(x = Day, y = mixed.mean, color = Group)) +
  geom_pointrange(aes(x = Day, y = mixed.mean, ymin = mixed.lower, ymax = mixed.upper, color = Group), alpha = 0.8) +
  annotate(geom = "text", x = prev.CI$Day, y = -0.05,
           label = c("Jul 09", "Nov 09", "Jan 10", "Mar 10", "May 10", "Jul 10", "Sep 10",
                     "Nov 10", "Mar 11", "Jul 11", "Jan 12", "", "", "Mar 12"),
           angle = 45, size = c(rep(2, 10), 2.25, rep(2, 3)), fontface = c(rep(1, 10), 2, rep(1, 3))) +
  scale_color_colorblind(name = "Group") +
  scale_x_continuous(name = "Day", limits = c(-5, 1005), breaks = seq(0, 1000, 250)) +
  scale_y_continuous(name = "Coinfection prevalence", limits = c(-0.05, 0.75), breaks = seq(0, 0.75, 0.25)) +
  theme_cowplot()
ggsave("Results/SuppFig_mixed.png", height=3.5, width=7, units="in", dpi=300, bg = "white")

### Sup Fig: Prevalence and infection load segmented regression ####

# Captive bats positive markers
one <- c(9, 13, 62, 57, 37, 44, 44, 40, 10, 18, 27, 33, 35, 39)
two <- c(6, 11, 55, 46, 27, 37, 36, 27, 3, 12, 14, 24, 22, 36)
three <- c(4, 8, 45, 38, 23, 27, 27, 23, 2, 5, 9, 17, 15, 30)
four <- c(3, 6, 18, 27, 17, 10, 12, 7, 0, 0, 2, 8, 9, 21)
total <- c(12, 22, 73, 61, 72, 63, 69, 67, 58, 81, 87, 80, 85, 81)
# Wilson score confidence intervals
prev.CI1 <- binom.confint(x=c(one, two, three, four), n=rep(total, 4), method="wilson")

# Read in captive prevalence data
prev.CI <- read.csv("Data/captive_colony_prevalence.csv", head=T)
# Reshape marker data into long format
markers <- melt(data=prev.CI, id.var="Day", measure.vars=c("one", "two", "three", "four"))

# Segmented regression analysis for all marker prevalence
# Read in dates and calculate days since start
prev.CI1$Date <- rep(c("2009-07-28", "2009-11-05", "2010-01-28", "2010-03-06", "2010-05-21", "2010-07-14",
                      "2010-09-23", "2010-11-05", "2011-03-04", "2011-07-13", "2012-01-17", "2012-01-31",
                      "2012-02-14", "2012-03-15"), 4)
prev.CI1$day <- as.numeric(as.Date(prev.CI1$Date)-as.Date(prev.CI1$Date[1]))
prev.CI1$f <- prev.CI1$n-prev.CI1$x
# Perform regression and segmented regression
mod1 <- glm(mean~day, prev.CI1, family=binomial(link="logit"), weights=n)
set.seed(20250430)
segmented.mod1 <- segmented(obj=mod1, seg.Z=~day, psi=c(200, 800),
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod1)
(sum.mod1 <- summary(segmented.mod1))
(slope.mod1 <- slope(segmented.mod1)$day)
colnames(slope.mod1) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod1
confint.lm(segmented.mod1)
# Extract breakpoints and predicted values
breaks1 <- data.frame(confint.segmented(segmented.mod1))
colnames(breaks1) <- c("estimate", "lower", "upper"); breaks1
xweight1 <- data.frame(day=c(0, breaks1$estimate, 961))
yweight1 <- predict.segmented(segmented.mod1, newdata=xweight1, type="response")
pred1 <- data.frame(cbind(xweight1, yweight1))

# Plot of segmented regression for all marker prevalence
(SuppFigSegInfa <- ggplot(data=markers, aes(x=Day, y=value)) +
  geom_line(data=pred1, aes(x=day, y=yweight1), linetype=2) +
  geom_point(data=markers, aes(x=Day, y=value, color=variable), alpha=0.8, size=3) +
  geom_line(data=data.frame(x=as.numeric(breaks1[1, 2:3]), y=-0.1), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_line(data=data.frame(x=as.numeric(breaks1[2, 2:3]), y=-0.1), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks1[,1], y=-0.1), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  scale_colour_manual(values=c("grey60", "grey40", "grey20", "grey0"),
                      name="Positive\nmarkers", labels=c("1", "2", "3", "4")) +
  xlab("") +
  scale_y_continuous(name = "Prevalence", limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Perform regression and segmented regression
# Join data and filter RT-PCR data to remove negative results
inf.dat.filter <- inf.dat %>%
  filter(!is.na(Ct), Ct<40)

mod2 <- glm(Ct~Day, inf.dat.filter, family=Gamma(link="inverse"))
set.seed(20250430)
segmented.mod2 <- segmented(obj=mod2, seg.Z=~Day, psi=c(200, 800),
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod2)
(sum.mod2 <- summary(segmented.mod2))
(slope.mod2 <- slope(segmented.mod2)$Day)
colnames(slope.mod2) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod2
confint.lm(segmented.mod2)
# Extract breakpoints and predicted values
breaks2 <- data.frame(confint.segmented(segmented.mod2))
colnames(breaks2) <- c("estimate", "lower", "upper"); breaks2
xweight2 <- data.frame(Day=c(0, breaks2$estimate, 961))
yweight2 <- predict.segmented(segmented.mod2, newdata=xweight2, type="response")
pred2 <- data.frame(cbind(xweight2, yweight2))

# Plot of segmented regression for real-time Ct values
(SuppFigSegInfb <- ggplot(data=inf.dat.filter, aes(x=Day, y=Ct)) +
    geom_line(data=pred2, aes(x=Day, y=yweight2), linetype=2) +
    geom_jitter(size = 0.75, shape = 21, width = 5, height = 0) +
    geom_line(data=data.frame(x=as.numeric(breaks2[1, 2:3]), y=45), aes(x=x, y=y), alpha=0.8, color="grey60") +
    geom_line(data=data.frame(x=as.numeric(breaks2[2, 2:3]), y=45), aes(x=x, y=y), alpha=0.8, color="grey60") +
    geom_point(data=data.frame(x=breaks2[,1], y=45), aes(x=x, y=y), shape=15, size=3, color="grey60") +
    theme_cowplot() +
    xlab("") +
    scale_y_reverse(name = "Ct value", breaks = seq(45, 20, -5)) +
    theme_cowplot() +
    coord_cartesian(xlim = c(-5, 1005)))

# Co-infections
mixed.pos <- c(4, 3, 24, 17, 7, 7, 5, 7, 0, 2, 7, 11, 8, 12)
mixed.total <- c(12, 22, 73, 61, 72, 63, 69, 67, 58, 81, 87, 80, 85, 81)
mixed.CI <- binom.confint(mixed.pos, mixed.total, method="wilson")

# Segmented regression analysis for mixed infection prevalence
# Read in dates and calculate days since start
mixed.CI$Date <- c("2009-07-28", "2009-11-05", "2010-01-28", "2010-03-06", "2010-05-21", "2010-07-14",
                  "2010-09-23", "2010-11-05", "2011-03-04", "2011-07-13", "2012-01-17", "2012-01-31",
                  "2012-02-14", "2012-03-15")
mixed.CI$day <- as.numeric(as.Date(mixed.CI$Date)-as.Date(mixed.CI$Date[1]))
# Perform regression and segmented regression
mod4 <- glm(mean~day, mixed.CI, family=binomial(link="logit"), weight=n)
set.seed(20250430)
segmented.mod4 <- segmented(obj=mod4, seg.Z=~day, psi=c(200, 800),
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod4)
(sum.mod4 <- summary(segmented.mod4))
(slope.mod4 <- slope(segmented.mod4)$day)
colnames(slope.mod4) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod4
confint.lm(segmented.mod4)
# Extract breakpoints and predicted values
breaks4 <- data.frame(confint.segmented(segmented.mod4))
colnames(breaks4) <- c("estimate", "lower", "upper"); breaks4
xweight4 <- data.frame(day=c(0, breaks4$estimate, 961))
yweight4 <- predict.segmented(segmented.mod4, newdata=xweight4, type="response")
pred4 <- data.frame(cbind(xweight4, yweight4))

# Plot of segmented regression for mixed infection prevalence
(SuppFigSegInfd <- ggplot(data=prev.CI, aes(x=Day, y=mixed.mean)) +
  geom_line(data=pred4, aes(x=day, y=yweight4), linetype=2) +
  geom_point(size=3, alpha=0.8) +
  geom_line(data=data.frame(x=as.numeric(breaks4[1, 2:3]), y=-0.1), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_line(data=data.frame(x=as.numeric(breaks4[2, 2:3]), y=-0.1), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks4[,1], y=-0.1), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  xlab("Day") +
  scale_y_continuous(name = "Coinfection prevalence", limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Combine plots
plot_grid(SuppFigSegInfa, SuppFigSegInfb, SuppFigSegInfd, labels=c("A", "B", "C"),
          nrow=3, align="hv")
ggsave("Results/SuppFig_segmented_prevalence+load.png", height=9, width=7, units="in", dpi=300, bg="white")

#### Sup Fig: Correlations between Bartonella species counts and abundance ####

# Test correlation between all timepoints
all.bat.count <- c(34, 33, 103, 64, 99, 220, 41, 7)
all.bat.abund <- c(46, 48, 157, 96, 192, 438, 46, 8)
cor.test(all.bat.count, all.bat.abund)

# Read in bat count data
bat.count.data <- as.data.frame(read.csv("Data/captive_colony_counts.csv", head=T))
# Read in bat abundance data
bat.abund.data <- as.data.frame(read.csv("Data/captive_colony_abundance.csv", head=T))
# Create a data frame to store results
bat.cor <- data.frame(array(0, c(14, 7)))
colnames(bat.cor) <- c("Date", "Day", "Estimate", "CI.lower", "CI.upper", "p.val")
bat.cor$Date <- bat.count.data$Date
bat.cor$Day <- bat.count.data$Day
for(i in 1:14){
  test <- cor.test(as.numeric(bat.count.data[i, 2:9]), as.numeric(bat.abund.data[i, 2:9]), method="pearson")
  bat.cor[i, 3:7] <- c(as.numeric(test$estimate), test$conf.int[[1]], test$conf.int[[2]], test$p.value, test$p.value<0.05)
}
bat.cor$R2 <- bat.cor$Estimate^2
bat.cor$R2lower <- bat.cor$CI.lower^2
bat.cor$R2upper <- bat.cor$CI.upper^2

SuppFigCorra <- ggplot(data=bat.cor, aes(x=Day, y=R2, group=1)) +
  geom_pointrange(data=bat.cor,
                  aes(x=Day, y=Estimate, ymin=R2lower, ymax=R2upper), shape = 20, alpha=0.8) +
  theme_cowplot() +
  ylim(0.7, 1) +
  xlab("Date") +
  ylab(bquote(R^2))

# Read in fly count data
fly.count.data <- as.data.frame(read.csv("Data/fly_comparisons_counts.csv", head=T))
# Read in fly abundance data
fly.abund.data <- as.data.frame(read.csv("Data/fly_comparisons_abundance.csv", head=T))
# Create a data frame to store results
fly.cor <- data.frame(array(0, c(3, 5)))
colnames(fly.cor) <- c("Group", "Estimate", "CI.lower", "CI.upper", "p.val")
fly.cor$Group <- fly.count.data$Group
for(i in 1:3){
  test <- cor.test(as.numeric(fly.count.data[i, 2:9]), as.numeric(fly.abund.data[i, 2:9]), method="pearson")
  fly.cor[i, 2:5] <- c(as.numeric(test$estimate), test$conf.int[[1]], test$conf.int[[2]], test$p.value)
}
fly.cor$R2 <- fly.cor$Estimate^2
fly.cor$R2lower <- fly.cor$CI.lower^2
fly.cor$R2upper <- fly.cor$CI.upper^2

SuppFigCorrb <- ggplot(data=fly.cor, aes(x=Group, y=Estimate, group=1)) +
  geom_pointrange(data=fly.cor,
                  aes(x=Group, y=Estimate, ymin=CI.lower, ymax=CI.upper), shape=18, alpha=0.8) +
  theme_cowplot() +
  ylim(0.7, 1) +
  xlab("Group") +
  ylab("")

# Combine plots
plot_grid(SuppFigCorra, SuppFigCorrb, labels=c("A", "B"), ncol=2, rel_widths=c(0.625, 0.375), align="hv")
ggsave("Results/SuppFig_correlations.png", height=3.5, width=10, units="in", dpi=300, bg="white")

#### Sup Fig: Bartonella diversity measures segmented regression I ####

# Segmented regression analysis for Bartonella species richness
# Perform regression and segmented regression
mod5 <- glm(Species.richness~Day, inf.div3, family=poisson(link="log"))
set.seed(20250430)
segmented.mod5 <- segmented(obj=mod5, seg.Z=~Day, psi=200,
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod5)
(sum.mod5 <- summary(segmented.mod5))
slope.mod5 <- slope(segmented.mod5)$Day
colnames(slope.mod5) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod5
confint.lm(segmented.mod5)
# Extract breakpoints and predicted values
breaks5 <- data.frame(confint.segmented(segmented.mod5))
colnames(breaks5) <- c("estimate", "lower", "upper"); breaks5
xweight5 <- data.frame(Day=c(0, breaks5$estimate, 961))
yweight5 <- predict.segmented(segmented.mod5, newdata=xweight5, type="response")
pred5 <- data.frame(cbind(xweight5, yweight5))

# Plot of segmented regression for species richness
(SuppFigSegDivIa <- ggplot(data=inf.div3, aes(x=Day, y=Species.richness)) +
  geom_line(data=pred5, aes(x=Day, y=yweight5), linetype=2) +
  geom_point(size=3, alpha=0.8) +
  geom_line(data=data.frame(x=as.numeric(breaks5[1, 2:3]), y=0), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks5[,1], y=0), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  xlab("") +
  scale_y_continuous(name = "Species richness", limits = c(0, 8), breaks = seq(0, 8, 2)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Segmented regression analysis for Shannon numbers
# Perform regression and segmented regression
mod6 <- glm(Shannon.number~Day, inf.div3, family=Gamma(link="identity"))
set.seed(20250430)
segmented.mod6 <- segmented(obj=mod6, seg.Z=~Day, psi=500,
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod6)
(sum.mod6 <- summary(segmented.mod6))
slope.mod6 <- slope(segmented.mod6)$Day
colnames(slope.mod6) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod6
confint.lm(segmented.mod6)
# Extract breakpoints and predicted values
breaks6 <- data.frame(confint.segmented(segmented.mod6))
colnames(breaks6) <- c("estimate", "lower", "upper"); breaks6
xweight6 <- data.frame(Day=c(0, breaks6$estimate, 961))
yweight6 <- predict.segmented(segmented.mod6, newdata=xweight6, type="response")
pred6 <- data.frame(cbind(xweight6, yweight6))

# Plot of segmented regression for Shannon numbers
(SuppFigSegDivIb <- ggplot(data=inf.div3, aes(x=Day, y=Shannon.number)) +
  geom_line(data=pred6, aes(x=Day, y=yweight6), linetype=2) +
  geom_point(size=3, alpha=0.8) +
  geom_line(data=data.frame(x=as.numeric(breaks6[1, 2:3]), y=0), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks6[,1], y=0), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  xlab("") +
  scale_y_continuous(name = "Shannon number", limits = c(0, 8), breaks = seq(0, 8, 2)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Segmented regression analysis for inverse Simpson index
# Perform regression and segmented regression
mod7 <- glm(InvSimpson~Day, inf.div3, family=Gamma(link="identity"))
set.seed(20250430)
segmented.mod7 <- segmented(obj=mod7, seg.Z=~Day, psi=500,
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod7)
(sum.mod7 <- summary(segmented.mod7))
slope.mod7 <- slope(segmented.mod7)$Day
colnames(slope.mod7) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod7
confint.lm(segmented.mod7)
# Extract breakpoints and predicted values
breaks7 <- data.frame(confint.segmented(segmented.mod7))
colnames(breaks7) <- c("estimate", "lower", "upper"); breaks7
xweight7 <- data.frame(Day=c(0, breaks7$estimate, 961))
yweight7 <- predict.segmented(segmented.mod7, newdata=xweight7, type="response")
pred7 <- data.frame(cbind(xweight7, yweight7))

# Plot of segmented regression for inverse Simpson index
(SuppFigSegDivIc <- ggplot(data=inf.div3, aes(x=Day, y=InvSimpson)) +
  geom_line(data=pred7, aes(x=Day, y=yweight7), linetype=2) +
  geom_point(size=3, alpha=0.8) +
  geom_line(data=data.frame(x=as.numeric(breaks7[1, 2:3]), y=0), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks7[,1], y=0), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  xlab("Day") +
  scale_y_continuous(name = "Inv. Simpson index", limits = c(0, 8), breaks = seq(0, 8, 2)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Combine plots
plot_grid(SuppFigSegDivIa, SuppFigSegDivIb, SuppFigSegDivIc, labels=c("A", "B", "C"),
          nrow=3, align="hv")
ggsave("Results/SuppFig_segmented_diversityI.png", height=9, width=7, units="in", dpi=300, bg="white")

#### Sup Fig: Bartonella diversity measures segmented regression II ####

# Read in captive RT-PCR data
inf.dat <- read.csv("Data/captive_colony_infection.csv", head=T)

# Segmented regression analysis for number of Bartonella species in each sample
# Perform regression and segmented regression
mod8 <- glm(NumSpecies~Day, inf.dat, family=poisson(link="log"))
set.seed(20250430)
segmented.mod8 <- segmented(obj=mod8, seg.Z=~Day, psi=500,
                            control=seg.control(n.boot=10, it.max=1000, maxit.glm=1000))
AICc(segmented.mod8)
(sum.mod8 <- summary(segmented.mod8))
slope.mod8 <- slope(segmented.mod8)$Day
colnames(slope.mod8) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod8
confint.lm(segmented.mod8)
# Extract breakpoints and predicted values
breaks8 <- data.frame(confint.segmented(segmented.mod8))
colnames(breaks8) <- c("estimate", "lower", "upper"); breaks8
xweight8 <- data.frame(Day=c(0, breaks8$estimate, 961))
yweight8 <- predict.segmented(segmented.mod8, newdata=xweight8, type="response")
pred8 <- data.frame(cbind(xweight8, yweight8))

# Plot of segmented regression for number of Bartonella species per sample
(SuppFigSegDivIIa <- ggplot(data=inf.dat, aes(x=Day, y=NumSpecies)) +
  geom_line(data=pred8, aes(x=Day, y=yweight8), linetype=2) +
  geom_count(shape=21) +
  geom_line(data=data.frame(x=as.numeric(breaks8[1, 2:3]), y=0), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks8[,1], y=0), aes(x=x, y=y), shape=15, size=3, alpha=0.8, color="grey60") +
  scale_size_area(name="Counts", breaks = c(1, 20, 40)) +
  xlab("") +
  scale_y_continuous(name = "Species in sample", breaks = seq(0, 4, 1)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Segmented regression analysis for Bartonella bat diversity
# Perform regression and segmented regression
mod9 <- glm(index1~days, dissim, family=binomial(link="logit"))
set.seed(20250430)
segmented.mod9 <- segmented(obj=mod9, seg.Z=~days, psi=c(200, 800),
                            control=seg.control(n.boot=10, it.max=1e6, maxit.glm=1e6))
AICc(segmented.mod9)
(sum.mod9 <- summary(segmented.mod9))
(slope.mod9 <- slope(segmented.mod9)$days)
colnames(slope.mod9) <- c("estimate", "st.err", "t.value", "lower", "upper"); slope.mod9
confint.lm(segmented.mod9)
# Extract breakpoints and predicted values
breaks9 <- data.frame(confint.segmented(segmented.mod9))
colnames(breaks9) <- c("estimate", "lower", "upper"); breaks9
xweight9 <- data.frame(days=c(0, breaks9$estimate, 961))
yweight9 <- predict.segmented(segmented.mod9, newdata=xweight9, type="response")
pred9 <- data.frame(cbind(xweight9, yweight9))

(SuppFigSegDivIIb <- ggplot(data=dissim, aes(x=days, y=index1, group=1)) +
  geom_line(data=pred9, aes(x=days, y=yweight9), linetype=2) +
  geom_count(shape=21) +
  geom_line(data=data.frame(x=as.numeric(breaks9[1, 2:3]), y=-0.15), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_line(data=data.frame(x=as.numeric(breaks9[2, 2:3]), y=-0.15), aes(x=x, y=y), alpha=0.8, color="grey60") +
  geom_point(data=data.frame(x=breaks9[,1], y=-0.15), aes(x=x, y=y), shape=15, size=3, color="grey60") +
  scale_size_area(name="Counts", breaks = c(1, 100, 500)) +
  xlab("Day") +
  scale_y_continuous(name = "Beta diversity", breaks = seq(0, 1, 0.2)) +
  theme_cowplot() +
  coord_cartesian(xlim = c(-5, 1005)))

# Combine plots
plot_grid(SuppFigSegDivIIa, SuppFigSegDivIIb, labels=c("A", "B"), nrow=2, align="hv")
ggsave("Results/SuppFig_segmented_diversityII.png", height=7, width=7, units="in", dpi=300, bg="white")

#### Testing correlations between measures of infection and diversity ####

#### Ct values ####

# Check which samples were real-time positive
inf.dat$RTpos <- is.na(inf.dat$Ct)
inf.dat$RTpos <- ifelse(inf.dat$RTpos==FALSE, 1, 0)
inf.dat[which(inf.dat$Ct>40),]$RTpos=0

# Determine correlation between Ct values and markers and between real-time positives and markers
cor.test(inf.dat$Ct, inf.dat$Markers)
cor.test(inf.dat$RTpos, inf.dat$Markers)
boxplot(inf.dat$RTpos, inf.dat$Markers)
RTpos.mark <- glm(RTpos~Markers, inf.dat, family=binomial(link="logit"))
summary(RTpos.mark)
RTpos.mark.xweight <- seq(1, 4, 0.01)
RTpos.mark.yweight <- predict(RTpos.mark, list(Markers=RTpos.mark.xweight), type="response")
plot(inf.dat$Markers, inf.dat$RTpos)
lines(RTpos.mark.xweight, RTpos.mark.yweight, col="red")

# Determine correlation between Ct values and markers
cor.test(inf.dat.filter$Ct, inf.dat.filter$Markers)
# Determine correlation between Ct values and number of species
cor.test(inf.dat$Ct, inf.dat$NumSpecies)
cor.test(inf.dat.filter$Ct, inf.dat.filter$NumSpecies)

# Test normality of Ct values
hist(inf.dat.filter$Ct)
shapiro.test(inf.dat.filter$Ct)

# Ct values between single and co-infections
boxplot(Ct~Coinfection, inf.dat.filter)
kruskal.test(Ct~Coinfection, inf.dat.filter)
Ct.coinf <- glm(Coinfection~Ct, inf.dat.filter, family=binomial(link="logit"))
summary(Ct.coinf)
Ct.coinf.xweight <- seq(20, 50, 0.1)
Ct.coinf.yweight <- predict(Ct.coinf, list(Ct=Ct.coinf.xweight), type="response")
plot(inf.dat.filter$Ct, inf.dat.filter$Coinfection)
lines(Ct.coinf.xweight, Ct.coinf.yweight, col="red")
data.frame(plyr::ddply(inf.dat.filter, ~Coinfection, summarise, mean=mean(Ct), median=median(Ct)))

# Ct values by markers
boxplot(Ct~Markers, inf.dat.filter)
kruskal.test(Ct~Markers, inf.dat.filter)
Ct.mark <- glm(Ct~Markers, inf.dat.filter, family=Gamma(link="identity"))
hist(Ct.mark$residuals)
qqnorm(Ct.mark$residuals)
shapiro.test(Ct.mark$residuals)
summary(Ct.mark)
plot(inf.dat.filter$Markers, inf.dat.filter$Ct)
abline(Ct.mark, col="red")
data.frame(plyr::ddply(inf.dat.filter, ~Markers, summarise, mean=mean(Ct), median=median(Ct)))

# Ct values by species for all infections
boxplot(Ct~Species, inf.dat.filter)
summary(aov(Ct~Species, inf.dat.filter))
kruskal.test(Ct~Species, inf.dat.filter)
data.frame(plyr::ddply(inf.dat.filter, ~Species, summarise, mean=mean(Ct), median=median(Ct)))

# Ct values by number of species
boxplot(Ct~NumSpecies, inf.dat.filter)
kruskal.test(Ct~NumSpecies, inf.dat.filter)
Ct.Numsp <- glm(Ct~NumSpecies, inf.dat.filter, family=Gamma(link="identity"))
hist(Ct.Numsp$residuals)
qqnorm(Ct.Numsp$residuals)
shapiro.test(Ct.Numsp$residuals)
summary(Ct.Numsp)
plot(inf.dat.filter$NumSpecies, inf.dat.filter$Ct)
abline(Ct.Numsp, col="red")
data.frame(plyr::ddply(inf.dat.filter, ~NumSpecies, summarise, mean=mean(Ct), median=median(Ct)))

#### Markers ####

# Test normality of mark values
hist(inf.dat$Marker)
shapiro.test(inf.dat$Marker)

# Marker values between single and co-infections
boxplot(Markers~Coinfection, inf.dat)
kruskal.test(Markers~Coinfection, inf.dat)
mark.coinf <- glm(Coinfection~Markers, inf.dat, family=binomial(link="logit"))
summary(mark.coinf)
mark.coinf.xweight <- seq(1, 4, 0.01)
mark.coinf.yweight <- predict(mark.coinf, list(Markers=mark.coinf.xweight), type="response")
plot(inf.dat$Markers, inf.dat$Coinfection)
lines(mark.coinf.xweight, mark.coinf.yweight, col="red")
data.frame(plyr::ddply(inf.dat, ~Coinfection, summarise, mean=mean(Markers), median=median(Markers)))

# Marker values by species for all infections
boxplot(Markers~Species, inf.dat)
kruskal.test(Markers~Species, inf.dat)
data.frame(plyr::ddply(inf.dat, ~Species, summarise, mean=mean(Markers), median=median(Markers)))

# Marker values by number of species
boxplot(Markers~NumSpecies, inf.dat)
kruskal.test(Markers~NumSpecies, inf.dat)
mark.Numsp <- glm(Markers~NumSpecies, inf.dat, family=Gamma(link="identity"))
summary(mark.Numsp)
mark.Numsp.xweight <- seq(1, 4, 0.01)
mark.Numsp.yweight <- predict(mark.Numsp, list(NumSpecies=mark.Numsp.xweight), type="response")
plot(inf.dat$NumSpecies, inf.dat$Markers)
lines(mark.Numsp.xweight, mark.Numsp.yweight, col="red")
data.frame(plyr::ddply(inf.dat, ~NumSpecies, summarise, mean=mean(Markers), median=median(Markers)))

