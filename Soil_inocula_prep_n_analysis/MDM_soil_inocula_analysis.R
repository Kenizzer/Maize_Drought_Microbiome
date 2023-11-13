# Processing and analysis of soil and inocula samples from across
# an E/W precipitation gradient in KS
# Code by: Natalie Ford & Joel Swift

#Packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(ALDEx2); packageVersion("ALDEx2")
library(lme4); packageVersion("lme4")
library(emmeans); packageVersion("emmeans")
library(ggpubr); packageVersion("ggpubr")

#Themes and color palette 
theme_set(theme_pubr())
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#7fd66f", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location

### Load datasets for use throughout ###
# Bacteria datasets
#Load cleaned, non-transformed dataset 
hmep.bact <- readRDS('./Intermediate_data/phyloseq_b_clean_soil_inocula.RDS')
#Load cleaned, CLR-transformed datasets 
hmep.bact.clr <- readRDS('./Intermediate_data/phyloseq_b_clean_soil_inocula_clr.RDS')
# Fungal datasets
#Load cleaned, non-transformed dataset 
hmep.fungi <- readRDS('./Intermediate_data/phyloseq_f_clean_soil_inocula.RDS')
#Load cleaned, CLR-transformed datasets 
hmep.fungi.clr <- readRDS('./Intermediate_data/phyloseq_f_clean_soil_inocula_clr.RDS')

#### Figure 2 ####
## RDA plots for bacterial and fungal soil and inocula samples

#bacteria
sample_data(hmep.bact.clr)$Location <- factor(sample_data(hmep.bact.clr)$Location, levels = c("SVR", "HAY", "TLI", "KNZ"))
#Ordinate (constrained PCOA)
RDA.ait.hmep.bact.clr <- ordinate(hmep.bact.clr, method='RDA',distance='euclidean',formula=~Location + Type + Condition(logObs))
#Color by Type (soil vs. inocula)
RDA_16S <- plot_ordination(hmep.bact.clr, RDA.ait.hmep.bact.clr)
RDA_16S[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- 'white'
RDA_16S <- RDA_16S + 
  geom_point(aes(fill = Location, shape = Type), size = 3, alpha = 0.80, color = "black") +
  scale_shape_manual(values=c(24, 21), name = "Substrate") +
  scale_fill_manual(values = location_pallete, name = "Soil Source")+
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))

#fungi
sample_data(hmep.fungi.clr)$Location <- factor(sample_data(hmep.fungi.clr)$Location, levels = c("SVR", "HAY", "TLI", "KNZ"))
#Ordinate 
RDA.ait.hmep.fungi.clr <- ordinate(hmep.fungi.clr, method = "RDA", distance = "euclidean", formula=~Location + Type + Condition(logObs))
#Color 
RDA_ITS <- plot_ordination(hmep.fungi.clr, RDA.ait.hmep.fungi.clr)
RDA_ITS[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- 'white'
RDA_ITS <- RDA_ITS + 
  geom_point(aes(fill = Location, shape = Type), size = 3, alpha = 0.80, color = "black") +
  scale_shape_manual(values=c(24, 21), name = "Substrate") +
  scale_fill_manual(values = location_pallete, name = "Soil Source")+
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))

#Combine plots
RDA_16S_ITS <- ggarrange(RDA_16S, RDA_ITS, labels = "AUTO", common.legend = TRUE, legend = "right")
RDA_16S_ITS
ggsave("./figures/Figure2_RDA_16S_ITS.svg", RDA_16S_ITS, width = 6, height = 4)
ggsave("./figures/Figure2_RDA_16S_ITS.png", RDA_16S_ITS, width = 6, height = 4)


#### Supplemental Figure 1/3 ####
## Supplemental Figure 1
## Taxonomic barplots for bacterial and fungal soil and inocula samples
phylo.16s_phylum <- phyloseq::tax_glom(hmep.bact, taxrank= "Phylum") # 11 taxa
# Make relative abundance 
phylo.16s_phylum_relab <- transform_sample_counts(phylo.16s_phylum, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_relab_top10_16S <- fantaxtic::get_top_taxa(physeq_obj = phylo.16s_phylum_relab, n = 11, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
phylum_top10_16S <- psmelt(phylum_relab_top10_16S)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
phylum_top10_16S$Phylum <- forcats::fct_relevel(as.factor(phylum_top10_16S$Phylum), "Other", after = Inf)
sapply(phylum_top10_16S, levels) #factor levels

palette_MOD_soil_16s <- c(
  "Armatimonadetes" = "#FF5733",
  "Chloroflexi" = "#AA4499",
  "Nitrospirae" = "#51C14B",
  "Planctomycetes" = "#6699CC",
  "Thaumarchaeota" = "#138D75",
  "Actinobacteria" = "#D55E00",
  "Firmicutes" = "#6948b8", 
  "Gemmatimonadetes" = "#CC6677",
  "Proteobacteria" = "#15632c",
  "Verrucomicrobia" = "#882255",
  "Bacteroidetes" = "#DDCC77",
  "Other" = "#888888")

# Fixing a couple things for cleaner plotting
# Order W -> E and rename the landuse terms to fit the paper
phylum_top10_16S$Location <- factor(phylum_top10_16S$Location, levels = c("SVR", "HAY", "TLI", "KNZ"))
# soil taxonomic barplot
top10_16S_soil <- ggplot(phylum_top10_16S[phylum_top10_16S$Type == "Soil", ], aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Phylum", values = palette_MOD_soil_16s, labels = c("Actinomycetota", "Armatimonadota", "Bacteroidota", "Chloroflexota", "Bacillota", "Gemmatimonadota", "Nitrospirota", "Planctomycetota", "Pseudomonadota", "Nitrososphaerota", "Verrucomicrobiota")) +
  facet_grid(~Location, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# inocula taxonomic barplot
top10_16S_inocula <- ggplot(phylum_top10_16S[phylum_top10_16S$Type == "Inocula", ], aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Phylum", values = palette_MOD_soil_16s, labels = c("Actinomycetota", "Armatimonadota", "Bacteroidota", "Chloroflexota", "Bacillota", "Gemmatimonadota", "Nitrospirota", "Planctomycetota", "Pseudomonadota", "Nitrososphaerota", "Verrucomicrobiota")) +
  facet_grid(~Location, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# combine and save
tax_barplot_16S <- ggarrange(top10_16S_soil, top10_16S_inocula, nrow = 2, common.legend = TRUE, labels = "AUTO", legend = "right")
ggsave("./figures/FigureS5_tax_barplot_16S.svg", tax_barplot_16S, width = 10, height = 6)
ggsave("./figures/FigureS5_tax_barplot_16S.png", tax_barplot_16S, width = 10, height = 6)

## Supplementary Figure 3
phylo.ITS_phylum <- phyloseq::tax_glom(hmep.fungi, taxrank= "Class") # 25 taxa
# Make relative abundance 
phylo.ITS_phylum_relab <- transform_sample_counts(phylo.ITS_phylum, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
phylum_relab_top10_ITS <- fantaxtic::get_top_taxa(physeq_obj = phylo.ITS_phylum_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
phylum_top10_ITS <- psmelt(phylum_relab_top10_ITS)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
phylum_top10_ITS$Class <- forcats::fct_relevel(as.factor(phylum_top10_ITS$Class), "Other", after = Inf)
sapply(phylum_top10_ITS, levels) #factor levels

palette_MOD_soil_ITS <- c(
  "Eurotiomycetes" = "#88CCEE",
  "Dothideomycetes" = "#CC6677",
  "Mortierellomycetes" = "#44AA99",
  "Sordariomycetes" = "#117733",
  "Malasseziomycetes" = "#AA4499", 
  "Leotiomycetes" = "#DDCC77",
  "Agaricomycetes" = "#843c35", 
  "Tremellomycetes" = "#6948b8",
  "unidentified" = "#6699CC",
  "Other" = "#882255")

# Fixing a couple things for cleaning plotting
# Order W -> E and rename the landuse terms to fit the paper
phylum_top10_ITS$Location <- factor(phylum_top10_ITS$Location, levels = c("SVR", "HAY", "TLI", "KNZ"))
# soil taxonomic barplot
top10_ITS_soil <- ggplot(phylum_top10_ITS[phylum_top10_ITS$Type == "Soil", ], aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Class", values = safe_colorblind_palette_MOD_FUN) +
  facet_grid(~Location, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# inocula taxonomic barplot
top10_ITS_inocula <- ggplot(phylum_top10_ITS[phylum_top10_ITS$Type == "Inocula", ], aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Class", values = safe_colorblind_palette_MOD_FUN) +
  facet_grid(~Location, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# combine and save
tax_barplot_ITS <- ggarrange(top10_ITS_soil, top10_ITS_inocula, nrow = 2, common.legend = TRUE, labels = "AUTO", legend = "right")
ggsave("./figures/FigureS6_tax_barplot_ITS.svg", tax_barplot_ITS, width = 10, height = 6)
ggsave("./figures/FigureS6_tax_barplot_ITS.png", tax_barplot_ITS, width = 10, height = 6)

#### Alpha Diversity Statistics ####
## Bacteria
#extract metadata from clr transformed phyloseq object
bact.clr.sampledata <- data.frame(sample_data(hmep.bact.clr))
#linear models
mdsi.16S.z.richness <- lm(log(Chao1) ~ Type+Location+logObs.z, data = bact.clr.sampledata)
mdsi.16S.shannon <- lm(log(Shannon) ~ Type+Location+logObs.z, data = bact.clr.sampledata)
mdsi.16S.InvSimp <- lm(log(InvSimpson) ~ Type+Location+logObs.z, data = bact.clr.sampledata)
#anova
anova(mdsi.16S.z.richness)
"            Df Sum Sq Mean Sq  F value    Pr(>F)    
  Type       1  1.118   1.118   2.9384  0.094642 .  
  Location   3  8.358   2.786   7.3240  0.000543 ***
  logObs.z   1 69.734  69.734 183.3238 4.057e-16 ***
  Residuals 38 14.455   0.380 "
anova(mdsi.16S.shannon)
"            Df Sum Sq Mean Sq F value    Pr(>F)    
  Type       1  2.170   2.170  4.0923   0.05016 .  
  Location   3  7.643   2.548  4.8035   0.00621 ** 
  logObs.z   1 50.986  50.986 96.1356 5.887e-12 ***
  Residuals 38 20.153   0.530"
anova(mdsi.16S.InvSimp)
"            Df Sum Sq Mean Sq F value    Pr(>F)    
  Type       1  2.012   2.012  3.8251   0.05787 .  
  Location   3  6.500   2.167  4.1179   0.01267 *  
  logObs.z   1 32.086  32.086 60.9866 2.031e-09 ***
  Residuals 38 19.992   0.526    "
#plot residuals
plot(resid(mdsi.16S.z.richness)~fitted(mdsi.16S.z.richness))
plot(resid(mdsi.16S.shannon)~fitted(mdsi.16S.shannon))
plot(resid(mdsi.16S.InvSimp)~fitted(mdsi.16S.InvSimp))

#posthocs
summary(emmeans(mdsi.16S.z.richness,~Location), type = "response") # back transformed
pairs(emmeans(mdsi.16S.z.richness,~Location))
summary(emmeans(mdsi.16S.z.richness,~Type), type = "response") # back transformed



summary(emmeans(mdsi.16S.shannon,~Location), type = "response") # back transformed
pairs(emmeans(mdsi.16S.shannon,~Location))
summary(emmeans(mdsi.16S.shannon,~Type), type = "response") # back transformed
pairs(emmeans(mdsi.16S.shannon,~Type))

summary(emmeans(mdsi.16S.InvSimp,~Location), type = "response") # back transformed
pairs(emmeans(mdsi.16S.InvSimp,~Location))
summary(emmeans(mdsi.16S.InvSimp,~Type), type = "response") # back transformed

# emmeans plots
a <- ggplot(plot(emmeans(mdsi.16S.z.richness,~Location, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Chao1") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = location_pallete, name = "Collection Site") +
  theme(legend.text.align = 0, axis.title.x=element_blank())
b <- ggplot(plot(emmeans(mdsi.16S.shannon,~Location, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Shannon Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = location_pallete, name = "Collection Site") +
  theme(legend.text.align = 0, axis.title.x=element_blank())
c <- ggplot(plot(emmeans(mdsi.16S.InvSimp,~Location, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Inverse Simpson's Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = location_pallete, name = "Collection Site") +
  theme(legend.text.align = 0, axis.title.x=element_blank())
emm_16s_combo <- ggarrange(a,b,c, nrow = 1, labels = "AUTO", common.legend = TRUE, legend = "right")
ggsave("./figures/FigureS3_alpha.div.emm.16s.svg", emm_16s_combo,  width = 12, height = 6)
ggsave("./figures/FigureS3_alpha.div.emm.16s.png", emm_16s_combo,  width = 12, height = 6)

## Fungi 
#extract metadata from clr transformed phyloseq object
fungi.clr.sampledata <- data.frame(sample_data(hmep.fungi.clr))
#linear models
mdsi.ITS.z.richness <- lm(log(Chao1) ~ Type+Location+logObs.z, data = fungi.clr.sampledata)
mdsi.ITS.shannon <- lm(log(Shannon) ~ Type+Location+logObs.z, data = fungi.clr.sampledata)
mdsi.ITS.InvSimp <- lm(log(InvSimpson) ~ Type+Location+logObs.z, data = fungi.clr.sampledata)
#anova
anova(mdsi.ITS.z.richness)
"          Df  Sum Sq Mean Sq F value    Pr(>F)    
Type       1  8.5825  8.5825 36.2722 1.479e-07 ***
Location   3  0.8439  0.2813  1.1889    0.3225    
logObs.z   1 10.8864 10.8864 46.0089 8.549e-09 ***
Residuals 55 13.0138  0.2366"   
anova(mdsi.ITS.shannon)
"          Df  Sum Sq Mean Sq F value Pr(>F)
Type       1  0.0230 0.02298  0.0638 0.8015
Location   3  2.2465 0.74883  2.0807 0.1133
logObs.z   1  0.0279 0.02789  0.0775 0.7818
Residuals 55 19.7942 0.35989"  
anova(mdsi.ITS.InvSimp)
"          Df  Sum Sq Mean Sq F value  Pr(>F)  
Type       1  0.8457 0.84573  1.8699 0.17705  
Location   3  3.4649 1.15496  2.5537 0.06474 .
logObs.z   1  2.4060 2.40596  5.3196 0.02488 *
Residuals 55 24.8753 0.45228" 

#plot residuals
plot(resid(mdsi.ITS.z.richness)~fitted(mdsi.ITS.z.richness))
plot(resid(mdsi.ITS.shannon)~fitted(mdsi.ITS.shannon))
plot(resid(mdsi.ITS.InvSimp)~fitted(mdsi.ITS.InvSimp))

#posthocs
summary(emmeans(mdsi.ITS.z.richness,~Type), type = "response") # back transformed
pairs(emmeans(mdsi.ITS.z.richness,~Type))

# emmeans plots
a <- ggplot(plot(emmeans(mdsi.ITS.z.richness,~Type, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Chao1") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = c("red", "blue"), name = "Substrate") +
  scale_x_discrete(limits = rev) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')
ggsave("./figures/FigureS4_alpha.div.emm.ITS.svg", a,  width = 4, height = 4)
ggsave("./figures/FigureS4_alpha.div.emm.ITS.png", a,  width = 4, height = 4)
#clean env


#### PERMANOVA ####
## Bacteria
set.seed(77764413)
perm_16S <- with(as(sample_data(hmep.bact.clr),'data.frame'),
             adonis2(t(as(otu_table(hmep.bact.clr),'matrix')) ~ Location + Type,
             data=as(sample_data(hmep.bact.clr),'data.frame'), method='euclidean'))
perm_16S
"         Df SumOfSqs      R2      F Pr(>F)    
Location  3    68034 0.21987 4.4117  0.001 ***
Type      1    40923 0.13225 7.9610  0.001 ***
Residual 39   200477 0.64788                  
Total    43   309433 1.00000 "
bact.dbRDA <- ordinate(hmep.bact.clr, distance = 'euclidean', method = 'RDA', formula=~Location + Type + Condition(logObs))
bact.dbRDA
anova.cca(bact.dbRDA)


## Fungi
set.seed(77764413)
perm_ITS <- with(as(sample_data(hmep.fungi.clr),'data.frame'),
             adonis2(t(as(otu_table(hmep.fungi.clr),'matrix')) ~ Location + Type,
             data=as(sample_data(hmep.fungi.clr),'data.frame'), method='euclidean'))
perm_ITS
"         Df SumOfSqs      R2      F Pr(>F)    
Location  3    71904 0.23687 6.6441  0.001 ***
Type      1    29645 0.09766 8.2178  0.001 ***
Residual 56   202016 0.66548                  
Total    60   303565 1.00000"
drt.fung.dbRDA <- ordinate(hmep.fungi.clr, distance = 'euclidean', method = 'RDA', formula=~Location+ Type+ Condition(logObs))
summary(drt.fung.dbRDA)