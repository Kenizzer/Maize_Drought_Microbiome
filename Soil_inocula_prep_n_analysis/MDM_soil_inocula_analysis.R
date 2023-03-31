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
soil_merged_pallete<- c("SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") # Combination of Soil habitat and location, colors from https://lospec.com/palette-list/zencillo14
theme_set(theme_pubr())

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
sample_data(hmep.bact.clr)$SoilSource <- factor(sample_data(hmep.bact.clr)$SoilSource, levels = c("SVR_Agriculture", "SVR_Native", "HAY_Native", "TLI_Agriculture", "TLI_Native", "KNZ_Native"))
#Ordinate (constrained PCOA)
RDA.ait.hmep.bact.clr <- ordinate(hmep.bact.clr, method='RDA',distance='euclidean',formula=~SoilSource + Type + Condition(logObs))
#Color by Type (soil vs. inocula)
RDA_16S <- plot_ordination(hmep.bact.clr, RDA.ait.hmep.bact.clr)
RDA_16S[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- 'white'
RDA_16S <- RDA_16S + 
  geom_point(aes(fill = SoilSource, shape = Type), size = 3, alpha = 0.80, color = "black") +
  scale_shape_manual(values=c(24, 21), name = "Substrate") +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Source", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P])))+
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))

#fungi
sample_data(hmep.fungi.clr)$SoilSource <- factor(sample_data(hmep.fungi.clr)$SoilSource, levels = c("SVR_Agriculture", "SVR_Native", "HAY_Native", "TLI_Agriculture", "TLI_Native", "KNZ_Native"))
#Ordinate 
RDA.ait.hmep.fungi.clr <- ordinate(hmep.fungi.clr, method = "RDA", distance = "euclidean", formula=~SoilSource + Type + Condition(logObs))
#Color 
RDA_ITS <- plot_ordination(hmep.fungi.clr, RDA.ait.hmep.fungi.clr)
RDA_ITS[["layers"]][[1]][["geom"]][["default_aes"]][["colour"]] <- 'white'
RDA_ITS <- RDA_ITS + 
  geom_point(aes(fill = SoilSource, shape = Type), size = 3, alpha = 0.80, color = "black") +
  scale_shape_manual(values=c(24, 21), name = "Substrate") +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Source", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P])))+
  theme(legend.text.align = 0) +
  guides(fill = guide_legend(override.aes = c(shape = 21, alpha = 1)), color = guide_legend(override.aes = c(alpha=1)))

#Combine plots
RDA_16S_ITS <- ggarrange(RDA_16S, RDA_ITS, labels = "AUTO", common.legend = TRUE, legend = "right")
RDA_16S_ITS
ggsave("./figures/RDA_16S_ITS_F2.svg", RDA_16S_ITS, width = 6, height = 4)
ggsave("./figures/RDA_16S_ITS_F2.png", RDA_16S_ITS, width = 6, height = 4)


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
phylum_top10_16S$Land_use <- gsub("Native", "Prairie", phylum_top10_16S$Land_use)

# soil taxonomic barplot
top10_16S_soil <- ggplot(phylum_top10_16S[phylum_top10_16S$Type == "Soil", ], aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Phylum", values = palette_MOD_soil_16s, labels = c("Actinomycetota", "Armatimonadota", "Bacteroidota", "Chloroflexota", "Bacillota", "Gemmatimonadota", "Nitrospirota", "Planctomycetota", "Pseudomonadota", "Nitrososphaerota", "Verrucomicrobiota")) +
  facet_grid(~Location+Land_use, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# inocula taxonomic barplot
top10_16S_inocula <- ggplot(phylum_top10_16S[phylum_top10_16S$Type == "Inocula", ], aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Phylum", values = palette_MOD_soil_16s, labels = c("Actinomycetota", "Armatimonadota", "Bacteroidota", "Chloroflexota", "Bacillota", "Gemmatimonadota", "Nitrospirota", "Planctomycetota", "Pseudomonadota", "Nitrososphaerota", "Verrucomicrobiota")) +
  facet_grid(~Location+Land_use, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# combine and save
tax_barplot_16S <- ggarrange(top10_16S_soil, top10_16S_inocula, nrow = 2, common.legend = TRUE, labels = "AUTO", legend = "right")
ggsave("./figures/tax_barplot_16S.svg", tax_barplot_16S, width = 10, height = 6)
ggsave("./figures/tax_barplot_16S.png", tax_barplot_16S, width = 10, height = 6)

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

safe_colorblind_palette_MOD_FUN <- c(
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
phylum_top10_ITS$Land_use <- gsub("Native", "Prairie", phylum_top10_ITS$Land_use)

# soil taxonomic barplot
top10_ITS_soil <- ggplot(phylum_top10_ITS[phylum_top10_ITS$Type == "Soil", ], aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Class", values = safe_colorblind_palette_MOD_FUN) +
  facet_grid(~Location+Land_use, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# inocula taxonomic barplot
top10_ITS_inocula <- ggplot(phylum_top10_ITS[phylum_top10_ITS$Type == "Inocula", ], aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  scale_fill_manual(name = "Class", values = safe_colorblind_palette_MOD_FUN) +
  facet_grid(~Location+Land_use, scale = "free") +
  theme(legend.position = 'right') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# combine and save
tax_barplot_ITS <- ggarrange(top10_ITS_soil, top10_ITS_inocula, nrow = 2, common.legend = TRUE, labels = "AUTO", legend = "right")
ggsave("./figures/tax_barplot_ITS.svg", tax_barplot_ITS, width = 10, height = 6)
ggsave("./figures/tax_barplot_ITS.png", tax_barplot_ITS, width = 10, height = 6)



#### Alpha Diversity Statistics ####

## Bacteria
#extract metadata from clr transformed phyloseq object
bact.clr.sampledata <- data.frame(sample_data(hmep.bact.clr))
#linear models
mdsi.16S.z.richness <- lm(log(Chao1) ~ Type+SoilSource+logObs.z, data = bact.clr.sampledata)
mdsi.16S.shannon <- lm(log(Shannon) ~ Type+SoilSource+logObs.z, data = bact.clr.sampledata)
mdsi.16S.InvSimp <- lm(log(InvSimpson) ~ Type+SoilSource+logObs.z, data = bact.clr.sampledata)
#anova
anova(mdsi.16S.z.richness)
# Df Sum Sq Mean Sq  F value    Pr(>F)    
# Type        1  0.569   0.569   1.1411 0.2903600    
# SoilSource  5 15.899   3.180   6.3715 0.0001101 ***
# logObs.z    1 81.019  81.019 162.3442 < 2.2e-16 ***
# Residuals  52 25.951   0.499  
anova(mdsi.16S.shannon)
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# Type        1  0.077   0.077  0.1125  0.738644    
# SoilSource  5 14.406   2.881  4.2367  0.002635 ** 
# logObs.z    1 58.198  58.198 85.5792 1.425e-12 ***
# Residuals  52 35.363   0.680 
anova(mdsi.16S.InvSimp)
#             Df Sum Sq Mean Sq F value    Pr(>F)    
# Type        1  0.006   0.006  0.0100   0.92070    
# SoilSource  5 10.666   2.133  3.2933   0.01171 *  
# logObs.z    1 36.957  36.957 57.0555 6.497e-10 ***
# Residuals  52 33.682   0.648

#plot residuals
plot(resid(mdsi.16S.z.richness)~fitted(mdsi.16S.z.richness))
plot(resid(mdsi.16S.shannon)~fitted(mdsi.16S.shannon))
plot(resid(mdsi.16S.InvSimp)~fitted(mdsi.16S.InvSimp))

#posthocs
summary(emmeans(mdsi.16S.z.richness,~SoilSource), type = "response") # back transformed
pairs(emmeans(mdsi.16S.z.richness,~SoilSource))

summary(emmeans(mdsi.16S.shannon,~SoilSource), type = "response") # back transformed
pairs(emmeans(mdsi.16S.shannon,~SoilSource))

summary(emmeans(mdsi.16S.InvSimp,~SoilSource), type = "response") # back transformed
pairs(emmeans(mdsi.16S.InvSimp,~SoilSource))

# emmeans plots
a <- ggplot(plot(emmeans(mdsi.16S.z.richness,~SoilSource, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Chao1") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = soil_merged_pallete, name = "Collection Site", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  scale_x_discrete(labels=c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())
b <- ggplot(plot(emmeans(mdsi.16S.shannon,~SoilSource, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Shannon Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = soil_merged_pallete, name = "Collection Site", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  scale_x_discrete(labels=c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())
c <- ggplot(plot(emmeans(mdsi.16S.InvSimp,~SoilSource, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Inverse Simpson's Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = soil_merged_pallete, name = "Collection Site", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  scale_x_discrete(labels=c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

emm_16s_combo <- ggarrange(a,b,c, nrow = 1, labels = "AUTO", common.legend = TRUE, legend = "right")
ggsave("./figures/alpha.div.emm.16s.svg", emm_16s_combo,  width = 12, height = 6)
ggsave("./figures/alpha.div.emm.16s.png", emm_16s_combo,  width = 12, height = 6)


## Fungi 
#extract metadata from clr transformed phyloseq object
fungi.clr.sampledata <- data.frame(sample_data(hmep.fungi.clr))
#linear models
mdsi.ITS.z.richness <- lm(log(Chao1) ~ Type+SoilSource+logObs.z, data = fungi.clr.sampledata)
mdsi.ITS.shannon <- lm(log(Shannon) ~ Type+SoilSource+logObs.z, data = fungi.clr.sampledata)
mdsi.ITS.InvSimp <- lm(log(InvSimpson) ~ Type+SoilSource+logObs.z, data = fungi.clr.sampledata)
#anova
anova(mdsi.ITS.z.richness)
#              Df  Sum Sq Mean Sq F value    Pr(>F)    
#  Type        1 10.4983 10.4983 45.0693 2.397e-09 ***
#  SoilSource  5  3.0722  0.6144  2.6378   0.02922 *  
#  logObs.z    1 14.6715 14.6715 62.9850 1.003e-11 ***
#  Residuals  81 18.8678  0.2329   
anova(mdsi.ITS.shannon)
#             Df  Sum Sq Mean Sq F value Pr(>F)
# Type        1  0.0247 0.02467  0.0809 0.7768
# SoilSource  5  2.6693 0.53386  1.7515 0.1323
# logObs.z    1  0.4488 0.44880  1.4725 0.2285
# Residuals  81 24.6884 0.30479  
anova(mdsi.ITS.InvSimp)
#             Df Sum Sq Mean Sq F value   Pr(>F)   
# Type        1  2.987  2.9875  7.1609 0.009015 **
# SoilSource  5  5.211  1.0423  2.4983 0.037218 * 
# logObs.z    1  1.282  1.2825  3.0741 0.083334 . 
# Residuals  81 33.793  0.4172 

#plot residuals
plot(resid(mdsi.ITS.z.richness)~fitted(mdsi.ITS.z.richness))
plot(resid(mdsi.ITS.shannon)~fitted(mdsi.ITS.shannon))
plot(resid(mdsi.ITS.InvSimp)~fitted(mdsi.ITS.InvSimp))

#posthocs
summary(emmeans(mdsi.ITS.z.richness,~Type), type = "response") # back transformed
pairs(emmeans(mdsi.ITS.z.richness,~Type))
summary(emmeans(mdsi.ITS.z.richness,~SoilSource), type = "response") # back transformed
pairs(emmeans(mdsi.ITS.z.richness,~SoilSource))

summary(emmeans(mdsi.ITS.InvSimp,~Type), type = "response") # back transformed
pairs(emmeans(mdsi.ITS.InvSimp,~Type))
summary(emmeans(mdsi.ITS.InvSimp,~SoilSource), type = "response") # back transformed
pairs(emmeans(mdsi.ITS.InvSimp,~SoilSource))

# emmeans plots
a <- ggplot(plot(emmeans(mdsi.ITS.z.richness,~SoilSource, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Chao1") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = soil_merged_pallete, name = "Collection Site", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  scale_x_discrete(labels=c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())
b <- ggplot(plot(emmeans(mdsi.ITS.z.richness,~Type, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Chao1") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = c("red", "blue"), name = "Substrate") +
  scale_x_discrete(limits = rev) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

c <- ggplot(plot(emmeans(mdsi.ITS.InvSimp,~SoilSource, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Inverse Simpson's Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = soil_merged_pallete, name = "Collection Site", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  scale_x_discrete(labels=c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())
d <- ggplot(plot(emmeans(mdsi.ITS.InvSimp,~Type, type = "response"), plotit = FALSE), aes(x = pri.fac, y = the.emmean, color = pri.fac)) +
  ylab("Inverse Simpson's Diversity") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), size = 0.75, linewidth = 2, fill = "black", shape = 22) +
  scale_color_manual(values = c("red", "blue"), name = "Substrate") +
  scale_x_discrete(limits = rev) +
  theme(legend.text.align = 0, axis.title.x=element_blank())


P1.1 <- ggarrange(a,c, common.legend = TRUE, legend = "right", labels = "AUTO")
P1.2 <- ggarrange(b,d, common.legend = TRUE, legend = "right", labels = c("C", "D"))
P1 <- ggarrange(P1.1, P1.2, nrow = 2, align = "hv")
ggsave("./figures/alpha.div.emm.ITS.svg", P1,  width = 12, height = 12)
ggsave("./figures/alpha.div.emm.ITS.png", P1,  width = 12, height = 12)
#clean env
rm(a,b,c,d,P1,P1.1,P1.2)


#### PERMANOVA ####
## Bacteria
set.seed(77764413)
perm_16S <- with(as(sample_data(hmep.bact.clr),'data.frame'),
             adonis2(t(as(otu_table(hmep.bact.clr),'matrix')) ~ SoilSource + Type,
             data=as(sample_data(hmep.bact.clr),'data.frame'), method='euclidean'))
perm_16S
#              Df SumOfSqs      R2      F Pr(>F)    
#  SoilSource  5   101664 0.22885 3.6164  0.001 ***
#  Type        1    44584 0.10036 7.9298  0.001 ***
#  Residual   53   297985 0.67078                  
#  Total      59   444234 1.00000 
bact.dbRDA <- ordinate(hmep.bact.clr, distance = 'euclidean', method = 'RDA', formula=~SoilSource + Type + Condition(logObs))
bact.dbRDA
anova.cca(bact.dbRDA)


## Fungi
set.seed(77764413)
perm_ITS <- with(as(sample_data(hmep.fungi.clr),'data.frame'),
             adonis2(t(as(otu_table(hmep.fungi.clr),'matrix')) ~ SoilSource + Type,
             data=as(sample_data(hmep.fungi.clr),'data.frame'), method='euclidean'))
perm_ITS
#              Df SumOfSqs      R2      F Pr(>F)    
#  SoilSource  5   101753 0.22065 5.0990  0.001 ***
#  Type        1    32134 0.06968 8.0515  0.001 ***
#  Residual   82   327269 0.70967                  
#  Total      88   461157 1.00000
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum+Condition(Plate + logObs))
summary(drt.bact.dbRDA)
