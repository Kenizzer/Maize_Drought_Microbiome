# Analysis of maize B73 and Mo17 inoculated plants under drought stress
# Samples were collected from a greenhouse study in Feb 2020
# Microbiomes are from the root compartment (rhizosphere/endosphere)
# Growth measurements were also recorded throughout the experiment (~50 days in length)
# Code by: Maggie Wagner, Matthew Kolp, Joel Swift

#### Design of experiment ####
# 2 genotypes: B73 / Mo17
# 2 treatments: Well watered / Drought
# 6 soil inoculates (Soil location:habitat):
## Smoky Valley Ranch - SVR (Grain field and Prairie)
## Hays Prairie - HAY (Prairie)
## The Land Institute - TLI (Grain field and Prairie)
## Konza Prairie - KZ (Prairie)
###

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(phyloseq); packageVersion('phyloseq')
library(vegan); packageVersion('vegan')
library(ALDEx2); packageVersion('ALDEx2')
library(lme4); packageVersion('lme4')
library(lmerTest); packageVersion('lmerTest')
library(emmeans); packageVersion('emmeans')
library(ggpubr); packageVersion('ggpubr')

# Theme set and Color Palettes
theme_set(theme_pubr())
genotype_pallete <- c("B73" = "#91ff26", "Mo17" = "#9426ff")# B73/Mo17 - Genotype
treatment_pallete <- c("W" = "#0000FF", "D" = "#DAA520") # Drought/WW - Treatment
habitat_pallete <- c("Agriculture" = "#332288", "Native" = "#44AA99") # Native/Ag.- Soil habitat
location_pallete <- c("SVR" = "#88CCEE", "HAY" = "#CC6677", "TLI" = "#DDCC77", "KNZ" = "#117733") # SVR/HAY/TLI/KNZ - Soil location
soil_merged_pallete<- c("SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") # Combination of Soil habitat and location, colors from https://lospec.com/palette-list/zencillo14

taxa_palette <- c("#88CCEE", "#CC6677",
                  "#DDCC77", "#AA4499",
                  "#999933", "#44AA99",
                  "#D55E00", "#117733",
                  "#6699CC", "#882255",
                  "#000000", "#6948b8",
                  "#888888" ) 

### Load datasets for use throughout ###
drt.fungi.late <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late.RDS')
drt.fungi.late.clr <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late_clr.RDS')


### Late Timepoint ###
#### Alpha diversity ####
# Split out a metadata df for use outside of phyloseq
sampledata <- sample_data(drt.fungi.late.clr)
# LMs w/ genotype included in the model
drt.fungi.late.z.richness <- lmerTest::lmer(Chao1  ~ Genotype + SoilInoculum * Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.fungi.late.Shannon    <- lmerTest::lmer(log(Shannon)  ~ Genotype + SoilInoculum * Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.fungi.late.InvSimpson <- lmerTest::lmer(log(InvSimpson)  ~ Genotype + SoilInoculum * Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
# Plot residuals vs fitted
plot(resid(drt.fungi.late.z.richness),data.frame(sampledata)$Chao1)
plot(resid(drt.fungi.late.Shannon),data.frame(sampledata)$Shannon)
plot(resid(drt.fungi.late.InvSimpson),data.frame(sampledata)$InvSimpson)
# Homogeneity of Variance (~equal above/below line)
plot(drt.fungi.late.z.richness)
plot(drt.fungi.late.Shannon)
plot(drt.fungi.late.InvSimpson)
# qq plots
lattice::qqmath(drt.fungi.late.z.richness, id=0.05)
lattice::qqmath(drt.fungi.late.Shannon, id=0.05)
lattice::qqmath(drt.fungi.late.InvSimpson, id=0.05)
#ANOVAs on richness, diversity, and evenness
anova(drt.fungi.late.z.richness)
anova(drt.fungi.late.Shannon) # interaction sig
anova(drt.fungi.late.InvSimpson) # interaction sig
# check the random effects
rand(drt.fungi.late.z.richness)
rand(drt.fungi.late.Shannon) # Plate significant
rand(drt.fungi.late.InvSimpson) # Plate significant 
# Posthocs
pairs(emmeans(drt.fungi.late.Shannon, ~ Drought.or.Watered|SoilInoculum))
pairs(emmeans(drt.fungi.late.InvSimpson, ~ Drought.or.Watered|SoilInoculum))

#### Examining taxonomy and changes across factors ####
# Stats of fungal taxonomy (phylum and class) 
table(phyloseq::tax_table(drt.fungi.late.clr)[, "Phylum"]) #total community
#Ascomycota     Basidiomycota Mortierellomycota      Mucoromycota      unidentified 
#69                 9                 4                 1                 5 

table(phyloseq::tax_table(drt.fungi.late.clr)[, "Class"]) #total community
# Dothideomycetes     Eurotiomycetes      Leotiomycetes  Malasseziomycetes Microbotryomycetes Mortierellomycetes 
# 27                 19                  3                  5                  3                  4 
# Mucoromycetes      Pezizomycetes    Sordariomycetes    Tremellomycetes       unidentified 
# 1                  1                 17                  1                  7 

# Which Genera are impacted by the main effects
drt.fungi.late.clr_gen <- psmelt(phyloseq::tax_glom(drt.fungi.late.clr, "Genus")) # 84 taxa
# LMs
fungi.clr.long.aovs_gen <- drt.fungi.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lmer(Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data=data)))) %>%
  summarize(broom::tidy(mod))
# get vector by term for fungal genera, correct pvalues with BH, extract the significant ones
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilInoculum:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilInoculum:Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilInoculum", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilInoculum", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype", ][c(temp), ]
# Conclusion 
# Genera significant for drought: None
# Genera significant for Interaction (SIxDT): None
# Genera significant for inoculum: Geomyces, Guehomyces, Mucor, Penicillium, Pseudogymnoascus
# Genera significant for Genotype: None


#### Taxonomic barplots ####
# Use non-transformed data for relative abundance taxonomic barplots
drt.fungi.late_phylum <- phyloseq::tax_glom(drt.fungi.late, "Phylum") # 5 taxa
drt.fungi.late_class <- phyloseq::tax_glom(drt.fungi.late, "Class") # 12 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
drt.fungi.late_phylum_relab <- transform_sample_counts(drt.fungi.late_phylum, function(x) x/sum(x))
drt.fungi.late_class_relab <- transform_sample_counts(drt.fungi.late_class, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
drt.fungi.late_phylum_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.fungi.late_phylum_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Low Abundance Taxa")
drt.fungi.late_class_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.fungi.late_class_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Low Abundance Taxa")
# Melt data frame with phyloseq function for plotting.
phylum_top10 <- psmelt(drt.fungi.late_phylum_relab_top10)
class_top10 <- psmelt(drt.fungi.late_class_relab_top10)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
phylum_top10$Phylum <- forcats::fct_relevel(as.factor(phylum_top10$Phylum), "Low Abundance Taxa", after = Inf)
class_top10$Class <- forcats::fct_relevel(as.factor(class_top10$Class), "Low Abundance Taxa", after = Inf)
# Relative abundance stats
# Phyla
levels(as.factor(phylum_top10$Phylum))
mean((phylum_top10[phylum_top10$Phylum == "Ascomycota",]$Abundance)*100)         #92.07%
mean((phylum_top10[phylum_top10$Phylum == "Basidiomycota",]$Abundance)*100)      #1.97%
mean((phylum_top10[phylum_top10$Phylum == "Mortierellomycota",]$Abundance)*100)  #0.21%
mean((phylum_top10[phylum_top10$Phylum == "Mucoromycota",]$Abundance)*100)       #0.01%
mean((phylum_top10[phylum_top10$Phylum == "unidentified",]$Abundance)*100)       #5.75%
# Class
levels(as.factor(class_top10$Class))
mean((class_top10[class_top10$Class == "Dothideomycetes",]$Abundance)*100)         #14.20%
mean((class_top10[class_top10$Class == "Eurotiomycetes",]$Abundance)*100)          #53.11%
mean((class_top10[class_top10$Class == "Leotiomycetes",]$Abundance)*100)           #0.21%
mean((class_top10[class_top10$Class == "Malasseziomycetes",]$Abundance)*100)       #1.22%
mean((class_top10[class_top10$Class == "Microbotryomycetes",]$Abundance)*100)      #0.66%
mean((class_top10[class_top10$Class == "Mortierellomycetes",]$Abundance)*100)      #0.21%
mean((class_top10[class_top10$Class == "Pezizomycetes",]$Abundance)*100)           #0.37%
mean((class_top10[class_top10$Class == "Sordariomycetes",]$Abundance)*100)         #20.65%
mean((class_top10[class_top10$Class == "unidentified",]$Abundance)*100)            #4.64%
mean((class_top10[class_top10$Class == "Other",]$Abundance)*100)                   #0.10%

## Plots
# phylum
# By genotype
a <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Genotype, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By treatment
b <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By soil location
c <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilLocation, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By soil habitat
d <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilHabitat, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By treatment by soil location
e <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilLocation + Drought.or.Watered, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By treatment by soil location
f <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilInoculum, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

a
ggsave("figures/barplot_genotype.svg", a, height = 8, width = 12)
b
ggsave("figures/barplot_treatment.svg", b, height = 8, width = 12)
c
ggsave("figures/barplot_location.svg", c, height = 8, width = 12)
d
ggsave("figures/barplot_habitat.svg", d, height = 8, width = 12)
e
ggsave("figures/barplot_loaction_treatment.svg", e, height = 8, width = 12)
f
ggsave("figures/barplot_inoculum.svg", e, height = 8, width = 12)

rm(a,b,c,d,e,f)

# Class
# By genotype
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Genotype, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By treatment
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By soil location
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilLocation, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By soil habitat
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilHabitat, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By treatment by soil location
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilLocation + Drought.or.Watered, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By Inoculum
ggplot(class_top10, aes(x=Sample, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilInoculum, scale = "free", ncol = 2) +
  scale_fill_manual(values=taxa_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

#### Adonis and Ordinations ####
####### Use perMANOVA to partition variance in community composition at the ASV level
set.seed(8945761)

# First assess the marginal effects of the terms w/ interaction 
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype + SoilInoculum*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
# interaction is non-significant, so remove from model and re-run to assess main effects.
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype + SoilInoculum + Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
# No difference in P values assessing by margin

# Present table with terms assessed sequential (okay, given the results above).
perm <- with(as(sample_data(drt.fungi.late.clr),'data.frame'),
             adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype + SoilInoculum*Drought.or.Watered,
                     strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "term"))

#                                                      Df SumOfSqs      R2      F Pr(>F)    
"
Genotype                          1      774 0.00749 1.0148  0.280    
SoilInoculum                      5     9152 0.08861 2.3999  0.001 ***
Drought.or.Watered                1      693 0.00670 0.9080  0.423    
SoilInoculum:Drought.or.Watered   5     4198 0.04064 1.1007  0.220    
Residual                        116    88478 0.85656                  
Total                           128   103295 1.00000"
#
#8.86% of variance in fungal community composition explained by soil inocula.
#2nd largest effect size (r2): SoilInoculum:Drought.or.Watered, 4.06% of variance in fungal community composition.


# Constrained ordinations | CAP euclidean
# Start by constraining on largest factors of the anova
# by Soil Inoculum
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum+Condition(Plate + logObs))
summary(drt.fungi.dbRDA)
#Color and shape of plotted points by Soil Inoculum
a <- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'SoilInoculum') + geom_point(size=3) + scale_color_manual(values = soil_merged_pallete) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilInoculum)) + scale_fill_manual(values = soil_merged_pallete)
anova.cca(drt.fungi.dbRDA) 
(64.77 / (64.77 + 604.02)) * 100 #9.68% of community composition variation explained by model

# by drought
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
summary(drt.fungi.dbRDA)
#Color and shape of plotted points by drought treatment
b<- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'Drought.or.Watered') + geom_point(size=3) + scale_color_manual(values = treatment_pallete, name = "Treatment")
anova.cca(drt.fungi.dbRDA) 
(6.45 / (6.45 + 662.33))*100 #0.96% of community composition variation explained by model

# by genotype
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~Genotype+Condition(Plate + logObs))
summary(drt.fungi.dbRDA)
#Color and shape of plotted points genotype
c <- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'Genotype') + geom_point(size=3) + scale_color_manual(values = genotype_pallete)
anova.cca(drt.fungi.dbRDA) 
(3.79 / (3.79 + 665.12))*100 #0.57% of community composition variation explained by model

# by Soil Inoculum and treatment
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum + Drought.or.Watered + Condition(Plate + logObs))
summary(drt.fungi.dbRDA)
#Color and shape of plotted points by soil inocula source location and drought treatment
d <- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'SoilInoculum', shape = "Drought.or.Watered") + geom_point(size=3) + scale_color_manual(values = soil_merged_pallete, name = "Soil Inoculum") + guides(shape=guide_legend(title="Treatment"))
anova.cca(drt.fungi.dbRDA) 
(43.91 / (43.91 + 624.99))*100 #6.56% of community composition variation explained by model

# plots controlled by plate and logObs and constrained by each factor
# main effects
a 
b # Used in figure 3
c 
plot1 <- ggarrange(c, align = 'h', legend = 'right')
ggsave("figures/main_effect_RDA_ordination.svg", plot1, height = 6, width = 6)
ggsave("figures/main_effect_RDA_ordination.png", plot1, height = 6, width = 6)
# 2 ways
d
plot2 <- ggarrange(d, align = 'h', labels = "AUTO")
ggsave("figures/two_way_RDA_ordination.svg", plot2, height = 8, width = 24)
ggsave("figures/two_way_RDA_ordination.png", plot2, height = 8, width = 24)


#### betadisper ####
# test for difference in variance among groups
# "homogeneity of dispersion among sample groups was assessed using the betadisper function" https://microbiomejournal.biomedcentral.com/track/pdf/10.1186/s40168-017-0370-7.pdf
Betadisp.fun <- function(phyloseq_obj, Grouping, plot = TRUE){
  set.seed(81231)
  ## Calculate Distance
  dis <- vegdist(t(as.matrix(phyloseq_obj@otu_table)), method = "euclidean")
  ## Groups
  groups <- factor(phyloseq_obj@sam_data[[Grouping]])
  ## Calculate multivariate dispersions
  mod <- betadisper(dis, groups)
  if(plot == TRUE) {
    return(plot(mod))
  } else{
    ## Perform test
    #return(anova(mod))
    perm<-how(nperm=999)
    setBlocks(perm)<-with(as(phyloseq_obj@sam_data, 'data.frame'), Plate)
    return(permutest(mod, pairwise = TRUE, permutations = perm))
  }
}

# Genotype
Betadisp.fun(drt.fungi.late.clr, 'Genotype', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'Genotype', plot = FALSE) #NS
# drought
Betadisp.fun(drt.fungi.late.clr, 'Drought.or.Watered', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'Drought.or.Watered', plot = FALSE) #NS
# soil location
Betadisp.fun(drt.fungi.late.clr, 'SoilInoculum', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'SoilInoculum', plot = FALSE) # TLI_Agriculture and TLI_Native sig 0.024 padj permutation
# Plate
Betadisp.fun(drt.fungi.late.clr, 'Plate', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'Plate', plot = FALSE)

#### Aldex (differential abundance analysis) ####
# extract OTU table matrix from phyloseq object
asvs <- as(otu_table(drt.fungi.late, taxa_are_rows = TRUE), 'matrix')
asvs <- t(asvs)
# Make vectors of explanatory variables to test <2 groups
Drought.or.Watered <- sample_data(drt.fungi.late)$Drought.or.Watered %>% as.character()
Genotype <- sample_data(drt.fungi.late)$Genotype %>% as.character()
# Transform counts with CLR
asvs.clr.drought <- aldex.clr(asvs, conds= Drought.or.Watered, mc.samples=1000)
asvs.clr.genotype <- aldex.clr(asvs, conds= Genotype, mc.samples=1000)

# T-test with factors with 2 levels
#Outputs a dataframe with the following information:
#wi.ep   a vector containing the expected P value of the Wilcoxon test for each feature
#wi.eBH  a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature

# Drought
ttest.drought <- aldex.ttest(asvs.clr.drought)
sum(ttest.drought$wi.eBH < 0.05) # 0 significant ASV drought treatment
# Calculate effect sizes
effect.drought <- aldex.effect(asvs.clr.drought)
ttest.effect.drought <- data.frame(ttest.drought, effect.drought)
aldex.plot(ttest.effect.drought, type='MA')
aldex.plot(ttest.effect.drought, type='MW') 

# Genotype
ttest.genotype <- aldex.ttest(asvs.clr.genotype)
sum(ttest.genotype$wi.eBH < 0.05) # 0 significant genotype
# Calculate effect sizes
effect.genotype <- aldex.effect(asvs.clr.genotype)
# plot results
ttest.effect.genotype <- data.frame(ttest.genotype, effect.genotype)
aldex.plot(ttest.effect.genotype, type='MA', test = 'wilcox')
aldex.plot(ttest.effect.genotype, type='MW', test = 'wilcox')

# more than two levels - Soil Location
SoilInoculum.df <- data.frame(SoilInoculum = sample_data(drt.fungi.late)$SoilInoculum %>% as.character())
# Setup model matrix (required for factors with >2 levels)
mm <- model.matrix(~ SoilInoculum, SoilInoculum.df)
SoilInoculum.clr <- aldex.clr(asvs, mm, mc.samples=128)
SoilInoculum.glm.test <- aldex.glm(SoilInoculum.clr)
# How many ASV differ between Hay vs others
colnames(SoilInoculum.glm.test)
sum(SoilInoculum.glm.test$`model.SoilInoculumKNZ_Native Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumSVR_Agriculture Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumSVR_Native Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumTLI_Agriculture Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumTLI_Native Pr(>|t|).BH` < 0.05)

which(SoilInoculum.glm.test$`model.SoilInoculumKNZ_Native Pr(>|t|).BH` < 0.05)
which(SoilInoculum.glm.test$`model.SoilInoculumSVR_Agriculture Pr(>|t|).BH` < 0.05)
which(SoilInoculum.glm.test$`model.SoilInoculumSVR_Native Pr(>|t|).BH` < 0.05)
which(SoilInoculum.glm.test$`model.SoilInoculumTLI_Agriculture Pr(>|t|).BH` < 0.05)
which(SoilInoculum.glm.test$`model.SoilInoculumTLI_Native Pr(>|t|).BH` < 0.05)

SoilInoculum.glm.test[9,] 
SoilInoculum.glm.test[6,]
SoilInoculum.glm.test[12,]
# no effect sizes for greater than 2 condition levels

#        Kingdom Phylum       Class            Order        Family           Genus        Species
#fASV_9  "Fungi" "Ascomycota" "Eurotiomycetes" "Eurotiales" "Aspergillaceae" "Penicillium" "Penicillium_subrubescens"
#fASV_6  "Fungi" "Ascomycota" "Eurotiomycetes" "Eurotiales" "Aspergillaceae" "Penicillium" "Penicillium_annulatum"
#fASV_12 "Fungi" "Ascomycota" "Eurotiomycetes" "Eurotiales" "Aspergillaceae" "Penicillium" "Penicillium_raperi"


#### Figure 5 panels A,B,C #####
# reorder the class dataframe in order to samples in order of the largest class
Class_top10_ordered <- class_top10[order(class_top10$Sample),]
# get a vector that will tell which samples have the highest abundance of Eurotiomycetes
Class_top10_ordered$plot_order <- rep(rank(-Class_top10_ordered$Abundance[Class_top10_ordered$Class == "Eurotiomycetes"], ties.method = 'first'), each = 11)
Class_top10_ordered$plot_order <- factor(Class_top10_ordered$plot_order)
Class_top10_ordered$Class <- factor(Class_top10_ordered$Class, levels = c("Eurotiomycetes", "Dothideomycetes", "Leotiomycetes", "Malasseziomycetes", "Microbotryomycetes", "Mortierellomycetes", "Pezizomycetes", "Sordariomycetes", "unidentified", "Other"))

# Fixing labeling and order issues for facets
Class_top10_ordered <- Class_top10_ordered %>%
  mutate(Drought.or.Watered = fct_recode(as.factor(Drought.or.Watered), 
                                         Drought = "D",
                                         "Well-Watered" = "W"))
Class_top10_ordered <- Class_top10_ordered %>%
  mutate(SoilInoculum = fct_recode(as.factor(SoilInoculum), 
                                   `HAY[P]` = "HAY_Native",
                                   `KNZ[P]` = "KNZ_Native",
                                   `SVR[Ag]` = "SVR_Agriculture",
                                   `SVR[P]` = "SVR_Native",
                                   `TLI[Ag]` = "TLI_Agriculture",
                                   `TLI[P]` = "TLI_Native"))
Class_top10_ordered$SoilInoculum <- factor(Class_top10_ordered$SoilInoculum, levels = c("SVR[Ag]", "SVR[P]", "HAY[P]", "TLI[Ag]", "TLI[P]", "KNZ[P]"))

a <- ggplot(Class_top10_ordered, aes(x=plot_order, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered + SoilInoculum, scale = "free_x", nrow = 2, labeller = label_parsed) +
  scale_fill_manual(values=taxa_palette) + guides(fill=guide_legend(nrow=5, byrow=TRUE))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top', 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=8))

# Pcoa treatment Panel E and RF CM treatment Panel F
# by treatment
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
anova.cca(drt.fungi.dbRDA)
#Color and shape of plotted points by soil habitat
b <- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'Drought.or.Watered') + geom_point(size=3) + scale_color_manual(values = treatment_pallete, name = "Treatment", labels = c("Drought (D)", "Well-Watered (W)"))
b <- b + annotate("text", x = 3, y = 5, label = "Treatment p = 0.18", fontface = 'italic')
b <- b + guides(color = guide_legend(nrow = 2))
# Random forest prediction for treatment
c <- readRDS("./Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
c <- c$CMatrixPLOT
part <- ggarrange(b,c, ncol = 1, heights = c(1,0.85), labels = c("B","C"))
fig5_def <- ggarrange(a ,part, ncol = 2, widths = c(1.3,0.55), labels = c("A"))
ggsave("figures/figure5_panels_abc.svg", fig3_def, height = 6, width = 10)
ggsave("figures/figure5_panels_abc.png", fig3_def, height = 6, width = 10)


# by soil inoculum, was cut from the figure above will now be a supplement
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum+Condition(Plate + logObs))
#Color and shape of plotted points by soil inocula source location and drought treatment
soil_inoc_RDA <- plot_ordination(drt.fungi.late.clr, drt.fungi.dbRDA, color = 'SoilInoculum') + geom_point(size=2) + scale_color_manual(values = soil_merged_pallete, name = "Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilInoculum)) + scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  guides(color=guide_legend(nrow=6, byrow=TRUE), fill=guide_legend(nrow=6, byrow=TRUE)) + theme(legend.text.align = 0, legend.position = 'right')

ggsave("figures/Soil_inoculum_RDA.svg", soil_inoc_RDA, height = 6, width = 6)
ggsave("figures/Soil_inoculum_RDA.png", soil_inoc_RDA, height = 6, width = 6)