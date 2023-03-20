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
treatment_pallete <- c("W" = "#0000FF", "D" = "#DAA520") # Drought/WW     - Treatment
habitat_pallete <- c("Agriculture" = "#332288", "Native" = "#44AA99") # Native/Ag.     - Soil habitat
location_pallete <- c("SVR" = "#88CCEE", "HAY" = "#CC6677", "TLI" = "#DDCC77", "KNZ" = "#117733") # SVR/HAY/TLI/KNZ - Soil location
soil_merged_pallete<- c("SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") # Combination of Soil habitat and location, colors from https://lospec.com/palette-list/zencillo14

phyla_palette <- c("Actinobacteria" = "#88CCEE", "Armatimonadetes"= "#CC6677",
                             "Bacteroidetes" = "#DDCC77",  "Chloroflexi" = "#AA4499",
                             "Firmicutes" = "#999933", "Gemmatimonadetes" = "#44AA99",
                            "Nitrospirae" = "#D55E00", "Proteobacteria" = "#117733",
                            "Planctomycetes" =  "#6699CC", "Verrucomicrobia" = "#882255",
                            "Deinococcus-Thermus" = "#000000", "Thaumarchaeota" = "#6948b8", 
                            "Other" = "#888888")

### Load datasets for use throughout ###
drt.bact.late <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS')
drt.bact.late.clr <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS')

### Late Timepoint ###
#### Alpha diversity ####
# Split out a metadata df for use outside of phyloseq
sampledata <- sample_data(drt.bact.late.clr)
# LMs w/ genotype included in the model
drt.bact.late.z.richness <- lmerTest::lmer(Chao1  ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.bact.late.Shannon    <- lmerTest::lmer(log(Shannon)  ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.bact.late.InvSimpson <- lmerTest::lmer(log(InvSimpson)  ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
# Plot residuals vs fitted
plot(resid(drt.bact.late.z.richness),data.frame(sampledata)$Chao1)
plot(resid(drt.bact.late.Shannon),log(data.frame(sampledata)$Shannon))
plot(resid(drt.bact.late.InvSimpson),log(data.frame(sampledata)$InvSimpson))
# Homogeneity of Variance (~equal above/below line)
plot(drt.bact.late.z.richness)
plot(drt.bact.late.Shannon)
plot(drt.bact.late.InvSimpson)
# qq plots
lattice::qqmath(drt.bact.late.z.richness, id=0.05)
lattice::qqmath(drt.bact.late.Shannon, id=0.05)
lattice::qqmath(drt.bact.late.InvSimpson, id=0.05)
#ANOVAs on richness, diversity, and evenness
anova(drt.bact.late.z.richness)
anova(drt.bact.late.Shannon)
anova(drt.bact.late.InvSimpson)
# check the random effects
rand(drt.bact.late.z.richness) # Plate significant
rand(drt.bact.late.Shannon) # Plate significant
rand(drt.bact.late.InvSimpson) # Plate significant 
# posthoc
emmeans(drt.bact.late.z.richness,~Drought.or.Watered)
pairs(emmeans(drt.bact.late.z.richness,~Drought.or.Watered))

#### Examining taxonomy and changes of phyla across factors ####
#Descriptive stats of bacterial taxonomy (phylum and class)
table(phyloseq::tax_table(drt.bact.late.clr)[, "Phylum"]) #total community
#Actinobacteria           Armatimonadetes             Bacteroidetes 
#129                      1                           57 
#Chloroflexi        Cyanobacteria/Chloroplast       Deinococcus-Thermus 
#2                  1                               1 
#Firmicutes          Gemmatimonadetes               Nitrospirae 
#17                  4                              2 
#Planctomycetes            Proteobacteria            Thaumarchaeota 
#6                         387                       1 
#Verrucomicrobia 
#4 

table(phyloseq::tax_table(drt.bact.late.clr)[, "Class"]) #total community
#Actinobacteria Alphaproteobacteria        Anaerolineae        Bacilli        Betaproteobacteria 
#129                 132                   1                  15                 200 
#Clostridia       Cyanobacteria          Cytophagia          Deinococci     Deltaproteobacteria 
#2                   1                   6                   1                   5 
#Fimbriimonadia      Flavobacteriia Gammaproteobacteria    Gemmatimonadetes   Nitrososphaerales 
#1                   4                  50                   4                   1 
#Nitrospira         Opitutae      Planctomycetia    Sphingobacteriia      Thermomicrobia 
#2                   2                   6                  47                   1 
#Verrucomicrobiae 
#2 


# Which phyla show largest changes by factors
drt.bact.late.clr_phylum <- phyloseq::tax_glom(drt.bact.late.clr, "Phylum")
phylum_df <- psmelt(drt.bact.late.clr_phylum)
levels(as.factor(phylum_df$Phylum))
mod.act <- lmerTest::lmer(filter(phylum_df, Phylum == "Actinobacteria")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Actinobacteria"))
mod.arm <- lmerTest::lmer(filter(phylum_df, Phylum == "Armatimonadetes")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Armatimonadetes"))
mod.bac <- lmerTest::lmer(filter(phylum_df, Phylum == "Bacteroidetes")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Bacteroidetes"))
mod.chl <- lmerTest::lmer(filter(phylum_df, Phylum == "Chloroflexi")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Chloroflexi"))
mod.cya <- lmerTest::lmer(filter(phylum_df, Phylum == "Cyanobacteria/Chloroplast")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Cyanobacteria/Chloroplast"))
mod.dei <- lmerTest::lmer(filter(phylum_df, Phylum == "Deinococcus-Thermus")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Deinococcus-Thermus"))
mod.fir <- lmerTest::lmer(filter(phylum_df, Phylum == "Firmicutes")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Firmicutes"))
mod.gem <- lmerTest::lmer(filter(phylum_df, Phylum == "Gemmatimonadetes")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Gemmatimonadetes"))
mod.nit <- lmerTest::lmer(filter(phylum_df, Phylum == "Nitrospirae")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Nitrospirae"))
mod.pla <- lmerTest::lmer(filter(phylum_df, Phylum == "Planctomycetes")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Planctomycetes"))
mod.pro <- lmerTest::lmer(filter(phylum_df, Phylum == "Proteobacteria")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Proteobacteria"))
mod.tha <- lmerTest::lmer(filter(phylum_df, Phylum == "Thaumarchaeota")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Thaumarchaeota"))
mod.ver <- lmerTest::lmer(filter(phylum_df, Phylum == "Verrucomicrobia")$Abundance ~ Genotype + SoilInoculum*Drought.or.Watered + (1|Block) + (1|Plate), data = filter(phylum_df, Phylum == "Verrucomicrobia"))

#### Correct P-vales of phylum models
Return_SS_proportion <- function(model){
  # **Model must be from Lmer (can be change to fit models from lm)**
  SS <- c(anova(model)$Sum, sum(resid(model)^2))
  output <- (SS / sum(SS)) * 100
  output <- as.data.frame(output) 
  rownames(output) <- c(rownames(anova(model)), "Residuals")
  output <- t(output)
  return(as.data.frame(output))
}
# Make list of phyla in order and run the function to generate the % variance by model by factor
phyla_list <- levels(as.factor(phylum_df$Phylum))
Percent_var <- bind_rows(Return_SS_proportion(mod.act),
                         Return_SS_proportion(mod.arm),
                         Return_SS_proportion(mod.bac),
                         Return_SS_proportion(mod.chl),
                         Return_SS_proportion(mod.cya),
                         Return_SS_proportion(mod.dei),
                         Return_SS_proportion(mod.fir),
                         Return_SS_proportion(mod.gem),
                         Return_SS_proportion(mod.nit),
                         Return_SS_proportion(mod.pla),
                         Return_SS_proportion(mod.pro),
                         Return_SS_proportion(mod.tha),
                         Return_SS_proportion(mod.ver))
rownames(Percent_var) <- phyla_list
# remove residual column
Percent_var_noint_nores <- Percent_var[, !names(Percent_var) %in% "Residuals"]
# Make dataframe of P values to connect to the % var df above
Get_Pval<- function(model){
  Pval<- as.data.frame(anova(model)$`Pr(>F)`)
  rownames(Pval) <- paste(rownames(anova(model)), "p", sep = "_")
  Pval<- t(Pval)
  return(as.data.frame(Pval))
}
# Run function above on all anova and bind output by row
Pval <- bind_rows(Get_Pval(mod.act),
                  Get_Pval(mod.arm),
                  Get_Pval(mod.bac),
                  Get_Pval(mod.chl),
                  Get_Pval(mod.cya),
                  Get_Pval(mod.dei),
                  Get_Pval(mod.fir),
                  Get_Pval(mod.gem),
                  Get_Pval(mod.nit),
                  Get_Pval(mod.pla),
                  Get_Pval(mod.pro),
                  Get_Pval(mod.tha),
                  Get_Pval(mod.ver))
rownames(Pval) <- phyla_list
# Combine variance and P val dataframes
SS_pval_combo_df <- cbind(Percent_var_noint_nores,Pval)
SS_pval_combo_df$Phylum <- rownames(SS_pval_combo_df)
# Get variance columns
total_var <- SS_pval_combo_df %>% dplyr::select(Phylum, "Genotype":"SoilInoculum:Drought.or.Watered")
# Reorganize 
total_var <- total_var %>% gather(key=factor, value=var, -Phylum)
# Get p-value columns
total_p_name <- SS_pval_combo_df %>% dplyr::select("Genotype_p":"SoilInoculum:Drought.or.Watered_p")
total_p <- data.frame(t(apply(total_p_name, 1, FUN=p.adjust, method='BH')))
colnames(total_p)<- colnames(total_p_name) #This is hacky but I need to preserve the colnames to correctly join them later on.
total_p$Phylum <- SS_pval_combo_df$Phylum
# Reorganize and rename to format for plotting
total_p <- total_p %>% gather(key=factor, value=p_value, -Phylum)
total_p <- total_p %>% mutate(factor=str_replace(factor,'_p$', ""))
# Join variance and p-value tables back together 
total_var_p <- full_join(total_var, total_p,by=c("Phylum", "factor"))
total_var_p_sig <- total_var_p %>% filter(p_value < 0.05)
total_var_p_sig
#  Phylum                                       factor      var    p_value
#  1 Deinococcus-Thermus                    SoilInoculum 6.928455 0.02815631
#  2      Proteobacteria              Drought.or.Watered 3.680463 0.03599070
#  3 Deinococcus-Thermus SoilInoculum:Drought.or.Watered 7.119791 0.02815631

# Conclusion: 
# Deinococcus-Thermus
# ~6.93% of the variance in CLR abundance explained by Soil inoculum
# ~7.11% of the variance in CLR abundance explained by Soil inoculum x treatment interaction
# Proteobacteria
# ~3.68% of the variance in CLR abundance explained by treatment

#### Taxonomic barplots ####
# Use non-transformed data for relative abundance taxonomic barplots
drt.bact.late_phylum <- phyloseq::tax_glom(drt.bact.late, "Phylum") # 13 taxa
drt.bact.late_class <- phyloseq::tax_glom(drt.bact.late, "Class") # 21 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
drt.bact.late_phylum_relab <- transform_sample_counts(drt.bact.late_phylum, function(x) x/sum(x))
drt.bact.late_class_relab <- transform_sample_counts(drt.bact.late_class, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
drt.bact.late_phylum_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.bact.late_phylum_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
drt.bact.late_class_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.bact.late_class_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
phylum_top10 <- psmelt(drt.bact.late_phylum_relab_top10)
class_top10 <- psmelt(drt.bact.late_class_relab_top10)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
phylum_top10$Phylum <- forcats::fct_relevel(as.factor(phylum_top10$Phylum), "Other", after = Inf)
class_top10$Class <- forcats::fct_relevel(as.factor(class_top10$Class), "Other", after = Inf)

# Stats on relative abundance
mean((phylum_top10[phylum_top10$Phylum == "Actinobacteria",]$Abundance)*100)    #13.13%
mean((phylum_top10[phylum_top10$Phylum == "Armatimonadetes",]$Abundance)*100)   #0.008%
mean((phylum_top10[phylum_top10$Phylum == "Bacteroidetes",]$Abundance)*100)     #4.24%
mean((phylum_top10[phylum_top10$Phylum == "Chloroflexi",]$Abundance)*100)       #0.007%
mean((phylum_top10[phylum_top10$Phylum == "Firmicutes",]$Abundance)*100)        #0.60%
mean((phylum_top10[phylum_top10$Phylum == "Gemmatimonadetes",]$Abundance)*100)  #0.11%
mean((phylum_top10[phylum_top10$Phylum == "Nitrospirae",]$Abundance)*100)       #0.007%
mean((phylum_top10[phylum_top10$Phylum == "Planctomycetes",]$Abundance)*100)    #0.07%
mean((phylum_top10[phylum_top10$Phylum == "Proteobacteria",]$Abundance)*100)    #81.69%
mean((phylum_top10[phylum_top10$Phylum == "Verrucomicrobia",]$Abundance)*100)   #0.12
mean((phylum_top10[phylum_top10$Phylum == "Other",]$Abundance)*100)             #0.01%
# Drought specific microbes D vs W plants
tapply(phylum_top10[phylum_top10$Phylum == "Actinobacteria",]$Abundance, phylum_top10[phylum_top10$Phylum == "Actinobacteria",]$Drought.or.Watered, summary)
tapply(phylum_top10[phylum_top10$Phylum == "Chloroflexi",]$Abundance, phylum_top10[phylum_top10$Phylum == "Chloroflexi",]$Drought.or.Watered, summary)
tapply(phylum_top10[phylum_top10$Phylum == "Proteobacteria",]$Abundance, phylum_top10[phylum_top10$Phylum == "Proteobacteria",]$Drought.or.Watered, summary)
tapply(phylum_top10[phylum_top10$Phylum == "Bacteroidetes",]$Abundance, phylum_top10[phylum_top10$Phylum == "Bacteroidetes",]$Drought.or.Watered, summary)

## Plots
# phylum
# By genotype
a <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Genotype, scale = "free", ncol = 2) +
  scale_fill_manual(values=phyla_palette) +
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
  scale_fill_manual(values=phyla_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

# By Soil Inoculum
c <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilInoculum, scale = "free", ncol = 2) +
  scale_fill_manual(values=phyla_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')


# By treatment by soil location
d <- ggplot(phylum_top10, aes(x=Sample, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~SoilInoculum + Drought.or.Watered, scale = "free", ncol = 2) +
  scale_fill_manual(values=phyla_palette) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

a
ggsave("figures/barplot_genotype.svg", a, height = 8, width = 12)
b
ggsave("figures/barplot_treatment.svg", b, height = 8, width = 12)
c
ggsave("figures/barplot_inoculum.svg", c, height = 8, width = 12)
d
ggsave("figures/barplot_inoculum_treatment.svg", d, height = 8, width = 12)
rm(a,b,c,d)

#### Adonis and Ordinations ####
####### Use perMANOVA to partition variance in community composition at the ASV level
set.seed(777)

# First assess the marginal effects of the terms w/ interaction 
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype + SoilInoculum*Drought.or.Watered,
     strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
# interaction is non-significant, so remove from model and re-run to assess main effects.
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype + SoilInoculum + Drought.or.Watered,
     strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
# No difference in P values assessing by margin

# Present table with terms assessed sequential (okay, given the results above).
perm <- with(as(sample_data(drt.bact.late.clr),'data.frame'),
             adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype + SoilInoculum*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "term"))
#                                 Df SumOfSqs      R2      F Pr(>F)    
"
Genotype                          1     4109 0.00605 1.2303  0.221    
SoilInoculum                      5    34055 0.05012 2.0395  0.001 ***
Drought.or.Watered                1     7483 0.01101 2.2409  0.001 ***
SoilInoculum:Drought.or.Watered   5    16011 0.02356 0.9589  0.500    
Residual                        185   617812 0.90926                  
Total                           197   679469 1.00000"
#
# two main effects are significant.
#5% of variance in bacterial community composition explained by soil inocula.
#1% of variance in bacterial community composition explained by treatment.

# Constrained ordinations | CAP euclidean
# Start by constraining on largest factors of the anova
# by Soil Inoculum
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum+Condition(Plate + logObs))
summary(drt.bact.dbRDA)
#Color and shape of plotted points by Soil Inoculum
a <- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'SoilInoculum') + geom_point(size=3) + scale_color_manual(values = soil_merged_pallete) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilInoculum)) + scale_fill_manual(values = soil_merged_pallete)
anova.cca(drt.bact.dbRDA) 
(156.89  / (156.89  + 2418.96)) * 100 #6.09% of community composition variation explained by model

# by drought
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
summary(drt.bact.dbRDA)
#Color and shape of plotted points by  drought treatment
b<- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'Drought.or.Watered') + geom_point(size=3) + scale_color_manual(values = treatment_pallete)
anova.cca(drt.bact.dbRDA) 
(40.39 / (40.39 + 2535.46))*100 #1.57% of community composition variation explained by model

# by  genotype
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~Genotype+Condition(Plate + logObs))
summary(drt.bact.dbRDA)
#Color and shape of plotted points by genotype
c <- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'Genotype') + geom_point(size=3) + scale_color_manual(values = genotype_pallete)
anova.cca(drt.bact.dbRDA) 
(14.87 / (14.87 + 2560.98))*100 #0.58% of community composition variation explained by model

# by Soil Inoculum and treatment
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum + Drought.or.Watered + Condition(Plate + logObs))
summary(drt.bact.dbRDA)
#Color and shape of plotted points by soil location and treatment
d <- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'SoilInoculum', shape = "Drought.or.Watered") + geom_point(size=3) + scale_color_manual(values = soil_merged_pallete, name = "Soil Inoculum") + guides(shape=guide_legend(title="Treatment"))
anova.cca(drt.bact.dbRDA) 
(194.95 / (194.95 + 2380.90))*100 #7.57% of community composition variation explained by model


# plots controlled by plate and logObs and constrained by each factor
# main effects
a 
b # Used in figure 3
c
plot1 <- ggarrange(c, ncol = 1, align = 'h', legend = 'right')
ggsave("figures/main_effect_RDA_ordination.svg", plot1, height = 6, width = 6)
ggsave("figures/main_effect_RDA_ordination.png", plot1, height = 6, width = 6)
# 2 ways
d
plot2 <- ggarrange(d, ncol = 1, nrow = 1, align = 'h', labels = "AUTO")
ggsave("figures/two_way_RDA_ordination.svg", plot2, height = 8, width = 24)
ggsave("figures/two_way_RDA_ordination.png", plot2, height = 8, width = 24)


#### betadisper ####
# test for difference in variance among groups
# "homogeneity of dispersion among sample groups was assessed using the 
# betadisper function" https://microbiomejournal.biomedcentral.com/track/pdf/10.1186/s40168-017-0370-7.pdf
Betadisp.fun <- function(phyloseq_obj, Grouping, plot = TRUE){
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
Betadisp.fun(drt.bact.late.clr, 'Genotype', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Genotype', plot = FALSE)
# drought
Betadisp.fun(drt.bact.late.clr, 'Drought.or.Watered', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Drought.or.Watered', plot = FALSE)
# soil Inoculum
Betadisp.fun(drt.bact.late.clr, 'SoilInoculum', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'SoilInoculum', plot = FALSE)
# block
Betadisp.fun(drt.bact.late.clr, 'Block', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Block', plot = FALSE)
# Plate
Betadisp.fun(drt.bact.late.clr, 'Plate', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Plate', plot = FALSE)


#### Aldex (differential abundance analysis) ####
#extract OTU table matrix from phyloseq object
asvs <- as(otu_table(drt.bact.late, taxa_are_rows = TRUE), 'matrix')
asvs <- t(asvs)
# Make vectors of explanatory variables to test <2 groups
Drought.or.Watered <- sample_data(drt.bact.late)$Drought.or.Watered %>% as.character()
Genotype <- sample_data(drt.bact.late)$Genotype %>% as.character()
# Transform counts with CLR
asvs.clr.drought <- aldex.clr(asvs, conds= Drought.or.Watered, mc.samples=1000)
asvs.clr.genotype <- aldex.clr(asvs, conds= Genotype, mc.samples=1000)

# T-test with factors with 2 levels
#Outputs a dataframe with the following information:
#wi.ep   a vector containing the expected P value of the Wilcoxon test for each feature
#wi.eBH  a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature

# Drought
ttest.drought <- aldex.ttest(asvs.clr.drought)
sum(ttest.drought$wi.eBH < 0.05) # 2 significant ASV drought treatment
which(ttest.drought$wi.eBH < 0.05)
ttest.drought[46,] # bASV_46 padj = 0.004
ttest.drought[55,] # bASV_55 padj = 0.042
# Calculate effect sizes
effect.drought <- aldex.effect(asvs.clr.drought)
effect.drought[46,]
effect.drought[55,]
#       rab.all rab.win.D  rab.win.W   diff.btw diff.win    effect overlap
"bASV_46 2.670767  1.056214  6.655424 3.727558 7.631719 0.4623896 0.3061388"
"bASV_55 1.890642 0.6704841  6.476532 3.236593 8.421709 0.3960539  0.3388"
ttest.effect.drought <- data.frame(ttest.drought, effect.drought)
aldex.plot(ttest.effect.drought, type='MA') #two asvs diff. abundant in watered plants: bASV_46 and bASV_55
aldex.plot(ttest.effect.drought, type='MW') #two asvs diff. abundant
# What are these ASVS
tax_table(drt.bact.late)[46]
tax_table(drt.bact.late)[55]
# Both are Sphingomonas sp. and there are only 3 nt changes between them.
# what do the boxplots of these look like?
plot_sphingo <- data.frame(sample_data(drt.bact.late.clr), t(otu_table(drt.bact.late.clr)[46]), t(otu_table(drt.bact.late.clr)[55]))
ggplot(plot_sphingo, aes(x = Drought.or.Watered, y = bASV_46, fill = Drought.or.Watered)) + geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete)
ggplot(plot_sphingo, aes(x = Drought.or.Watered, y = bASV_55, fill = Drought.or.Watered)) + geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete)
# Are they correlated?
cor.test(plot_sphingo$bASV_46, plot_sphingo$bASV_55, method = "pearson") # yes

# Genotype
ttest.genotype <- aldex.ttest(asvs.clr.genotype)
sum(ttest.genotype$wi.eBH < 0.05) # 0 significant genotype
# Calculate effect sizes
effect.genotype <- aldex.effect(asvs.clr.genotype)
# plot results
ttest.effect.genotype <- data.frame(ttest.genotype, effect.genotype)
aldex.plot(ttest.effect.genotype, type='MA', test = 'wilcox')
aldex.plot(ttest.effect.genotype, type='MW', test = 'wilcox')

# more than two levels - Soil inoculum
SoilInoculum.df <- data.frame(SoilInoculum = sample_data(drt.bact.late)$SoilInoculum %>% as.character())
# Setup model matrix (required for factors with >2 levels)
mm <- model.matrix(~ SoilInoculum, SoilInoculum.df)
SoilInoculum.clr <- aldex.clr(asvs, mm, mc.samples=128)
SoilInoculum.glm.test <- aldex.glm(SoilInoculum.clr)
# How many ASV differ between Hay vs others
sum(SoilInoculum.glm.test$`model.SoilInoculumKNZ_Native Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumSVR_Agriculture Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumSVR_Native Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumTLI_Agriculture Pr(>|t|).BH` < 0.05)
sum(SoilInoculum.glm.test$`model.SoilInoculumTLI_Native Pr(>|t|).BH` < 0.05)
which(SoilInoculum.glm.test$`model.SoilInoculumKNZ_Native Pr(>|t|).BH` < 0.05)
SoilInoculum.glm.test[10,] # bASV_10 padj = 0.006065548
# no effect sizes for greater than 2 condition levels

#### Figure 1 Panel D: Graph of shoot height by treatment ####
####### Load cleaned, non-transformed dataset (late timepoint only, no uninoculated controls)
drt.bact.late <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS')
####### Load cleaned, CLR-transformed dataset (late timepoint only, no uninoculated controls)
drt.bact.late.clr <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS')
# will use this for experimental design figure
sample.df <- as(sample_data(drt.bact.late.clr), "data.frame")
sample.df$Drought.or.Watered <- as.factor(sample.df$Drought.or.Watered)

# mutate data to correct dim for plotting
sample.df.long <- dplyr::select(sample.df,Drought.or.Watered,Genotype,starts_with("Height")) %>% # height data only
  mutate(SampleID=factor(rownames(.)))
sample.df.long <- gather(sample.df.long,key='Height',value='Height_cm',starts_with('Height')) %>%
  mutate(Day=as.factor(str_remove(Height,'Height_')))

sample.df.long$Day<- recode(sample.df.long$Day, Jan29  = "13",
                  Feb5   = "20", 
                  Feb12  = "27", 
                  Feb19  = "34")

# Stats for legend
mod.df <- lmer(Height_cm ~ Drought.or.Watered + Day + Drought.or.Watered*Day + (1|Genotype), data = sample.df.long)
anova(mod.df)
pairs(emmeans(mod.df, ~ Drought.or.Watered|Day))
# significant effect of treatment by day interaction.

# convert day to numeric for plotting
sample.df.long$Day <- as.numeric(as.character((sample.df.long$Day)))
# Make a dataframe to hold height measure averages by treatment/day
treatment_mean <- sample.df.long %>% group_by(Drought.or.Watered, Day) %>% 
                           summarise(Height_cm  = mean(Height_cm, na.rm=TRUE))
ht_plot <- ggplot(sample.df.long, aes(x=Day,y=Height_cm,color=Drought.or.Watered)) +
  geom_vline(xintercept=13,linetype='dotted') + geom_vline(xintercept=20,linetype='dotted') +
  geom_vline(xintercept=27,linetype='dotted') + geom_vline(xintercept=34,linetype='dotted') +
  geom_point(size=2, alpha = 0.3) +
  geom_line(aes(group=SampleID), alpha = 0.3) +
  geom_line(data = treatment_mean, size = 3) +
  scale_color_manual(values= treatment_pallete) +
  labs(x = "Day", y = "Height (cm)") +
  theme(axis.title=element_text(size=18,face='bold'), axis.text=element_text(size=18), strip.text=element_text(size=18,face='bold'))+
  theme(legend.text=element_text(size=18),legend.title=element_blank(),legend.background=element_rect(color='grey77'),legend.position='right') +
  scale_x_continuous(breaks=seq(13, 34, 7))
ggsave("figures/Shoot_height_by_treatment_timeseries.svg", ht_plot, height = 6,width = 8)


#### Figure 4 panels A,B,C: Taxonomic barplot, RDA and CM by treatment ####
# reorder the class dataframe in order to samples in order of the largest class
class_top10_ordered <- class_top10[order(class_top10$Sample),]
class_top10_ordered$Phylum <- forcats::fct_relevel(as.factor(class_top10_ordered$Phylum), "Other", after = Inf)
# get a vector that will tell which samples have the highest abundance of Proteobacteria
class_top10_ordered$plot_order <- rep(rank(-class_top10_ordered$Abundance[class_top10_ordered$Phylum == "Actinobacteria"], ties.method = 'first'), each = 11)
class_top10_ordered$plot_order <- factor(class_top10_ordered$plot_order)
# convert proteobacteria to classes
# Alphaproteobacteria  Betaproteobacteria  Deltaproteobacteria  Gammaproteobacteria
class_top10_ordered$Phylum <- as.character(class_top10_ordered$Phylum)
class_top10_ordered$Class <- as.character(class_top10_ordered$Class)
class_top10_ordered[class_top10_ordered$Class == "Alphaproteobacteria",]$Phylum <- class_top10_ordered[class_top10_ordered$Class == "Alphaproteobacteria",]$Class
class_top10_ordered[class_top10_ordered$Class == "Betaproteobacteria",]$Phylum <- class_top10_ordered[class_top10_ordered$Class == "Betaproteobacteria",]$Class
class_top10_ordered[class_top10_ordered$Class == "Deltaproteobacteria",]$Phylum <- class_top10_ordered[class_top10_ordered$Class == "Deltaproteobacteria",]$Class
class_top10_ordered[class_top10_ordered$Class == "Gammaproteobacteria",]$Phylum <- class_top10_ordered[class_top10_ordered$Class == "Gammaproteobacteria",]$Class

class_top10_ordered$Phylum <- factor(class_top10_ordered$Phylum, levels = c("Actinobacteria",
                                                                              "Armatimonadetes",
                                                                              "Bacteroidetes",
                                                                              "Chloroflexi",
                                                                              "Firmicutes",
                                                                              "Gemmatimonadetes",
                                                                              "Nitrospirae",
                                                                              "Planctomycetes",
                                                                              "Alphaproteobacteria",
                                                                              "Betaproteobacteria",
                                                                              "Deltaproteobacteria",
                                                                              "Gammaproteobacteria",
                                                                              "Proteobacteria",
                                                                              "Verrucomicrobia",
                                                                              "Other"))

phyla_palette_MOD <- c("Actinobacteria" = "#D55E00",
                             "Bacteroidetes" = "#999933",
                             "Firmicutes" = "#6948b8", 
                             "Gemmatimonadetes" = "#CC6677",
                             #"Proteobacteria" = "#117733",
                             "Alphaproteobacteria" = "#90ba94",
                             "Betaproteobacteria" = "#15632c",
                             "Deltaproteobacteria" = "#458d53",
                             "Gammaproteobacteria" = "#153c1d",
                             "Verrucomicrobia" = "#882255",
                             "Other" = "#888888")

a <- ggplot(class_top10_ordered, aes(x=plot_order, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered, scale = "free", nrow = 2) +
  scale_fill_manual(values=phyla_palette_MOD, labels = c("Actinomycetota", 
                                                         "Bacteroidota",
                                                         "Bacillota",
                                                         "Gemmatimonadota",
                                                         "Alphaproteobacteria",
                                                         "Betaproteobacteria",
                                                         "Deltaproteobacteria",
                                                         "Gammaproteobacteria",
                                                         "Verrucomicrobiota",
                                                         "Other")) +
  guides(fill=guide_legend(nrow=5, byrow=TRUE)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top', 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=8))

# Pcoa treatment Panel (B) RF CM for drought treatment (C)
# by drought treatment
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
anova.cca(drt.bact.dbRDA) 
#Color and shape of plotted points drought treatment
b <- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'Drought.or.Watered') + geom_point(size=2) + scale_color_manual(values = treatment_pallete, name = "Treatment")
b <- b + annotate("text", x = 7, y = 5, label = "Treatment p <0.01", fontface = 'italic')
# Random forest prediction for treatment
c <- readRDS("./Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
c <- c$CMatrixPLOT
part <- ggarrange(b,c, ncol = 1, heights = c(1,0.85), align = 'v', labels = c("B","C"))
fig3_abc <- ggarrange(a,part, ncol = 2, widths = c(1,0.75), labels = c("A"))
ggsave("figures/figure4_panels_abc.svg", fig3_abc, height = 6, width = 10)
ggsave("figures/figure4_panels_abc.png", fig3_abc, height = 6, width = 10)


#### Figure S2: RDA by soil inoculum ####
# by soil inoculum, was cut from the figure above will now be a supplement
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilInoculum+Condition(Plate + logObs))
#Color and shape of plotted points by soil inocula source location and drought treatment
soil_inoc_RDA <- plot_ordination(drt.bact.late.clr, drt.bact.dbRDA,color = 'SoilInoculum') + geom_point(size=2) + scale_color_manual(values = soil_merged_pallete, name = "Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilInoculum)) + scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]), expression(HAY[P]),expression(TLI[Ag]), expression(TLI[P]),expression(KNZ[P]))) +
  guides(color=guide_legend(nrow=6, byrow=TRUE), fill=guide_legend(nrow=6, byrow=TRUE)) + theme(legend.text.align = 0, legend.position = 'right')

ggsave("figures/Soil_inoculum_RDA.svg", soil_inoc_RDA, height = 6, width = 6)
ggsave("figures/Soil_inoculum_RDA.png", soil_inoc_RDA, height = 6, width = 6)