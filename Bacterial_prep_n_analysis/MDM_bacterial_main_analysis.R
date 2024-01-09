# Analysis of maize B73 and Mo17 inoculated plants under drought stress
# Samples were collected from a greenhouse study in Feb 2020
# Microbiomes are from the root compartment (endosphere)
# Growth measurements were also recorded throughout the experiment (~50 days in length)
# Code by: Maggie Wagner, Matthew Kolp, Joel Swift

#### Design of experiment ####
# 2 genotypes: B73 / Mo17
# 2 treatments: Well watered / Drought
# 4 soil inoculates
## Smoky Valley Ranch - SVR
## Hays Prairie - HAY
## The Land Institute - TLI
## Konza Prairie - KNZ
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
library(fantaxtic); packageVersion('fantaxtic')

# Theme set and Color Palettes
theme_set(theme_pubr())
genotype_pallete <- c("B73" = "#91ff26", "Mo17" = "#9426ff")# B73/Mo17 - Genotype
treatment_pallete <- c("W" = "#0000FF", "D" = "#DAA520") # Drought/WW     - Treatment
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#ffe785", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location


### Load datasets for use throughout ###
drt.bact.late <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS')
drt.bact.late.clr <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS')

### Late Timepoint ###
#### Alpha diversity ####
# Split out a metadata df for use outside of phyloseq
sampledata <- sample_data(drt.bact.late.clr)
# LMs w/ genotype included in the model
drt.bact.late.z.richness <- lmerTest::lmer(Chao1  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.bact.late.Shannon    <- lmerTest::lmer(log(Shannon)  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.bact.late.InvSimpson <- lmerTest::lmer(log(InvSimpson)  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
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
rand(drt.bact.late.z.richness)
rand(drt.bact.late.Shannon) # Plate significant
rand(drt.bact.late.InvSimpson) # Plate significant 
# posthoc
emmeans(drt.bact.late.z.richness,~SoilLocation)
pairs(emmeans(drt.bact.late.z.richness,~SoilLocation))

####  Which Families are impacted by the main effects ####
drt.bact.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 70 taxa
# LMs
bact.clr.long.aovs_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lmer(Abundance ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data=data)))) %>%
  summarize(broom::tidy(mod))
# get vector by term for bacterial genera, correct pvalues with BH, extract the significant ones
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype", ][c(temp), ]
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "SoilLocation", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "SoilLocation", ][c(temp), ]
# 2 way interactions
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:SoilLocation", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:SoilLocation", ][c(temp), ]
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "SoilLocation:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "SoilLocation:Drought.or.Watered", ][c(temp), ]
# 3 way
temp <- which(p.adjust(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:SoilLocation:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$term == "Genotype:SoilLocation:Drought.or.Watered", ][c(temp), ]

# Conclusion 
# Families significant for drought: Sphingomonadaceae
# Families significant for Genotype: None
# Families significant for inoculum: Nocardiaceae
# Families significant for Interaction (GxSI): None
# Families significant for Interaction (GxDT): None
# Families significant for Interaction (SIxDT): None
# Families significant for Interaction (GxSIxDT): None

### Percent variance explained 
(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$Family == "Sphingomonadaceae", ]$sumsq[3] /
    sum(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$Family == "Sphingomonadaceae", ]$sumsq)) * 100 # 49.0% variance explained
(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$Family == "Nocardiaceae", ]$sumsq[2] /
    sum(bact.clr.long.aovs_fam[bact.clr.long.aovs_fam$Family == "Nocardiaceae", ]$sumsq)) * 100 # 52.6% variance explained

### Relative abundance stats
drt.bact.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 70 taxa
drt.bact.late_family <- phyloseq::tax_glom(drt.bact.late, "Family")
drt.bact.late_family_relab <- transform_sample_counts(drt.bact.late_family, function(x) x/sum(x))
family_relab <- psmelt(drt.bact.late_family_relab)
tapply(family_relab[family_relab$Family == "Sphingomonadaceae",]$Abundance, family_relab[family_relab$Family == "Sphingomonadaceae",]$Drought.or.Watered, summary)
tapply(family_relab[family_relab$Family == "Nocardiaceae",]$Abundance, family_relab[family_relab$Family == "Nocardiaceae",]$SoilLocation, summary)

# Stats for three ML ASVs
# ASV 46, 47, 55
ML_relab <- transform_sample_counts(drt.bact.late, function(x) x/sum(x))
ML_relab_df <- psmelt(ML_relab)
# Sphingomonadaceae 
tapply(ML_relab_df[ML_relab_df$OTU == "bASV_46",]$Abundance, ML_relab_df[ML_relab_df$OTU == "bASV_46",]$Drought.or.Watered, summary)
tapply(ML_relab_df[ML_relab_df$OTU == "bASV_55",]$Abundance, ML_relab_df[ML_relab_df$OTU == "bASV_55",]$Drought.or.Watered, summary)
# Sphingobacteriaceae 
tapply(ML_relab_df[ML_relab_df$OTU == "bASV_47",]$Abundance, ML_relab_df[ML_relab_df$OTU == "bASV_47",]$Drought.or.Watered, summary)


#### Taxonomic barplots ####
# Use non-transformed data for relative abundance taxonomic barplots
drt.bact.late_class <- phyloseq::tax_glom(drt.bact.late, "Class") # 21 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
drt.bact.late_class_relab <- transform_sample_counts(drt.bact.late_class, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
drt.bact.late_class_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.bact.late_class_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Other")
# Melt data frame with phyloseq function for plotting.
class_top10 <- psmelt(drt.bact.late_class_relab_top10)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
class_top10$Class <- forcats::fct_relevel(as.factor(class_top10$Class), "Other", after = Inf)
class_top10_ordered <- class_top10[order(class_top10$Sample),]
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

# Set the order for the barplot to use
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

# Fixing labeling and order issues for facets + defining palette to use
class_top10_ordered <- class_top10_ordered %>%
  mutate(Drought.or.Watered = fct_recode(as.factor(Drought.or.Watered), 
                                         Drought = "D",
                                         "Well-Watered" = "W"))
class_top10_ordered$SoilLocation <- factor(class_top10_ordered$SoilLocation, levels = c("SVR", "HAY", "TLI", "KNZ"))

palette_TB <- c("Actinobacteria" = "#D55E00",
                "Bacteroidetes" = "#999933",
                "Firmicutes" = "#6948b8", 
                "Gemmatimonadetes" = "#CC6677",
                "Proteobacteria" = "#117733",
                "Alphaproteobacteria" = "#90ba94",
                "Betaproteobacteria" = "#15632c",
                "Deltaproteobacteria" = "#458d53",
                "Gammaproteobacteria" = "#153c1d",
                "Verrucomicrobia" = "#882255",
                "Other" = "#888888")

# plot barplot
Figure4A <- ggplot(class_top10_ordered, aes(x=plot_order, y=Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered + SoilLocation, scale = "free_x", nrow = 2, labeller = label_parsed) +
  scale_fill_manual(name = "Phylum/Class", values=palette_TB, labels = c("Actinomycetota",
                                                         "Bacteroidota",
                                                         "Bacillota",
                                                         "Gemmatimonadota",
                                                         "Alphaproteobacteria",
                                                         "Betaproteobacteria",
                                                         "Deltaproteobacteria",
                                                         "Gammaproteobacteria",
                                                         "Verrucomicrobiota",
                                                         "Low Abundance Taxa")) +
  guides(fill=guide_legend(nrow=5, byrow=TRUE)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top', 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=8))

#### Permanova and Ordination ####
####### Use perMANOVA to partition variance in community composition at the ASV level
set.seed(777)
# First assess the marginal effects of the terms w/ 3 way interaction 
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype*SoilLocation*Drought.or.Watered,
     strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
# 2 way interactions, none significant 
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype*SoilLocation + Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype + SoilLocation*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype*Drought.or.Watered + SoilLocation,
             strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))

# interactions are non-significant, so remove from model and re-run to assess main effects.
with(as(sample_data(drt.bact.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype + SoilLocation + Drought.or.Watered,
     strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "margin"))
# No difference in P values assessing by margin

# Present table with terms assessed sequential to include interaction terms.
perm <- with(as(sample_data(drt.bact.late.clr),'data.frame'),
             adonis2(t(as(otu_table(drt.bact.late.clr),'matrix')) ~ Genotype*SoilLocation*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.bact.late.clr),'data.frame'), method='euclidean', by = "term"))
"                                          Df SumOfSqs      R2      F Pr(>F)
Genotype                                   1     2687 0.00664 0.8849  0.466    
SoilLocation                               3    19091 0.04714 2.0958  0.001 ***
Drought.or.Watered                         1     5484 0.01354 1.8062  0.014 *  
Genotype:SoilLocation                      3     9849 0.02432 1.0812  0.471    
Genotype:Drought.or.Watered                1     2289 0.00565 0.7537  0.715    
SoilLocation:Drought.or.Watered            3    10026 0.02476 1.1007  0.434    
Genotype:SoilLocation:Drought.or.Watered   3     6352 0.01568 0.6973  0.949    
Residual                                 115   349180 0.86226                  
Total                                    130   404957 1.0000"

# two main effects are significant.
#5% of variance in bacterial community composition explained by soil inocula.
#1% of variance in bacterial community composition explained by treatment.
# Pcoa treatment Panel (Figure4B)
# by drought treatment
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
anova.cca(drt.bact.dbRDA) 
#Color and shape of plotted points drought treatment
b <- plot_ordination(drt.bact.late.clr,drt.bact.dbRDA,color = 'Drought.or.Watered') + geom_point(size=2) + scale_color_manual(values = treatment_pallete, name = "Treatment", labels = c("Drought (D)", "Well-Watered (W)"))
b <- b + annotate("text", x = 4, y = 8, label = "Treatment p <0.01", fontface = 'italic')
Figure4B <- b + guides(color = guide_legend(nrow = 2))

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
Betadisp.fun(drt.bact.late.clr, 'SoilLocation', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'SoilLocation', plot = FALSE)
# block
Betadisp.fun(drt.bact.late.clr, 'Block', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Block', plot = FALSE)
# Plate
Betadisp.fun(drt.bact.late.clr, 'Plate', plot = TRUE)
Betadisp.fun(drt.bact.late.clr, 'Plate', plot = FALSE)

#### Figure 4 panels A,B,C: Taxonomic barplot, RDA and CM by treatment ####
Figure4A
Figure4B
# Panel C Random forest prediction for treatment
ML_obj <- readRDS("./Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
Figure4C <- ML_obj$CMatrixPLOT
part <- ggarrange(Figure4B, Figure4C, ncol = 1, heights = c(1,0.85), align = 'v', labels = c("B","C"))
Figure4 <- ggarrange(Figure4A, part, ncol = 2, widths = c(1.0,0.75), labels = c("A"))
ggsave("figures/figure4_panels_ABC.svg", Figure4, height = 6, width = 10)
ggsave("figures/figure4_panels_ABC.png", Figure4, height = 6, width = 10)

#### Figure S8: RDA by soil inoculum ####
# by soil inoculum, was cut from the figure above will now be a supplement
drt.bact.dbRDA <- ordinate(drt.bact.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilLocation+Condition(Plate + logObs))
#Color and shape of plotted points by soil inocula source location and drought treatment
drt.bact.late.clr@sam_data$SoilLocation <- factor(drt.bact.late.clr@sam_data$SoilLocation , levels = c("SVR", "HAY", "TLI", "KNZ"))
soil_inoc_RDA <- plot_ordination(drt.bact.late.clr, drt.bact.dbRDA,color = 'SoilLocation') + geom_point(size=2) + scale_color_manual(values = location_pallete, name = "Soil Inoculum") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilLocation)) + scale_fill_manual(values = location_pallete, name = "Soil Inoculum") +
  guides(color=guide_legend(nrow=6, byrow=TRUE), fill=guide_legend(nrow=6, byrow=TRUE)) + theme(legend.text.align = 0, legend.position = 'right')
ggsave("figures/Soil_inoculum_RDA.svg", soil_inoc_RDA, height = 6, width = 6)
ggsave("figures/Soil_inoculum_RDA.png", soil_inoc_RDA, height = 6, width = 6)

