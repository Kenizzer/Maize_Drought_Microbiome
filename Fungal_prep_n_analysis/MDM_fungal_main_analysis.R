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
drt.fungi.late <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late.RDS')
drt.fungi.late.clr <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late_clr.RDS')

### Late Timepoint ###
#### Alpha diversity ####
# Split out a metadata df for use outside of phyloseq
sampledata <- sample_data(drt.fungi.late.clr)
# LMs w/ genotype included in the model
drt.fungi.late.z.richness <- lmerTest::lmer(Chao1  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.fungi.late.Shannon    <- lmerTest::lmer(log(Shannon)  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
drt.fungi.late.InvSimpson <- lmerTest::lmer(log(InvSimpson)  ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data = data.frame(sampledata))
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
emmeans(drt.fungi.late.Shannon, ~ Drought.or.Watered)
emmeans(drt.fungi.late.InvSimpson, ~ Drought.or.Watered)

pairs(emmeans(drt.fungi.late.Shannon, ~ Drought.or.Watered))
pairs(emmeans(drt.fungi.late.InvSimpson, ~ Drought.or.Watered))

#### Which Families are impacted by the main effects ####
drt.fungi.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.fungi.late.clr, "Family")) # 84 taxa
# LMs
fungi.clr.long.aovs_gen <- drt.fungi.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lmer(Abundance ~ Genotype*SoilLocation*Drought.or.Watered + (1|Block) + (1|Plate), data=data)))) %>%
  summarize(broom::tidy(mod))
# get vector by term for fungal genera, correct pvalues with BH, extract the significant ones
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilLocation", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilLocation", ][c(temp), ]
# 2 way interactions
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:SoilLocation", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:SoilLocation", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilLocation:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "SoilLocation:Drought.or.Watered", ][c(temp), ]
# 3 way
temp <- which(p.adjust(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:SoilLocation:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$term == "Genotype:SoilLocation:Drought.or.Watered", ][c(temp), ]
# Conclusion 
# Families significant for drought: None
# Families significant for Genotype: None
# Families significant for inoculum: Aspergillaceae, Mucoraceae
# Families significant for Interaction (GxSI): None
# Families significant for Interaction (GxDT): None
# Families significant for Interaction (SIxDT): None
# Families significant for Interaction (GxSIxDT): None

### Percent variance explained 
(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$Family == "Aspergillaceae", ]$sumsq[2] /
    sum(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$Family == "Aspergillaceae", ]$sumsq)) * 100 # 62.6% variance explained
(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$Family == "Mucoraceae", ]$sumsq[2] /
    sum(fungi.clr.long.aovs_gen[fungi.clr.long.aovs_gen$Family == "Mucoraceae", ]$sumsq)) * 100 # 47.8% variance explained

# Stats for two ML ASVs
# ASV 4 and 5
ML_relab <- transform_sample_counts(drt.fungi.late, function(x) x/sum(x))
ML_relab_df <- psmelt(ML_relab)
# Aspergillaceae 
tapply(ML_relab_df[ML_relab_df$OTU == "fASV_4",]$Abundance, ML_relab_df[ML_relab_df$OTU == "fASV_4",]$Drought.or.Watered, summary)
# Sarocladiaceae 
tapply(ML_relab_df[ML_relab_df$OTU == "fASV_5",]$Abundance, ML_relab_df[ML_relab_df$OTU == "fASV_5",]$Drought.or.Watered, summary)



#### Taxonomic barplots ####
# Use non-transformed data for relative abundance taxonomic barplots
drt.fungi.late_class <- phyloseq::tax_glom(drt.fungi.late, "Class") # 12 taxa
# Make relative abundance AKA transform sample counts by diving by the ASV total.
drt.fungi.late_class_relab <- transform_sample_counts(drt.fungi.late_class, function(x) x/sum(x))
# Take only the n number of taxa per taxonomic rank based on relative abundance.
# Taxa >n will be added to a other label.
drt.fungi.late_class_relab_top10 <- fantaxtic::get_top_taxa(physeq_obj = drt.fungi.late_class_relab, n = 10, relative = TRUE, discard_other = FALSE, other_label = "Low Abundance Taxa")
# Melt data frame with phyloseq function for plotting.
class_top10 <- psmelt(drt.fungi.late_class_relab_top10)
# Reorder levels to put other to the end, otherwise taxa are in alphabetical order.
class_top10$Class <- forcats::fct_relevel(as.factor(class_top10$Class), "Low Abundance Taxa", after = Inf)
# reorder the class dataframe in order to samples in order of the largest class
class_top10_ordered <- class_top10[order(class_top10$Sample),]
# get a vector that will tell which samples have the highest abundance of Eurotiomycetes
class_top10_ordered$plot_order <- rep(rank(-class_top10_ordered$Abundance[class_top10_ordered$Class == "Eurotiomycetes"], ties.method = 'first'), each = 11)
class_top10_ordered$plot_order <- factor(class_top10_ordered$plot_order)
class_top10_ordered$Class <- factor(class_top10_ordered$Class, levels = c("Eurotiomycetes", "Dothideomycetes", "Leotiomycetes", "Malasseziomycetes", "Microbotryomycetes", "Mortierellomycetes", "Pezizomycetes", "Sordariomycetes", "unidentified", "Low Abundance Taxa"))

# Fixing labeling and order issues for facets + defining palette to use
class_top10_ordered <- class_top10_ordered %>%
  mutate(Drought.or.Watered = fct_recode(as.factor(Drought.or.Watered), 
                                         Drought = "D",
                                         "Well-Watered" = "W"))
class_top10_ordered$SoilLocation <- factor(class_top10_ordered$SoilLocation, levels = c("SVR", "HAY", "TLI", "KNZ"))

taxa_palette <- c("Eurotiomycetes" = "#88CCEE", "Dothideomycetes" = "#CC6677",
                  "Leotiomycetes"= "#DDCC77", "Malasseziomycetes" = "#AA4499",
                  "Microbotryomycetes" = "#999933", "Mortierellomycetes"= "#44AA99",
                  "Pezizomycetes" = "#D55E00", "Sordariomycetes" = "#117733",
                  "unidentified" = "#6699CC", "Low Abundance Taxa" = "#882255") 

Figure5A<- ggplot(class_top10_ordered, aes(x=plot_order, y=Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  facet_wrap(~Drought.or.Watered + SoilLocation, scale = "free_x", nrow = 2, labeller = label_parsed) +
  scale_fill_manual(values=taxa_palette) + 
  guides(fill=guide_legend(nrow=5, byrow=TRUE))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top', 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=8))


#### Adonis and Ordinations ####
####### Use perMANOVA to partition variance in community composition at the ASV level
set.seed(8945761)

# First assess the marginal effects of the terms w/ 3 way interaction 
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype*SoilLocation*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
# 2 way interactions, none significant 
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype*SoilLocation + Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype + SoilLocation*Drought.or.Watered,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype*Drought.or.Watered + SoilLocation,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))

# interactions are non-significant, so remove from model and re-run to assess main effects.
with(as(sample_data(drt.fungi.late.clr),'data.frame'),
     adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype + Drought.or.Watered + SoilLocation,
             strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "margin"))
# No difference in P values assessing by margin


# Present table with terms assessed sequential (okay, given the results above).
perm <- with(as(sample_data(drt.fungi.late.clr),'data.frame'),
             adonis2(t(as(otu_table(drt.fungi.late.clr),'matrix')) ~ Genotype*SoilLocation*Drought.or.Watered,
                     strata = Block, data=as(sample_data(drt.fungi.late.clr),'data.frame'), method='euclidean', by = "term"))
"                                         Df SumOfSqs      R2      F Pr(>F)    
Genotype                                  1      586 0.00914 0.8933  0.486    
SoilLocation                              3     5896 0.09203 2.9982  0.001 ***
Drought.or.Watered                        1      541 0.00844 0.8250  0.633    
Genotype:SoilLocation                     3     2181 0.03404 1.1089  0.602    
Genotype:Drought.or.Watered               1      740 0.01156 1.1295  0.175    
SoilLocation:Drought.or.Watered           3     2335 0.03645 1.1875  0.167    
Genotype:SoilLocation:Drought.or.Watered  3     1969 0.03073 1.0010  0.609    
Residual                                 76    49822 0.77762                  
Total                                    91    64070 1.00000"
#
#9.2% of variance in fungal community composition explained by soil inocula.
#2nd largest effect size (r2): SoilInoculum:Drought.or.Watered, 3.6% of variance in fungal community composition.

# Pcoa treatment Panel (Figure5B)
# by treatment
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~Drought.or.Watered+Condition(Plate + logObs))
anova.cca(drt.fungi.dbRDA)
#Color and shape of plotted points by soil habitat
b <- plot_ordination(drt.fungi.late.clr,drt.fungi.dbRDA,color = 'Drought.or.Watered') + geom_point(size=3) + scale_color_manual(values = treatment_pallete, name = "Treatment", labels = c("Drought (D)", "Well-Watered (W)"))
b <- b + annotate("text", x = 4, y = 4, label = "Treatment p = 0.42", fontface = 'italic')
Figure5B <- b + guides(color = guide_legend(nrow = 2))

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
Betadisp.fun(drt.fungi.late.clr, 'SoilLocation', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'SoilLocation', plot = FALSE)
# block
Betadisp.fun(drt.fungi.late.clr, 'Block', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'Block', plot = FALSE)
# Plate
Betadisp.fun(drt.fungi.late.clr, 'Plate', plot = TRUE)
Betadisp.fun(drt.fungi.late.clr, 'Plate', plot = FALSE)

#### Figure 5 panels A,B,C #####
Figure5A
Figure5B
# Panel C Random forest prediction for treatment
ML_obj <- readRDS("./Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
Figure5C <- ML_obj$CMatrixPLOT
part <- ggarrange(Figure5B, Figure5C, ncol = 1, heights = c(1,0.85), labels = c("B","C"))
Figure5 <- ggarrange(Figure5A, part, ncol = 2, widths = c(1.0,0.75), labels = c("A"))
ggsave("figures/figure5_panels_ABC.svg", Figure5, height = 6, width = 10)
ggsave("figures/figure5_panels_ABC.png", Figure5, height = 6, width = 10)

#### Figure S10: RDA by soil inoculum ####
# by soil inoculum, was cut from the figure above will now be a supplement
drt.fungi.dbRDA <- ordinate(drt.fungi.late.clr, distance = 'euclidean', method = 'RDA', formula=~SoilLocation+Condition(Plate + logObs))
#Color and shape of plotted points by soil inocula source location and drought treatment
drt.fungi.late.clr@sam_data$SoilLocation <- factor(drt.fungi.late.clr@sam_data$SoilLocation , levels = c("SVR", "HAY", "TLI", "KNZ"))
soil_inoc_RDA <- plot_ordination(drt.fungi.late.clr, drt.fungi.dbRDA,color = 'SoilLocation') + geom_point(size=2) + scale_color_manual(values = location_pallete, name = "Soil Inoculum") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.1, aes(fill=SoilLocation)) + scale_fill_manual(values = location_pallete, name = "Soil Inoculum") +
  guides(color=guide_legend(nrow=6, byrow=TRUE), fill=guide_legend(nrow=6, byrow=TRUE)) + theme(legend.text.align = 0, legend.position = 'right')
ggsave("figures/Soil_inoculum_RDA.svg", soil_inoc_RDA, height = 6, width = 6)
ggsave("figures/Soil_inoculum_RDA.png", soil_inoc_RDA, height = 6, width = 6)

