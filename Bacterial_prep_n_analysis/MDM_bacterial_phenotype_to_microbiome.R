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

#=======================================================#
#==== Connecting plant traits to microbe abundance =====#
#=======================================================#

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
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#7fd66f", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location

### Load data set for use throughout ###
drt.bact.late.clr <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS') # CLR normalized
drt.bact.late     <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS') # non-normalized
# Calculate RSR
sample_data(drt.bact.late.clr)$RootShootRatio <- sample_data(drt.bact.late.clr)$ThirdDryRootWeight / sample_data(drt.bact.late.clr)$ThirdDryShootWeight
sample_data(drt.bact.late)$RootShootRatio <- sample_data(drt.bact.late)$ThirdDryRootWeight / sample_data(drt.bact.late)$ThirdDryShootWeight
# Remove one sample that lacks root weight measures
drt.bact.late.clr <- subset_samples(drt.bact.late.clr, !(SampleID %in% c("JJ-B-11")))
drt.bact.late     <- subset_samples(drt.bact.late, !(SampleID %in% c("JJ-B-11")))
# Extract data frame from phyloseq sample data
bact.samp.df.clr <- data.frame(sample_data(drt.bact.late.clr))

##### ShootMassRate #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = bact.samp.df.clr)
#add residuals of SMR to dataframe and phyloseq obj
bact.samp.df.clr$ShootMassRateResid <- resid(mod)
sample_data(drt.bact.late.clr)$ShootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualSMR 
drt.bact.late.clr_phy <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Phylum")) # 13 taxa
drt.bact.late.clr_cla <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Class")) # 21 taxa
drt.bact.late.clr_ord <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Order")) # 37 taxa
drt.bact.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 84 taxa
drt.bact.late.clr_gen <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Genus")) # 161 taxa

# LMs
bact.clr.long.SMR.aovs_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.SMR.aovs_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.SMR.aovs_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.SMR.aovs_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.SMR.aovs_gen <- drt.bact.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))

# Adjust P values across taxonomic level but separately per term
# Get corrected P values from here

Padj.term.func <- function(factor, modelDF){
  # Check if any are significant, save those their indexes to temp
  temp <- which(p.adjust(modelDF[modelDF$term == factor, ]$p.value, method = "BH") < 0.05)
    temp2 <- filter(modelDF, term == factor)
    # adjust p values by term and add this to a term specific tibble save to temp2
    temp2$p.adj <- p.adjust(modelDF[modelDF$term == factor, ]$p.value, method = "BH")
    return(temp2[c(temp), ]) # index temp2 with temp
    # Returns an empty tibble if none are significant after p.adjust
}

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.SMR.aovs_phy)
Padj.term.func("Abundance", bact.clr.long.SMR.aovs_phy)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.SMR.aovs_cla)
Padj.term.func("Abundance", bact.clr.long.SMR.aovs_cla)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.SMR.aovs_ord)
Padj.term.func("Abundance", bact.clr.long.SMR.aovs_ord)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.SMR.aovs_fam)
Padj.term.func("Abundance", bact.clr.long.SMR.aovs_fam)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.SMR.aovs_gen)
Padj.term.func("Abundance", bact.clr.long.SMR.aovs_gen)
  
# phylum 
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.phy <-  phyloseq::tax_glom(drt.bact.late, "Phylum") # 9 taxa
drt.bact.late.phy.relab <- psmelt(transform_sample_counts(drt.bact.late.phy, function(x) x/sum(x)))

bact.clr.long.SMR.mod_phy[bact.clr.long.SMR.mod_phy$Phylum == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_phy %>% group_by(Phylum) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Phylum == "Actinobacteria")
mean((drt.bact.late.phy.relab[drt.bact.late.phy.relab$Phylum == "Actinobacteria",]$Abundance)*100)


# class
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.cla <-  phyloseq::tax_glom(drt.bact.late, "Class") # 15 taxa
drt.bact.late.cla.relab <- psmelt(transform_sample_counts(drt.bact.late.cla, function(x) x/sum(x)))

bact.clr.long.SMR.mod_cla[bact.clr.long.SMR.mod_cla$Class == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Actinobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Actinobacteria",]$Abundance)*100) 


# Order
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.ord <-  phyloseq::tax_glom(drt.bact.late, "Order") # 29 taxa
drt.bact.late.ord.relab <- psmelt(transform_sample_counts(drt.bact.late.ord, function(x) x/sum(x)))

bact.clr.long.SMR.mod_ord[bact.clr.long.SMR.mod_ord$Order == "Actinomycetales",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Actinomycetales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Actinomycetales",]$Abundance)*100) 

bact.clr.long.SMR.mod_ord[bact.clr.long.SMR.mod_ord$Order == "Rhodocyclales",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Rhodocyclales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Rhodocyclales",]$Abundance)*100) 

bact.clr.long.SMR.mod_ord[bact.clr.long.SMR.mod_ord$Order == "Sphingomonadales",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Sphingomonadales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Sphingomonadales",]$Abundance)*100) 


# Family - w/ interaction
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.fam <-  phyloseq::tax_glom(drt.bact.late, "Family") # 84 taxa
drt.bact.late.fam.relab <- psmelt(transform_sample_counts(drt.bact.late.fam, function(x) x/sum(x)))

#interaction significant 
bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Streptococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Streptococcaceae",]$Abundance)*100) 

# STATS
car::Anova(bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]], type = "III")
mod <- bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]
emtrends(mod, ~Drought.or.Watered, var = "Abundance")

Strep_interaction_plot <- filter(drt.bact.late.clr_fam, Family == "Streptococcaceae") %>%
  ggplot(aes(x = Abundance, y = ShootMassRateResid, color = Drought.or.Watered)) +
  xlab("CLR Transformed Abundance") + ylab("Residual Shoot Mass Rate (g/day)") +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment", values = treatment_pallete, labels = c("Drought", "Well-Watered")) +
  annotate("text", x = 8, y = 0.040, label = "W = -0.0039", fontface = 'italic') +
  annotate("text", x = 8, y = 0.035, label = "D =  0.0001", fontface = 'italic') +
  annotate("text", x = 8, y = 0.030, label = "PVEres = 5.05", fontface = 'italic') +
  annotate("text", x = 8, y = 0.025, label = "Abund.×DT q = 0.02", fontface = 'italic') +
  theme(legend.position = 'right')

ggsave("figures/Streptococcaceae_interaction_plot.svg", Strep_interaction_plot, height = 6, width = 6)
ggsave("figures/Streptococcaceae_interaction_plot.png", Strep_interaction_plot, height = 6, width = 6)


# Main effect significant
bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Acetobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Acetobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Acetobacteraceae",]$Abundance)*100)

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Corynebacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Corynebacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Corynebacteriaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Dietziaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Dietziaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Dietziaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Lactobacillaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Lactobacillaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Lactobacillaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Microbacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Microbacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Microbacteriaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Micrococcaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Micrococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Micrococcaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Nocardiaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Nocardiaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Nocardiaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Rhodocyclaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Rhodocyclaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Rhodocyclaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Sinobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sinobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sinobacteraceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Sphingomonadaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sphingomonadaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sphingomonadaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Streptococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Streptococcaceae",]$Abundance)*100) 

# Genus
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_gen <- drt.bact.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.gen <-  phyloseq::tax_glom(drt.bact.late, "Genus") # 84 taxa
drt.bact.late.gen.relab <- psmelt(transform_sample_counts(drt.bact.late.gen, function(x) x/sum(x)))

# Main Effects
bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Corynebacterium",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Corynebacterium")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Corynebacterium",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Cupriavidus",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Cupriavidus")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Cupriavidus",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Dietzia",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Dietzia")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Dietzia",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Kribbella",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Kribbella")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Kribbella",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Lactobacillus",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Lactobacillus")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Lactobacillus",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Minicystis",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Minicystis")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Minicystis",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Modestobacter",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Modestobacter")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Modestobacter",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Rhodococcus",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Rhodococcus")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Rhodococcus",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Siphonobacter",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Siphonobacter")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Siphonobacter",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Sphingobium",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Sphingobium")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Sphingobium",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Sphingomonas",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Sphingomonas")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Sphingomonas",]$Abundance)*100) 

bact.clr.long.SMR.mod_gen[bact.clr.long.SMR.mod_gen$Genus == "Streptococcus",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Streptococcus")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Streptococcus",]$Abundance)*100) 









##### RootMassRate #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(sqrt(RootMassRate) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = bact.samp.df.clr)
# add residuals of RMR to data frame
bact.samp.df.clr$RootMassRateResid <- resid(mod)
sample_data(drt.bact.late.clr)$RootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualRMR
drt.bact.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Phylum")) # 13 taxa
drt.bact.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Class")) # 21 taxa
drt.bact.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Order")) # 37 taxa
drt.bact.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 84 taxa
drt.bact.late.clr_gen <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Genus")) # 161 taxa
# LMs
bact.clr.long.RMR.aovs_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RMR.aovs_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RMR.aovs_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RMR.aovs_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RMR.aovs_gen <- drt.bact.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))

#Adjust P values across term
Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RMR.aovs_phy)
Padj.term.func("Abundance", bact.clr.long.RMR.aovs_phy)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RMR.aovs_cla)
Padj.term.func("Abundance", bact.clr.long.RMR.aovs_cla)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RMR.aovs_ord)
Padj.term.func("Abundance", bact.clr.long.RMR.aovs_ord)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RMR.aovs_fam)
Padj.term.func("Abundance", bact.clr.long.RMR.aovs_fam)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RMR.aovs_gen)
Padj.term.func("Abundance", bact.clr.long.RMR.aovs_gen)

# Re-run models but dont send through anova and broom
bact.clr.long.RMR.mod_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(RootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.fam <-  phyloseq::tax_glom(drt.bact.late, "Family") # 84 taxa
drt.bact.late.fam.relab <- psmelt(transform_sample_counts(drt.bact.late.fam, function(x) x/sum(x)))


bact.clr.long.RMR.mod_fam[bact.clr.long.RMR.mod_fam$Family == "Acetobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.RMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Acetobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Acetobacteraceae",]$Abundance)*100) 

bact.clr.long.RMR.mod_fam[bact.clr.long.RMR.mod_fam$Family == "Corynebacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.RMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Corynebacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Corynebacteriaceae",]$Abundance)*100) 





##### Root/Shoot Ratio #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(log(RootShootRatio) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = bact.samp.df.clr)
# add residuals of RSR to data frame
bact.samp.df.clr$RootShootRatioResid <- resid(mod)
sample_data(drt.bact.late.clr)$RootShootRatioResid <- resid(mod)
# Correlate taxon abundance to residualRSR
drt.bact.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Phylum")) # 13 taxa
drt.bact.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Class")) # 21 taxa
drt.bact.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Order")) # 37 taxa
drt.bact.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 84 taxa
drt.bact.late.clr_gen <-  psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Genus")) # 161 taxa
# LMs
bact.clr.long.RSR.aovs_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RSR.aovs_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RSR.aovs_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RSR.aovs_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
bact.clr.long.RSR.aovs_gen <- drt.bact.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))

#Adjust P values across taxonomic level
Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RSR.aovs_phy)
Padj.term.func("Abundance", bact.clr.long.RSR.aovs_phy)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RSR.aovs_cla)
Padj.term.func("Abundance", bact.clr.long.RSR.aovs_cla)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RSR.aovs_ord)
Padj.term.func("Abundance", bact.clr.long.RSR.aovs_ord)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RSR.aovs_fam)
Padj.term.func("Abundance", bact.clr.long.RSR.aovs_fam)

Padj.term.func("Abundance:Drought.or.Watered", bact.clr.long.RSR.aovs_gen)
Padj.term.func("Abundance", bact.clr.long.RSR.aovs_gen)


# phylum 
# Re-run models but dont send through anova and broom
bact.clr.long.RSR.mod_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.phy <-  phyloseq::tax_glom(drt.bact.late, "Phylum") # 9 taxa
drt.bact.late.phy.relab <- psmelt(transform_sample_counts(drt.bact.late.phy, function(x) x/sum(x)))

bact.clr.long.RSR.mod_phy[bact.clr.long.RSR.mod_phy$Phylum == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_phy %>% group_by(Phylum) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Phylum == "Actinobacteria")
mean((drt.bact.late.phy.relab[drt.bact.late.phy.relab$Phylum == "Actinobacteria",]$Abundance)*100)

bact.clr.long.RSR.mod_phy[bact.clr.long.RSR.mod_phy$Phylum == "Nitrospirae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_phy %>% group_by(Phylum) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Phylum == "Nitrospirae")
mean((drt.bact.late.phy.relab[drt.bact.late.phy.relab$Phylum == "Nitrospirae",]$Abundance)*100)

# class
# Re-run models but dont send through anova and broom
bact.clr.long.RSR.mod_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.cla <-  phyloseq::tax_glom(drt.bact.late, "Class") # 15 taxa
drt.bact.late.cla.relab <- psmelt(transform_sample_counts(drt.bact.late.cla, function(x) x/sum(x)))

bact.clr.long.RSR.mod_cla[bact.clr.long.RSR.mod_cla$Class == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Actinobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Actinobacteria",]$Abundance)*100) 

bact.clr.long.RSR.mod_cla[bact.clr.long.RSR.mod_cla$Class == "Deltaproteobacteria",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Deltaproteobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Deltaproteobacteria",]$Abundance)*100) 

bact.clr.long.RSR.mod_cla[bact.clr.long.RSR.mod_cla$Class == "Nitrospira",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Nitrospira")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Nitrospira",]$Abundance)*100) 


# Order
# Re-run models but dont send through anova and broom
bact.clr.long.RSR.mod_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.ord <-  phyloseq::tax_glom(drt.bact.late, "Order") # 29 taxa
drt.bact.late.ord.relab <- psmelt(transform_sample_counts(drt.bact.late.ord, function(x) x/sum(x)))

bact.clr.long.RSR.mod_ord[bact.clr.long.RSR.mod_ord$Order == "Actinomycetales",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Actinomycetales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Actinomycetales",]$Abundance)*100)

bact.clr.long.RSR.mod_ord[bact.clr.long.RSR.mod_ord$Order == "Sphingomonadales",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Sphingomonadales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Sphingomonadales",]$Abundance)*100)

# Family
# Re-run models but dont send through anova and broom
bact.clr.long.RSR.mod_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.fam <-  phyloseq::tax_glom(drt.bact.late, "Family") # 84 taxa
drt.bact.late.fam.relab <- psmelt(transform_sample_counts(drt.bact.late.fam, function(x) x/sum(x)))

# Main effect significant
bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Polyangiaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Polyangiaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Polyangiaceae",]$Abundance)*100)

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Sinobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sinobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sinobacteraceae",]$Abundance)*100)

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Sphingomonadaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sphingomonadaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sphingomonadaceae",]$Abundance)*100)

# Genus
# Re-run models but dont send through anova and broom
bact.clr.long.RSR.mod_gen <- drt.bact.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.gen <-  phyloseq::tax_glom(drt.bact.late, "Genus") # 84 taxa
drt.bact.late.gen.relab <- psmelt(transform_sample_counts(drt.bact.late.gen, function(x) x/sum(x)))

# Main Effects
bact.clr.long.RSR.mod_gen[bact.clr.long.RSR.mod_gen$Genus == "Minicystis",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_gen %>% group_by(Genus) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Genus == "Minicystis")
mean((drt.bact.late.gen.relab[drt.bact.late.gen.relab$Genus == "Minicystis",]$Abundance)*100) 
