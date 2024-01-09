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

### Load dataset for use throughout ###
drt.fung.late.clr <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late_clr.RDS') # CLR normalized
drt.fung.late     <- readRDS('Intermediate_data/phyloseq_f_asv_clean_inoculated_late.RDS') # non-normalized
# Calculate RSR
sample_data(drt.fung.late.clr)$RootShootRatio <- sample_data(drt.fung.late.clr)$ThirdDryRootWeight / sample_data(drt.fung.late.clr)$ThirdDryShootWeight
sample_data(drt.fung.late)$RootShootRatio <- sample_data(drt.fung.late)$ThirdDryRootWeight / sample_data(drt.fung.late)$ThirdDryShootWeight
# Extract dataframe from phyloseq sample data
fung.samp.df.clr <- data.frame(sample_data(drt.fung.late.clr))


##### ShootMassRate #####
# Lmm with genotype + soilInoculum
mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = fung.samp.df.clr)
#add residuals of SMR to dataframe and phyloseq obj
fung.samp.df.clr$ShootMassRateResid <- resid(mod)
sample_data(drt.fung.late.clr)$ShootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualSMR 
drt.fung.late.clr_phy <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum")) # 5 taxa
drt.fung.late.clr_cla <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class")) # 12 taxa
drt.fung.late.clr_ord <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order")) # 14 taxa
drt.fung.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family")) # 24 taxa
# LMs
fung.clr.long.SMR.aovs_phy <- drt.fung.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.SMR.aovs_cla <- drt.fung.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.SMR.aovs_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.SMR.aovs_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
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

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.SMR.aovs_phy)
Padj.term.func("Abundance", fung.clr.long.SMR.aovs_phy)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.SMR.aovs_cla)
Padj.term.func("Abundance", fung.clr.long.SMR.aovs_cla)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.SMR.aovs_ord)
Padj.term.func("Abundance", fung.clr.long.SMR.aovs_ord)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.SMR.aovs_fam)
Padj.term.func("Abundance", fung.clr.long.SMR.aovs_fam)

# Class 
# Re-run models but dont send through anova and broom
fung.clr.long.SMR.mod_cla <- drt.fung.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.fung.late.cla <-  phyloseq::tax_glom(drt.fung.late, "Class")
drt.fung.late.cla.relab <- psmelt(transform_sample_counts(drt.fung.late.cla, function(x) x/sum(x)))
#STATS
fung.clr.long.SMR.mod_cla[fung.clr.long.SMR.mod_cla$Class == "Microbotryomycetes",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Microbotryomycetes")
mean((drt.fung.late.cla.relab[drt.fung.late.cla.relab$Class == "Microbotryomycetes",]$Abundance)*100) 

# Order
# Re-run models but dont send through anova and broom
fung.clr.long.SMR.mod_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.fung.late.ord <-  phyloseq::tax_glom(drt.fung.late, "Order")
drt.fung.late.ord.relab <- psmelt(transform_sample_counts(drt.fung.late.ord, function(x) x/sum(x)))
#STATS
fung.clr.long.SMR.mod_ord[fung.clr.long.SMR.mod_ord$Order == "Sporidiobolales",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Sporidiobolales")
mean((drt.fung.late.ord.relab[drt.fung.late.ord.relab$Order == "Sporidiobolales",]$Abundance)*100) 

# Family 
# Re-run models but dont send through anova and broom
fung.clr.long.SMR.mod_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.fung.late.fam <-  phyloseq::tax_glom(drt.fung.late, "Family")
drt.fung.late.fam.relab <- psmelt(transform_sample_counts(drt.fung.late.fam, function(x) x/sum(x)))
#STATS
fung.clr.long.SMR.mod_fam[fung.clr.long.SMR.mod_fam$Family == "Sporidiobolaceae",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sporidiobolaceae")
mean((drt.fung.late.fam.relab[drt.fung.late.fam.relab$Family == "Sporidiobolaceae",]$Abundance)*100)

##### RootMassRate #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(sqrt(RootMassRate) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = fung.samp.df.clr)
# add residuals of RMR to data frame
fung.samp.df.clr$RootMassRateResid <- resid(mod)
sample_data(drt.fung.late.clr)$RootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualRMR
# All
drt.fung.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum"))
drt.fung.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class"))
drt.fung.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order"))
drt.fung.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family"))
# LMs
fung.clr.long.RMR.aovs_phy <- drt.fung.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RMR.aovs_cla <- drt.fung.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RMR.aovs_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RMR.aovs_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(RootMassRateResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))

#Adjust P values across taxonomic level
Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RMR.aovs_phy)
Padj.term.func("Abundance", fung.clr.long.RMR.aovs_phy)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RMR.aovs_cla)
Padj.term.func("Abundance", fung.clr.long.RMR.aovs_cla)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RMR.aovs_ord)
Padj.term.func("Abundance", fung.clr.long.RMR.aovs_ord)

Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RMR.aovs_fam)
Padj.term.func("Abundance", fung.clr.long.RMR.aovs_fam)

# Family 
# Re-run models but dont send through anova and broom
fung.clr.long.RMR.mod_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(RootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.fung.late.fam <-  phyloseq::tax_glom(drt.fung.late, "Family")
drt.fung.late.fam.relab <- psmelt(transform_sample_counts(drt.fung.late.fam, function(x) x/sum(x)))
#STATS
fung.clr.long.RMR.mod_fam[fung.clr.long.RMR.mod_fam$Family == "Sporormiaceae",]$mod[[1]]$coefficients
fung.clr.long.RMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sporormiaceae")
mean((drt.fung.late.fam.relab[drt.fung.late.fam.relab$Family == "Sporormiaceae",]$Abundance)*100)

# STATS
car::Anova(fung.clr.long.RMR.mod_fam[fung.clr.long.RMR.mod_fam$Family == "Sporormiaceae",]$mod[[1]], type = "III")
mod <- fung.clr.long.RMR.mod_fam[fung.clr.long.RMR.mod_fam$Family == "Sporormiaceae",]$mod[[1]]
emtrends(mod, ~Drought.or.Watered, var = "Abundance")

Sporormiaceae_interaction_plot <- filter(drt.fung.late.clr_fam, Family == "Sporormiaceae") %>%
  ggplot(aes(x = Abundance, y = RootMassRateResid, color = Drought.or.Watered)) +
  xlab("CLR Transformed Abundance") +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment", values = treatment_pallete) +
  annotate("text", x = 3, y = 0.040, label = "W = -0.0023", fontface = 'italic') +
  annotate("text", x = 3, y = 0.035, label = "D =  0.0010", fontface = 'italic') +
  annotate("text", x = 3, y = 0.030, label = "PVEres = 9.66", fontface = 'italic') +
  annotate("text", x = 3, y = 0.025, label = "Abund.×DT q = 0.02", fontface = 'italic')

ggsave("figures/Sporormiaceae_interaction_plot.svg", Sporormiaceae_interaction_plot, height = 6, width = 6)
ggsave("figures/Sporormiaceae_interaction_plot.png", Sporormiaceae_interaction_plot, height = 6, width = 6)






# None Significant 
# ##### Root/Shoot Ratio #####
# # Lmm with genotype + soilInoculum 
# mod <- lmerTest::lmer(log(RootShootRatio) ~ Genotype + SoilLocation + logObs.z + (1|Block), data = fung.samp.df.clr)
# # add residuals of RSR to data frame
# fung.samp.df.clr$RootShootRatioResid <- resid(mod)
# sample_data(drt.fung.late.clr)$RootShootRatioResid <- resid(mod)
# # Correlate Taxon abundance to residualRSR
# drt.fung.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum"))
# drt.fung.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class"))
# drt.fung.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order"))
# drt.fung.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family"))
# # LMs
# fung.clr.long.RSR.aovs_phy <- drt.fung.late.clr_phy %>% nest_by(Phylum) %>%
#   mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
#   summarize(broom::tidy(mod))
# fung.clr.long.RSR.aovs_cla <- drt.fung.late.clr_cla %>% nest_by(Class) %>%
#   mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
#   summarize(broom::tidy(mod))
# fung.clr.long.RSR.aovs_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
#   mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
#   summarize(broom::tidy(mod))
# fung.clr.long.RSR.aovs_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
#   mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
#   summarize(broom::tidy(mod))
# 
# #Adjust P values across taxonomic level
# Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RSR.aovs_phy)
# Padj.term.func("Abundance", fung.clr.long.RSR.aovs_phy)
# 
# Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RSR.aovs_cla)
# Padj.term.func("Abundance", fung.clr.long.RSR.aovs_cla)
# 
# Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RSR.aovs_ord)
# Padj.term.func("Abundance", fung.clr.long.RSR.aovs_ord)
# 
# Padj.term.func("Abundance:Drought.or.Watered", fung.clr.long.RSR.aovs_fam)
# Padj.term.func("Abundance", fung.clr.long.RSR.aovs_fam)
