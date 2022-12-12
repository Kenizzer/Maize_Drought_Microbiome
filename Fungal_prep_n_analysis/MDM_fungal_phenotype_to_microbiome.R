# Analysis of maize B73 and Mo17 inoculated plants under drought stress
# Samples were collected from a greenhouse study in Feb 2020
# Microbiomes are from the root compartment (rhizosphere/endosphere)
# Growth measurements were also recorded throughout the experiment (~50 days in length)
# Code by: Joel Swift

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
habitat_pallete <- c("Agriculture" = "#332288", "Native" = "#44AA99") # Native/Ag.     - Soil habitat
location_pallete <- c("SVR" = "#88CCEE", "HAY" = "#CC6677", "TLI" = "#DDCC77", "KNZ" = "#117733") # SVR/HAY/TLI/KNZ - Soil location
soil_merged_pallete<- c("SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") # Combination of Soil habitat and location, colors from https://lospec.com/palette-list/zencillo14

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
mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = fung.samp.df.clr)
#add residuals of SMR to dataframe and phyloseq obj
fung.samp.df.clr$ShootMassRateResid <- resid(mod)
sample_data(drt.fung.late.clr)$ShootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualSMR 
drt.fung.late.clr_phy <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum")) # 5 taxa
drt.fung.late.clr_cla <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class")) # 12 taxa
drt.fung.late.clr_ord <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order")) # 16 taxa
drt.fung.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family")) # 30 taxa
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
sum(p.adjust(fung.clr.long.SMR.aovs_phy[fung.clr.long.SMR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.SMR.aovs_phy[fung.clr.long.SMR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.SMR.aovs_cla[fung.clr.long.SMR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(fung.clr.long.SMR.aovs_cla[fung.clr.long.SMR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)
fung.clr.long.SMR.aovs_cla[fung.clr.long.SMR.aovs_cla$term == "Abundance", ][c(temp), ]

sum(p.adjust(fung.clr.long.SMR.aovs_ord[fung.clr.long.SMR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(fung.clr.long.SMR.aovs_ord[fung.clr.long.SMR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)
fung.clr.long.SMR.aovs_ord[fung.clr.long.SMR.aovs_ord$term == "Abundance", ][c(temp), ]
# 3


# Order
# Re-run models but dont send through anova and broom
fung.clr.long.SMR.mod_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.fung.late.ord <-  phyloseq::tax_glom(drt.fung.late, "Order") # 84 taxa
drt.fung.late.ord.relab <- psmelt(transform_sample_counts(drt.fung.late.ord, function(x) x/sum(x)))
#STATS
fung.clr.long.SMR.mod_ord[fung.clr.long.SMR.mod_ord$Order == "Hypocreales",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Hypocreales")
mean((drt.fung.late.ord.relab[drt.fung.late.ord.relab$Order == "Hypocreales",]$Abundance)*100) 

fung.clr.long.SMR.mod_ord[fung.clr.long.SMR.mod_ord$Order == "Malasseziales",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Malasseziales")
mean((drt.fung.late.ord.relab[drt.fung.late.ord.relab$Order == "Malasseziales",]$Abundance)*100) 

fung.clr.long.SMR.mod_ord[fung.clr.long.SMR.mod_ord$Order == "Sporidiobolales",]$mod[[1]]$coefficients
fung.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Sporidiobolales")
mean((drt.fung.late.ord.relab[drt.fung.late.ord.relab$Order == "Sporidiobolales",]$Abundance)*100) 


##### RootMassRate #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(sqrt(RootMassRate) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = fung.samp.df.clr)
# add residuals of RMR to data frame
fung.samp.df.clr$RootMassRateResid <- resid(mod)
sample_data(drt.fung.late.clr)$RootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualRMR
# All
drt.fung.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum"))
drt.fung.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class"))
drt.fung.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order"))
drt.fung.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family"))
drt.fung.late.clr_gen <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Genus")) 
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
sum(p.adjust(fung.clr.long.RMR.aovs_phy[fung.clr.long.RMR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RMR.aovs_phy[fung.clr.long.RMR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RMR.aovs_cla[fung.clr.long.RMR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RMR.aovs_cla[fung.clr.long.RMR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RMR.aovs_ord[fung.clr.long.RMR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RMR.aovs_ord[fung.clr.long.RMR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RMR.aovs_fam[fung.clr.long.RMR.aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RMR.aovs_fam[fung.clr.long.RMR.aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)


##### Root/Shoot Ratio #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(log(RootShootRatio) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = fung.samp.df.clr)
# add residuals of RSR to data frame
fung.samp.df.clr$RootShootRatioResid <- resid(mod)
sample_data(drt.fung.late.clr)$RootShootRatioResid <- resid(mod)
# Correlate Taxon abundance to residualRSR
drt.fung.late.clr_phy <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Phylum"))
drt.fung.late.clr_cla <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Class"))
drt.fung.late.clr_ord <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Order"))
drt.fung.late.clr_fam <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Family"))
drt.fung.late.clr_gen <-  psmelt(phyloseq::tax_glom(drt.fung.late.clr, "Genus")) 
# LMs
fung.clr.long.RSR.aovs_phy <- drt.fung.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RSR.aovs_cla <- drt.fung.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RSR.aovs_ord <- drt.fung.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RSR.aovs_fam <- drt.fung.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
fung.clr.long.RSR.aovs_gen <- drt.fung.late.clr_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered , data=data)))) %>%
  summarize(broom::tidy(mod))
#Adjust P values across taxonomic level
sum(p.adjust(fung.clr.long.RSR.aovs_phy[fung.clr.long.RSR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RSR.aovs_phy[fung.clr.long.RSR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RSR.aovs_cla[fung.clr.long.RSR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RSR.aovs_cla[fung.clr.long.RSR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RSR.aovs_ord[fung.clr.long.RSR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RSR.aovs_ord[fung.clr.long.RSR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(fung.clr.long.RSR.aovs_fam[fung.clr.long.RSR.aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(fung.clr.long.RSR.aovs_fam[fung.clr.long.RSR.aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)