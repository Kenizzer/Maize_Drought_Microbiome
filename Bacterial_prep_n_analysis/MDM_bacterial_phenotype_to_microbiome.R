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
mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = bact.samp.df.clr)
#add residuals of SMR to dataframe and phyloseq obj
bact.samp.df.clr$ShootMassRateResid <- resid(mod)
sample_data(drt.bact.late.clr)$ShootMassRateResid <- resid(mod)
# Correlate Taxon abundance to residualSMR 
drt.bact.late.clr_phy <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Phylum")) # 13 taxa
drt.bact.late.clr_cla <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Class")) # 21 taxa
drt.bact.late.clr_ord <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Order")) # 37 taxa
drt.bact.late.clr_fam <- psmelt(phyloseq::tax_glom(drt.bact.late.clr, "Family")) # 84 taxa
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


# Adjust P values across taxonomic level but separately per term
# Get corrected P values from here
sum(p.adjust(bact.clr.long.SMR.aovs_phy[bact.clr.long.SMR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.SMR.aovs_phy[bact.clr.long.SMR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.SMR.aovs_phy[bact.clr.long.SMR.aovs_phy$term == "Abundance", ][c(temp), ]

sum(p.adjust(bact.clr.long.SMR.aovs_cla[bact.clr.long.SMR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.SMR.aovs_cla[bact.clr.long.SMR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.SMR.aovs_cla[bact.clr.long.SMR.aovs_cla$term == "Abundance", ][c(temp), ]

sum(p.adjust(bact.clr.long.SMR.aovs_ord[bact.clr.long.SMR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.SMR.aovs_ord[bact.clr.long.SMR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.SMR.aovs_ord[bact.clr.long.SMR.aovs_ord$term == "Abundance", ][c(temp), ]

temp <- which(p.adjust(bact.clr.long.SMR.aovs_fam[bact.clr.long.SMR.aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
  bact.clr.long.SMR.aovs_fam[bact.clr.long.SMR.aovs_fam$term == "Abundance:Drought.or.Watered", ][c(temp), ]
temp <- which(p.adjust(bact.clr.long.SMR.aovs_fam[bact.clr.long.SMR.aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)
  bact.clr.long.SMR.aovs_fam[bact.clr.long.SMR.aovs_fam$term == "Abundance", ][c(temp), ]


# phylum 
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.phy <-  phyloseq::tax_glom(drt.bact.late, "Phylum") # 13 taxa
drt.bact.late.phy.relab <- psmelt(transform_sample_counts(drt.bact.late.phy, function(x) x/sum(x)))

bact.clr.long.SMR.mod_phy[bact.clr.long.SMR.mod_phy$Phylum == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_phy %>% group_by(Phylum) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Phylum == "Actinobacteria")
mean((drt.bact.late.phy.relab[drt.bact.late.phy.relab$Phylum == "Actinobacteria",]$Abundance)*100)


# class
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.cla <-  phyloseq::tax_glom(drt.bact.late, "Class") # 21 taxa
drt.bact.late.cla.relab <- psmelt(transform_sample_counts(drt.bact.late.cla, function(x) x/sum(x)))

bact.clr.long.SMR.mod_cla[bact.clr.long.SMR.mod_cla$Class == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Actinobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Actinobacteria",]$Abundance)*100) 

bact.clr.long.SMR.mod_cla[bact.clr.long.SMR.mod_cla$Class == "Deltaproteobacteria",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Deltaproteobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Deltaproteobacteria",]$Abundance)*100) 


# Order
# Re-run models but dont send through anova and broom
bact.clr.long.SMR.mod_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.ord <-  phyloseq::tax_glom(drt.bact.late, "Order") # 37 taxa
drt.bact.late.ord.relab <- psmelt(transform_sample_counts(drt.bact.late.ord, function(x) x/sum(x)))

bact.clr.long.SMR.mod_ord[bact.clr.long.SMR.mod_ord$Order == "Actinomycetales",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Actinomycetales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Actinomycetales",]$Abundance)*100) 

bact.clr.long.SMR.mod_ord[bact.clr.long.SMR.mod_ord$Order == "Myxococcales",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Myxococcales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Myxococcales",]$Abundance)*100) 

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

#interactions significant 
bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Corynebacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Corynebacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Corynebacteriaceae",]$Abundance)*100) 

# STATS
car::Anova(bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Corynebacteriaceae",]$mod[[1]], type = "III")
mod <- bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Corynebacteriaceae",]$mod[[1]]
emtrends(mod, ~Drought.or.Watered, var = "Abundance")

Cory_interaction_plot <- filter(drt.bact.late.clr_fam, Family == "Corynebacteriaceae") %>%
  ggplot(aes(x = Abundance, y = ShootMassRateResid, fill = Drought.or.Watered, color = Drought.or.Watered)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_fill_manual(name = "Treatment", values = treatment_pallete) +
  scale_color_manual(name = "Treatment", values = treatment_pallete) +
  annotate("text", x = 10, y = 0.040, label = "W = -0.0032", fontface = 'italic') +
  annotate("text", x = 10, y = 0.035, label = "D = -0.0001", fontface = 'italic') +
  annotate("text", x = 10, y = 0.030, label = "PVEres = 4.79", fontface = 'italic') +
  annotate("text", x = 10, y = 0.025, label = "Abund.×DT q <0.01", fontface = 'italic')
  

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Streptococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Streptococcaceae",]$Abundance)*100) 

# STATS
car::Anova(bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]], type = "III")
mod <- bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]
emtrends(mod, ~Drought.or.Watered, var = "Abundance")

Strep_interaction_plot <- filter(drt.bact.late.clr_fam, Family == "Streptococcaceae") %>%
  ggplot(aes(x = Abundance, y = ShootMassRateResid, fill = Drought.or.Watered, color = Drought.or.Watered)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_fill_manual(name = "Treatment", values = treatment_pallete) +
  scale_color_manual(name = "Treatment", values = treatment_pallete) +
  annotate("text", x = 10, y = 0.040, label = "W = -0.0031", fontface = 'italic') +
  annotate("text", x = 10, y = 0.035, label = "D = -0.0002", fontface = 'italic') +
  annotate("text", x = 10, y = 0.030, label = "PVEres = 3.29", fontface = 'italic') +
  annotate("text", x = 10, y = 0.025, label = "Abund.×DT q <0.01", fontface = 'italic')
  
p1 <- ggarrange(Cory_interaction_plot, Strep_interaction_plot, labels = "AUTO", nrow = 1, align = "h", common.legend = TRUE, legend = 'right')
ggsave("figures/Cory_Strep_interaction_plot.svg", height = 6, width = 12)
ggsave("figures/Cory_Strep_interaction_plot.png", height = 6, width = 12)


# Main effect significant
bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Actinomycetaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Actinomycetaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Actinomycetaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Alcaligenaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Alcaligenaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Alcaligenaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Beijerinckiaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Beijerinckiaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Beijerinckiaceae",]$Abundance)*100) 

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

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Rhodospirillaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Rhodospirillaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Rhodospirillaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Sinobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sinobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sinobacteraceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Sphingomonadaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sphingomonadaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sphingomonadaceae",]$Abundance)*100) 

bact.clr.long.SMR.mod_fam[bact.clr.long.SMR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]$coefficients
bact.clr.long.SMR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Streptococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Streptococcaceae",]$Abundance)*100) 



##### RootMassRate #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(sqrt(RootMassRate) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = bact.samp.df.clr)
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

#Adjust P values across taxonomic level
sum(p.adjust(bact.clr.long.RMR.aovs_phy[bact.clr.long.RMR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(bact.clr.long.RMR.aovs_phy[bact.clr.long.RMR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(bact.clr.long.RMR.aovs_cla[bact.clr.long.RMR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(bact.clr.long.RMR.aovs_cla[bact.clr.long.RMR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(bact.clr.long.RMR.aovs_ord[bact.clr.long.RMR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(bact.clr.long.RMR.aovs_ord[bact.clr.long.RMR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(bact.clr.long.RMR.aovs_fam[bact.clr.long.RMR.aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(bact.clr.long.RMR.aovs_fam[bact.clr.long.RMR.aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)

sum(p.adjust(bact.clr.long.RMR.aovs_gen[bact.clr.long.RMR.aovs_gen$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(bact.clr.long.RMR.aovs_gen[bact.clr.long.RMR.aovs_gen$term == "Abundance", ]$p.value, method = "BH") < 0.05)



##### Root/Shoot Ratio #####
# Lmm with genotype + soilInoculum 
mod <- lmerTest::lmer(log(RootShootRatio) ~ Genotype + SoilInoculum + logObs.z + (1|Block), data = bact.samp.df.clr)
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
sum(p.adjust(bact.clr.long.RSR.aovs_phy[bact.clr.long.RSR.aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.RSR.aovs_phy[bact.clr.long.RSR.aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.RSR.aovs_phy[bact.clr.long.RSR.aovs_phy$term == "Abundance", ][c(temp), ]

sum(p.adjust(bact.clr.long.RSR.aovs_cla[bact.clr.long.RSR.aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.RSR.aovs_cla[bact.clr.long.RSR.aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.RSR.aovs_cla[bact.clr.long.RSR.aovs_cla$term == "Abundance", ][c(temp), ]

sum(p.adjust(bact.clr.long.RSR.aovs_ord[bact.clr.long.RSR.aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.RSR.aovs_ord[bact.clr.long.RSR.aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.RSR.aovs_ord[bact.clr.long.RSR.aovs_ord$term == "Abundance", ][c(temp), ]

sum(p.adjust(bact.clr.long.RSR.aovs_fam[bact.clr.long.RSR.aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
temp <- which(p.adjust(bact.clr.long.RSR.aovs_fam[bact.clr.long.RSR.aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)
bact.clr.long.RSR.aovs_fam[bact.clr.long.RSR.aovs_fam$term == "Abundance", ][c(temp),]


#phylum
# Re-run models but don't send through anova and broom
bact.clr.long.RSR.mod_phy <- drt.bact.late.clr_phy %>% nest_by(Phylum) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.phy <-  phyloseq::tax_glom(drt.bact.late, "Phylum") # 13 taxa
drt.bact.late.phy.relab <- psmelt(transform_sample_counts(drt.bact.late.phy, function(x) x/sum(x)))

bact.clr.long.RSR.mod_phy[bact.clr.long.RSR.mod_phy$Phylum == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_phy %>% group_by(Phylum) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Phylum == "Actinobacteria")
mean((drt.bact.late.phy.relab[drt.bact.late.phy.relab$Phylum == "Actinobacteria",]$Abundance)*100) 


#class
# Re-run models but don't send through anova and broom
bact.clr.long.RSR.mod_cla <- drt.bact.late.clr_cla %>% nest_by(Class) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.cla <-  phyloseq::tax_glom(drt.bact.late, "Class") # 21 taxa
drt.bact.late.cla.relab <- psmelt(transform_sample_counts(drt.bact.late.cla, function(x) x/sum(x)))

bact.clr.long.RSR.mod_cla[bact.clr.long.RSR.mod_cla$Class == "Actinobacteria",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_cla %>% group_by(Class) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Class == "Actinobacteria")
mean((drt.bact.late.cla.relab[drt.bact.late.cla.relab$Class == "Actinobacteria",]$Abundance)*100) 


#order
# Re-run models but don't send through anova and broom
bact.clr.long.RSR.mod_ord <- drt.bact.late.clr_ord %>% nest_by(Order) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.ord <-  phyloseq::tax_glom(drt.bact.late, "Order") # 37 taxa
drt.bact.late.ord.relab <- psmelt(transform_sample_counts(drt.bact.late.ord, function(x) x/sum(x)))

bact.clr.long.RSR.mod_ord[bact.clr.long.RSR.mod_ord$Order == "Actinomycetales",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Actinomycetales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Actinomycetales",]$Abundance)*100) 

bact.clr.long.RSR.mod_ord[bact.clr.long.RSR.mod_ord$Order == "Sphingomonadales",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_ord %>% group_by(Order) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Order == "Sphingomonadales")
mean((drt.bact.late.ord.relab[drt.bact.late.ord.relab$Order == "Sphingomonadales",]$Abundance)*100) 


#Family
# Re-run models but don't send through anova and broom
bact.clr.long.RSR.mod_fam <- drt.bact.late.clr_fam %>% nest_by(Family) %>%
  mutate(mod = list(lm(RootShootRatioResid ~ Abundance*Drought.or.Watered, data=data)))
# Relative abundance calculation and extraction
drt.bact.late.fam <-  phyloseq::tax_glom(drt.bact.late, "Family") # 37 taxa
drt.bact.late.fam.relab <- psmelt(transform_sample_counts(drt.bact.late.fam, function(x) x/sum(x)))

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Dietziaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Dietziaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Dietziaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Microbacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Microbacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Microbacteriaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Micrococcaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Micrococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Micrococcaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Mycobacteriaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Mycobacteriaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Mycobacteriaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Polyangiaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Polyangiaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Polyangiaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Rhodospirillaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Rhodospirillaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Rhodospirillaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Sinobacteraceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sinobacteraceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sinobacteraceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Sphingomonadaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Sphingomonadaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Sphingomonadaceae",]$Abundance)*100) 

bact.clr.long.RSR.mod_fam[bact.clr.long.RSR.mod_fam$Family == "Streptococcaceae",]$mod[[1]]$coefficients
bact.clr.long.RSR.aovs_fam %>% group_by(Family) %>% mutate(PVE = (sumsq / sum(sumsq)) * 100) %>% filter(Family == "Streptococcaceae")
mean((drt.bact.late.fam.relab[drt.bact.late.fam.relab$Family == "Streptococcaceae",]$Abundance)*100) 