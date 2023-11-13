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
library(lme4); packageVersion('lme4')
library(lmerTest); packageVersion('lmerTest')
library(emmeans); packageVersion('emmeans')
library(ggpubr); packageVersion('ggpubr')
library(car); packageVersion('car')

# Theme set and Color Palettes
theme_set(theme_pubr())
genotype_pallete <- c("B73" = "#91ff26", "Mo17" = "#9426ff")# B73/Mo17 - Genotype
treatment_pallete <- c("Well-Watered" = "#0000FF", "Drought" = "#DAA520") # Drought/WW     - Treatment
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#7fd66f", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location


### Load dataset for use throughout
sampledata <- read.csv("../Raw_data_n_metadata/Drought_Experiment_Data_FULL.csv", header = TRUE)
### Remove Ag samples
sampledata <- sampledata[sampledata$SoilHabitat != "Agriculture",]
### Remove controls
sampledata_w_controls <- sampledata
sampledata <- sampledata[sampledata$SoilLocation != "Control",]

#### Emergence Rate ####
# No effect of drought since the drought treatment wasn't applied to the plants prior to germination
sampledata$Germination <- as.numeric(as.factor(sampledata$Germination))
Germ_mod <- glmer(Germination ~ Genotype * SoilLocation + (1|Block), data = sampledata)
car::Anova(Germ_mod, type = "III")
# Count # germinates by genotype (coded numerically, No = 1, yes = 2)
summary(as.factor(paste(sampledata$Genotype, sampledata$Germination, sep = "_")))
(200/(46 + 200)) * 100 #B73 81.3%
(150/(92 + 150)) * 100 #Mo17 62.0%

# Now remove non-germinants (coded numerically, No = 1, yes = 2)
sampledata_germinants <- sampledata[sampledata$Germination == 2,]
### Calculating Root shoot ratio
# Dry weight Roots / Dry weight Shoots
a <- sampledata_germinants$FirstDryRootWeight / sampledata_germinants$FirstDryShootWeight
b <- sampledata_germinants$SecondDryRootWeight / sampledata_germinants$SecondDryShootWeight
c <- sampledata_germinants$ThirdDryRootWeight / sampledata_germinants$ThirdDryShootWeight
df <- data.frame(cbind(a,b,c))
sampledata_germinants$RootShootRatio <- rowSums(df, na.rm = TRUE)
sampledata_germinants$RootShootRatio[sampledata_germinants$RootShootRatio == 0] <- NA
rm(a,b,c,df)

#### Re-work the factor levels to make plots fully write out treatment levels
sampledata_germinants <- sampledata_germinants %>%
  mutate(Drought.or.Watered = fct_recode(as.factor(Drought.or.Watered), 
                                     Drought = "D",
                                     `Well-Watered` = "W"))

#### Linear modeling on phenotypic data ####
RMR_mod <- lmer(sqrt(RootMassRate) ~ Drought.or.Watered * Genotype * SoilLocation + (1|Block), data = sampledata_germinants, na.action = na.exclude)
SMR_mod <- lmer(sqrt(ShootMassRate) ~ Drought.or.Watered * Genotype * SoilLocation + (1|Block), data = sampledata_germinants, na.action = na.exclude)
# control for timepoint of harvest here as we do not account for days post-germination in raw measure to get the ratio
RSR_mod <- lmer(log(RootShootRatio) ~ Drought.or.Watered * Genotype * SoilLocation + (1|Block) + (1|Timepoint), data = sampledata_germinants, na.action = na.exclude)
# anovas
anova(SMR_mod)
anova(RMR_mod)
anova(RSR_mod)
# Posthocs
# Drought main effects 
pairs(emmeans(SMR_mod, ~ Drought.or.Watered))
pairs(emmeans(RMR_mod, ~ Drought.or.Watered))
pairs(emmeans(RSR_mod, ~ Drought.or.Watered))
# Genotype
pairs(emmeans(SMR_mod, ~ Genotype))
pairs(emmeans(RMR_mod, ~ Genotype))
pairs(emmeans(RSR_mod, ~ Genotype))
# SoilInoculum
pairs(emmeans(SMR_mod, ~ SoilLocation))
pairs(emmeans(RSR_mod, ~ SoilLocation))
# Homogeneity of Variance (~equal above/below line)
plot(RMR_mod)
plot(SMR_mod)
plot(RSR_mod)
# qq plots
lattice::qqmath(RMR_mod, id=0.05)
lattice::qqmath(SMR_mod, id=0.05)
lattice::qqmath(RSR_mod, id=0.05)

# This tibble is used to place the significance letters at a consist height over each barplot
labels_df <- tibble(Drought.or.Watered=levels(as.factor(sampledata_germinants$Drought.or.Watered)),
                    Genotype=levels(as.factor(sampledata_germinants$Genotype)),
                    RMR=max(sampledata_germinants$RootMassRate, na.rm = TRUE) * 1.1,
                    SMR=max(sampledata_germinants$ShootMassRate, na.rm = TRUE) * 1.1,
                    RSR=max(sampledata_germinants$RootShootRatio, na.rm = TRUE) * 1.1)

#### Figure 3 and S7 ####
# Drought effects
a <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = ShootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Shoot Mass Rate (g/day)") + xlab("Treatment") + 
  geom_text(data=labels_df, aes(Drought.or.Watered, SMR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
b <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = RootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Root Mass Rate (g/day)") + xlab("Treatment")+ 
  geom_text(data=labels_df, aes(Drought.or.Watered, RMR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
c  <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = RootShootRatio, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Root/Shoot Ratio") + xlab("Treatment") + 
  geom_text(data=labels_df, aes(Drought.or.Watered, RSR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
part1 <- ggarrange(a,b,c, align = "hv", common.legend = TRUE, legend = 'right',labels = c("A","B","C"), nrow=1)
# Genotype effects
d <- ggplot(sampledata_germinants, aes(x = Genotype, y = ShootMassRate, fill = Genotype)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Shoot Mass Rate (g/day)") + xlab("Genotype")  + 
  geom_text(data=labels_df, aes(Genotype, SMR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
e <- ggplot(sampledata_germinants, aes(x = Genotype, y = RootMassRate, fill = Genotype)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Root Mass Rate (g/day)") + xlab("Genotype")  + 
  geom_text(data=labels_df, aes(Genotype, RMR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
f <- ggplot(sampledata_germinants, aes(x = Genotype, y = RootShootRatio, fill = Genotype)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Root/Shoot Ratio") + xlab("Genotype") + 
  geom_text(data=labels_df, aes(Genotype, RSR, label=c("a","b")), size = 6, color = '#3F3F3F') +
  theme(axis.title.x = element_blank())
part2 <- ggarrange(d,e,f, align = "hv", common.legend = TRUE, legend = 'right',labels = c("D","E","F"), nrow=1)
### Figure 3 ###
part1
ggsave("figures/figure3_phenotypic_treatment.svg", part1, width = 10, height = 5)
ggsave("figures/figure3_phenotypic_treatment.png", part1, width = 10, height = 5)
### Figure S7 ###
part2 <- ggarrange(d,e,f, align = "hv", common.legend = TRUE, legend = 'right',labels = "AUTO", nrow=1)
ggsave("figures/figureS6_phenotypic_genotype.svg", part2, width = 10, height = 5)
ggsave("figures/figureS6_phenotypic_genotype.png", part2, width = 10, height = 5)



#### Figure 6 and S16 ####
## Plots by soil inocula w/ controls
sampledata_w_controls
# Now remove non-germinants (coded numerically, No = 1, yes = 2)
sampledata_w_controls$Germination <- as.numeric(as.factor(sampledata_w_controls$Germination))
sampledata_w_controls_germinants <- sampledata_w_controls[sampledata_w_controls$Germination == 2,]
### Recalculating Root shoot ratio
# Dry weight Roots / Dry weight Shoots
a <- sampledata_w_controls_germinants$FirstDryRootWeight / sampledata_w_controls_germinants$FirstDryShootWeight
b <- sampledata_w_controls_germinants$SecondDryRootWeight / sampledata_w_controls_germinants$SecondDryShootWeight
c <- sampledata_w_controls_germinants$ThirdDryRootWeight / sampledata_w_controls_germinants$ThirdDryShootWeight
df <- data.frame(cbind(a,b,c))
sampledata_w_controls_germinants$RootShootRatio <- rowSums(df, na.rm = TRUE)
sampledata_w_controls_germinants$RootShootRatio[sampledata_w_controls_germinants$RootShootRatio == 0] <- NA
rm(a,b,c,df)
#### Re-work the factor levels to make plots fully write out treatment levels
sampledata_w_controls_germinants <- sampledata_w_controls_germinants %>%
  mutate(Drought.or.Watered = fct_recode(as.factor(Drought.or.Watered), 
                                         Drought = "D",
                                         `Well-Watered` = "W"))

# Reorder dry to wet sites
sampledata_w_controls_germinants$SoilLocation <- factor(sampledata_w_controls_germinants$SoilLocation, levels = c("Control","SVR", "HAY", "TLI", "KNZ"))
soil_merged_pallete_w_control <- c("Control" = "#FFFFFF","SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#7fd66f", "KNZ" = "#3ba150")
# Both plots together FACET 
# ShootMassRate
Figure6_A <- ggplot(sampledata_w_controls_germinants, aes(x = SoilLocation, y = ShootMassRate, fill = SoilLocation)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)") +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum") +
  theme(legend.text.align = 0, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Drought.or.Watered)
# RootMassRate
b <- ggplot(sampledata_w_controls_germinants, aes(x = SoilLocation, y = RootMassRate, fill = SoilLocation)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)") +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum") +
  theme(legend.text.align = 0, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Drought.or.Watered)
# RootShootRatio
c <- ggplot(sampledata_w_controls_germinants, aes(x = SoilLocation, y = RootShootRatio, fill = SoilLocation)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root/Shoot Ratio") +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum") +
  theme(legend.text.align = 0, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Drought.or.Watered) + ylim(0,12) # This removes 3 points, noted in legend
#### Figure S16 ####
b
c
RMR_RSR_plot <- ggarrange(b,c, nrow = 2, common.legend = TRUE, labels = "AUTO", align = 'hv', legend = "right")
# save
ggsave("figures/FigureS16_Controls_phenotypic_responses_RMR_RSR.svg", RMR_RSR_plot, height = 8, width = 10)
ggsave("figures/FigureS16_Controls_phenotypic_responses_RMR_RSR.png", RMR_RSR_plot,  height = 8, width = 10)
#### Figure 6BC Association of plant growth with Normal Annual Precipitation of the soil inocula sites####
# 30 year normal MAP 1991-2021
sampledata_with_MAP <- sampledata_germinants %>%
  group_by(SoilLocation) %>%
  mutate( MAP = case_when(
    SoilLocation == "SVR" ~ 478.9250,
    SoilLocation == "HAY" ~ 601.7375,
    SoilLocation == "TLI" ~ 783.8156,
    SoilLocation == "KNZ" ~ 912.3438))
# Quick plots
P1 <- ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank())
P2 <- ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Mass Rate (g/day)") + theme(axis.title.x = element_blank())
P3 <- ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Shoot Ratio") + xlab("Normal Annual Precipitation (mm)")
P4 <- ggarrange(P1,P2,P3, nrow = 3, align = 'hv', common.legend = TRUE, legend = 'right', labels = "AUTO")
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_MAP))
anova(lmerTest::lmer(sqrt(RootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_MAP))
anova(lmerTest::lmer(log(RootShootRatio) ~ MAP*Drought.or.Watered*Genotype + (1|Block) + (1|Timepoint), data = sampledata_with_MAP))
# SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_MAP)
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "MAP")
emtrends(SMR.native.mod, ~ Genotype, var = "MAP")
# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_MAP)
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "MAP")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "MAP")
# Panels B/C
P1 <- ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + annotate("text", x = 775, y= 0.022, label = "W = -2.88 e-6", fontface = 'italic') +
  annotate("text", x = 775, y= 0.020, label = "D =  4.53 e-7", fontface = 'italic') +
  annotate("text", x = 775, y= 0.018, label = "MAP×DT p = 0.04", fontface = 'italic')
P2 <- ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Normal Annual Precipitation (mm)") +
  ylim(0.000, 0.022) + annotate("text", x = 775, y=0.022, label = "B73 = -2.94 e-6", fontface = 'italic') +
  annotate("text", x = 775, y=0.020, label = "Mo17 = 5.10 e-7", fontface = 'italic') +
  annotate("text", x = 775, y=0.018, label = "MAP×G p = 0.02", fontface = 'italic') 
Figure6_BC <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'v', legend = 'right', labels = c('B', 'C'))
# add line plots
Controls_phenotypic_responses_with_line_plot <- ggarrange(Figure6_A, Figure6_BC, ncol = 2, widths = c(0.9, 0.5), legend = 'right', labels = "AUTO")
# save
ggsave("figures/Figure6_Controls_phenotypic_responses_facet_lineplots.svg", Controls_phenotypic_responses_with_line_plot, height = 8, width = 14)
ggsave("figures/Figure6_Controls_phenotypic_responses_facet_lineplots.png", Controls_phenotypic_responses_with_line_plot,  height = 8, width = 14)


#### Figure S17 Alternative Metrics for association with plant phenotypic measures####
##### Evapotranspiration of sites #####
sampledata_with_evap <- sampledata_germinants %>%
  group_by(SoilLocation) %>%
  mutate(EVAP = case_when(
    SoilLocation == "SVR" ~ 1408.700,
    SoilLocation == "HAY" ~ 1333.322,
    SoilLocation == "TLI" ~ 1266.691,
    SoilLocation == "KNZ" ~ 1170.869))
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ EVAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_evap))
# SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ EVAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_evap)
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "EVAP")
emtrends(SMR.native.mod, ~ Genotype, var = "EVAP")
# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ EVAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_evap)
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "EVAP")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "EVAP")
# Both interactions with genotype and treatment
P1 <- ggplot(sampledata_with_evap, aes(x = EVAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + 
  annotate("text", x = 1350, y=0.022, label = "W =  4.37 e-4", fontface = 'italic') +
  annotate("text", x = 1350, y=0.020, label = "D = -1.13 e-6", fontface = 'italic') +
  annotate("text", x = 1350, y=0.018, label = "ETo×DT p = 0.08", fontface = 'italic')
P2 <- ggplot(sampledata_with_evap, aes(x = EVAP, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Evapotranspiration") +
  ylim(0.000, 0.022) + 
  annotate("text", x = 1350, y=0.022, label = "B73 = 4.73 e-6", fontface = 'italic') +
  annotate("text", x = 1350, y=0.020, label = "Mo17 = -1.49 e-6", fontface = 'italic') +
  annotate("text", x = 1350, y=0.018, label = "ETo×G p = 0.03", fontface = 'italic') 
EVAP_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'hv', legend = "none", labels = c('A', 'B'))

# ##### Longitude of sites #####
# sampledata_with_long <- sampledata_germinants %>%
#   group_by(SoilLocation) %>%
#   mutate(LONG = case_when(
#     SoilLocation == "SVR" ~ -100.9951,
#     SoilLocation == "HAY" ~ -99.3033,
#     SoilLocation == "TLI" ~ -97.4690,
#     SoilLocation == "KNZ" ~ -96.6099))
# #STATS
# anova(lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_long))
# # SMR across just native soils
# SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_long)
# emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "LONG")
# emtrends(SMR.native.mod, ~ Genotype, var = "LONG")
# # Get non-transformed slopes for easier interpretation of results.
# SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_long)
# emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "LONG")
# emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "LONG")
# # Both interactions with genotype and treatment
# P1 <- ggplot(sampledata_with_long, aes(x = LONG, y = ShootMassRate, color = Drought.or.Watered)) +
#   geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
#   scale_color_manual(name = "Treatment" , values = treatment_pallete) +
#   ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
#   ylim(0.000, 0.022) + annotate("text", x = -98, y=0.022, label = "W = -3.29 e-4", fontface = 'italic') +
#   annotate("text", x = -98, y=0.020, label = "D = 3.04 e-5", fontface = 'italic') +
#   annotate("text", x = -98, y=0.018, label = "Long×DT p = 0.03", fontface = 'italic')
# P2 <- ggplot(sampledata_with_long, aes(x = LONG, y = ShootMassRate, color = Genotype)) +
#   geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
#   scale_color_manual(name = "Genotype" , values = genotype_pallete) +
#   ylab("Shoot Mass Rate (g/day)") + xlab("Longitude (Lower MAP to Higher MAP)") +
#   ylim(0.000, 0.022) + annotate("text", x = -98, y=0.022, label = "B73 = -3.07 e-4", fontface = 'italic') +
#   annotate("text", x = -98, y=0.020, label = "Mo17 = 9.12 e-6", fontface = 'italic') +
#   annotate("text", x = -98, y=0.018, label = "Long×G p = 0.03", fontface = 'italic') 
# # combine
# Long_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'hv', legend = "none", labels = c('A', 'B'))

##### Aridity index #####
sampledata_with_ARID <- sampledata_germinants %>%
  group_by(SoilLocation) %>%
  mutate(Aridity_index = case_when(
    SoilLocation == "SVR" ~ 0.3399766,
    SoilLocation == "HAY" ~ 0.4513070,
    SoilLocation == "TLI" ~ 0.6187901,
    SoilLocation == "KNZ" ~ 0.7792024))
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ Aridity_index*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_ARID))
# SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Aridity_index*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_ARID)
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "Aridity_index")
emtrends(SMR.native.mod, ~ Genotype, var = "Aridity_index")
# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ Aridity_index*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_ARID)
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "Aridity_index")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "Aridity_index")
# Both interactions with genotype and treatment
P1 <- ggplot(sampledata_with_ARID, aes(x = Aridity_index, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + annotate("text", x = 0.65, y=0.022, label = "W = -2.51 e-3", fontface = 'italic') +
  annotate("text", x = 0.65, y=0.020, label = "D = 5.44 e-4", fontface = 'italic') +
  annotate("text", x = 0.65, y=0.018, label = "Arid×DT p = 0.06", fontface = 'italic')
P2 <- ggplot(sampledata_with_ARID, aes(x = Aridity_index, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Aridity Index") +
  ylim(0.000, 0.022) + annotate("text", x = 0.65, y=0.022, label = "B73 = -2.72 e-3", fontface = 'italic') +
  annotate("text", x = 0.65, y=0.020, label = "Mo17 = 7.49 e-4", fontface = 'italic') +
  annotate("text", x = 0.65, y=0.018, label = "Arid×G p = 0.02", fontface = 'italic') 
# combine
Arid_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'hv', legend = "none", labels = c('C', 'D'))

##### Soil Moisture #####
sampledata_with_soilmoisture <- sampledata_germinants %>%
  group_by(SoilLocation) %>%
  mutate(SoilMoisture = case_when(
    SoilLocation == "SVR" ~ 37.38437,
    SoilLocation == "HAY" ~ 92.87812,
    SoilLocation == "TLI" ~ 177.09688,
    SoilLocation == "KNZ" ~ 348.25937))
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ SoilMoisture*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_soilmoisture))
# SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ SoilMoisture*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_soilmoisture)
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "SoilMoisture")
emtrends(SMR.native.mod, ~ Genotype, var = "SoilMoisture")
# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ SoilMoisture*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_soilmoisture)
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "SoilMoisture")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "SoilMoisture")
# Both interactions with genotype and treatment
P1 <- ggplot(sampledata_with_soilmoisture, aes(x = SoilMoisture, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + annotate("text", x = 275, y=0.022, label = "W = -1.92 e-6", fontface = 'italic') +
  annotate("text", x = 275, y=0.020, label = "D = 1.11 e-6", fontface = 'italic') +
  annotate("text", x = 275, y=0.018, label = "SM×DT p = 0.18", fontface = 'italic')
P2 <- ggplot(sampledata_with_soilmoisture, aes(x = SoilMoisture, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Normal Annual Soil Moisture (mm)") +
  ylim(0.000, 0.022) + annotate("text", x = 275, y=0.022, label = "B73 = -2.82 e-6", fontface = 'italic') +
  annotate("text", x = 275, y=0.020, label = "Mo17 = 2.01 e-6", fontface = 'italic') +
  annotate("text", x = 275, y=0.018, label = "SM×G p = 0.03", fontface = 'italic') 
# combine
Soil_moisture_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'hv', legend = "none", labels = c('E', 'F'))

##### PDSI #####
sampledata_with_PDSI <- sampledata_germinants %>%
  group_by(SoilLocation) %>%
  mutate(PDSI = case_when(
    SoilLocation == "SVR" ~ 0.2431771,
    SoilLocation == "HAY" ~ 0.6842187,
    SoilLocation == "TLI" ~ 0.8786198,
    SoilLocation == "KNZ" ~ 0.3597917))
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ PDSI*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_PDSI))
# SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ PDSI*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_PDSI)
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "PDSI")
emtrends(SMR.native.mod, ~ Genotype, var = "PDSI")
# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ PDSI*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_PDSI)
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "PDSI")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "PDSI")
# Both interactions with genotype and treatment
P1 <- ggplot(sampledata_with_PDSI, aes(x = PDSI, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + annotate("text", x = 0.7, y=0.022, label = "W = -3.31 e-3", fontface = 'italic') +
  annotate("text", x = 0.7, y=0.020, label = "D = -7.2 e-4", fontface = 'italic') +
  annotate("text", x = 0.7, y=0.018, label = "PDSI×DT p = 0.03", fontface = 'italic') +
  annotate(geom = "point", x = 0.2431771, y = 0.016357143, color = "#de3a68", size = 5, alpha = 0.8) + 
  annotate(geom = "point", x = 0.2431771, y = 0.016357143, color = "#0000FF") + 
  annotate(geom = "text", x = 0.2531771, y = 0.016357143, label = "SVR", hjust = "left") +
  annotate(geom = "point", x = 0.3597917, y = 0.01434211, color = "#3ba150", size = 5, alpha = 0.8) + 
  annotate(geom = "point", x = 0.3597917, y = 0.01434211, color = "#0000FF") + 
  annotate(geom = "text", x = 0.3697917, y = 0.01434211, label = "KNZ", hjust = "left") +
  annotate(geom = "point", x = 0.6842187, y = 0.01464865, color = "#f5805d", size = 5, alpha = 0.8) + 
  annotate(geom = "point", x = 0.6842187, y = 0.01464865, color = "#0000FF") + 
  annotate(geom = "text", x = 0.6742187, y = 0.01464865, label = "HAY", hjust = "right") +
  annotate(geom = "point", x = 0.8786198, y = 0.015375, color = "#7fd66f", size = 5, alpha = 0.8) + 
  annotate(geom = "point", x = 0.8786198, y = 0.015375, color = "#0000FF") + 
  annotate(geom = "text", x = 0.8686198, y = 0.015375, label = "TLI", hjust = "right")
P2 <- ggplot(sampledata_with_PDSI, aes(x = PDSI, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Normal Monthly PDSI") +
  ylim(0.000, 0.022) + annotate("text", x = 0.7, y=0.022, label = "B73 = -2.09 e-3", fontface = 'italic') +
  annotate("text", x = 0.7, y=0.020, label = "Mo17 = -1.94 e-3", fontface = 'italic') +
  annotate("text", x = 0.7, y=0.018, label = "PDSI×G p = 0.92", fontface = 'italic')
# combine
PDSI_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'hv', legend = 'none', labels = c('G', 'H'))

## Put together all the plots side by side
Alt_metrics_smooth <- ggarrange(EVAP_smooth, Arid_smooth, Soil_moisture_smooth, PDSI_smooth, ncol = 4, nrow = 1, common.legend = TRUE)
# Save them
ggsave("./figures/Figure_S17_Alt_metrics_SMR_Treatment_Genotype.svg", Alt_metrics_smooth, height = 6, width = 16)
ggsave("./figures/Figure_S17_Alt_metrics_SMR_Treatment_Genotype.png", Alt_metrics_smooth, height = 6, width = 16)
rm(P1, P2, Alt_metrics_smooth, PDSI_smooth, Arid_smooth, Soil_moisture_smooth)

#### Figure S19 Root/shoot mass rate - boxplots ####
sampledata_germinants$RootMassRate # third timepoint
sampledata_removed_NAs <- filter(sampledata_germinants, !(Timepoint == ""))
sampledata_removed_NAs$Timepoint <- factor(sampledata_removed_NAs$Timepoint, levels = c("early", "middle", "late"))

P1 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = RootMassRate)) +
  geom_jitter(width = 0.3) + geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  ylab("Root Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late")) 

P2 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = ShootMassRate)) +
  geom_jitter(width = 0.3) + geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  ylab("Shoot Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late"))

Plantmassrates_by_time <- ggarrange(P1, P2, labels = "AUTO")
#stats
mod.root  <- lmerTest::lmer(sqrt(RootMassRate) ~ Timepoint*Drought.or.Watered + (1|Block), data = sampledata_removed_NAs)
mod.shoot <- lmerTest::lmer(sqrt(ShootMassRate) ~ Timepoint*Drought.or.Watered + (1|Block), data = sampledata_removed_NAs)
anova(mod.root)
anova(mod.shoot)
emmeans(mod.root,~Timepoint|Drought.or.Watered)
emmeans(mod.shoot,~Timepoint)
ggsave("./figures/S19_Plant_mass_rates_by_time.svg", Plantmassrates_by_time, height = 4, width = 8)
ggsave("./figures/S19_Plant_mass_rates_by_time.png", Plantmassrates_by_time, height = 4, width = 8)
