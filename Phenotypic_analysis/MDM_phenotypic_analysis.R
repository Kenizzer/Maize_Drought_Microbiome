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
library(lmerTest); packageVersion('lmerTest')
library(emmeans); packageVersion('emmeans')
library(ggpubr); packageVersion('ggpubr')
library(car); packageVersion('car')

# Theme set and Color Palettes
theme_set(theme_pubr())
genotype_pallete <- c("B73" = "#91ff26", "Mo17" = "#9426ff")# B73/Mo17 - Genotype
treatment_pallete <- c("Well-Watered" = "#0000FF", "Drought" = "#DAA520") # Drought/WW     - Treatment
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#ffe785", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location


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
Germ_mod <- lmer(Germination ~ Genotype * SoilLocation + (1|Block), data = sampledata)
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

#### Within-inoculum Coefficient of variation for SMR in DRT vs in WW ####
cv <- function(x, na.rm = FALSE) (sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm))*100

sampledata_germinants %>%
  select(Drought.or.Watered, SoilLocation, ShootMassRate) %>%
  group_by(Drought.or.Watered, SoilLocation) %>%
  summarise(ShootMassRate_cv = cv(ShootMassRate, na.rm = TRUE))

sampledata_germinants %>%
  select(Drought.or.Watered, SoilLocation, ShootMassRate) %>%
  group_by(Drought.or.Watered, SoilLocation) %>%
  summarise(ShootMassRate_mean = mean(ShootMassRate, na.rm = TRUE))


#### Figure 3 and S6 ####
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
### Figure S6 ###
part2 <- ggarrange(d,e,f, align = "hv", common.legend = TRUE, legend = 'right',labels = "AUTO", nrow=1)
ggsave("figures/figureS6_phenotypic_genotype.svg", part2, width = 10, height = 5)
ggsave("figures/figureS6_phenotypic_genotype.png", part2, width = 10, height = 5)



#### Figure 6 and S7 ####
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
soil_merged_pallete_w_control <- c("Control" = "#FFFFFF","SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#ffe785", "KNZ" = "#3ba150")
# Both plots together FACET 
# ShootMassRate
Figure6a <- ggplot(sampledata_w_controls_germinants, aes(x = SoilLocation, y = ShootMassRate, fill = SoilLocation)) +
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
#### Figure S7 ####
SMR_RMR_RSR_plot <- ggarrange(Figure6a,b,c, nrow = 3, common.legend = TRUE, labels = "AUTO", align = 'hv', legend = "right")
# save
ggsave("figures/FigureS7_Controls_phenotypic_responses_SMR_RMR_RSR.svg", SMR_RMR_RSR_plot, height = 8, width = 10)
ggsave("figures/FigureS7_Controls_phenotypic_responses_SMR_RMR_RSR.png", SMR_RMR_RSR_plot,  height = 8, width = 10)
#### Figure 6A/B Association of plant growth with aridity index and PDSI of the soil inocula sites####
# Base models
SMR.native.mod <- lmer(sqrt(ShootMassRate) ~ SoilLocation*Drought.or.Watered*Genotype + (1|Block), data = sampledata_germinants)
RMR.native.mod <- lmer(sqrt(RootMassRate) ~ SoilLocation*Drought.or.Watered*Genotype + (1|Block), data = sampledata_germinants)
RSR.native.mod <- lmer(log(RootShootRatio) ~ SoilLocation*Drought.or.Watered*Genotype + (1|Block) + (1|Timepoint), data = sampledata_germinants)

# Function to extract estimated marginal means from the model and add climate data
Model_to_DF_w_climate <- function(model){
  # Extract emmeans and back transform.
  mod.df <- summary(emmeans(model, ~SoilLocation*Drought.or.Watered*Genotype, type = "response"))
  # Add climate data to the DF, matching to soil locations
  clim.mod.df <- mod.df %>%
    group_by(SoilLocation) %>%
    mutate(MAP = case_when( SoilLocation == "SVR" ~ 478.9250, SoilLocation == "HAY" ~ 601.7375,
                            SoilLocation == "TLI" ~ 783.8156, SoilLocation == "KNZ" ~ 912.3438),
           EVAP = case_when(SoilLocation == "SVR" ~ 1408.700, SoilLocation == "HAY" ~ 1333.322,
                            SoilLocation == "TLI" ~ 1266.691, SoilLocation == "KNZ" ~ 1170.869),
           Aridity_index = case_when(SoilLocation == "SVR" ~ 0.3399766, SoilLocation == "HAY" ~ 0.4513070,
                                     SoilLocation == "TLI" ~ 0.6187901, SoilLocation == "KNZ" ~ 0.7792024),
           SoilMoisture = case_when(SoilLocation == "SVR" ~ 37.38437, SoilLocation == "HAY" ~ 92.87812,
                                    SoilLocation == "TLI" ~ 177.09688, SoilLocation == "KNZ" ~ 348.25937),
           PDSI = case_when(SoilLocation == "SVR" ~ 0.2431771, SoilLocation == "HAY" ~ 0.6842187,
                            SoilLocation == "TLI" ~ 0.8786198, SoilLocation == "KNZ" ~ 0.3597917))
  return(clim.mod.df)
}

SMR.emmean.df <- Model_to_DF_w_climate(SMR.native.mod)
RSR.emmean.df <- Model_to_DF_w_climate(RSR.native.mod)
RMR.emmean.df <- Model_to_DF_w_climate(RMR.native.mod)

### Aridity Index
# New model with climate metric, non-transformed for slopes
SMR.AI.mod <- lm(response ~ Aridity_index*Drought.or.Watered*Genotype, data = SMR.emmean.df)
RSR.AI.mod <- lm(response ~ Aridity_index*Drought.or.Watered*Genotype, data = RSR.emmean.df)
RMR.AI.mod <- lm(response ~ Aridity_index*Drought.or.Watered*Genotype, data = RMR.emmean.df)
# Anovas
anova(SMR.AI.mod)
anova(RSR.AI.mod)
anova(RMR.AI.mod)
# Slope estimates
emtrends(SMR.AI.mod, ~ Drought.or.Watered, var = "Aridity_index")
emtrends(RSR.AI.mod, ~ Drought.or.Watered, var = "Aridity_index")
emtrends(RMR.AI.mod, ~ Drought.or.Watered, var = "Aridity_index")
# Plots
Figure6b <- ggplot(SMR.emmean.df, aes(x = Aridity_index, y = response, color = Drought.or.Watered)) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Aridity Index | Dry to Wet") +
  geom_vline(xintercept = 0.3399766, linetype="solid", color = "#780c72", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.3399766, y = 0.0130, label = "SVR", hjust = "middle", fontface = 2, color = "black") +
  geom_vline(xintercept = 0.4513070, linetype="solid", color = "#f5805d", size=12, alpha = 0.4) +
  annotate(geom = "text", x = 0.4513070, y = 0.0130, label = "HAY", hjust = "middle", fontface = 2, color = "black" ) +
  geom_vline(xintercept = 0.6187901, linetype="solid", color = "#ffe785", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.6187901, y = 0.0130, label = "TLI", hjust = "middle", fontface = 2, color = "black") +
  geom_vline(xintercept = 0.7792024, linetype="solid", color = "#3ba150", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.7792024, y = 0.0130, label = "KNZ", hjust = "middle", fontface = 2, color = "black" ) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, shape = Genotype), size = 0.75, linewidth = 1.75,
                  fill = "white", position = position_dodge(width = 0.06)) +
  scale_shape_manual(values = c(22,24)) +
  scale_color_manual(values = treatment_pallete, name = "Treatment") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  guides(color = guide_legend(title.position="top"), 
         shape = guide_legend(title.position="top"))

# New model with climate metric, non-transformed for slopes
SMR.PD.mod <- lm(response ~ PDSI*Drought.or.Watered*Genotype, data = SMR.emmean.df)
RSR.PD.mod <- lm(response ~ PDSI*Drought.or.Watered*Genotype, data = RSR.emmean.df)
RMR.PD.mod <- lm(response ~ PDSI*Drought.or.Watered*Genotype, data = RMR.emmean.df)
# Anovas
anova(SMR.PD.mod)
anova(RSR.PD.mod)
anova(RMR.PD.mod)
# Slope estimates
emtrends(SMR.PD.mod, ~ Drought.or.Watered, var = "PDSI")
emtrends(RSR.PD.mod, ~ Drought.or.Watered, var = "PDSI")
emtrends(RMR.PD.mod, ~ Drought.or.Watered, var = "PDSI")

Figure_6c <- ggplot(SMR.emmean.df, aes(x = PDSI, y = response, color = Drought.or.Watered)) +
  ylab("Shoot Mass Rate (g/day)") + xlab("PDSI | Dry to Wet") +
  # Annotations to denote sites (PDSI is not a linear measure)
  geom_vline(xintercept = 0.2431771, linetype="solid", color = "#780c72", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.2431771, y = 0.0130, label = "SVR", hjust = "middle", fontface = 2, color = "black") +
  geom_vline(xintercept = 0.3597917, linetype="solid", color = "#3ba150", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.3597917, y = 0.0130, label = "KNZ", hjust = "middle", fontface = 2, color = "black" ) +
  geom_vline(xintercept = 0.6842187, linetype="solid", color = "#f5805d", size=12, alpha = 0.4) +
  annotate(geom = "text", x = 0.6842187, y = 0.0130, label = "HAY", hjust = "middle", fontface = 2, color = "black" ) +
  geom_vline(xintercept = 0.8786198, linetype="solid", color = "#ffe785", size=12, alpha = 0.8) +
  annotate(geom = "text", x = 0.8786198, y = 0.0130, label = "TLI", hjust = "middle", fontface = 2, color = "black") +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL, shape = Genotype), size = 0.75, linewidth = 1.75,
                  fill = "white", position = position_dodge(width = 0.06)) +
  scale_shape_manual(values = c(22,24)) +
  scale_color_manual(values = treatment_pallete, name = "Treatment") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  guides(color = guide_legend(title.position="top"), 
         shape = guide_legend(title.position="top"))


library(patchwork)
Figure6_final <- Figure6a + (Figure6b / Figure_6c + plot_layout(guides = "collect") & theme(legend.box = "vertical", legend.spacing.y = unit(0, "cm"), legend.box.just = "left", legend.margin=margin(0,0,0,0))) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1, 0.4)) 

ggsave("figures/Figure6.svg", Figure6_final, width = 300, height = 175, units = "mm")
ggsave("figures/Figure6.png", Figure6_final, width = 300, height = 175, units = "mm")


#### Figure S17 Root/shoot mass rate - boxplots ####
sampledata_germinants$RootMassRate # third timepoint
sampledata_removed_NAs <- filter(sampledata_germinants, !(Timepoint == ""))
sampledata_removed_NAs$Timepoint <- factor(sampledata_removed_NAs$Timepoint, levels = c("early", "middle", "late"))

P1 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = RootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  ylab("Root Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late"))  +
  scale_fill_manual(values = treatment_pallete) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))

P2 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = ShootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  ylab("Shoot Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late")) +
  scale_fill_manual(values = treatment_pallete) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))

P3 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = RootShootRatio, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  ylab("Root/Shoot Ratio") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late")) +
  scale_fill_manual(values = treatment_pallete) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  

Plantmassrates_by_time <- ggarrange(P1, P2, P3, labels = "AUTO", legend = "none", nrow = 1)
#stats
mod.root  <- lmerTest::lmer(sqrt(RootMassRate) ~ Timepoint*Drought.or.Watered + (1|Block), data = sampledata_removed_NAs)
mod.shoot <- lmerTest::lmer(sqrt(ShootMassRate) ~ Timepoint*Drought.or.Watered + (1|Block), data = sampledata_removed_NAs)
mod.RS <- lmerTest::lmer(log(RootShootRatio) ~ Timepoint*Drought.or.Watered + (1|Block), data = sampledata_removed_NAs)
anova(mod.root)
anova(mod.shoot)
anova(mod.RS)
emmeans(mod.root,~Timepoint|Drought.or.Watered)
emmeans(mod.shoot,~Timepoint)
emmeans(mod.RS, ~Timepoint|Drought.or.Watered)
ggsave("./figures/S17_Plant_mass_rates_by_time.svg", Plantmassrates_by_time, height = 4, width = 8)
ggsave("./figures/S17_Plant_mass_rates_by_time.png", Plantmassrates_by_time, height = 4, width = 8)
