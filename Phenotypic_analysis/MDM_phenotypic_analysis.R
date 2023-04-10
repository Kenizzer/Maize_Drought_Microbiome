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
library(car); packageVersion('car')

# Theme set and Color Palettes
theme_set(theme_pubr())
genotype_pallete <- c("B73" = "#91ff26", "Mo17" = "#9426ff")# B73/Mo17 - Genotype
treatment_pallete <- c("W" = "#0000FF", "D" = "#DAA520") # Drought/WW     - Treatment
habitat_pallete <- c("Agriculture" = "#332288", "Native" = "#44AA99") # Native/Ag.     - Soil habitat
location_pallete <- c("SVR" = "#88CCEE", "HAY" = "#CC6677", "TLI" = "#DDCC77", "KNZ" = "#117733") # SVR/HAY/TLI/KNZ - Soil location
soil_merged_pallete<- c("SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") # Combination of Soil habitat and location, colors from https://lospec.com/palette-list/zencillo14


### Phenotypic analysis ###
### Load dataset for use throughout
sampledata <- read.csv("../Raw_data_n_metadata/Drought_Experiment_Data_FULL.csv", header = TRUE)
### Remove controls
sampledata_w_controls <- sampledata
sampledata <- sampledata[sampledata$SoilLocation != "Control",]
# Combine soil local and habitat into a merged soil factor
sampledata$SoilInoculum <- as.factor(paste(sampledata$SoilLocation, sampledata$SoilHabitat, sep = "_"))

#### Emergence Rate ####
# No effect of drought since the drought treatment wasn't applied to the plants prior to germination
sampledata$Germination <- as.numeric(as.factor(sampledata$Germination))
Germ_mod <- glmer(Germination ~ Genotype * SoilInoculum + (1|Block), data = sampledata)
car::Anova(Germ_mod, type = "III")
# Count # germinates by genotype (coded numerically, No = 1, yes = 2)
summary(as.factor(paste(sampledata$Genotype, sampledata$Germination, sep = "_")))
(302/(64 + 302)) * 100 #B73 82.5%
(225/(138 + 225)) * 100 #Mo17 62.0%

# Now remove non-germinants (coded numerically, No = 1, yes = 2)
sampledata_germinants <- sampledata[sampledata$Germination == 2,]
### Recalculating Root shoot ratio
# Dry weight Roots / Dry weight Shoots
a <- sampledata_germinants$FirstDryRootWeight / sampledata_germinants$FirstDryShootWeight
b <- sampledata_germinants$SecondDryRootWeight / sampledata_germinants$SecondDryShootWeight
c <- sampledata_germinants$ThirdDryRootWeight / sampledata_germinants$ThirdDryShootWeight
df <- data.frame(cbind(a,b,c))
sampledata_germinants$RootShootRatio <- rowSums(df, na.rm = TRUE)
sampledata_germinants$RootShootRatio[sampledata_germinants$RootShootRatio == 0] <- NA
rm(a,b,c,df)

#### Figure S12 Root/shoot mass rate - boxplots ####
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

# Do we see patterns with soil inocula for the root/shoot mass rate?
sampledata_removed_NAs$SoilInoculum <- factor(sampledata_removed_NAs$SoilInoculum, levels = c("SVR_Agriculture", "SVR_Native", "HAY_Native", "TLI_Agriculture", "TLI_Native", "KNZ_Native"))
                                                                                            
P3 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = RootMassRate, fill = SoilInoculum)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  ylab("Root Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late")) +
  scale_fill_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0)

P4 <- ggplot(sampledata_removed_NAs, aes(x = Timepoint, y = ShootMassRate, fill = SoilInoculum)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA, color = "black") +
  ylab("Shoot Mass Rate (g/day)") + facet_wrap(~Drought.or.Watered) +
  scale_x_discrete(labels=c("early" = "Early", "middle" = "Middle", "late" = "Late")) +
  scale_fill_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                                  expression(HAY[P]),expression(TLI[Ag]),
                                                  expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0)


Plantmassrates_by_time_soilinoculum <- ggarrange(P3, P4, common.legend = TRUE, labels = "AUTO", legend = 'right')

ggsave("./figures/Plant_mass_rates_by_time.svg", Plantmassrates_by_time, height = 4, width = 8)
ggsave("./figures/Plant_mass_rates_by_time.png", Plantmassrates_by_time, height = 4, width = 8)
ggsave("./figures/Plant_mass_rates_by_time_soilinoculum.svg", Plantmassrates_by_time_soilinoculum, height = 6, width = 12)
ggsave("./figures/Plant_mass_rates_by_time_soilinoculum.png", Plantmassrates_by_time_soilinoculum, height = 6, width = 12)

# Cleanup environment
rm(P1, P2, P3, P4, Plantmassrates_by_time, Plantmassrates_by_time_soilinoculum, 
   sampledata_removed_NAs, mod.root, mod.shoot)

#### Testing for an association of plant growth with longitude of the soil inocula ####
#SVR_Native	-100.9951
#SVR_Agriculture	-100.9828
#HAY_Native	-99.3033
#TLI_Agriculture	-97.5912
#TLI_Native	-97.4690
#KZ_Native	-96.6099

sampledata_with_long <- sampledata_germinants %>%
  group_by(SoilInoculum) %>%
  mutate( LONG = case_when(
    SoilInoculum == "SVR_Native" ~ -100.9951,
    SoilInoculum == "SVR_Agriculture" ~ -100.9828,
    SoilInoculum == "HAY_Native" ~ -99.3033,
    SoilInoculum == "TLI_Agriculture" ~ -97.5912,
    SoilInoculum == "TLI_Native" ~ -97.4690,
    SoilInoculum == "KNZ_Native" ~ -96.6099))

P1 <- ggplot(sampledata_with_long, aes(x = LONG, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank())

P2 <- ggplot(sampledata_with_long, aes(x = LONG, y = RootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Mass Rate (g/day)") + theme(axis.title.x = element_blank())

P3 <- ggplot(sampledata_with_long, aes(x = LONG, y = RootShootRatio, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root/Shoot Ratio") + xlab("Longitude (Lower MAP to Higher MAP)")


P4 <- ggarrange(P1,P2,P3, nrow = 3, align = 'hv', common.legend = TRUE, legend = 'right', labels = "AUTO")

ggsave("./figures/Long_by_plant_growth.svg", P4, height = 8, width = 4)
ggsave("./figures/Long_by_plant_growth.png", P4, height = 8, width = 4)

#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = data.frame(sampledata_with_long)))
anova(lmerTest::lmer(sqrt(RootMassRate) ~ LONG*Drought.or.Watered*LONG*Genotype + (1|Block), data = data.frame(sampledata_with_long)))
anova(lmerTest::lmer(log(RootShootRatio) ~ LONG*Drought.or.Watered*LONG*Genotype + (1|Block) + (1|Timepoint), data = data.frame(sampledata_with_long)))

# just native soils 
filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
ggplot(aes(x = LONG, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank())

filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = LONG, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank())

filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
ggplot(aes(x = LONG, y = RootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Mass Rate (g/day)") + theme(axis.title.x = element_blank())

filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = LONG, y = RootShootRatio, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Shoot Ratio") + xlab("Longitude (Lower MAP to Higher MAP)")
#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_long))
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
anova(lmerTest::lmer(sqrt(RootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
anova(lmerTest::lmer(log(RootShootRatio) ~ LONG*Drought.or.Watered*Genotype + (1|Block) + (1|Timepoint), data = filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
# Something interesting here with SMR across just native soils
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")))
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "LONG")
emtrends(SMR.native.mod, ~ Genotype, var = "LONG")

# Get non-transformed slopes for easier interpretation of results.
SMR.native.mod.non.transformed <- lmerTest::lmer(ShootMassRate ~ LONG*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")))
emtrends(SMR.native.mod.non.transformed, ~ Drought.or.Watered, var = "LONG")
emtrends(SMR.native.mod.non.transformed, ~ Genotype, var = "LONG")


# Both long interactions with genotype and treatment are significant
# PLOT them
P1 <- filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = LONG, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Treatment" , values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + theme(axis.title.x = element_blank()) +
  ylim(0.000, 0.022) + annotate("text", x = -98, y=0.022, label = "W = -3.29 e-4", fontface = 'italic') +
  annotate("text", x = -98, y=0.020, label = "D = 3.04 e-5", fontface = 'italic') +
  annotate("text", x = -98, y=0.018, label = "Long×DT p = 0.03", fontface = 'italic')

P2 <- filter(sampledata_with_long, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = LONG, y = ShootMassRate, color = Genotype)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(name = "Genotype" , values = genotype_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("Longitude (Lower MAP to Higher MAP)") +
  ylim(0.000, 0.022) + annotate("text", x = -98, y=0.022, label = "B73 = -3.07 e-4", fontface = 'italic') +
  annotate("text", x = -98, y=0.020, label = "Mo17 = 9.12 e-6", fontface = 'italic') +
  annotate("text", x = -98, y=0.018, label = "Long×G p = 0.03", fontface = 'italic') 



Native_only_smooth <- ggarrange(P1,P2, nrow = 2, ncol = 1, align = 'v', legend = 'right', labels = c('B', 'C'))

ggsave("./figures/NATIVE_only_Long_SMR_Genotype.svg", Native_only_smooth, height = 5.3, width = 4)
ggsave("./figures/NATIVE_only_Long_SMR_Genotype.png", Native_only_smooth, height = 5.3, width = 4)

# With Mean Annual Precipitation (county level) in-place of longitude #
# 30 year normal MAP 1991-2021
#SVR 19.44
#HAY 24.99
#KNZ 33.18
#TLI 30.21
sampledata_with_MAP <- sampledata_germinants %>%
  group_by(SoilInoculum) %>%
  mutate( MAP = case_when(
    SoilInoculum == "SVR_Native" ~ 19.44,
    SoilInoculum == "SVR_Agriculture" ~ 19.44,
    SoilInoculum == "HAY_Native" ~ 24.99,
    SoilInoculum == "TLI_Agriculture" ~ 30.21,
    SoilInoculum == "TLI_Native" ~ 30.21,
    SoilInoculum == "KNZ_Native" ~ 33.18))

filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("MAP 1991-2021")
ggplot(sampledata_with_MAP, aes(x = MAP, y = ShootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Shoot Mass Rate (g/day)") + xlab("MAP 1991-2021")

filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = MAP, y = RootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Mass Rate (g/day)") + xlab("MAP 1991-2021")
ggplot(sampledata_with_MAP, aes(x = MAP, y = RootMassRate, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root Mass Rate (g/day)") + xlab("MAP 1991-2021")

filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")) %>%
  ggplot(aes(x = MAP, y = RootShootRatio, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root/Shoot ratio") + xlab("MAP 1991-2021")
ggplot(sampledata_with_MAP, aes(x = MAP, y = RootShootRatio, color = Drought.or.Watered)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "lm") +
  scale_color_manual(values = treatment_pallete) +
  ylab("Root/Shoot ratio") + xlab("MAP 1991-2021")

#STATS
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = sampledata_with_MAP))
anova(lmerTest::lmer(sqrt(ShootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
anova(lmerTest::lmer(sqrt(RootMassRate)  ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
anova(lmerTest::lmer(log(RootShootRatio) ~ MAP*Drought.or.Watered*Genotype + (1|Block) + (1|Timepoint), data = filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native"))))
# W/o Ag. we see that MAP:treatment and MAP:Genotype
SMR.native.mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ MAP*Drought.or.Watered*Genotype + (1|Block), data = filter(sampledata_with_MAP, SoilInoculum %in% c("SVR_Native", "HAY_Native", "TLI_Native", "KNZ_Native")))
emtrends(SMR.native.mod, ~ Genotype, var = "MAP")
emtrends(SMR.native.mod, ~ Drought.or.Watered, var = "MAP")

#### Figure removed from manuscript: Timeseries of shoot height by treatment ####
sampledata$Drought.or.Watered <- as.factor(sampledata$Drought.or.Watered)

# mutate data to correct dim for plotting
sampledata.long <- dplyr::select(sampledata,Drought.or.Watered,Genotype,starts_with("Height")) %>% # height data only
  mutate(SampleID=factor(rownames(.)))
sampledata.long <- gather(sampledata.long,key='Height',value='Height_cm',starts_with('Height')) %>%
  mutate(Day=as.factor(str_remove(Height,'Height_')))

sampledata.long$Day<- recode_factor(sampledata.long$Day, Jan29  = "13",
                  Feb5   = "20", 
                  Feb12  = "27", 
                  Feb19  = "34")

# Stats for legend
mod.df <- lmer(Height_cm ~ Drought.or.Watered + Day + Drought.or.Watered*Day + (1|Genotype), data = sampledata.long)
anova(mod.df)
pairs(emmeans(mod.df, ~ Drought.or.Watered|Day))
# significant effect of treatment by day interaction.

# convert day to numeric for plotting
sampledata.long$Day <- as.numeric(as.character((sampledata.long$Day)))
# Make a dataframe to hold height measure averages by treatment/day
treatment_mean <- sampledata.long %>% group_by(Drought.or.Watered, Day) %>% 
                           summarise(Height_cm  = mean(Height_cm, na.rm=TRUE))

ht_plot <- ggplot(sampledata.long, aes(x=Day,y=Height_cm,color=Drought.or.Watered)) +
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

ggsave("./figures/Shoot_height_by_treatment_timeseries.svg", ht_plot, height = 6,width = 8)




#### Linear modeling on phenotypic data ####
AVR_mod <- lmer(AvgGrowthRate ~ Drought.or.Watered * Genotype * SoilInoculum + (1|Block), data = sampledata_germinants, na.action = na.exclude)
RMR_mod <- lmer(sqrt(RootMassRate) ~ Drought.or.Watered * Genotype * SoilInoculum + (1|Block), data = sampledata_germinants, na.action = na.exclude)
SMR_mod <- lmer(sqrt(ShootMassRate) ~ Drought.or.Watered * Genotype * SoilInoculum + (1|Block), data = sampledata_germinants, na.action = na.exclude)
TMR_mod <- lmer(sqrt(TotalMassRate) ~ Drought.or.Watered * Genotype * SoilInoculum + (1|Block), data = sampledata_germinants, na.action = na.exclude)
# control for timepoint of harvest here as we do not account for days post-germination in raw measure to get the ratio
RSR_mod <- lmer(log(RootShootRatio) ~ Drought.or.Watered * Genotype * SoilInoculum + (1|Block) + (1|Timepoint), data = sampledata_germinants, na.action = na.exclude)
# anovas
anova(AVR_mod)
anova(TMR_mod)
anova(RMR_mod)
anova(SMR_mod)
anova(RSR_mod)
# Posthocs and backtransformation XXX #
pairs(emmeans(RMR_mod, ~ Drought.or.Watered))
pairs(emmeans(RMR_mod, ~ Genotype))

pairs(emmeans(SMR_mod, ~ Drought.or.Watered))
pairs(emmeans(SMR_mod, ~ Genotype))
pairs(emmeans(SMR_mod, ~ SoilInoculum))

pairs(emmeans(RSR_mod, ~ Drought.or.Watered))
pairs(emmeans(RSR_mod, ~ Genotype))
pairs(emmeans(RSR_mod, ~ SoilInoculum:Genotype))

# Plot residuals vs fitted
plot(resid(AVR_mod),data.frame(sampledata_germinants)$AvgGrowthRate)
plot(resid(RMR_mod),data.frame(sampledata_germinants)$RootMassRate)
plot(resid(SMR_mod),data.frame(sampledata_germinants)$ShootMassRate)
plot(resid(TMR_mod),data.frame(sampledata_germinants)$TotalMassRate)
plot(resid(RSR_mod),data.frame(sampledata_germinants)$RootShootRatio)
# Homogeneity of Variance (~equal above/below line)
plot(AVR_mod)
plot(RMR_mod)
plot(SMR_mod)
plot(TMR_mod)
plot(RSR_mod)
# qq plots
lattice::qqmath(AVR_mod, id=0.05)
lattice::qqmath(RMR_mod, id=0.05)
lattice::qqmath(SMR_mod, id=0.05)
lattice::qqmath(TMR_mod, id=0.05)
lattice::qqmath(RSR_mod, id=0.05)
# correlation tests, mostly all correlated
cor.test(sampledata_germinants$AvgGrowthRate, sampledata_germinants$ShootMassRate, method = "pearson")
cor.test(sampledata_germinants$AvgGrowthRate, sampledata_germinants$RootMassRate, method = "pearson")
cor.test(sampledata_germinants$AvgGrowthRate, sampledata_germinants$RootShootRatio, method = "pearson")
cor.test(sampledata_germinants$ShootMassRate, sampledata_germinants$RootMassRate, method = "pearson")
cor.test(sampledata_germinants$ShootMassRate, sampledata_germinants$TotalMassRate, method = "pearson")
cor.test(sampledata_germinants$ShootMassRate, sampledata_germinants$RootShootRatio, method = "pearson")
cor.test(sampledata_germinants$TotalMassRate, sampledata_germinants$RootMassRate, method = "pearson")
cor.test(sampledata_germinants$TotalMassRate, sampledata_germinants$RootShootRatio, method = "pearson")
cor.test(sampledata_germinants$RootMassRate, sampledata_germinants$RootShootRatio, method = "pearson")

### Plotting of results
# This tibble is used to place the significance letters at a consist height over each barplot
labels_df <- tibble(Drought.or.Watered=levels(as.factor(sampledata_germinants$Drought.or.Watered)),
                    Genotype=levels(as.factor(sampledata_germinants$Genotype)),
                    #SoilLocation=levels(as.factor(sampledata_germinants$SoilLocation)),
                    RMR=max(sampledata_germinants$RootMassRate, na.rm = TRUE) * 1.1,
                    SMR=max(sampledata_germinants$ShootMassRate, na.rm = TRUE) * 1.1,
                    TMR=max(sampledata_germinants$TotalMassRate, na.rm = TRUE) * 1.1,
                    RSR=max(sampledata_germinants$RootShootRatio, na.rm = TRUE) * 1.1,
                    AGR=max(sampledata_germinants$AvgGrowthRate, na.rm = TRUE) * 1.1)




# Average growth rate
a <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = AvgGrowthRate, fill = Drought.or.Watered)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete, name = "Treatment") +
     ylab("Average Growth Rate(g/day)") + xlab("Treatment") #+ 
     #geom_text(data=labels_df, aes(Drought.or.Watered, AGR, label=c("a","b")), size = 6)

# Root Mass Rate
b <- ggplot(sampledata_germinants, aes(x = Genotype, y = RootMassRate, fill = Drought.or.Watered)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter( position = position_jitterdodge(0.2)) +
     scale_fill_manual(values = treatment_pallete, name = "Treatment") +
     ylab("Root Mass Rate (g/day)") + xlab("Genotype")

# Shoot mass rate
c <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = ShootMassRate, fill = Drought.or.Watered)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete, name = "Treatment") +
     ylab("Shoot Mass Rate (g/day)") + xlab("Treatment") #+ 
     #geom_text(data=labels_df, aes(Drought.or.Watered, SMR, label=c("a","b")), size = 6)
d <- ggplot(sampledata_germinants, aes(x = SoilLocation, y = ShootMassRate, fill = SoilLocation)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = location_pallete, name = "Soil Location") +
     ylab("Shoot Mass Rate (g/day)") + xlab("Soil Location")
e <- ggplot(sampledata_germinants, aes(x = Genotype, y = ShootMassRate, fill = Genotype)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = genotype_pallete, name = "Genotype") +
     ylab("Shoot Mass Rate (g/day)") + xlab("Genotype") #+
     #geom_text(data=labels_df, aes(Genotype, SMR, label=c("a","b")), size = 6)

# Total Mass Rate
f <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = TotalMassRate, fill = Drought.or.Watered)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete, name = "Treatment") +
     ylab("Total Mass Rate (g/day)") + xlab("Treatment")
g <- ggplot(sampledata_germinants, aes(x = Genotype, y = TotalMassRate, fill = Genotype)) +
     geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = genotype_pallete, name = "Genotype") +
     ylab("Total Mass Rate (g/day)") + xlab("Genotype")

# Log root shoot ratio
h <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = RootShootRatio, fill = Drought.or.Watered)) +
  geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Root/Shoot Ratio") + xlab("Treatment")

i <- ggplot(sampledata_germinants, aes(x = SoilHabitat, y = RootShootRatio, fill = Genotype)) +
  geom_boxplot(outlier.colour = NA) + geom_jitter(position = position_jitterdodge(0.2)) + scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Root/Shoot Ratio") + xlab("Genotype")

ggplot(sampledata_germinants, aes(x = Genotype, y = RootShootRatio, fill = Genotype)) +
  geom_boxplot(outlier.colour = NA) + geom_jitter(width = 0.2) + scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Root/Shoot Ratio") + xlab("Genotype")

# Drought effects
fig1 <- ggarrange(a,c,f,h,b, common.legend = TRUE, labels = "AUTO", legend = "right", align = "hv")
# Genotypic effects
fig2 <- ggarrange(e,g,i, common.legend = TRUE, labels = "AUTO", legend = "right", align = "hv", ncol = 3)
# Location effects
fig3 <- ggarrange(d, legend = "right")


ggsave("figures/drought_phenotypic_effects.svg",  fig1, width = 12, height = 8)
ggsave("figures/genotype_phenotypic_effects.svg", fig2, width = 12, height = 8)
ggsave("figures/location_phenotypic_effects.svg", fig3, width = 12, height = 8)
rm(a,c,f,h,b,e,g,i,d)
# combo figure for paper
# Drought effects
a <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = ShootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Shoot Mass Rate (g/day)") + xlab("Treatment") + 
  geom_text(data=labels_df, aes(Drought.or.Watered, SMR, label=c("a","b")), size = 6, color = '#3F3F3F')

b <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = RootMassRate, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Root Mass Rate (g/day)") + xlab("Treatment")+ 
  geom_text(data=labels_df, aes(Drought.or.Watered, RMR, label=c("a","b")), size = 6, color = '#3F3F3F')

c  <- ggplot(sampledata_germinants, aes(x = Drought.or.Watered, y = RootShootRatio, fill = Drought.or.Watered)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  scale_fill_manual(values = treatment_pallete, name = "Treatment") +
  ylab("Root/Shoot Ratio") + xlab("Genotype") + 
  geom_text(data=labels_df, aes(Drought.or.Watered, RSR, label=c("a","b")), size = 6, color = '#3F3F3F')

part1 <- ggarrange(a,b,c, align = "hv", common.legend = TRUE, legend = 'right',labels = c("A","B","C"), nrow=1)
# Genotype effects
d <- ggplot(sampledata_germinants, aes(x = Genotype, y = ShootMassRate, fill = Genotype)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Shoot Mass Rate (g/day)") + xlab("Genotype")  + 
  geom_text(data=labels_df, aes(Genotype, SMR, label=c("a","b")), size = 6, color = '#3F3F3F')
  
  
e <- ggplot(sampledata_germinants, aes(x = Genotype, y = RootMassRate, fill = Genotype)) +
   geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") +
  scale_fill_manual(values = genotype_pallete, name = "Genotype") +
  ylab("Root Mass Rate (g/day)") + xlab("Genotype")  + 
  geom_text(data=labels_df, aes(Genotype, RMR, label=c("a","b")), size = 6, color = '#3F3F3F')

f <- ggplot(sampledata_germinants, aes(x = Genotype, y = RootShootRatio, fill = Genotype)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  scale_fill_manual(values = genotype_pallete, name = "Treatment") +
  ylab("Root/Shoot Ratio") + xlab("Genotype") + 
  geom_text(data=labels_df, aes(Genotype, RSR, label=c("a","b")), size = 6, color = '#3F3F3F')

part2 <- ggarrange(d,e,f, align = "hv", common.legend = TRUE, legend = 'right',labels = c("D","E","F"), nrow=1)

Combo_fig <- ggarrange(part1, part2, nrow = 2)
ggsave("figures/phenotypic_figure_final.svg", Combo_fig, width = 12, height = 9)
ggsave("figures/phenotypic_figure_final.png", Combo_fig, width = 12, height = 9)

#### Figure 3 #####
part1
ggsave("figures/figure3_phenotypic_treatment.svg", part1, width = 10, height = 5)
ggsave("figures/figure3_phenotypic_treatment.png", part1, width = 10, height = 5)

#### Figure S6 #####
part2 <- ggarrange(d,e,f, align = "hv", common.legend = TRUE, legend = 'right',labels = "AUTO", nrow=1)
part2
ggsave("figures/figureS6_phenotypic_genotype.svg", part2, width = 10, height = 5)
ggsave("figures/figureS6_phenotypic_genotype.png", part2, width = 10, height = 5)

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
pairs(emmeans(SMR_mod, ~ SoilInoculum))
pairs(emmeans(RSR_mod, ~ SoilInoculum))

# Full posthocs
pairs(emmeans(SMR_mod, ~ Genotype | SoilInoculum ))
pairs(emmeans(RSR_mod, ~ Genotype | SoilInoculum ))


#### Plots by soil inocula ####
# Reorder dry to wet sites
levels(sampledata_germinants$SoilInoculum)
sampledata_germinants$SoilInoculum <- factor(sampledata_germinants$SoilInoculum, levels = c("SVR_Agriculture", "SVR_Native", "HAY_Native",
                                              "TLI_Agriculture", "TLI_Native", "KNZ_Native"))

# Droughted plants

a <- filter(sampledata_germinants, Drought.or.Watered == "D") %>%
  ggplot(aes(x = SoilInoculum, y = ShootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum",
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                    expression(HAY[P]),expression(TLI[Ag]),
                    expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

b <- filter(sampledata_germinants, Drought.or.Watered == "D") %>%
ggplot(aes(x = SoilInoculum, y = RootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                    expression(HAY[P]),expression(TLI[Ag]),
                    expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())


c <- filter(sampledata_germinants, Drought.or.Watered == "D") %>%
ggplot(aes(x = SoilInoculum, y = RootShootRatio, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root/Shoot Ratio") +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                    expression(HAY[P]),expression(TLI[Ag]),
                    expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0,  axis.title.x=element_blank())

Drought_pheno_soilinoculum <- ggarrange(a,b,c, labels = "AUTO", nrow = 1, common.legend = TRUE, legend = 'right')

ggsave("figures/Drought_phenotype_response_to_soilinoculum.svg", Drought_pheno_soilinoculum, width = 12, height = 6)
ggsave("figures/Drought_phenotype_response_to_soilinoculum.png", Drought_pheno_soilinoculum, width = 12, height = 6)

# Stats
anova(lmer(sqrt(ShootMassRate) ~ SoilInoculum + Genotype + (1|Block), data = sampledata_germinants[sampledata_germinants$Drought.or.Watered == "D",]))
anova(lmer(sqrt(RootMassRate) ~ SoilInoculum +  Genotype + (1|Block), data = sampledata_germinants[sampledata_germinants$Drought.or.Watered == "D",]))
anova(lmer(log(RootShootRatio) ~ SoilInoculum + Genotype + (1|Block), data = sampledata_germinants[sampledata_germinants$Drought.or.Watered == "D",]))
# All non-significant with regards to soil inoculum source
rm(a,b,c)

# WW
a<- filter(sampledata_germinants, Drought.or.Watered == "W") %>%
  ggplot(aes(x = SoilInoculum, y = ShootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)")  +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0,  axis.title.x=element_blank())

b<- filter(sampledata_germinants, Drought.or.Watered == "W") %>%
  ggplot(aes(x = SoilInoculum, y = RootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)")  +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0,  axis.title.x=element_blank())



c<- filter(sampledata_germinants, Drought.or.Watered == "W") %>%
  ggplot( aes(x = SoilInoculum, y = RootShootRatio, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root/Shoot Ratio")  +
  scale_x_discrete("Soil Inoculum", labels = c(expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete, name = "Soil Inoculum", 
                    labels = c(expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) + 
  theme(legend.text.align = 0,  axis.title.x=element_blank())

WellWatered_pheno_soilinoculum <- ggarrange(a,b,c, labels = "AUTO", nrow = 1, common.legend = TRUE, legend = 'right')

ggsave("figures/WellWatered_phenotype_response_to_soilinoculum.svg", WellWatered_pheno_soilinoculum, width = 12, height = 6)
ggsave("figures/WellWatered_phenotype_response_to_soilinoculum.png", WellWatered_pheno_soilinoculum, width = 12, height = 6)


#### Plots by soil inocula w/ controls ####
sampledata_w_controls
# Combine soil local and habitat into a merged soil factor
sampledata_w_controls$SoilInoculum <- as.factor(paste(sampledata_w_controls$SoilLocation, sampledata_w_controls$SoilHabitat, sep = "_"))
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


# Reorder dry to wet sites
levels(sampledata_w_controls_germinants$SoilInoculum)
sampledata_w_controls_germinants$SoilInoculum <- factor(sampledata_w_controls_germinants$SoilInoculum, levels = c("Control_Control","SVR_Agriculture", "SVR_Native", "HAY_Native",
                                                                                            "TLI_Agriculture", "TLI_Native", "KNZ_Native"))

soil_merged_pallete_w_control <- c("Control_Control" = "#FFFFFF","SVR_Agriculture" = "#780c72","SVR_Native" = "#de3a68", "HAY_Native" = "#f5805d",
                        "TLI_Agriculture" = "#ffe785", "TLI_Native" = "#7fd66f", "KNZ_Native" = "#3ba150") 


#### Figure 6 ####
# Both plots together FACET 
# ShootMassRate
a <- ggplot(sampledata_w_controls_germinants, aes(x = SoilInoculum, y = ShootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Drought.or.Watered)


# RootMassRate
b <- ggplot(sampledata_w_controls_germinants, aes(x = SoilInoculum, y = RootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~Drought.or.Watered)

# RootShootRatio
c <- ggplot(sampledata_w_controls_germinants, aes(x = SoilInoculum, y = RootShootRatio, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root/Shoot Ratio") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank()) +
  facet_wrap(~Drought.or.Watered) + ylim(0,12) # This removes 4 points, noted in legend

# add line plots
Controls_phenotypic_responses_with_line_plot <- ggarrange(a, Native_only_smooth, ncol = 2, widths = c(1, 0.4), legend = 'right', labels = "AUTO")
# save
ggsave("figures/Figure6_Controls_phenotypic_responses_facet_lineplots.svg", Controls_phenotypic_responses_with_line_plot, height = 8, width = 14)
ggsave("figures/Figure6_Controls_phenotypic_responses_facet_lineplots.png", Controls_phenotypic_responses_with_line_plot,  height = 8, width = 14)

#### Figure S15 ####
b
c
RMR_RSR_plot <- ggarrange(b,c, nrow = 2, common.legend = TRUE, labels = "AUTO", align = 'hv', legend = "right")
# save
ggsave("figures/FigureS15_Controls_phenotypic_responses_RMR_RSR.svg", RMR_RSR_plot, height = 8, width = 10)
ggsave("figures/FigureS15_Controls_phenotypic_responses_RMR_RSR.png", RMR_RSR_plot,  height = 8, width = 10)

#Seperate plots for D vs WW

a <- filter(sampledata_w_controls_germinants, Drought.or.Watered == "D") %>%
  ggplot(aes(x = SoilInoculum, y = ShootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())


b <- filter(sampledata_w_controls_germinants, Drought.or.Watered == "D") %>%
  ggplot(aes(x = SoilInoculum, y = RootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())


c <- filter(sampledata_w_controls_germinants, Drought.or.Watered == "D") %>%
  ggplot(aes(x = SoilInoculum, y = RootShootRatio, fill = SoilInoculum)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root/Shoot Ratio") +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

Drought_pheno_soilinoculum_w_controls <- ggarrange(a,b,c, labels = "AUTO", nrow = 1, common.legend = TRUE, legend = 'right')

ggsave("figures/Drought_phenotype_response_to_soilinoculum_w_controls.svg", Drought_pheno_soilinoculum_w_controls, width = 16, height = 6)
ggsave("figures/Drought_phenotype_response_to_soilinoculum_w_controls.png", Drought_pheno_soilinoculum_w_controls, width = 16, height = 6)

rm(a,b,c)

# WW
a<- filter(sampledata_w_controls_germinants, Drought.or.Watered == "W") %>%
  ggplot(aes(x = SoilInoculum, y = ShootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Shoot Mass Rate (g/day)")  +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

b<- filter(sampledata_w_controls_germinants, Drought.or.Watered == "W") %>%
  ggplot(aes(x = SoilInoculum, y = RootMassRate, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Mass Rate (g/day)")  +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())



c<- filter(sampledata_w_controls_germinants, Drought.or.Watered == "W") %>%
  ggplot( aes(x = SoilInoculum, y = RootShootRatio, fill = SoilInoculum)) +
  geom_jitter(width = 0.3) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Root Shoot Ratio")  +
  scale_x_discrete("Soil Inoculum", labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                                               expression(HAY[P]),expression(TLI[Ag]),
                                               expression(TLI[P]),expression(KNZ[P]))) +
  scale_fill_manual(values = soil_merged_pallete_w_control, name = "Soil Inoculum", 
                    labels = c("Control", expression(SVR[Ag]),expression(SVR[P]),
                               expression(HAY[P]),expression(TLI[Ag]),
                               expression(TLI[P]),expression(KNZ[P]))) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

WellWatered_pheno_soilinoculum_w_controls <- ggarrange(a,b,c, labels = "AUTO", nrow = 1, common.legend = TRUE, legend = 'right')

ggsave("figures/WellWatered_phenotype_response_to_soilinoculum_w_controls.svg", WellWatered_pheno_soilinoculum_w_controls, width = 16, height = 6)
ggsave("figures/WellWatered_phenotype_response_to_soilinoculum_w_controls.png", WellWatered_pheno_soilinoculum_w_controls, width = 16, height = 6)