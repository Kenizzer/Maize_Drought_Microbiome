# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(ggpubr); packageVersion('ggpubr')

# Theme set 
theme_set(theme_pubr())
location_pallete <- c("SVR" = "#780c72", "HAY" = "#f5805d", "TLI" = "#ffe785", "KNZ" = "#3ba150") # SVR/HAY/TLI/KNZ - Soil location



# Take terraclimate point data and get mean annual precipitation, Evapotranspiration and PDSI
# then calculate an aridity index.
Get_aridity <- function(terraclimdf){
  temp  <- terraclimdf %>% group_by(Year) %>% summarize(Annual_precip_mm = sum(ppt.mm.), Annual_ETo_mm = sum(pet.mm.), Annual_PDSI = mean(PDSI.unitless.), Annual_soilmoisture_mm = sum(soil.mm.))
  # get the last normal AKA 1990 to 2021 
  temp2 <- temp %>% filter(between(Year, 1990, 2021)) %>% summarize(MA_Prec = mean(Annual_precip_mm), MA_ETo = mean(Annual_ETo_mm), M_monthly_PDSI = mean(Annual_PDSI), MA_Soil_moisture = mean(Annual_soilmoisture_mm))
  temp2$Aridity_index <- temp2$MA_Prec / temp2$MA_ETo
  return(temp2)
}

#Aridity Index Value chart https://doi.org/10.1038/s41597-022-01493-1
# Climate Class 
# <0.03 Hyper Arid
# 0.03–0.2 Arid
# 0.2–0.5 Semi-Arid
# 0.5–0.65 Dry sub-humid
# >0.65 Humid

##### Coords #####
#Smoky Valley Ranch prairie [latitude 38.8665, longitude -100.9951]
#Hayes prairie [38.8355, -99.3033]
#The Land Institute prairie [38.9698, -97.4690]
#Konza native prairie [39.1056, -96.6099])

# Load point data
SVR_Native <- read.csv("terraclimate_38.8665N_100.9951W.csv", skip = 14, header = TRUE)
HAY_Native <- read.csv("terraclimate_38.8355N_99.3033W.csv", skip = 14, header = TRUE)
TLI_Native <- read.csv("terraclimate_38.9698N_97.4690W.csv", skip = 14, header = TRUE)
KNZ_Native <- read.csv("terraclimate_39.1056N_96.6099W.csv", skip = 14, header = TRUE)


# Get aridity indexes
arid.df <- data.frame(rbind(Get_aridity(SVR_Native),
Get_aridity(HAY_Native),
Get_aridity(TLI_Native),
Get_aridity(KNZ_Native)))
# rename rows and add site factor for plotting
rownames(arid.df) <- c("SVR", "HAY", "TLI", "KNZ")
arid.df$site <- as.factor(row.names(arid.df))
arid.df$site <- factor(arid.df$site, levels = c("SVR", "HAY","TLI", "KNZ"))
# These values will be utilized in other scripts for linear modeling
arid.df

# MA_Prec   MA_ETo M_monthly_PDSI MA_Soil_moisture Aridity_index site
# SVR 478.9250 1408.700      0.2431771         37.38437     0.3399766  SVR
# HAY 601.7375 1333.322      0.6842187         92.87812     0.4513070  HAY
# TLI 783.8156 1266.691      0.8786198        177.09688     0.6187901  TLI
# KNZ 912.3438 1170.869      0.3597917        348.25937     0.7792024  KNZ



# Barplots for the normal measures (1991-2020) by collection site
a <- ggplot(arid.df, aes(x = site, y = MA_Prec, fill = site)) +
  geom_col(color = "black") +
  ylab("Normal Annual Precipitation 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')
  
b <- ggplot(arid.df, aes(x = site, y = MA_ETo, fill = site)) +
  geom_col(color = "black") +
  ylab("Normal Annual Evapotranspiration 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')

c <- ggplot(arid.df, aes(x = site, y = M_monthly_PDSI, fill = site)) +
  geom_col(color = "black") +
  ylab("Normal Monthly PDSI 1990-2021") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')


d <- ggplot(arid.df, aes(x = site, y = MA_Soil_moisture, fill = site)) +
  geom_col(color = "black") +
  ylab("Normal Annual Soil Moisture 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')

e <- ggplot(arid.df, aes(x = site, y = Aridity_index, fill = site)) +
  geom_col(color = "black") +
  ylab("Aridity Index") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank(), legend.position = 'right')

f <- ggarrange(a,b,c,d,e, common.legend = TRUE, legend = 'right', align = 'hv')

ggsave("Normals_by_site.svg", f, height = 9, width = 12)
ggsave("Normals_by_site.png", f, height = 9, width = 12)
rm(a,b,c,d,e,f) # cleanup

# Plots with the values over 30 years represented
# first reformat data to closer to what I need
Get_data_30_years <- function(terraclimdf){
  temp  <- terraclimdf %>% group_by(Year) %>% summarize(Annual_precip_mm = sum(ppt.mm.), Annual_ETo_mm = sum(pet.mm.), Annual_PDSI = mean(PDSI.unitless.), Annual_soilmoisture_mm = sum(soil.mm.))
  # get the last normal AKA 1990 to 2021 
  temp2 <- temp %>% filter(between(Year, 1990, 2021)) 
  temp2$Aridity_index <- temp2$Annual_precip_mm / temp2$Annual_ETo_mm
  return(temp2)
}

long.df <- data.frame(rbind(Get_data_30_years(SVR_Native),
                            Get_data_30_years(HAY_Native),
                            Get_data_30_years(TLI_Native),
                            Get_data_30_years(KNZ_Native)))



# rename rows and add site factor for plotting
long.df$site <- as.factor(rep(c("SVR", "HAY","TLI", "KNZ"), each = 32))
long.df$site <- factor(long.df$site, levels = c("SVR", "HAY", "TLI", "KNZ"))
colnames(long.df)

a <- ggplot(long.df, aes(x = site, y = Annual_precip_mm, fill = site)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Annual Precipitation 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

b <- ggplot(long.df, aes(x = site, y = Annual_ETo_mm, fill = site)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Annual Evapotranspiration 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

c <- ggplot(long.df, aes(x = site, y = Aridity_index, fill = site)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Aridity Index 1990-2021") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

d <- ggplot(long.df, aes(x = site, y = Annual_PDSI, fill = site)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Average Monthly PDSI 1990-2021") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank())

e <- ggplot(long.df, aes(x = site, y = Annual_soilmoisture_mm, fill = site)) +
  geom_jitter(width = 0.2) + geom_boxplot(outlier.colour = NA, alpha=0.8, color = "black") + 
  ylab("Annual Soil Moisture 1990-2021 (mm)") +
  scale_fill_manual(name = "Site", values = location_pallete) +
  theme(legend.text.align = 0, axis.title.x=element_blank())


f <- ggarrange(a,b,c,d,e, common.legend = TRUE, legend = "right", align = 'hv', labels = "AUTO")

ggsave("Annual_boxplots_by_site.svg", f, height = 9, width = 12)
ggsave("Annual_boxplots_by_site.png", f, height = 9, width = 12)

rm(a,b,c,d,e,f) # cleanup