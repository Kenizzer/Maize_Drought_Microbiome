# Maize Drought Microbiome

## bioRxiv Preprint
<https://doi.org/10.1101/2023.04.11.536405>

## Abstract

**Background and Aims:** As the climate changes, plants and their associated microbiomes face greater water limitation and increased frequency of drought. Soil- and root-associated microbes perform important functions that affect plant drought resilience, but the dynamics of their interactions with the host are poorly understood. Plant growth responses to natural soil microbiomes, as well as patterns of root colonization by microbes, are challenging to predict.
**Methods:** We collected soil microbiomes from four native prairies across a steep precipitation gradient in Kansas, USA. Seedlings of two Zea mays genotypes were inoculated with each soil microbiome in a factorial drought experiment. We investigated plant phenotypic and root microbiome responses to drought and compared the effects of natural microbiome variation on plant growth under water-limited and well-watered conditions. 
**Results:** Drought caused plants to accumulate shoot mass more slowly and achieve greater root/shoot mass ratios. Drought restructured the bacterial root-associated microbiome via depletion of Pseudomonadota and enrichment of Actinomycetota, whereas the fungal microbiome was largely unaffected. Taxonomically distinct soil microbiomes from the four contrasting environments affected plant growth under well-watered but not drought conditions. 
**Conclusion:** We demonstrated that the functional consequences of naturally-occurring soil microbiome variation are dependent on water availability, suggesting that future drying climates may dampen plants’ responsiveness to beneficial and/or pathogenic microbes.


Table of markers used within the study.

| Marker   | F Primer| F Primer Sequence       | R Primer| R Primer Sequence    | Size (bp) | Citation                                       |
|:--------:|:-------:|:-----------------------:|:-------:|:--------------------:|:---------:|:----------------------------------------------:|
| 16S      | 515F    | GTGYCAGCMGCCGCGGTAA     | 806R    | GGACTACNVGGGTWTCTAAT | 390       | Parada *et al.* 2016 and Apprill *et al.* 2015 |
| ITS      | ITS1f   | CTTGGTCATTTAGAGGAAGTAA  | ITS2    | GCTGCGTTCTTCATCGATGC | variable  | Smith and Peay 2014                            |


## Experimental Design

**A)** Soil samples were collected from native prairies at four locations in Kansas along precipitation gradient. **B-E** Boxplots to the right of the map depict annual climatic metrics extracted from the TerraClimate database for each collection site from 1990-2021: B) precipitation, C) reference evapotranspiration, D) aridity index, and E) soil moisture. Boxplot hinges represent the 1st and 3rd quartiles; whiskers represent 1.5 times the interquartile range.

![Image of experimental design 1](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Experimental_design_1.png)

 **A)** A complete randomized block design was used; depicted is a single representative block (number of blocks=27). Each pot was filled with sterile calcined clay soil mix and received a seed of either Zea mays genotype B73 or Mo17, was subject to either well-watered or drought conditions, and was inoculated with one of the four inocula. One uninoculated control (denoted with white stripe) and three additional randomized replicates were included per randomized block. **B)** Plant biomass was collected in a time series: for time point one (day 25) and two (day 39), 68 plants each were destructively harvested, then dried for root and shoot biomass. For the third time point (day 50), all remaining plants (n=379 – non-germinates) were harvested. Prior to drying for biomass measures, one nodal root from each plant was collected for amplicon sequencing for bacterial and fungal microbiome composition.

![Image of experimental design 2](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Experimental_design_2.png)

## Phenotypic processing and analysis

[***See Phenotypic_analysis folder***](https://github.com/Kenizzer/Maize_Drought_Microbiome/tree/main/Phenotypic_analysis)

- Main analysis was conducted in [***MDM_phenotypic_analysis.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Phenotypic_analysis/MDM_phenotypic_analysis.R).
- Historical enviromental paramemters were processing in [***TerraClimate_pointdata_analysis.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Phenotypic_analysis/TerraClimate_environmental_parameters/TerraClimate_pointdata_analysis.R).

## Soil and inocula processing and analysis

[***See Soil_inocula_prep_n_analysis folder***](https://github.com/Kenizzer/Maize_Drought_Microbiome/tree/main/Soil_inocula_prep_n_analysis)

- Bacterial soil and inocula sequences were processed, filtered, and CLR transformed in [***MDM_soil_inocula_bacterial_dataprep.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Soil_inocula_prep_n_analysis/MDM_soil_inocula_bacterial_dataprep.R).
- Fungal soil and inocula sequences were processed, filtered, and CLR transformed in [***MDM_soil_inocula_fungal_dataprep.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Soil_inocula_prep_n_analysis/MDM_soil_inocula_fungal_dataprep.R).
- Main soil and inocula analysis was conducted in [***MDM_soil_inocula_analysis.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Soil_inocula_prep_n_analysis/MDM_soil_inocula_analysis.R).

## Bacterial processing and analysis

[***See Bacterial_prep_n_analysis folder***](https://github.com/Kenizzer/Maize_Drought_Microbiome/tree/main/Bacterial_prep_n_analysis)

- Bacterial sequences were processed, filtered, and CLR transformed in [***MDM_bacterial_dataprep.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Bacterial_prep_n_analysis/MDM_bacterial_dataprep.R).
- Bacterial main analysis was conducted in [***MDM_bacterial_main_analysis.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Bacterial_prep_n_analysis/MDM_bacterial_main_analysis.R).
- Bacterial machine learning was conducted in [***MDM_bacterial_machine_learning.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Bacterial_prep_n_analysis/MDM_bacterial_machine_learning.R).
- Bacterial associations between plant phenotype and taxon abundance was conducted in [***MDM_bacterial_phenotype_to_microbiome.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Bacterial_prep_n_analysis/MDM_bacterial_phenotype_to_microbiome.R).

## Fungal processing and analysis

[***See Fungal_prep_n_analysis folder***](https://github.com/Kenizzer/Maize_Drought_Microbiome/tree/main/Fungal_prep_n_analysis)

- Fungal sequences were processed, filtered, and CLR transformed in [***MDM_fungal_dataprep.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Fungal_prep_n_analysis/MDM_fungal_dataprep.R).
- Fungal main analysis was conducted in [***MDM_fungal_main_analysis.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Fungal_prep_n_analysis/MDM_fungal_main_analysis.R).
- Fungal machine learning was conducted in [***MDM_fungal_machine_learning.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Fungal_prep_n_analysis/MDM_fungal_machine_learning.R).
- Fungal associations between plant phenotype and taxon abundance was conducted in [***MDM_fungal_phenotype_to_microbiome.R***](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Fungal_prep_n_analysis/MDM_fungal_phenotype_to_microbiome.R).