# Maize Drought Microbiome

## bioRxiv Preprint
<https://doi.org/10.1101/2023.04.11.536405>

## Abstract

As the climate changes, plants and their associated microbiomes face greater water limitation and increased frequency of drought. Historical precipitation patterns can leave a legacy effect on soil and root-associated microbiomes, but the impact of this conditioning on future drought performance is poorly understood. Moreover, agricultural practices, like tilling and fertilization, may override precipitation legacies in natural systems. Precipitation gradients provide a means to assess these legacy effects. We collected six soil microbiomes across a Kansas precipitation gradient from both agricultural fields and native prairies. Seedlings of two *Zea mays* genotypes were inoculated with each soil microbiome in a factorial drought experiment. Droughted plants exhibited decreased shoot mass accumulation rates and greater root mass relative to shoot mass. Restructuring of the bacterial root-associated microbiome was apparent, with depletion observed in Pseudomonadota and enrichment in Actinomycetota, while the fungal microbiome was largely unaffected by drought. A historical precipitation legacy effect on soil microbiomes interacted with plants during drought treatment, but only among prairie soils. Prairie soils from historically wetter locations increased maize shoot biomass under drought more so than agricultural or historically drier prairie soils. We demonstrate links between legacy effects and drought performance, suggesting that future drying climates may condition soils to negatively impact plant performance.

Table of markers used within the study.

| Marker   | F Primer| F Primer Sequence       | R Primer| R Primer Sequence    | Size (bp) | Citation                                       |
|:--------:|:-------:|:-----------------------:|:-------:|:--------------------:|:---------:|:----------------------------------------------:|
| 16S      | 515F    | GTGYCAGCMGCCGCGGTAA     | 806R    | GGACTACNVGGGTWTCTAAT | 390       | Parada *et al.* 2016 and Apprill *et al.* 2015 |
| ITS      | ITS1f   | CTTGGTCATTTAGAGGAAGTAA  | ITS2    | GCTGCGTTCTTCATCGATGC | variable  | Smith and Peay 2014                            |


## Experimental Design

**A)** Soil samples were collected from either agricultural grain fields or native prairies at four locations in Kansas along an increasing precipitation gradient. **B)** A complete randomized block design (blocks = 27) was used. Each pot was filled with sterile calcined clay soil mix and received a seed of either Zea mays genotype B73 or Mo17, was subject to either well-watered or drought conditions, and was inoculated with one of the six inocula. One uninoculated control was used per randomized block (denoted with white stripe). **C)** Plant biomass was collected in a time series: for time point one (day 25) and two (day 39), 100 plants each were destructively harvested, then dried for root and shoot biomass. For the third time point (day 50), all remaining plants (n = 529 – non-germinates) were harvested. Prior to drying for biomass measures, one nodal root from each plant was collected for amplicon sequencing for bacterial and fungal microbiome composition. **D)** Plant height was measured every seven days in weeks 2-5 of the experiment.

![Image of experimental design](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Experimental_design_600dpi.png)


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