# Maize Drought Microbiome

## Abridged abstract

Coming soon!

Table of markers used within the study.

| Marker   | F Primer| F Primer Sequence       | R Primer| R Primer Sequence    | Size (bp) | Citation                                       |
|:--------:|:-------:|:-----------------------:|:-------:|:--------------------:|:---------:|:----------------------------------------------:|
| 16S      | 515F    | GTGYCAGCMGCCGCGGTAA     | 806R    | GGACTACNVGGGTWTCTAAT | 390       | Parada *et al.* 2016 and Apprill *et al.* 2015 |
| ITS      | ITS1f   | CTTGGTCATTTAGAGGAAGTAA  | ITS2    | GCTGCGTTCTTCATCGATGC | variable       | Smith and Peay 2014                       |


## Experimental Design

**A)** Soil samples were collected from either agricultural grain fields or native prairies at four locations in Kansas along an increasing precipitation gradient. **B)** A complete randomized block design (blocks = 27) was used. Each pot was filled with sterile calcined clay soil mix and received a seed of either Zea mays genotype B73 or Mo17, was subject to either well-watered or drought conditions, and was inoculated with one of the six inocula. One uninoculated control was used per randomized block (denoted with white stripe). **C)** Plant biomass was collected in a time series: for time point one (day 25) and two (day 39), 100 plants each were destructively harvested, then dried for root and shoot biomass. For the third time point (day 50), all remaining plants (n = 529 – non-germinates) were harvested. Prior to drying for biomass measures, one nodal root from each plant was collected for amplicon sequencing for bacterial and fungal microbiome composition. **D)** Plant height was measured every seven days in weeks 2-5 of the experiment.

![Image of experimental design](https://github.com/Kenizzer/Maize_Drought_Microbiome/blob/main/Experimental_design_600dpi.png)


## Phenotypic processing and analysis

Coming Soon!

## Soil processing and analysis

Coming Soon!

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