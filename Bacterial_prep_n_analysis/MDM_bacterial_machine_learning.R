# Analysis of maize B73 and Mo17 inoculated plants under drought stress
# Samples were collected from a greenhouse study in Feb 2020
# Microbiomes are from the root compartment (rhizosphere/endosphere)
# Growth measurements were also recorded throughout the experiment (~50 days in length)
# Code by: Joel Swift

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
library('tidyverse'); packageVersion('tidyverse')
library("phyloseq"); packageVersion('phyloseq')
library("ggpubr"); packageVersion('ggpubr')
library("MASS"); packageVersion('MASS')
library("scales"); packageVersion('scales')
library("caret"); packageVersion('caret')
library("AppliedPredictiveModeling"); packageVersion('AppliedPredictiveModeling')
library("ranger"); packageVersion('ranger')
library("e1071"); packageVersion('e1071')
library("randomForest"); packageVersion('randomForest')
library("alluvial"); packageVersion('alluvial')
library("matrixStats"); packageVersion("matrixStats")

# Theme set and Color Palettes
theme_set(theme_pubr())
phyla_palette <- c("Actinobacteria" = "#D55E00", 
                   "Bacteroidetes" = "#999933",
                   "Gemmatimonadetes" = "#CC6677",
                   "Nitrospirae" = "#51C14B", 
                   "Proteobacteria" = "#15632c") # reduced so all do not print on important ASV charts

# Functions
# Function to return metadata df from phyloseq object
pssd2veg <- function(physeq) {
  # From a phyloseq object return a dataframe of the sample metadata for use in vegan
  # From: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Function to plot confusion matrix using ggtile plot from a confusion matrix object
# By user: Enrique Perez Herrero 
# on https://stackoverflow.com/questions/46063234/how-to-produce-a-confusion-matrix-and-find-the-misclassification-rate-of-the-na%C3%AF
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]))
                   #,"Kappa", percent_format()(m$overall[2])) no kappa on plot
  
  d <- as.data.frame.matrix(m$table)
  drn <- colnames(d)
  drr <- rownames(d)
  d <- as.data.frame(t(t(d)/colSums(d))) # Standardize the shading in the ggplot by # samples for a given factor level in the reference. 
  d <- d %>% gather(x, value)
  Y <- cbind(as.data.frame(m$table), Proportion = d$value)
  Y$Reference <- fct_rev(Y$Reference) # Added this line to get a downward diagonal 
  p <-
    ggplot(data = Y, aes(x = Reference, y = Prediction, fill= Proportion)) +
    geom_tile( colour = "white") +
    scale_fill_gradient(low = "white", high = "#14A02E", na.value = "white", limits=c(0,1)) +
    ggtitle(mytitle) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 60, hjust = 1)) +
    guides(fill = guide_colorbar(frame.colour = "black", ticks = FALSE))
  return(p)
}

# Define a custom Ranger RF model that saves the per fold info to a seperate folder.
ranger_funcs <- getModelInfo("ranger", regex = FALSE)[[1]]
ranger_funcs$fit <- function (x, y, wts, param, lev, last, classProbs, ...) 
{
  if ((!is.data.frame(x)) || dplyr::is.tbl(x)) 
    x <- as.data.frame(x, stringsAsFactors = TRUE)
  x$.outcome <- y
  if (!is.null(wts)) {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = min(param$mtry, ncol(x)), min.node.size = param$min.node.size, 
                          splitrule = as.character(param$splitrule), write.forest = TRUE, 
                          probability = classProbs, case.weights = wts, ...)
  }
  else {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = min(param$mtry, ncol(x)), min.node.size = param$min.node.size, 
                          splitrule = as.character(param$splitrule), write.forest = TRUE, 
                          probability = classProbs, ...)
  }
  if (!last) 
    out$y <- y
  save(out, file = paste("./rangerML_stats/ranger", param$mtry, param$splitrule, param$min.node.size, format(Sys.time(), "%H_%M_%S.RData"), sep = "_"))
  out
}

# Machine learning main function
MachineLearning_RF_ranger <- function(PHYSEQ_OBJ_1, GROUPING, TREES) {
  # Remove ASV Table and meta data from phyloseq objects
  ASV.df <- as.data.frame(otu_table(PHYSEQ_OBJ_1))
  ASV_metadata.df <- as.data.frame(sample_data(PHYSEQ_OBJ_1))
  # Format ASV table to be used for machine learning applications and make metadata df
  ASV.df <- t(ASV.df)
  # Shorten metadata table
  ASV_meta.df <- data.frame(Sample = rownames(ASV_metadata.df), SoilLocation = ASV_metadata.df$SoilLocation, SoilHabitat = ASV_metadata.df$SoilHabitat, 
                            Drought.or.Watered = ASV_metadata.df$Drought.or.Watered, Genotype = ASV_metadata.df$Genotype, SoilInoculum = ASV_metadata.df$SoilInoculum)
  ASV_prefiltered.df <- cbind(ASV.df, ASV_meta.df)
  # 80% for train set
  set.seed(1095678254) # seed for sampling of samples
  train_index <- as.data.frame(ASV_prefiltered.df %>% sample_n(158))
  rownames(train_index) <- train_index$Sample
  train_index <- match(rownames(train_index), rownames(ASV_prefiltered.df))
  train_x <- as.data.frame(ASV.df[train_index, ])
  test_y <- as.data.frame(ASV.df[-train_index, ])
  dim(train_x) #158
  dim(test_y) #40
  # Train set, 158 samples
  train_x$Sample <- rownames(train_x)
  Training_meta.df <- merge(train_x, ASV_meta.df, by = 'Sample')
  train_x <- subset(Training_meta.df, select = -c(SoilLocation, SoilHabitat, Drought.or.Watered, Genotype, SoilInoculum))
  rownames(train_x) <- train_x$Sample
  train_x <- subset(train_x, select = -c(Sample))
  Training_meta.df <- subset(Training_meta.df, select = c(Sample, SoilLocation, SoilHabitat, Drought.or.Watered, Genotype, SoilInoculum))
  rownames(Training_meta.df) <- Training_meta.df$Sample 
  # Test set, 40 samples
  test_y$Sample <- rownames(test_y)
  Testing_meta.df <- merge(test_y, ASV_meta.df, by = "Sample")
  test_y <- subset(Testing_meta.df, select = -c(SoilLocation, SoilHabitat, Drought.or.Watered, Genotype, SoilInoculum))
  rownames(test_y) <- test_y$Sample
  test_y <- subset(test_y, select = -c(Sample))
  Testing_meta.df <- subset(Testing_meta.df, select = c(Sample, SoilLocation, SoilHabitat, Drought.or.Watered, Genotype, SoilInoculum))
  rownames(Testing_meta.df) <- Testing_meta.df$Sample 
  # Training model
  Training_grid <- expand.grid(.mtry = seq(10, length(train_x), round(length(train_x)*0.1)), .splitrule= "gini",
                               .min.node.size = c(1, 5, 10))
  train_control <- trainControl(method="cv", number=10)
  RF_CM <- list()
  RF_CM[["RF_model"]] <- train(x = train_x, y = Training_meta.df[[GROUPING]], method = ranger_funcs, importance = "permutation",
                               tuneGrid = Training_grid, trControl = train_control, num.trees = TREES)
  RF_prediction <- predict(RF_CM[["RF_model"]], test_y)
  RF_CM[["CMatrix"]] <- confusionMatrix(RF_prediction, as.factor(Testing_meta.df[[GROUPING]]), mode = "everything")
  RF_CM[["CMatrixPLOT"]] <- ggplotConfusionMatrix(RF_CM[["CMatrix"]])
  RF_CM[["VarImporance"]] <- varImp(RF_CM[["RF_model"]])
  return(RF_CM)
}

# Function extract taxonomy for ASVs with importance data from either a optimal model or per fold model
Get_taxa_importance <- function(RF_MODEL = NULL, FOLD_MODEL = NULL, TOP_N = 100){
  if (!is.null(RF_MODEL)){
    # Make vector of ASV numbers with high importance (specified or limited by TOP_N)
    X <- c(order(varImp(RF_model[["RF_model"]], scale = FALSE)$importance, decreasing=TRUE)[1:TOP_N])
    my_vect <- c()
    # For loop to go through ASV numbers and get ASV_hashs
    for (i in X) {
      Z <- rownames(varImp(RF_model[["RF_model"]], scale = FALSE)$importance)[i]
      my_vect <- append(my_vect, c(Z))  
    }
    #return(my_vect) # for testing XXX
    # For loop to get taxonomic assignments of the ASV_Hashs
    taxa_df <- data.frame(Kingdom=character(), Phylum=character(), Class=character(), Order=character(), Family=character(), Genus=character(), Species=character(), stringsAsFactors=FALSE)
    for (i in my_vect){
      taxa_df[i, ] <- c(drt.bact.late.clr@tax_table[i,])
    }
    #return(taxa_df) # for testing XXX
    # Connect these back to their importance values
    final_df <- merge(taxa_df,varImp(RF_model[["RF_model"]], scale = FALSE)$importance, by="row.names")
    colnames(final_df) <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Importance")
    # Fix first column name, coerce all columns but importance to factor()
    cols <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    final_df[,cols] <- lapply(final_df[cols], factor) 
    return(final_df)
  } else if (!is.null(FOLD_MODEL)){
    X <- as.data.frame(FOLD_MODEL$variable.importance)
    X <- c(order(X, decreasing=TRUE)[1:TOP_N])
    my_vect <- c()
    # For loop to go through ASV numbers and get ASV_hashs
    for (i in X) {
      Z <- rownames(as.data.frame(FOLD_MODEL$variable.importance))[i]
      my_vect <- append(my_vect, c(Z))  
      #return(my_vect) # for testing XXX
    }
    taxa_df <- data.frame(Kingdom=character(), Phylum=character(), Class=character(), Order=character(), Family=character(), Genus=character(), Species=character(), stringsAsFactors=FALSE)
    for (i in my_vect){
      taxa_df[i, ] <- c(drt.bact.late.clr@tax_table[i,])
    }
    #return(taxa_df) # for testing XXX
    # Connect these back to their importance values
    final_df <- merge(taxa_df, as.data.frame(FOLD_MODEL$variable.importance, to = c(0,100)), by="row.names")
    colnames(final_df) <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Importance")
    # Fix first column name, coerce all columns but importance to factor()
    cols <- c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    final_df[,cols] <- lapply(final_df[cols], factor) 
    return(final_df)
  }
}

# Function to obtain importance values with taxonomy for ASVs from a list of Ranger random forest fold outputs from caret.
# This should take ~3 mins per list (its getting almost 8k values from 10 folds).
Extract_importance <- function(list_of_imp = NULL){
  VarI_list<- c()
  for (i in seq(1,10,1)){
    X <- list_of_imp[[i]]$out
    X <- Get_taxa_importance(FOLD_MODEL = X, TOP_N = 612)
    VarI_list[[i]]<- X
  }
  return(VarI_list)
}

# Merge importance columns from all entries
merge_imp <- function(obj){
  X <- cbind(obj[[1]],
             "Importance_02" = obj[[2]][,c(9)],
             "Importance_03" = obj[[3]][,c(9)],
             "Importance_04" = obj[[4]][,c(9)],
             "Importance_05" = obj[[5]][,c(9)],
             "Importance_06" = obj[[6]][,c(9)],
             "Importance_07" = obj[[7]][,c(9)],
             "Importance_08" = obj[[8]][,c(9)],
             "Importance_09" = obj[[9]][,c(9)],
             "Importance_10" = obj[[10]][,c(9)])
  return(X)
}


#### Load datasets for use throughout ####
drt.bact.late <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS')
drt.bact.late.clr <- readRDS('Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS')

# Running with 10 CV on 80:20 dataspilt
## Already Ran
set.seed(1095678254) # also set above in the ML function
# RF_CM_Treatment <- MachineLearning_RF_ranger(drt.bact.late.clr, "Drought.or.Watered", 300)
# saveRDS(RF_CM_Treatment, "Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
# RF_CM_SoilInoculum <- MachineLearning_RF_ranger(drt.bact.late.clr, "SoilInoculum", 300)
# saveRDS(RF_CM_SoilInoculum, "Intermediate_data/RF_CM_SoilInoculum_CV10_300trees_rand_split.rds")
# RF_CM_Genotype <- MachineLearning_RF_ranger(drt.bact.late.clr, "Genotype", 300)
# saveRDS(RF_CM_Genotype, "Intermediate_data/RF_CM_Genotype_CV10_300trees_rand_split.rds")

# Loading completed ML runs from previous run
RF_CM_Treatment     <- readRDS("Intermediate_data/RF_CM_Treatment_CV10_300trees_rand_split.rds")
RF_CM_SoilInoculum  <- readRDS("Intermediate_data/RF_CM_SoilInoculum_CV10_300trees_rand_split.rds")
RF_CM_Genotype      <- readRDS("Intermediate_data/RF_CM_Genotype_CV10_300trees_rand_split.rds")
# Confusion Matrix plot for supplement
a <- RF_CM_Treatment$CMatrixPLOT
b <- RF_CM_SoilInoculum$CMatrixPLOT
# fix labels for soil inoculum
b <- b + scale_x_discrete(name = "Reference", labels = c(expression(TLI[P]), expression(TLI[Ag]), expression(SVR[P]), expression(SVR[Ag]), expression(KNZ[P]), expression(HAY[P])))
b <- b + scale_y_discrete(name = "Prediction", labels = c(expression(HAY[P]), expression(KNZ[P]),  expression(SVR[Ag]), expression(SVR[P]), expression(TLI[Ag]), expression(TLI[P])))
c <- RF_CM_Genotype$CMatrixPLOT
CM_PLOT <- ggarrange(b,c, align = "hv", common.legend = TRUE, nrow = 1, legend = 'right', labels = "AUTO")
ggsave("figures/Confusion_matrix_plots.svg", CM_PLOT, width = 8, height = 4)
ggsave("figures/Confusion_matrix_plots.png", CM_PLOT, width = 8, height = 4)

# Optimal hyperparameters
RF_CM_Treatment$RF_model
RF_CM_SoilInoculum$RF_model
RF_CM_Genotype$RF_model

# Optimal hyperparameters 
"factor mtry  min.node.size
Treatment 559 10
Soil Inoculum 437 5
Genotype  599  5
"

# Plotting ASVs that aid in predictions of the main effects
# Load a list of files to get mean and standard error for the OOB samples
# Treatment
files <- list.files(path = "./rangerML_stats/treatment/", pattern = ".*559_gini_10.*")
files <- paste("./rangerML_stats/treatment/", files, sep = "")
RF_list_treatment <- lapply(files, function(x) mget(load(x)))
RF_list_treatment <- Extract_importance(RF_list_treatment)
RF_list_treatment <- merge_imp(RF_list_treatment)
RF_list_treatment[19] <- as.data.frame(rowSds(as.matrix(RF_list_treatment[,9:18])))
RF_list_treatment[20] <- as.data.frame(rowMeans((RF_list_treatment[,9:18])))
colnames(RF_list_treatment)[19] <- "SD"
colnames(RF_list_treatment)[20] <- "mean"
RF_list_treatment <- plyr::join(data.frame('ASV' = rownames(drt.bact.late.clr@tax_table)), RF_list_treatment) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_treatment$ASV_numb <- paste("ASV", rownames(RF_list_treatment))
# soilInoculum
files <- list.files(path = "./rangerML_stats/soil_inoculum/", pattern = ".*437_gini_5.*")
files <- paste("./rangerML_stats/soil_inoculum/", files, sep = "")
RF_list_soil_inoculum <- lapply(files, function(x) mget(load(x)))
RF_list_soil_inoculum <- Extract_importance(RF_list_soil_inoculum)
RF_list_soil_inoculum <- merge_imp(RF_list_soil_inoculum)
RF_list_soil_inoculum[19] <- as.data.frame(rowSds(as.matrix(RF_list_soil_inoculum[,9:18])))
RF_list_soil_inoculum[20] <- as.data.frame(rowMeans((RF_list_soil_inoculum[,9:18])))
colnames(RF_list_soil_inoculum)[19] <- "SD"
colnames(RF_list_soil_inoculum)[20] <- "mean"
RF_list_soil_inoculum <- plyr::join(data.frame('ASV' = rownames(drt.bact.late.clr@tax_table)), RF_list_soil_inoculum) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_soil_inoculum$ASV_numb <- paste("ASV", rownames(RF_list_soil_inoculum))
# Genotype
files <- list.files(path = "./rangerML_stats/genotype/", pattern = ".*559_gini_5.*")
files <- paste("./rangerML_stats/genotype/", files, sep = "")
RF_list_genotype <- lapply(files, function(x) mget(load(x)))
RF_list_genotype <- Extract_importance(RF_list_genotype)
RF_list_genotype <- merge_imp(RF_list_genotype)
RF_list_genotype[19] <- as.data.frame(rowSds(as.matrix(RF_list_genotype[,9:18])))
RF_list_genotype[20] <- as.data.frame(rowMeans((RF_list_genotype[,9:18])))
colnames(RF_list_genotype)[19] <- "SD"
colnames(RF_list_genotype)[20] <- "mean"
RF_list_genotype <- plyr::join(data.frame('ASV' = rownames(drt.bact.late.clr@tax_table)), RF_list_genotype) # Fix to adjust order of dataframe to match phyloseq_taxtable
RF_list_genotype$ASV_numb <- paste("ASV", rownames(RF_list_genotype))

# Subplots for each main effect
a <- ggplot(RF_list_treatment[RF_list_treatment$mean > 0.001, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = phyla_palette, labels = c("Actinomycetota", "Bacteroidota", "Gemmatimonadota", "Pseudomonadota")) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Mean Decrease in Accuracy ") + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

b <- ggplot(RF_list_soil_inoculum[RF_list_soil_inoculum$mean > 0.001, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = phyla_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Mean Decrease in Accuracy ") + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

c <- ggplot(RF_list_genotype[RF_list_genotype$mean > 0.001, ], aes(x= mean, y = reorder(ASV_numb, mean), fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = phyla_palette) +
  geom_errorbar(aes(xmin=mean - SD/sqrt(length(mean)), xmax=mean + SD/sqrt(length(mean))), width=.2) +
  xlab("Mean Decrease in Accuracy ") + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

# Combination plot
ggarrange(a,b,c, common.legend = TRUE, labels = c("Treatment","Soil Inoculum", "Genotype"), align = "hv", legend = "right")
ASVs_deciding_ML <- ggarrange(a,b,c, nrow = 1, common.legend = TRUE, labels = c("AUTO"), align = "hv", legend = "right")

ggsave("figures/Important_ASV_by_factor.svg", ASVs_deciding_ML, width =12, height = 6)
ggsave("figures/Important_ASV_by_factor.png", ASVs_deciding_ML, width =12, height = 6)
