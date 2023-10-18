# Analysis of maize B73 and Mo17 inoculated plants under drought stress
# Samples were collected from a greenhouse study in Feb 2020
# Microbiomes are from the root compartment (endosphere)
# Growth measurements were also recorded throughout the experiment (~50 days in length)
# Code by: Maggie Wagner & Joel Swift

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(readxl); packageVersion('readxl')
library(phyloseq); packageVersion('phyloseq')
library(genefilter); packageVersion('genefilter')
library(ALDEx2); packageVersion('ALDEx2')
library(vegan); packageVersion('vegan')
source('./microfiltR_source_code.R') # source code for within-sample thresholding by reference to positive controls

####### Fungal  sequence prefix #######
ASVprefix <- 'f' # use this if these are ITS sequences

####### Make sample metadata file #######
smd <- read_xlsx('../Raw_data_n_metadata/maize_drought_sampleDNA_list.xlsx') # load key linking samples to DNA libraries
plants <- read.delim('../Raw_data_n_metadata/Drought_Experiment_Data_FULL.csv',sep=',',header=TRUE) # load key with plant information
# Convert PlantIDs to factors
smd$Plant_ID <- as.factor(smd$Plant_ID)
plants$RepID <- as.factor(plants$RepID)

# Merge together:
smd <- right_join(plants,smd,by=c('RepID'='Plant_ID'))
smd$SampleID <- str_replace_all(smd$SampleID,'_','-') # replace _ with - for compatibility with sequencing output
rownames(smd) <- smd$SampleID # make the sample IDs the row names for import into Phyloseq

####### Load ASV table and taxa info #######
ASV <- readRDS('../Raw_data_n_metadata/drtITS-merged.seqtab.bothLanes.RDS')
tax <- read.delim('../Raw_data_n_metadata/taxa_fungi.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
fdrt <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
# remove agricultural soil sites
fdrt <- subset_samples(fdrt, SoilHabitat != "Agriculture")
# Stats
nsamples(fdrt)# 282 samples
ntaxa(fdrt) # 1699 ASVs
sum(taxa_sums(fdrt))# 16696929 reads

# Save original phyloseq object:
saveRDS(fdrt,'Intermediate_data/phyloseq_f_original.RDS')


####### Fix issue with the first column (kingdom) of the taxa_names #######
#  Current: Talaromyces_purpureogenus|JN899372|SH209381.07FU|refs|k__Fungi
# Expected: k__Fungi
# Solution: Grab first column of taxa table, search for K__Fungi, remove anything before it with a regex
tax_table(fdrt)[,1] <- gsub(as.data.frame(tax_table(fdrt))[,1], pattern = ".*k__Fungi", replacement = "k__Fungi")

###### Fix issue with all taxa_names columns being preceded by A__ ########
tax_table(fdrt) <-  gsub(tax_table(fdrt)[, colnames(tax_table(fdrt))], pattern = "[a-z]__", replacement = "")

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level 
fdrt.nobadASVs<- subset_taxa(fdrt, Kingdom == 'Fungi')
ntaxa(fdrt.nobadASVs) # 980 ASVs remain
sample_data(fdrt.nobadASVs)$UsableReads <- sample_sums(fdrt.nobadASVs)

####### Calculate alpha diversity before thresholding: #######
# Calculate alpha diversity and add to data frame:
alphadiv <- estimate_richness(fdrt.nobadASVs,measures=c('Chao1', 'Shannon','InvSimpson')) %>%
  as.data.frame() 
# Convert to true diversity https://twitter.com/derekseveri/status/1106359475756240898?lang=en
alphadiv$Shannon <- exp(alphadiv$Shannon)
#alphadiv$test <- 1/(1-alphadiv$Simpson) # This is the exact same as InvSimspon from vegan
rownames(alphadiv) <- str_replace_all(rownames(alphadiv),'\\.','-') # change row names . to - for compatibility with phyloseq
merge(alphadiv,sample_data(fdrt.nobadASVs),by='row.names') %>%
  column_to_rownames('SampleID') %>% dplyr::rename(SampleID=Row.names) -> smd
sample_data(fdrt.nobadASVs) <- smd  # update in phyloseq object

####### Across-sample thresholding: 5x25 (throw out "non-reproducible" ASVs) #######
threshold<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
fdrt.nobadASVs.thresholded<-filter_taxa(fdrt.nobadASVs,threshold,TRUE)
ntaxa(fdrt.nobadASVs.thresholded) # 69 ASVs remain
ntaxa(fdrt.nobadASVs) # 980 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(fdrt.nobadASVs.thresholded))/(sum(taxa_sums(fdrt.nobadASVs)))
# 94.5% of reads remain after thresholding

####### Remove samples with <400 usable reads #######
fdrt.nobadASVs.thresholded.highcoverage <- subset_samples(fdrt.nobadASVs.thresholded,UsableReads>=400) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(fdrt.nobadASVs.thresholded.highcoverage) # 133 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs <- sample_sums(fdrt.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(fdrt.nobadASVs.thresholded.highcoverage)$logObs <- log(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
sum(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # total = 5963181 reads
mean(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # mean = 44836 observations
sd(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # sd = 118238
ntaxa(fdrt.nobadASVs.thresholded.highcoverage) # 69 ASVs 
sqrt(var(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) /
       length(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs)) # SE = 10253 sqrt(var(X) / length(X))

# how many observations in late dataset?
temp_late <- prune_samples(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Timepoint == "late", fdrt.nobadASVs.thresholded.highcoverage)
mean(sample_data(temp_late)$Obs) # mean = 41081 observations
sd(sample_data(temp_late)$Obs) # sd = 99241
sqrt(var(sample_data(temp_late)$Obs) /
       length(sample_data(temp_late)$Obs)) # SE = 10129 sqrt(var(X) / length(X))

####### Quick unconstrained ordination to check for outliers #######
# transform ASV counts to proportions / relative abundance and ordinate:
drt.fungi.prop <- transform_sample_counts(fdrt.nobadASVs.thresholded.highcoverage, function(ASV) ASV/sum(ASV))
# Ordinate:
pcoa.drt.fungi.prop <- ordinate(drt.fungi.prop,method='CAP',distance='bray',formula=~SoilLocation*Drought.or.Watered+Condition(log(UsableReads)))
# Visualize quick ordination:
plot_ordination(drt.fungi.prop, pcoa.drt.fungi.prop,color='SoilLocation',shape='SoilHabitat') +
  facet_wrap(~Drought.or.Watered)+
  scale_color_manual(values=c('black','sky blue','forest green','magenta','goldenrod','darkorchid2','cadetblue'))+
  theme_classic()# no extreme outliers

with(as(sample_data(drt.fungi.prop),'data.frame'),
     adonis2(as(otu_table(drt.fungi.prop),'matrix') ~ SoilLocation*Drought.or.Watered*Genotype,
             strata = Plate, data=as(sample_data(drt.fungi.prop),'data.frame')))
"                                          Df SumOfSqs      R2      F Pr(>F)    
SoilLocation                               4    5.456 0.11952 4.4425  0.001 ***
Drought.or.Watered                         1    0.496 0.01086 1.6148  0.115    
Genotype                                   1    0.483 0.01059 1.5744  0.081 .  
SoilLocation:Drought.or.Watered            4    1.667 0.03651 1.3570  0.050 *  
SoilLocation:Genotype                      4    1.587 0.03476 1.2920  0.106    
Drought.or.Watered:Genotype                1    0.304 0.00667 0.9911  0.575    
SoilLocation:Drought.or.Watered:Genotype   3    0.655 0.01435 0.7111  0.883    
Residual                                 114   35.002 0.76675                  
Total                                    132   45.650 1.00000 "

# Clean up
rm(pcoa.drt.fungi.prop,drt.fungi.prop)

####### Extract sample data for "final" dataset #######
drt.fungi <- fdrt.nobadASVs.thresholded.highcoverage # rename final phyloseq object
smd.fungi <- as(sample_data(drt.fungi),'data.frame')
smd.fungi$SampleID <- row.names(smd.fungi) # Store Sample IDs as a column
summary(smd.fungi) # Make sure everything looks right
mutate(smd.fungi,Plate=factor(Plate),Location=factor(Location), Timepoint=factor(Timepoint)) -> smd.fungi #recode variables as factors
####### Update sample data with Z-transformation of sequencing depth #######
smd.fungi <- mutate(smd.fungi, logObs.z = (logObs-mean(logObs))/sd(logObs))
row.names(smd.fungi) <- smd.fungi$SampleID
# Add back into phyloseq object:
sample_data(drt.fungi) <- smd.fungi

####### Relabel ASVs for convenience #######
tax <- as(tax_table(drt.fungi),'matrix')
ASV <- as(otu_table(drt.fungi),'matrix')

# Make sure ASV sequences match in taxonomy & ASV tables
setdiff(rownames(tax),colnames(ASV)) # no differences
setdiff(colnames(ASV),rownames(tax)) # no differences

# For security, make a named list of sequences and ASV_IDs to use in relabeling:
ASVlist <- paste0(ASVprefix,'ASV_',1:ncol(ASV)) # list of ASV_IDs 
names(ASVlist) <- colnames(ASV) # name ASV_IDs by their sequences

# Tax table: 
# add column to hold sequences:
tax <- transform(tax, Sequence = row.names(tax))
# then, relabel the rows:
for (i in 1:nrow(tax)) {
  seq <- tax[['Sequence']][i] # get the ASV sequence for this row
  row.names(tax)[i] <- ASVlist[[seq]] # then match it to the named list to get the short ASV_ID, and store that as the row name
  print(paste0('relabeling ',ASVprefix,'ASV_',i,'...')) }

# ASV table: relabel columns
for (i in 1:ncol(ASV)) { 
  seq <- colnames(ASV)[i] # Currently the columns of the ASV table are named as the full ASV sequences; get the sequence for each column
  colnames(ASV)[i] <- ASVlist[[seq]] # then match it to the named list to get the short ASV_ID, and store that as the column name
  print(paste0('relabeling ',ASVprefix,'ASV_',i,' out of ',ncol(ASV),'...'))
}

####### Construct and save final Phyloseq object: #######
# Add relabeled tax & OTU tables back into phyloseq object:
drt.fungi <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as(tax,'matrix')))
# Add back into phyloseq object:
sample_data(drt.fungi) <- smd.fungi
# Save cleaned dataset to file:
saveRDS(drt.fungi,'Intermediate_data/phyloseq_f_asv_clean_allRoots.RDS')
####### Separate out uninoculated controls #######
drt.fungi.uninoculated <- subset_samples(drt.fungi, SoilLocation=='Control')
saveRDS(drt.fungi.uninoculated,'Intermediate_data/phyloseq_f_asv_clean_uninoculated.RDS')
drt.fungi.inoculated <- subset_samples(drt.fungi, SoilLocation !='Control')
saveRDS(drt.fungi.inoculated,'Intermediate_data/phyloseq_f_asv_clean_inoculated.RDS')

####### Define function to extract median CLR-transformed value for each cell in ASV table #######
extractMedianCLRs <- function(clr.obj){
  for (sample in getSampleIDs(clr.obj)) {
    MCinstances <- getMonteCarloInstances(clr.obj)[[sample]]
    medianCLRs <- apply(MCinstances,MARGIN=1,function(x) median(x))
    if (sample == getSampleIDs(clr.obj)[1]) {
      medianCLRs.all <- data.frame(medianCLRs)
      colnames(medianCLRs.all) <- sample }
    else {medianCLRs.all[[sample]] <- medianCLRs} }
  return(medianCLRs.all) }

####### CLR transformation for inoculated root data #######
### Always double check whether taxa_are_rows should be TRUE or FALSE! 
asvs.root <- as(otu_table(drt.fungi.inoculated,taxa_are_rows=FALSE),'matrix')
conds.root <- sample_data(drt.fungi.inoculated)$Type %>% as.character() # define "conditions" vector as the sample Type (in this case it's the same for all samples)

# Create CLR-transformed object:
asvs.root.clr <-aldex.clr(t(asvs.root),conds=conds.root,mc.samples=128,useMC=TRUE) %>% extractMedianCLRs(.)

drt.fungi.inoculated.clr <- drt.fungi.inoculated # copy phyloseq object
### Always double check whether taxa_are_rows should be TRUE or FALSE! 
colnames(asvs.root.clr) # note that ALDEx2 converted the - in sample names to . 
colnames(asvs.root.clr) <- str_replace_all(colnames(asvs.root.clr),'\\.','-')
rownames(asvs.root.clr) # confirms that taxa_are_rows = TRUE
otu_table(drt.fungi.inoculated.clr) <- otu_table(asvs.root.clr,taxa_are_rows = TRUE) # replace ASV table with CLR-transformed version

# Save transformed phyloseq object to file:
saveRDS(drt.fungi.inoculated.clr, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_clr.RDS')

####### Separate dataset into parts by timepoint #######
# Non-transformed ASV tables:
drt.fungi.early <- subset_samples(drt.fungi.inoculated, Timepoint=='early')
drt.fungi.middle <- subset_samples(drt.fungi.inoculated, Timepoint=='middle')
drt.fungi.late <- subset_samples(drt.fungi.inoculated, Timepoint=='late')
# CLR-transformed ASV tables:
drt.fungi.early.clr <- subset_samples(drt.fungi.inoculated.clr, Timepoint=='early')
drt.fungi.middle.clr <- subset_samples(drt.fungi.inoculated.clr, Timepoint=='middle')
drt.fungi.late.clr <- subset_samples(drt.fungi.inoculated.clr, Timepoint=='late')
# Save into RDS objects
saveRDS(drt.fungi.early.clr, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_early_clr.RDS')
saveRDS(drt.fungi.early, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_early.RDS')
saveRDS(drt.fungi.middle.clr, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_middle_clr.RDS')
saveRDS(drt.fungi.middle, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_middle.RDS')
saveRDS(drt.fungi.late.clr, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_late_clr.RDS')
saveRDS(drt.fungi.late, 'Intermediate_data/phyloseq_f_asv_clean_inoculated_late.RDS')
