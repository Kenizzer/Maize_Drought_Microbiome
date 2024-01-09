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

####### Bacterial sequence prefix #######
ASVprefix <- 'b' # use this if these are 16S sequences

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
ASV <- readRDS('../Raw_data_n_metadata/merged.seqtab.nochim.16S.RDS')
tax <- read.delim('../Raw_data_n_metadata/taxa_bacteria.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
bdrt <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
# remove agricultural soil sites
bdrt <- subset_samples(bdrt, SoilHabitat != "Agriculture")

# Save original phyloseq object:
saveRDS(bdrt,'Intermediate_data/phyloseq_b_original.RDS')

####### Look for plant contamination: #######
plantASVs <- subset_taxa(bdrt,Family=='mitochondria' | str_detect(Class,'Chloroplast')) 
sum(taxa_sums(plantASVs)) # 100955 observations
sum(taxa_sums(bdrt)) # out of 29196937 = 0.35% plant contamination
ntaxa(bdrt) # 17278 ASVs
nsamples(bdrt) # 289 samples

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bdrt.nobadASVs <- subset_taxa(bdrt,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bdrt.nobadASVs) # 10888 ASVs remain
sample_data(bdrt.nobadASVs)$UsableReads <- sample_sums(bdrt.nobadASVs)

####### Calculate alpha diversity before thresholding: #######
# Calculate alpha diversity and add to data frame:
alphadiv <- estimate_richness(bdrt.nobadASVs,measures=c('Chao1','Shannon', 'InvSimpson')) %>%
  as.data.frame() 
# Convert to true diversity https://twitter.com/derekseveri/status/1106359475756240898?lang=en
alphadiv$Shannon <- exp(alphadiv$Shannon)
rownames(alphadiv) <- str_replace_all(rownames(alphadiv),'\\.','-') # change row names . to - for compatibility with phyloseq
merge(alphadiv,sample_data(bdrt.nobadASVs),by='row.names') %>%
  column_to_rownames('SampleID') %>% dplyr::rename(SampleID=Row.names) -> smd
sample_data(bdrt.nobadASVs) <- smd  # update in phyloseq object

####### Across-sample thresholding: 5x25 (throw out "non-reproducible" ASVs) #######
threshold<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
bdrt.nobadASVs.thresholded<-filter_taxa(bdrt.nobadASVs,threshold,TRUE)
ntaxa(bdrt.nobadASVs.thresholded) # 468 ASVs remain
ntaxa(bdrt.nobadASVs) # 10888 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bdrt.nobadASVs.thresholded))/(sum(taxa_sums(bdrt.nobadASVs)))
# 98% of reads remain after thresholding

####### Remove samples with <400 usable reads #######
bdrt.nobadASVs.thresholded.highcoverage<-subset_samples(bdrt.nobadASVs.thresholded,UsableReads>=400) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bdrt.nobadASVs.thresholded.highcoverage) # 200 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs<-sample_sums(bdrt.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bdrt.nobadASVs.thresholded.highcoverage)$logObs<-log(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
sum(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # total = 28098900 reads
mean(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # mean = 140494 observations
sd(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # sd = 235675
ntaxa(bdrt.nobadASVs.thresholded.highcoverage) # 468 ASVs 
sqrt(var(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) /
       length(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs)) # SE = 16665 sqrt(var(X) / length(X))

# how many observations in late dataset?
temp_late <- prune_samples(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Timepoint == "late", bdrt.nobadASVs.thresholded.highcoverage)
mean(sample_data(temp_late)$Obs) # mean = 123554 observations
sd(sample_data(temp_late)$Obs) # sd = 213919
sqrt(var(sample_data(temp_late)$Obs) /
       length(sample_data(temp_late)$Obs)) # SE = 18276 sqrt(var(X) / length(X))
# how many ASVs?
ntaxa(bdrt.nobadASVs.thresholded.highcoverage) # 468 bASVs

####### Quick unconstrained ordination to check for outliers #######
# transform ASV counts to proportions / relative abundance and ordinate:
drt.bact.prop <- transform_sample_counts(bdrt.nobadASVs.thresholded.highcoverage, function(ASV) ASV/sum(ASV))
# Ordinate:
pcoa.drt.bact.prop <- ordinate(drt.bact.prop,method='CAP',distance='bray',formula=~SoilLocation*Drought.or.Watered+Condition(log(UsableReads)))
# Visualize quick ordination:
plot_ordination(drt.bact.prop,pcoa.drt.bact.prop,color='SoilLocation') +
  facet_wrap(~Drought.or.Watered)+
  scale_color_manual(values=c('black','sky blue','forest green','magenta','goldenrod','darkorchid2','cadetblue'))+
  theme_classic()# no extreme outliers

with(as(sample_data(drt.bact.prop),'data.frame'),
     adonis2(as(otu_table(drt.bact.prop),'matrix') ~ SoilLocation*Drought.or.Watered*Genotype,
             strata = Plate, data=as(sample_data(drt.bact.prop),'data.frame')))
"                                          Df SumOfSqs      R2      F Pr(>F)    
SoilLocation                               4    2.641 0.06448 3.4336  0.001 ***
Drought.or.Watered                         1    0.925 0.02257 4.8083  0.001 ***
Genotype                                   1    0.081 0.00198 0.4224  0.980    
SoilLocation:Drought.or.Watered            4    0.959 0.02340 1.2461  0.091 .  
SoilLocation:Genotype                      4    0.890 0.02173 1.1573  0.250    
Drought.or.Watered:Genotype                1    0.155 0.00377 0.8036  0.605    
SoilLocation:Drought.or.Watered:Genotype   4    0.699 0.01707 0.9090  0.628    
Residual                                 180   34.618 0.84500                  
Total                                    199   40.968 1.00000"

# Clean up environment
rm(pcoa.drt.bact.prop,drt.bact.prop)

####### Extract sample data for "final" dataset #######
drt.bact <- bdrt.nobadASVs.thresholded.highcoverage # rename final phyloseq object
smd.bact <- as(sample_data(drt.bact),'data.frame')
smd.bact$SampleID <- row.names(smd.bact) # Store Sample IDs as a column
summary(smd.bact) # Make sure everything looks right
mutate(smd.bact,Plate=factor(Plate),Location=factor(Location), Timepoint=factor(Timepoint)) -> smd.bact #recode variables as factors
####### Update sample data with Z-transformation of sequencing depth #######
smd.bact <- mutate(smd.bact, logObs.z = (logObs-mean(logObs))/sd(logObs))
row.names(smd.bact) <- smd.bact$SampleID
# Add back into phyloseq object:
sample_data(drt.bact) <- smd.bact

####### Relabel ASVs for convenience #######
tax <- as(tax_table(drt.bact),'matrix')
ASV <- as(otu_table(drt.bact),'matrix')

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
drt.bact <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as(tax,'matrix')))
# Add back into phyloseq object:
sample_data(drt.bact) <- smd.bact
# Save cleaned dataset to file:
saveRDS(drt.bact,'Intermediate_data/phyloseq_b_asv_clean_allRoots.RDS')
####### Separate out uninoculated controls #######
drt.bact.uninoculated <- subset_samples(drt.bact, SoilLocation=='Control')
saveRDS(drt.bact.uninoculated,'Intermediate_data/phyloseq_b_asv_clean_uninoculated.RDS')
drt.bact.inoculated <- subset_samples(drt.bact, SoilLocation !='Control')
saveRDS(drt.bact.inoculated,'Intermediate_data/phyloseq_b_asv_clean_inoculated.RDS')

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
asvs.root <- as(otu_table(drt.bact.inoculated,taxa_are_rows=FALSE),'matrix')
conds.root <- sample_data(drt.bact.inoculated)$Type %>% as.character() # define "conditions" vector as the sample Type (in this case it's the same for all samples)

# Create CLR-transformed object:
asvs.root.clr <- aldex.clr(t(asvs.root),conds=conds.root,mc.samples=128,useMC=TRUE) %>% extractMedianCLRs(.)
drt.bact.inoculated.clr <- drt.bact.inoculated # copy phyloseq object

### Always double check whether taxa_are_rows should be TRUE or FALSE! 
colnames(asvs.root.clr) # note that ALDEx2 converted the - in sample names to . 
colnames(asvs.root.clr) <- str_replace_all(colnames(asvs.root.clr),'\\.','-')
rownames(asvs.root.clr) # confirms that taxa_are_rows = TRUE
otu_table(drt.bact.inoculated.clr) <- otu_table(asvs.root.clr,taxa_are_rows = TRUE) # replace ASV table with CLR-transformed version

# Save transformed phyloseq object to file:
saveRDS(drt.bact.inoculated.clr, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_clr.RDS')

####### Separate data set into parts by time point #######
# Non-transformed ASV tables:
drt.bact.early <- subset_samples(drt.bact.inoculated, Timepoint=='early')
drt.bact.middle <- subset_samples(drt.bact.inoculated, Timepoint=='middle')
drt.bact.late <- subset_samples(drt.bact.inoculated, Timepoint=='late')
# CLR-transformed ASV tables:
drt.bact.early.clr <- subset_samples(drt.bact.inoculated.clr, Timepoint=='early')
drt.bact.middle.clr <- subset_samples(drt.bact.inoculated.clr, Timepoint=='middle')
drt.bact.late.clr <- subset_samples(drt.bact.inoculated.clr, Timepoint=='late')
# Save into RDS objects
saveRDS(drt.bact.early.clr, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_early_clr.RDS')
saveRDS(drt.bact.early, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_early.RDS')
saveRDS(drt.bact.middle.clr, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_middle_clr.RDS')
saveRDS(drt.bact.middle, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_middle.RDS')
saveRDS(drt.bact.late.clr, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_late_clr.RDS')
saveRDS(drt.bact.late, 'Intermediate_data/phyloseq_b_asv_clean_inoculated_late.RDS')