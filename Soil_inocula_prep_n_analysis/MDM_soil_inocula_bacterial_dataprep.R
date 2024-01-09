# Pre-processing of soil and soil-derived inocula bacterial communities
# Soils collected across an E/W precipitation gradient across KS
# Code by: Natalie Ford, Maggie Wagner, and Joel Swift

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(readxl); packageVersion('readxl')
library(phyloseq); packageVersion('phyloseq')
library(genefilter); packageVersion('genefilter')
library(ALDEx2); packageVersion('ALDEx2')
library(vegan); packageVersion('vegan')
source('./microfiltR_source_code.R') # source code for within-sample thresholding by reference to positive controls

####### Bacterial sequence prefix #######
ASVprefix <- 'b' 

####### Make sample metadata file #######
smd <- read_xlsx('./Intermediate_data/mdsi_sample_list_16S.xlsx') # load key linking samples to DNA libraries
soils <- read_xlsx('./Intermediate_data/mdsi_sample_metadata.xlsx') # load key with plant information
# Convert PlantIDs to factors
smd$Soil_ID <- as.factor(smd$Soil_ID)
soils$RepID <- as.factor(soils$RepID)

# Merge together:
smd <- right_join(soils,smd,by=c('RepID'='Soil_ID'))
smd$Sample_ID <- str_replace_all(smd$Sample_ID,'_','-') # replace _ with - for compatibility with sequencing output
# Make new column to combine soil Location and Type
smd$SoilSource <- factor(paste0(smd$Location, "_",smd$Type))
smd <- as.data.frame(smd)
rownames(smd) <- smd$Sample_ID # make the sample IDs the row names for import into Phyloseq

####### Load ASV table and taxa info #######
ASV <- readRDS('./Intermediate_data/16Smerged.seqtab.nochim.RDS')
tax <- read.delim('./Intermediate_data/16S_taxa_bacteria.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
bdrt <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
# Remove ag samples
bdrt <- subset_samples(bdrt, Land_use != "Agriculture")
ntaxa(bdrt) # 23059 ASVs

# Save original phyloseq object:
saveRDS(bdrt,'./Intermediate_data/phyloseq_b_original.RDS')

####### Look for plant contamination: #######
plantASVs <- subset_taxa(bdrt,Family=='mitochondria' | str_detect(Class,'Chloroplast')) 
sum(taxa_sums(plantASVs)) # 1850 observations
sum(taxa_sums(bdrt)) # out of 3033502 = 0.06% plant contamination
ntaxa(bdrt) # 23059 ASVs
nsamples(bdrt) # 47 samples

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bdrt.nobadASVs<-subset_taxa(bdrt,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bdrt.nobadASVs) # 11146 ASVs remain
sample_data(bdrt.nobadASVs)$UsableReads<-sample_sums(bdrt.nobadASVs)

####### Calculate alpha diversity before thresholding: #######
# Calculate alpha diversity and add to data frame:
alphadiv <- estimate_richness(bdrt.nobadASVs,measures=c('Chao1','Shannon', 'InvSimpson')) %>%
  as.data.frame() 
# Convert to true diversity https://twitter.com/derekseveri/status/1106359475756240898?lang=en
alphadiv$Shannon <- exp(alphadiv$Shannon)
rownames(alphadiv) <- str_replace_all(rownames(alphadiv),'\\.','-') # change row names . to - for compatibility with phyloseq
merge(alphadiv,sample_data(bdrt.nobadASVs),by='row.names') %>%
  column_to_rownames('Sample_ID') %>% dplyr::rename(Sample_ID=Row.names) -> smd
sample_data(bdrt.nobadASVs) <- smd  # update in phyloseq object

####### Across-sample thresholding: 5x25 (throw out "non-reproducible" ASVs) #######
threshold<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
bdrt.nobadASVs.thresholded<-filter_taxa(bdrt.nobadASVs,threshold,TRUE)
ntaxa(bdrt.nobadASVs.thresholded) # 750 ASVs remain
ntaxa(bdrt.nobadASVs) # 11146 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bdrt.nobadASVs.thresholded))/(sum(taxa_sums(bdrt.nobadASVs)))
# 71.6% of reads remain after thresholding

####### Remove samples with <400 usable reads #######
bdrt.nobadASVs.thresholded.highcoverage<-subset_samples(bdrt.nobadASVs.thresholded,UsableReads>=400) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bdrt.nobadASVs.thresholded.highcoverage) # 44 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs<-sample_sums(bdrt.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bdrt.nobadASVs.thresholded.highcoverage)$logObs<-log(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
sum(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # total = 1105088 reads
mean(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # mean = 25115.64 observations
sd(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) # sd = 24699.16
ntaxa(bdrt.nobadASVs.thresholded.highcoverage) # 750 ASVs 
sqrt(var(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs) /
       length(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Obs)) # SE = 3723.538 sqrt(var(X) / length(X))

# how many observations in Inocula dataset? 
temp_Inocula <- prune_samples(sample_data(bdrt.nobadASVs.thresholded.highcoverage)$Type == "Inocula", bdrt.nobadASVs.thresholded.highcoverage)
mean(sample_data(temp_Inocula)$Obs) # mean = 24326.43 observations
sd(sample_data(temp_Inocula)$Obs) # sd = 28035.39
sqrt(var(sample_data(temp_Inocula)$Obs) /
       length(sample_data(temp_Inocula)$Obs)) # SE = 5845.784 sqrt(var(X) / length(X))

# how many ASVs?
ntaxa(bdrt.nobadASVs.thresholded.highcoverage) # 750 bASVs
####### Quick unconstrained ordination to check for outliers #######
# transform ASV counts to proportions / relative abundance and ordinate: 
drt.bact.prop <- transform_sample_counts(bdrt.nobadASVs.thresholded.highcoverage, function(ASV) ASV/sum(ASV))
# Ordinate: 
pcoa.drt.bact.prop <- ordinate(drt.bact.prop,method='CAP',distance='bray',formula=~Type*Location+Condition(log(UsableReads)))
# Visualize quick ordination: --> didn't do this NF
plot_ordination(drt.bact.prop,pcoa.drt.bact.prop, color ='Location', shape='Type') +
  facet_wrap(~Type)+
  scale_color_manual(values=c('black','sky blue','forest green','magenta','goldenrod','darkorchid2','cadetblue'))+
  theme_classic()# no extreme outliers

with(as(sample_data(drt.bact.prop),'data.frame'),
     adonis2(as(otu_table(drt.bact.prop),'matrix') ~Location + Type,
            data=as(sample_data(drt.bact.prop),'data.frame')))
#             Df SumOfSqs      R2      F Pr(>F)    
#   Location  3   2.6211 0.21350 4.2280  0.001 ***
#   Type      1   1.5964 0.13004 7.7253  0.001 ***
#   Residual 39   8.0592 0.65646                  
#   Total    43  12.2767 1.00000

# Clean up environment
rm(pcoa.drt.bact.prop,drt.bact.prop)

####### Extract sample data for "final" dataset #######
drt.bact <- bdrt.nobadASVs.thresholded.highcoverage # rename final phyloseq object
smd.bact <- as(sample_data(drt.bact),'data.frame')
smd.bact$SampleID <- row.names(smd.bact) # Store Sample IDs as a column
summary(smd.bact) # Make sure everything looks right
mutate(smd.bact,Location=factor(Location), Land_use=factor(Land_use), Type=factor(Type)) -> smd.bact #recode variables as factors
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
saveRDS(drt.bact,'./Intermediate_data/phyloseq_b_clean_soil_inocula.RDS')

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
asvs <- as(otu_table(drt.bact,taxa_are_rows=FALSE),'matrix')
conds <- sample_data(drt.bact)$Type %>% as.character() # define "conditions" vector as the sample Type (in this case it's the same for all samples)

# Create CLR-transformed object:
asvs.clr <- aldex.clr(t(asvs),conds=conds,mc.samples=128,useMC=TRUE) %>% extractMedianCLRs(.)
drt.bact.clr <- drt.bact # copy phyloseq object

### Always double check whether taxa_are_rows should be TRUE or FALSE! 
colnames(asvs.clr) # note that ALDEx2 converted the - in sample names to . 
colnames(asvs.clr) <- str_replace_all(colnames(asvs.clr),'\\.','-')
rownames(asvs.clr) # confirms that taxa_are_rows = TRUE
otu_table(drt.bact.clr) <- otu_table(asvs.clr,taxa_are_rows = TRUE) # replace ASV table with CLR-transformed version

# Save transformed phyloseq object to file:
saveRDS(drt.bact.clr, './Intermediate_data/phyloseq_b_clean_soil_inocula_clr.RDS')
