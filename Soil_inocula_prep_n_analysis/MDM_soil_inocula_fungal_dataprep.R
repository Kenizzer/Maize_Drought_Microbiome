# Pre-processing of soil and soil-derived inocula fungal communities
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

####### Fungal  sequence prefix #######
ASVprefix <- 'f' 

####### Make sample metadata file #######
smd <- read_xlsx('./Intermediate_data/mdsi_sample_list_ITS.xlsx') # load key linking samples to DNA libraries
soils <- read_xlsx('./Intermediate_data/mdsi_sample_metadata.xlsx') # load key with plant information
# Convert PlantIDs to factors
smd$Soil_ID <- as.factor(smd$Soil_ID)
soils$RepID <- as.factor(soils$RepID)

# Merge together:
smd <- right_join(soils,smd,by=c('RepID'='Soil_ID'))
smd$Sample_ID <- str_replace_all(smd$Sample_ID,'_','-') # replace _ with - for compatibility with sequencing output
# Make new column to combine SoilLocation and SoilHabitat
smd$SoilHabitat <- factor(paste0(smd$Location,"_",smd$Land_use,"_",smd$Type))
smd$SoilSource <- factor(paste0(smd$Location, "_",smd$Land_use))
smd <- as.data.frame(smd)
rownames(smd) <- smd$Sample_ID # make the sample IDs the row names for import into Phyloseq

####### Load ASV table and taxa info #######
ASV <- readRDS('./Intermediate_data/ITS-merged.seqtab.bothLanes.RDS')
tax <- read.delim('./Intermediate_data/ITS_taxa_fungi.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
fdrt <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
nsamples(fdrt)# 95 samples
ntaxa(fdrt) # 9843 ASVs
sum(taxa_sums(fdrt))# 8110722 reads

# Save original phyloseq object:
saveRDS(fdrt,'./Intermediate_data/phyloseq_f_original.RDS')


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
ntaxa(fdrt.nobadASVs) # 5085 ASVs remain
sample_data(fdrt.nobadASVs)$UsableReads<-sample_sums(fdrt.nobadASVs)

####### Calculate alpha diversity before thresholding: #######
# Calculate alpha diversity and add to data frame:
alphadiv <- estimate_richness(fdrt.nobadASVs,measures=c('Chao1', 'Shannon','InvSimpson')) %>%
  as.data.frame() 
# Convert to true diversity https://twitter.com/derekseveri/status/1106359475756240898?lang=en
alphadiv$Shannon <- exp(alphadiv$Shannon)
#alphadiv$test <- 1/(1-alphadiv$Simpson) # This is the exact same as InvSimspon from vegan
rownames(alphadiv) <- str_replace_all(rownames(alphadiv),'\\.','-') # change row names . to - for compatibility with phyloseq
merge(alphadiv,sample_data(fdrt.nobadASVs),by='row.names') %>%
  column_to_rownames('Sample_ID') %>% dplyr::rename(Sample_ID=Row.names) -> smd
sample_data(fdrt.nobadASVs) <- smd  # update in phyloseq object

####### Across-sample thresholding: 5x25 (throw out "non-reproducible" ASVs) #######
threshold<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
fdrt.nobadASVs.thresholded<-filter_taxa(fdrt.nobadASVs,threshold,TRUE)
ntaxa(fdrt.nobadASVs.thresholded) # 564 ASVs remain
ntaxa(fdrt.nobadASVs) # 5085 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(fdrt.nobadASVs.thresholded))/(sum(taxa_sums(fdrt.nobadASVs)))
# 76.05% of reads remain after thresholding

####### Remove samples with <400 usable reads #######
fdrt.nobadASVs.thresholded.highcoverage<-subset_samples(fdrt.nobadASVs.thresholded,UsableReads>=400) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(fdrt.nobadASVs.thresholded.highcoverage) # 89 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs<-sample_sums(fdrt.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(fdrt.nobadASVs.thresholded.highcoverage)$logObs<-log(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
sum(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # total = 3394118 reads
mean(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # mean = 38136.16 observations
sd(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) # sd = 36935.86
ntaxa(fdrt.nobadASVs.thresholded.highcoverage) # 564 ASVs 
sqrt(var(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs) /
       length(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Obs)) # SE = 3915.194 sqrt(var(X) / length(X))

# how many observations in late dataset?
temp_late <- prune_samples(sample_data(fdrt.nobadASVs.thresholded.highcoverage)$Type == "Inocula", fdrt.nobadASVs.thresholded.highcoverage)
mean(sample_data(temp_late)$Obs) # mean = 22978.11 observations
sd(sample_data(temp_late)$Obs) # sd = 24878.05
sqrt(var(sample_data(temp_late)$Obs) /
       length(sample_data(temp_late)$Obs)) # SE = 3668.069 sqrt(var(X) / length(X))


####### Quick unconstrained ordination to check for outliers ####### 
# transform ASV counts to proportions / relative abundance and ordinate:
drt.fungi.prop <- transform_sample_counts(subset_samples(fdrt.nobadASVs.thresholded.highcoverage,Type!='NA_NA'),function(ASV) ASV/sum(ASV))
# Ordinate:
pcoa.drt.fungi.prop <- ordinate(drt.fungi.prop,method='CAP',distance='bray',formula=~Type*Location+Condition(log(UsableReads)))
# Visualize quick ordination:
plot_ordination(drt.fungi.prop, pcoa.drt.fungi.prop,color='Type',shape='Location') +
  facet_wrap(~Land_use)+
  scale_color_manual(values=c('black','sky blue','forest green','magenta','goldenrod','darkorchid2','cadetblue'))+
  theme_classic()# no extreme outliers

with(as(sample_data(drt.fungi.prop),'data.frame'),
     adonis2(as(otu_table(drt.fungi.prop),'matrix') ~ Land_use:Location+Type,
            data=as(sample_data(drt.fungi.prop),'data.frame')))

#                       Df   SumOfSqs      R2       F Pr(>F)    
# Type               1    3.782 0.11268 13.6577  0.001 ***
# Land_use:Location  5    7.074 0.21077  5.1094  0.001 ***
# Residual          82   22.706 0.67654                   
# Total             88   33.562 1.00000 

# Clean up 
rm(pcoa.drt.fungi.prop,drt.fungi.prop)

####### Extract sample data for "final" dataset #######
drt.fungi <- fdrt.nobadASVs.thresholded.highcoverage # rename final phyloseq object
smd.fungi <- as(sample_data(drt.fungi),'data.frame')
smd.fungi$Sample_ID <- row.names(smd.fungi) # Store Sample IDs as a column
summary(smd.fungi) # Make sure everything looks right
mutate(smd.fungi,Land_use=factor(Land_use),Location=factor(Location), Type=factor(Type)) -> smd.fungi #recode variables as factors
####### Update sample data with Z-transformation of sequencing depth #######
smd.fungi <- mutate(smd.fungi, logObs.z = (logObs-mean(logObs))/sd(logObs))
row.names(smd.fungi) <- smd.fungi$Sample_ID
# Add back into phyloseq object:
sample_data(drt.fungi) <- smd.fungi

####### Get list of ASVs for Venn diagram ########
soil_asvs.vec <- colnames(otu_table(prune_taxa(taxa_sums(subset_samples(drt.fungi, Type=="Soil")) > 1, subset_samples(drt.fungi, Type=="Soil"))))
inoc_asvs.vec <- colnames(otu_table(prune_taxa(taxa_sums(subset_samples(drt.fungi, Type=="Inocula")) > 1, subset_samples(drt.fungi, Type=="Inocula"))))

write.csv(soil_asvs.vec, "Intermediate_data/ASV_list_for_venn_diagram_soil_ITS.csv")
write.csv(inoc_asvs.vec, "Intermediate_data/ASV_list_for_venn_diagram_inoc_ITS.csv")

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
saveRDS(drt.fungi,'./Intermediate_data/phyloseq_f_clean_soil_inocula.RDS')

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
asvs <- as(otu_table(drt.fungi,taxa_are_rows=FALSE),'matrix')
conds <- sample_data(drt.fungi)$Type %>% as.character() # define "conditions" vector as the sample Type (in this case it's the same for all samples)

# Create CLR-transformed object:
asvs.clr <-aldex.clr(t(asvs),conds=conds,mc.samples=128,useMC=TRUE) %>% extractMedianCLRs(.)

drt.fungi.clr <- drt.fungi # copy phyloseq object
### Always double check whether taxa_are_rows should be TRUE or FALSE! 
colnames(asvs.clr) # note that ALDEx2 converted the - in sample names to . 
colnames(asvs.clr) <- str_replace_all(colnames(asvs.clr),'\\.','-')
rownames(asvs.clr) # confirms that taxa_are_rows = TRUE
otu_table(drt.fungi.clr) <- otu_table(asvs.clr,taxa_are_rows = TRUE) # replace ASV table with CLR-transformed version

# Save transformed phyloseq object to file:
saveRDS(drt.fungi.clr, './Intermediate_data/phyloseq_f_clean_soil_inocula_clr.RDS')
