library(DESeq2)
# make the metadata
# # ensure the count files are where you think they are
list.files(directory)
# 
sampleFiles <- list.files(directory, pattern = ".*counts")
length(sampleFiles)

# create a vector of sample names. ensure these are in the same order as the "sampleFiles" object!----
sampleNames <- str_match(sampleFiles, "trim_\\s*(.*?)\\s*_S")[,2]
last10names <-c( "gv-01-10-sp","gv-03-05-sp","gv-03-06-sp","gv-04-06-sp","gv-06-12-sp","gv-01-05-vs","gv-03-03-vs",
                 "gv-03-06-vs","gv-04-06-vs","gv-03-05-vs")
sampleNames[42:51]<-last10names

## read in the metadata ----

metadata <- read.csv("yourpath/Data Repo/metadata_for_woodfrog_winter_RNAseq.csv")



# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
  sampleName = metadata$sampleNames,
  fileName = sampleFiles,
  period =metadata$sample_period,
  sex = metadata$sex,
  tissue = metadata$tissue,
  svl = metadata$svl,
  mass =metadata$mass,
  parasitized = metadata$parasitized,
  rhabdias = metadata$rhabdias,
  oswaldocruzia = metadata$oswaldocruzia,
  enclosure = metadata$enclosure_id, 
  Fall = metadata$Fall,
  Winter = metadata$Winter,
  PostWinter1 = metadata$PostWinter1,
  PostWinter2= metadata$PostWinter2
)





# look at the data frame to ensure it is what you expect:
str(sampleTable)


sampleTable <- sampleTable %>%mutate(names=sub("(.+)_S.*", "\\1", fileName))%>%
  mutate_all(~ ifelse(is.na(.), "0", .))

sampleTable$rhabdias<-as.numeric(sampleTable$rhabdias)
sampleTable$oswaldocruzia<-as.numeric(sampleTable$oswaldocruzia)


#log transform (after +1 worm data)
sampleTable$rhabdias_log<- log(sampleTable$rhabdias+1)
sampleTable$oswaldocruzia_log<- log(sampleTable$oswaldocruzia+1)
sampleTable$period<- as.factor(sampleTable$period)

# I want to analyze ventral skin and spleen separately so I'm going to split up the sampleTable (metadata)
vs_sampleTable = filter(sampleTable, tissue == "vs") %>% mutate(sample_period = as.factor(period))
sp_sampleTable = filter(sampleTable, tissue == "sp")%>% mutate(sample_period = as.factor(period))
sampleTable_all<-sampleTable %>%mutate(names=sub("(.+)_S.*", "\\1", fileName))




# Analyze our data using DESeq2----


# we've already created the necessary component objects
# use them to create a "DESeqDataSet" object

#all samples to compare tissues----
m_all<- m[,colnames(m) %in%sampleTable$names]
# Remove rows where row names start with "NA"
m_all <- m_all[!grepl("^NA", rownames(m_all)), , drop = FALSE]

all_ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = m_all,
  colData = sampleTable_all,
  design = ~ tissue * period + rhabdias_log + sex+ oswaldocruzia_log 
)
colnames(all_ddsHTSeq)<-sampleTable_all$sampleName
# sum counts for each gene across samples
all_sumcounts <- rowSums(counts(all_ddsHTSeq))


#sample counts
hist(colSums(counts(all_ddsHTSeq)))

# take the log
all_logsumcounts <- log(all_sumcounts,base=10)

# get genes with summed counts greater than 50
all_keep <- all_sumcounts > 50

# keep only the genes for which the vector "keep" is TRUE
all_ddsHTSeq <- all_ddsHTSeq[all_keep,]

#keep only samples with counts over 5
all_sampcounts <-colSums(counts(all_ddsHTSeq))
all_keepSamps <- all_sampcounts > 5

#filter the dataset
all_ddsHTSeq <- all_ddsHTSeq[,all_keepSamps]

#do the analysis#run deseq2

all_dds <- DESeq(all_ddsHTSeq)
all_ddsClean<- all_dds[which(mcols(all_dds)$betaConv),]

# list coefficient
resultsNames(all_ddsClean)


##### just ventral_skin samples----

vs_sampleTable <- vs_sampleTable %>%mutate(names=sub("(.+)_S.*", "\\1", fileName))
m_vs<- m[,colnames(m) %in%vs_sampleTable$names]

# Remove rows where row names start with "NA"
m_vs <- m_vs[!grepl("^NA", rownames(m_vs)), , drop = FALSE]


vs_ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = m_vs,
  colData = vs_sampleTable,
  design = ~  period + rhabdias_log + sex+ oswaldocruzia_log
)

colnames(vs_ddsHTSeq)<-vs_sampleTable$sampleName

# what does expression look like across genes?

# sum counts for each gene across samples
vs_sumcounts <- rowSums(counts(vs_ddsHTSeq))

# take the log
vs_logsumcounts <- log(vs_sumcounts,base=10)

# get genes with summed counts greater than 50
vs_keep <- vs_sumcounts > 50

# keep only the genes for which the vector "keep" is TRUE
vs_ddsHTSeq <- vs_ddsHTSeq[vs_keep,]

#keep only samples with counts over 5
vs_sampcounts <-colSums(counts(vs_ddsHTSeq))
vs_keepSamps <- vs_sampcounts > 5

#filter the dataset
vs_ddsHTSeq <- vs_ddsHTSeq[,vs_keepSamps]

#do the analysis#run deseq2
vs_dds <- DESeq(vs_ddsHTSeq)
vs_ddsClean<- vs_dds[which(mcols(vs_dds)$betaConv),]

# list coefficients
resultsNames(vs_ddsClean)



##### just spleen samples----
m_sp<- m[,colnames(m) %in%sp_sampleTable$names]
dim(m_sp)
# Remove rows where row names start with "NA"
m_sp <- m_sp[!grepl("^NA", rownames(m_sp)), , drop = FALSE]


sp_ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = m_sp,
  colData = sp_sampleTable,
  design = ~  period + rhabdias_log + sex+ oswaldocruzia_log
)

colnames(sp_ddsHTSeq)<-sp_sampleTable$sampleName

# sum counts for each gene across samples
sp_sumcounts <- rowSums(counts(sp_ddsHTSeq))

# take the log
sp_logsumcounts <- log(sp_sumcounts,base=10)

# get genes with summed counts greater than 50
sp_keep <- sp_sumcounts > 50

# keep only the genes for which the vector "keep" is TRUE
sp_ddsHTSeq <- sp_ddsHTSeq[sp_keep,]

#keep only samples with counts over 5
sp_sampcounts <-colSums(counts(sp_ddsHTSeq))
sp_keepSamps <- sp_sampcounts > 5

#filter the dataset
sp_ddsHTSeq <- sp_ddsHTSeq[,sp_keepSamps]

#do the analysis#run deseq2

sp_dds <- DESeq(sp_ddsHTSeq)
sp_ddsClean <- sp_dds[which(mcols(sp_dds)$betaConv),]

# list coefficients
resultsNames(sp_ddsClean)

