library(purrr)
library(dplyr)
# create a vector containing transcript lengths
# read in length
tab_gtf <- read.delim("yourpath/Data Repo/merged.gtf",header=FALSE, skip = 2)

# filter down to exons only, create a column with the gene ID and the transcript ID
exons <- filter(tab_gtf,V3=="exon") %>%
  mutate(exlen=V5-V4) %>%
  mutate(tsid=str_extract(V9,"(?<=transcript_id )[^;]*")) %>%
  mutate(gid=str_extract(V9,"(?<=gene_id )[^;]*"))

# get transcript lengths by grouping by transcript/gene ID and summing over exons
(tslen <- group_by(exons,tsid,gid) %>%
    dplyr::summarize(tslen=sum(exlen)))

# make sure each transcript just has one length
tslen<-as.data.frame(tslen)

#################################################
# Read in annotation data
#################################################

# Here we read in the table containing our gene annotations from EnTAP
# the sep, quote, and comment.char options ensure this complex table is parsed correctly
ann <- read.table("yourpath/Data Repo/final_annotations_lvl0.tsv",header=TRUE, sep="\t", quote="", comment.char="")
dim(ann)    


#let's look at how many variants there are for each query sequence
ann$geneid<- str_match(ann$Query.Sequence, "^(.*?\\..*?)\\.")[, 2] 
table(ann$geneid) #this tells me how many transcripts per gene id

# Use str_extract with a regex to get all text before the third period (excluding the third period)
ann<-ann %>% mutate(tsid=str_extract(Query.Sequence,"^(?:[^.]+\\.){2}[^.]+"),
                    tdscore=as.numeric(str_extract(Query.Sequence, "(?<=score=)[0-9.-]+")))

#get a list of the unique gene ids in the annotation 
(length(unique(ann$geneid))) #27586

# so from above, I can see that there are obviously multiple transcripts assoiated with a single gene.Each of these transcripts was 
# retained from transdecoder output because it was the isoform with the highest score. 
# Now what I want to do is to find and retain one single transcript for each gene that has the longest transcript length

#select transcript ID with longest exon

(longts<- group_by(tslen, gid) %>%
    dplyr::summarize(longts=tsid[which.max(tslen)],len=max(tslen)))


mylength<-longts


# filter the annotation table to just keep 1 transcript per gene (the longest transcript)
ann2<- ann %>% filter(tsid %in% longts$longts)

# Get my list of transcripts with > 1 annotation
newtable = group_by(ann2,tsid) %>% 
  dplyr::summarize(anncounts=sum(Predicted.Gene != "" | Description!=""| GO.Biological != "" |
                                   GO.Molecular != "" | GO.Cellular != ""))

# At this point, I should just have one line in my entap file that corresponds to a geneid in my count file

# merge the counts data frame with the annotations
tab<- full_join(co,ann2, join_by("gene_id"=="geneid"))

#i want to add in the transcript length of each transcript to my tab data table
tab = tab %>% mutate(tslen=longts$len[match(tab$gene_id, longts$gid)])

map(tab, ~ sum(is.na(.)))


# add a column to my table indicating whether a probe (MSTRG.XXXXX) is associated with a peptide ID'd from transdecoder, 
# and if so, whether it has an annotation or not
tab = tab[-c(1:5),] %>% mutate(peptide_status = case_when(is.na(Query.Sequence) == TRUE ~ "no peptide",
                                                          Query.Sequence !="<NA>" & Predicted.Gene == "" & Description=="" & GO.Biological == "" &
                                                            GO.Molecular == "" & GO.Cellular == ""~ "peptide but no annotation",
                                                          Query.Sequence !="<NA>" & Predicted.Gene != "" | Description!=""| GO.Biological != "" |
                                                            GO.Molecular != "" | GO.Cellular != ""~ "peptide with annotation"))


# make the missing data convention consistent across all columns (entap output uses "" whereas full_join gave a lot of NAs)

map(tab, ~ sum(is.na(.)))


tab <-tab %>%
  mutate_all(~ ifelse(is.na(.), "", .))


# Three modifications:
# 1. EnTAP flags possible contaminant transcripts, but leaves the field blank if
# it does not find a hit for the transcript, change the blank to "Unknown" for "Unknown"
# 2. Do the same for the Taxonomic Scope variable. 
# 3. Add a Prevalence column

tab <- mutate(tab,
              Contaminant = str_replace(Contaminant, "^$","Unknown"),
              Tax.Scope = str_replace(Tax.Scope, "^$","Unknown"),
              Prevalence = rowSums(tab[,2:52]>0)
)

