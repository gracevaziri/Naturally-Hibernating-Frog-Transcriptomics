#################################################
# create data frame 'm' where columns are counts, rownames are transcript IDs
#################################################
m <- tab %>% mutate_at(vars(2:52), as.integer)
rownames(m) <- m[,1]

#add in probe/gene IDs
m = m %>% mutate(gid = str_match(m$Query.Sequence, "^(.*?\\..*?)\\.")[, 2])

# also subset transcript lengths
mylength <- mylength %>% filter(gid %in% m$gid)

#################################################
# Preliminary data exploration and filtering
#################################################

# how many transcripts did we quantify expression for?
dim(m)

# how many (non-contaminant) transcripts have expression for n samples?
prevalence<- data.frame(tab %>% dplyr::select(Prevalence))
table(prevalence)


# what proportion of expression data maps to low prevalence transcripts?
colSums(m[prevalence<2,c(2:52)])/colSums(m[,c(2:52)]) 

# what is the distribution of expression across transcripts?
rowSums(m[,c(2:52)]) %>% log(.,10) %>% hist(.,100)

# Before moving forward, we're going to exclude genes that are most likely to be useless in the analysis
# let's throw out all the transcripts with
# a) less than 20 fragments mapping in total across all samples
# b) that have non-zero expression for only a single sample

subg <- prevalence > 1 & rowSums(m[,c(2:52)]) > 20
table(subg)

m <- m[subg,]
mylength <- mylength$gid[subg]

m <- m[!grepl("^NA", rownames(m)), , drop = FALSE]

m


