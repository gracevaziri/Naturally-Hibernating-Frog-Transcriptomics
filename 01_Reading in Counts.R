# # Load the libraries we'll need in the following code:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("readxl")
require(stringr)


######################################################
# Point R to count data, set an output file prefix ----
######################################################
# create an object with the directory containing your counts:
# !!edit this to point to your own count file directory!!


getwd()
setwd("yourpath/Data Repo/Counts/")
directory <- "yourpath/Data Repo/Counts/"
# create a empty dataframe called co to merge the data into
co <- data.frame()
dim(co)
# using for loop read all the count files in the count_dir path
# iterate over each "abundance.tsv" file
# merge each file, i, into a single data frame

count_files <- list.files(path = directory,recursive=FALSE,pattern = ".*counts",full.names=FALSE)

for (i in count_files) {
  
  # print the file that is being loaded
  print(paste0("reading file: ", i))
  
  # read file i as a data frame
  f <- read.table(i, sep = "\t", header = TRUE)[,c(1,2)]
  
  # extract the sample ID from the file name:
  sam <- str_extract(i,"trim_(.*?)(?=_S)")
  
  # rename the columns
  colnames(f) <- c("gene_id", sam)
  
  #if the counts object is empty just copy the f to m
  if(length(co) == 0){
    co <- f
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    co <- merge(co, f, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f)
}

#grab the rows from the 1st column and use it as the row-names in the dataframe
co
rownames(co) <- co[,1]
dim(co)

