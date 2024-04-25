# Automate the WGCNA to GOterm analysis pipeline----
require(goseq)
library(GO.db)
# Get inputs
#   1. List of genes in the module A  
#   2. Annotation table B (GOmap)
#   3. List of genes in the experiment C 
#   4. Transcript lengths D

# A - List of genes in module - this will be tissue-specific



# Ventral Skin
ventralskin_module_genelists

mods_list_vs <- list()

exp_genes_vs <- rownames(dat_vs)
str(exp_genes_vs)
length(exp_genes_vs)


for (name in names(ventralskin_module_genelists)) {
  i <- ventralskin_module_genelists[[name]]
  
  # Apply the condition to filter rows
  de <-  exp_genes_vs %in% i
  
  # Rename the rows based on the 'gene_id' column
  names(de) <- exp_genes_vs
  
  table(de)
  # Store the 'de' dataframe in the 'de_list' with the name of the original dataframe
  mods_list_vs[[paste("module", name, sep = "")]] <- de
  
}

# Spleen
spleen_module_genelists

exp_genes_sp <- rownames(dat_sp)
str(exp_genes_sp)
length(exp_genes_sp)


mods_list_sp <- list()

for (name in names(spleen_module_genelists)) {
  i <- spleen_module_genelists[[name]]
  
  # Apply the condition to filter rows
  de <-  exp_genes_sp %in% i
  
  # Rename the rows based on the 'gene_id' column
  names(de) <- exp_genes_sp
  
  table(de)
  # Store the 'de' dataframe in the 'de_list' with the name of the original dataframe
  mods_list_sp[[paste("module", name, sep = "")]] <- de
  
}

# B - Annotation table
Gomap <- dplyr::select(m,gene_id,GO.Biological,GO.Cellular,GO.Molecular) %>%
  pivot_longer(cols=!gene_id,names_to="GO.type",values_to="GO.term") %>%
  mutate(GO.term = str_split(GO.term,regex(".(?=GO:)"))) %>%   
  unnest(GO.term) %>% 
  separate(col=GO.term,into=c("GO.term","GO.description"),sep="(?<=GO:.......)-") %>%
  mutate(
    GO.term = na_if(GO.term, ""),
    GO.description = str_replace(GO.description,regex(",$"),"")
  )
head(Gomap)

# make the data frame that maps transcript IDs to GO terms. 
go_map <- data.frame(Gomap[,c(1,3)])



# C - List of genes in the experiment with nonzero read counts - this will be tissue specific
exp_genes_vs <- rownames(dat_vs)
exp_genes_sp <- rownames(dat_sp)


# D - Transcript lengths- this will be tissue specific
transcript_lengths_vs<- longts[,c(1,3)] %>% filter(gid %in% exp_genes_vs)
transcript_lengths_sp<- longts[,c(1,3)] %>% filter(gid %in% exp_genes_sp)


# Write function that does the following things
#  1. Look up the annotation from B associated with genes that are present in the module A > x
#  2. Assign each gene in input C a 1/0 based on whether it is present or absent in input A > y
#  3. Put y and transcript length (D) into the nullp function > z
#  4. Take output z and the annotations for "present in module genes" (x), use as input in the goseq function > desired output


# WGCNA.2.GO is a function that takes a list of genes (mod_gene) from a module output by WGCNA analysis, and 3 other inputs
# a list of background genes, a list of transcript lengths, and an annotation table, and runs the nullp and the goseq functions 
# with the goal of outputting a result showing which go terms are enriched or depleted in a particular module



# Now, 'de_list' contains output dataframes with de (T/F) data for each dataframe in 'vs_res2_WvF'.
WGCNA.2.GO<- function (de_gene,annotation_table, transcript_lengths ) {
  
  x <- go_map %>% filter(gene_id %in% rownames(data.frame(de_gene)))
  
  y <- de_gene
  table(y)
  # 
  len<-longts[,c(1,3)] %>% filter(gid %in% rownames(data.frame(de_gene)))
  
  z <- nullp(DEgenes = de_gene , bias.data = len$len )
  
  rownames(z) <- len$gid
  
  GO.output<- goseq(pwf=z,gene2cat=x)
  
  # do FDR correction on p-values using Benjamini-Hochberg, add to output object
  GO.output <- cbind(
    GO.output,
    padj_overrepresented=p.adjust(GO.output$over_represented_pvalue, method="BH"),
    padj_underrepresented=p.adjust(GO.output$under_represented_pvalue, method="BH")
  )
}



# make some lists of GO terms of interest
adap_offspring_bp = GOBPOFFSPRING[["GO:0002250"]]
inna_offspring_bp = GOBPOFFSPRING[["GO:0045087"]]
immun_offspring_bp = GOBPOFFSPRING[["GO:0002376"]]

adap<-which(immun_offspring_bp %in% adap_offspring_bp)
inna<-which(immun_offspring_bp %in% inna_offspring_bp)
adap_inna<- c(adap,inna)

immun_other_bp = immun_offspring_bp[-adap_inna]

# The trick is for me to get mod_gene to update by creating it in a loop that looks through the files I made with WGCNA
# So I want This loop will look at the list of files, for one file, assign it to the variable mod_gene, then pipe it into 
# the function WGCNA.2.GO. Then the output of the function will have the name of the module appended to it. ----
mod_out_vs<-list()


for (i in seq_along(mods_list_vs)) {
  # Get the dataframe from de_list
  de_gene <- mods_list_vs[[i]]
  
  # Get the name of the dataframe
  name <- names(mods_list_vs)[i]
  
  # Apply the DESeq.2.GO function and store the result
  mod_out_vs[[i]] <- WGCNA.2.GO(de_gene, annotation_table, transcript_lengths)
  names(mod_out_vs)[i] <- name
}


#look at the GO enrichment results for ventral skin----
# "MElightsteelblue1" 0.72 fall up
mod_out_vs$modulelightsteelblue %>% filter(padj_overrepresented <0.05)
mod_out_vs$modulelightsteelblue %>% filter(padj_overrepresented <0.05 & ontology == "BP") %>%filter(category %in%adap_offspring_bp |
                                                                                                   category %in% inna_offspring_bp|
                                                                                                   category %in% immun_other_bp )


# "MEturquoise"(0.92 fall up) 
mod_out_vs$moduleturquoise %>% filter(padj_overrepresented <0.05)
mod_out_vs$moduleturquoise %>% filter(padj_overrepresented <0.05 & ontology == "BP") %>% filter(category %in%adap_offspring_bp |
                                                                                                  category %in% inna_offspring_bp|
                                                                                                  category %in% immun_other_bp )

#"MEgreen"    0.71 down fall 0.73 up spring  
mod_out_vs$modulegreen %>% filter(padj_overrepresented <0.05)
mod_out_vs$modulegreen %>% filter(padj_overrepresented <0.05 & ontology == "BP") %>% filter(category %in%adap_offspring_bp |
                                                                                              category %in% inna_offspring_bp|
                                                                                              category %in% immun_other_bp )
# MEpink
mod_out_vs$modulepink %>% filter(padj_overrepresented <0.05)
mod_out_vs$modulepink%>% filter(padj_overrepresented <0.05 & ontology == "BP")%>% filter(category %in%adap_offspring_bp |
                                                                                           category %in% inna_offspring_bp|
                                                                                           category %in% immun_other_bp )
#   MEdarkolivegreen" 0.7 up spring
mod_out_vs$moduledarkolivegreen %>% filter(padj_overrepresented <0.05  & ontology == "BP") 
mod_out_vs$moduledarkolivegreen %>% filter(padj_overrepresented <0.05 & ontology == "BP") %>% filter(ontology == "BP") %>% filter(category %in%adap_offspring_bp |
                                                                                                                                    category %in% inna_offspring_bp|
                                                                                                                                    category %in% immun_other_bp )

#     "MEblue"  (-0.9 fall)  
mod_out_vs$moduleblue %>% filter(padj_overrepresented <0.05  & ontology == "BP") 
mod_out_vs$moduleblue %>% filter(padj_overrepresented <0.05 & ontology == "BP") %>% filter(ontology == "BP") %>% filter(category %in%adap_offspring_bp |
                                                                                                                          category %in% inna_offspring_bp|
                                                                                                                          category %in% immun_other_bp )

#  put WGCNA.2.GO into a loop for spleen----
mod_out_sp<-list()


for (i in seq_along(mods_list_sp)) {
  # Get the dataframe from de_list
  de_gene <- mods_list_sp[[i]]
  
  # Get the name of the dataframe
  name <- names(mods_list_sp)[i]
  
  # Apply the DESeq.2.GO function and store the result
  mod_out_sp[[i]] <- WGCNA.2.GO(de_gene, annotation_table, transcript_lengths)
  names(mod_out_sp)[i] <- name
}




#look at the GO enrichment results for spleen----
# yellow
mod_out_sp$moduleyellow  %>% filter(padj_overrepresented<0.05)
mod_out_sp$moduleyellow %>% filter(padj_overrepresented<0.05 & ontology == "BP") %>%filter(category %in%adap_offspring_bp |
                                                                                             category %in% inna_offspring_bp|
                                                                                             category %in% immun_other_bp )

# "MEturquoise" (0.73 fall, -0.71 winter)
mod_out_sp$moduleturquoise %>% filter(padj_overrepresented<0.05) 
mod_out_sp$moduleturquoise %>% filter(padj_overrepresented<0.05& ontology == "BP" ) %>%filter(category %in%adap_offspring_bp |
                                                                                                category %in% inna_offspring_bp|
                                                                                                category %in% immun_other_bp )

# "darkturquoise"(-0.83 winter)  
mod_out_sp$moduledarkturquoise  %>% filter(padj_overrepresented<0.05 )
mod_out_sp$moduledarkturquoise %>% filter(padj_overrepresented<0.05 & ontology == "BP")%>% filter(ontology == "BP")%>% filter(category %in%adap_offspring_bp |
                                                                                                                                category %in% inna_offspring_bp|
                                                                                                                                category %in% immun_other_bp )
# "pink" 
mod_out_sp$modulepink %>% filter(padj_overrepresented<0.05 )
mod_out_sp$modulepink %>% filter(padj_overrepresented<0.05 & ontology == "BP")%>% filter(ontology == "BP")%>% filter(category %in%adap_offspring_bp |
                                                                                                                       category %in% inna_offspring_bp|
                                                                                                                       category %in% immun_other_bp )

