library(goseq)

###########GOterm analysis----
#need a vector thats 1 or 0 for whether each gene diff exp or not , add gene names to that
#a vector of transcript lengths (set all to 1 for now)
#get GO terms out of entap file 

###############################################
# gene ontology enrichment analysis on DE genes ----
###############################################
# Now, 'de_list' contains output dataframes with de (T/F) data for each dataframe in 'vs_res2_WvF'.
DESeq.2.GO<- function (de_gene,annotation_table, transcript_lengths ) {
  
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


# All


de_list_up_all <- list(); de_list_down_all <- list()

# Assuming you have a list of dataframes named 'vs_results'

for (name in names(all_results)) {
  i <- all_results[[name]]
  
  # Apply the condition to filter rows
  de_up <- i$padj < 0.05 & i$log2FoldChange > 0
  de_down <- i$padj < 0.05 & i$log2FoldChange < 0
  
  # Rename the rows based on the 'gene_id' column
  names(de_up) <- i$gene_id
  names(de_down) <- i$gene_id
  
  # Store the 'de' dataframe in the 'de_list' with the name of the original dataframe
  de_list_up_all[[paste("de_up", name, sep = "")]] <- de_up
  de_list_down_all[[paste("de_down", name, sep = "")]] <- de_down
}

# Initialize the list to store the results
up_go_de_out_all <- list(); down_go_de_out_all <- list()

# Loop through the index of de_list
for (i in seq_along(de_list_up_all)) {
  # Get the dataframe from de_list
  de_gene <- de_list_up_all[[i]]
  
  # Get the name of the dataframe
  name <- names(de_list_up_all)[i]
  
  # Apply the DESeq.2.GO function and store the result
  up_go_de_out_all[[i]] <- DESeq.2.GO(de_gene, go_map, transcript_lengths)
  names(up_go_de_out_all)[i] <- name
}

# Loop through the index of de_list
for (i in seq_along(de_list_down_all)) {
  # Get the dataframe from de_list
  de_gene <- de_list_down_all[[i]]
  table(de_gene)
  # Get the name of the dataframe
  name <- names(de_list_down_all)[i]
  
  # Apply the DESeq.2.GO function and store the result
  down_go_de_out_all[[i]] <- DESeq.2.GO(de_gene, go_map, transcript_lengths)
  names(down_go_de_out_all)[i] <- name
}


#enriched GO terms in Upregulated gene sets 
up_go_de_out_all$de_upvs_v_sp_res %>% filter(padj_overrepresented<0.05   &
                                               (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
#there are a few immune prcesses upregulated in the ventral skin vs. the spleen

up_go_de_out_all$de_uptiss_int_winter_v_fall_res %>% filter(padj_overrepresented<0.05   &
                                                              (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
up_go_de_out_all$de_uptiss_int_emergence_v_winter_res %>% filter(padj_overrepresented<0.05   &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
up_go_de_out_all$de_uptiss_int_spring_v_emergence_res %>% filter(padj_overrepresented<0.05   &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 


#Depleted GO terms in upregulated gene sets
up_go_de_out_all$de_upvs_v_sp_res %>% filter(padj_underrepresented<0.05   &
                                               (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 

#way more depleted GO terms in the genes upregulated in vs vs spleen (111)

up_go_de_out_all$de_uptiss_int_winter_v_fall_res %>% filter(padj_underrepresented<0.05   &
                                                              (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
up_go_de_out_all$de_uptiss_int_emergence_v_winter_res %>% filter(padj_underrepresented<0.05   &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
up_go_de_out_all$de_uptiss_int_spring_v_emergence_res %>% filter(padj_underrepresented<0.05   &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 



#Enriched GO terms in downregulated gene sets 
down_go_de_out_all$de_downvs_v_sp_res %>% filter(padj_overrepresented<0.05   &
                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 

#among DOWN regulated genes, there are a lot of enriched Immune go terms in the skin compared to the spleen

up_go_de_out_all$de_uptiss_int_spring_v_emergence_res %>% filter(padj_overrepresented<0.05   &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 


down_go_de_out_all$de_downtiss_int_winter_v_fall_res %>% filter(padj_overrepresented<0.05   &
                                                                  (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) 
# In the downregulated genes in the ventralskin vs spleen, a good number of interaction winter v fall. Map these!


# Ventral skin ----
#make lists of the up and down regulated genes from each analysis
de_list_up_vs <- list(); de_list_down_vs <- list()

# Assuming you have a list of dataframes named 'vs_results'

for (name in names(vs_results)) {
  i <- vs_results[[name]]
  
  # Apply the condition to filter rows
  de_up <- i$padj < 0.05 & i$log2FoldChange > 0
  de_down <- i$padj < 0.05 & i$log2FoldChange < 0
  
  # Rename the rows based on the 'gene_id' column
  names(de_up) <- i$gene_id
  names(de_down) <- i$gene_id
  
  # Store the 'de' dataframe in the 'de_list' with the name of the original dataframe
  de_list_up_vs[[paste("de_up", name, sep = "")]] <- de_up
  de_list_down_vs[[paste("de_down", name, sep = "")]] <- de_down
}

# # D - Transcript lengths
# 
# exp_genes <- rownames(dat_vs)
# str(exp_genes)
# length(exp_genes)
# 
# transcript_lengths<- longts[,c(1,3)] %>% filter(gid %in% exp_genes)
# # dim(transcript_lengths)



# Initialize the list to store the results
up_go_de_out_vs <- list(); down_go_de_out_vs <- list()

# Loop through the index of de_list
for (i in seq_along(de_list_up_vs)) {
  # Get the dataframe from de_list
  de_gene <- de_list_up_vs[[i]]
  
  # Get the name of the dataframe
  name <- names(de_list_up_vs)[i]
  
  # Apply the DESeq.2.GO function and store the result
  up_go_de_out_vs[[i]] <- DESeq.2.GO(de_gene, annotation_table, transcript_lengths)
  names(up_go_de_out_vs)[i] <- name
}

# Loop through the index of de_list
for (i in seq_along(de_list_down_vs)) {
  # Get the dataframe from de_list
  de_gene <- de_list_down_vs[[i]]
  
  # Get the name of the dataframe
  name <- names(de_list_down_vs)[i]
  
  # Apply the DESeq.2.GO function and store the result
  down_go_de_out_vs[[i]] <- DESeq.2.GO(de_gene, annotation_table, transcript_lengths)
  names(down_go_de_out_vs)[i] <- name
}


names(down_go_de_out_vs)

#Enriched and depeleted Immune go terms from upregulated gene sets
(de_upvs_res_WvF<-up_go_de_out_vs$de_upvs_res_WvF %>% filter(padj_overrepresented<0.05 | padj_underrepresented<0.05) %>%
  filter(category %in% immun_other_bp |
         category %in%adap_offspring_bp |
         category %in% inna_offspring_bp))


(de_upvs_res_PH1vW<-up_go_de_out_vs$de_upvs_res_PH1vW %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_upvs_res_PH1vPH2<-up_go_de_out_vs$de_upvs_res_PH1vPH2 %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                       (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_upvs_res_rhab<-up_go_de_out_vs$de_upvs_res_rhab %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                 (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_upvs_res_oswa<-up_go_de_out_vs$de_upvs_res_oswa %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                 (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp))) 

(de_upvs_res_sex<-up_go_de_out_vs$de_upvs_res_MvF%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                              (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )


#Enriched and depeleted Immune go terms from downregulated gene sets
(de_downvs_res_WvF<-down_go_de_out_vs$de_downvs_res_WvF %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                     (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )


(de_downvs_res_PH1vW<-down_go_de_out_vs$de_downvs_res_PH1vW %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                         (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )


(de_downvs_res_PH1vPH2<-down_go_de_out_vs$de_downvs_res_PH1vPH2%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                            (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_downvs_res_rhab<-down_go_de_out_vs$de_downvs_res_rhab%>% filter(padj_overrepresented<0.05 | padj_underrepresented<0.05))

(de_downvs_res_oswa<-down_go_de_out_vs$de_downvs_res_oswa%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                      (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_downvs_res_sex<-down_go_de_out_vs$de_downvs_res_MvF%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                    (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

#Put all the Enriched GO term results into a list
vs_go_sig_res_over<- list(de_upvs_res_WvF = de_upvs_res_WvF %>% filter(padj_overrepresented<0.05),
                          de_upvs_res_PH1vW=de_upvs_res_PH1vW%>% filter(padj_overrepresented<0.05),
                          de_upvs_res_PH1vPH2=de_upvs_res_PH1vPH2%>% filter(padj_overrepresented<0.05),
                          de_upvs_res_rhab=de_upvs_res_rhab%>% filter(padj_overrepresented<0.05),
                          de_upvs_res_oswa=de_upvs_res_oswa%>% filter(padj_overrepresented<0.05),
                          de_upvs_res_sex=de_upvs_res_sex%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_WvF=de_downvs_res_WvF%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_PH1vW=de_downvs_res_PH1vW%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_PH1vPH2=de_downvs_res_PH1vPH2%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_rhab=de_downvs_res_rhab%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_oswa=de_downvs_res_oswa%>% filter(padj_overrepresented<0.05),
                          de_downvs_res_sex=de_downvs_res_sex%>% filter(padj_overrepresented<0.05) )

# saveRDS(vs_go_sig_res_over, file = "yourpath/ventralskin_go_overrepresented.rds")

#Put all the depleted GO term results into a list
vs_go_sig_res_under<- list(de_upvs_res_WvF = de_upvs_res_WvF %>% filter(padj_underrepresented<0.05),
                           de_upvs_res_PH1vW=de_upvs_res_PH1vW%>% filter(padj_underrepresented<0.05),
                           de_upvs_res_PH1vPH2=de_upvs_res_PH1vPH2%>% filter(padj_underrepresented<0.05),
                           de_upvs_res_rhab=de_upvs_res_rhab%>% filter(padj_underrepresented<0.05),
                           de_upvs_res_oswa=de_upvs_res_oswa%>% filter(padj_underrepresented<0.05),
                           de_upvs_res_sex=de_upvs_res_sex%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_WvF=de_downvs_res_WvF%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_PH1vW=de_downvs_res_PH1vW%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_PH1vPH2=de_downvs_res_PH1vPH2%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_rhab=de_downvs_res_rhab%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_oswa=de_downvs_res_oswa%>% filter(padj_underrepresented<0.05),
                           de_downvs_res_sex=de_downvs_res_sex%>% filter(padj_underrepresented<0.05) )

# saveRDS(vs_go_sig_res_under, file = "yourpath/ventralskin_go_underrepresented.rds")

# Spleen----
de_list_up <- list(); de_list_down <- list()

# Assuming you have a list of dataframes named 'vs_res2_WvF'
for (name in names(sp_results)) {
  i <- sp_results[[name]]
  
  # Apply the condition to filter rows
  de_up <- i$padj < 0.05 & i$log2FoldChange > 0
  de_down <- i$padj < 0.05 & i$log2FoldChange < 0
  
  # Rename the rows based on the 'gene_id' column
  names(de_up) <- i$gene_id
  names(de_down) <- i$gene_id
  # Store the 'de' dataframe in the 'de_list' with the name of the original dataframe
  de_list_up[[paste("de_up", name, sep = "")]] <- de_up
  de_list_down[[paste("de_down", name, sep = "")]] <- de_down
}

# Initialize the list to store the results
up_go_de_out <- list(); down_go_de_out <- list()

# Loop through the index of de_list
for (i in seq_along(de_list_up)) {
  # Get the dataframe from de_list
  de_gene <- de_list_up[[i]]
  
  # Get the name of the dataframe
  name <- names(de_list_up)[i]
  
  # Apply the DESeq.2.GO function and store the result
  up_go_de_out[[i]] <- DESeq.2.GO(de_gene, annotation_table, transcript_lengths)
  names(up_go_de_out)[i] <- name
}

# Loop through the index of de_list
for (i in seq_along(de_list_down)) {
  # Get the dataframe from de_list
  de_gene <- de_list_down[[i]]
  
  # Get the name of the dataframe
  name <- names(de_list_down)[i]
  
  # Apply the DESeq.2.GO function and store the result
  down_go_de_out[[i]] <- DESeq.2.GO(de_gene, annotation_table, transcript_lengths)
  names(down_go_de_out)[i] <- name
}

#Enriched and depeleted Immune go terms from upregulated gene sets

(de_upsp_res_WvF<-up_go_de_out$de_upsp_res_WvF %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                            (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_upsp_res_PH1vW<-up_go_de_out$de_upsp_res_PH1vW %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05) & 
                                                                (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_upsp_res_PH1vPH2<-up_go_de_out$de_upsp_res_PH1vPH2 %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                    (category %in% immun_other_bp | category %in%adap_offspring_bp |category %in% inna_offspring_bp)))

(de_upsp_res_rhab<-up_go_de_out$de_upsp_res_rhab%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                             (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)))

(de_upsp_res_oswa<-up_go_de_out$de_upsp_res_oswa%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                             (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)))

(de_upsp_res_sex<-up_go_de_out$de_upsp_res_sex%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05) &
                                                           (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)))



#Enriched and depeleted Immune go terms from downregulated gene sets

(de_downsp_res_WvF<-down_go_de_out$de_downsp_res_WvF %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                  (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)))

(down_go_de_out$de_downsp_res_WvF %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                               (category %in% inna_offspring_bp)))

(de_downsp_res_PH1vW<-down_go_de_out$de_downsp_res_PH1vW %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                      (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )# go_de_out$de_sp_res_PH2vW %>% filter(padj_overrepresented<0.05 | padj_underrepresented<0.05)

(de_downsp_res_PH1vPH2<-down_go_de_out$de_downsp_res_PH1vPH2 %>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05)  &
                                                                          (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_downsp_res_rhab<-down_go_de_out$de_downsp_res_rhab%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05) &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_downsp_res_oswa<-down_go_de_out$de_downsp_res_oswa%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05) &
                                                                   (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )

(de_downsp_res_sex<-down_go_de_out$de_downsp_res_sex%>% filter((padj_overrepresented<0.05 | padj_underrepresented<0.05) &
                                                                 (category %in% immun_other_bp | category %in%adap_offspring_bp | category %in% inna_offspring_bp)) )


#Put all the enriched GO term results into a list
sp_go_sig_res_over<- list(de_upsp_res_WvF = de_upsp_res_WvF %>% filter(padj_overrepresented<0.05),
                          de_upsp_res_PH1vW=de_upsp_res_PH1vW%>% filter(padj_overrepresented<0.05),
                          de_upsp_res_PH1vPH2=de_upsp_res_PH1vPH2%>% filter(padj_overrepresented<0.05),
                          de_upsp_res_rhab=de_upsp_res_rhab%>% filter(padj_overrepresented<0.05),
                          de_upsp_res_oswa=de_upsp_res_oswa%>% filter(padj_overrepresented<0.05),
                          de_upsp_res_sex=de_upsp_res_sex%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_WvF=de_downsp_res_WvF%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_PH1vW=de_downsp_res_PH1vW%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_PH1vPH2=de_downsp_res_PH1vPH2%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_rhab=de_downsp_res_rhab%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_oswa=de_downsp_res_oswa%>% filter(padj_overrepresented<0.05),
                          de_downsp_res_sex=de_downsp_res_sex%>% filter(padj_overrepresented<0.05) )

# saveRDS(sp_go_sig_res_over, file = "yourpath/spleen_go_overrepresented.rds")

#Put all the depleted GO term results into a list
sp_go_sig_res_under<- list(de_upsp_res_WvF = de_upsp_res_WvF %>% filter(padj_underrepresented<0.05),
                           de_upsp_res_PH1vW=de_upsp_res_PH1vW%>% filter(padj_underrepresented<0.05),
                           de_upsp_res_PH1vPH2=de_upsp_res_PH1vPH2%>% filter(padj_underrepresented<0.05),
                           de_upsp_res_rhab=de_upsp_res_rhab%>% filter(padj_underrepresented<0.05),
                           de_upsp_res_oswa=de_upsp_res_oswa%>% filter(padj_underrepresented<0.05),
                           de_upsp_res_sex=de_upsp_res_sex%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_WvF=de_downsp_res_WvF%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_PH1vW=de_downsp_res_PH1vW%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_PH1vPH2=de_downsp_res_PH1vPH2%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_rhab=de_downsp_res_rhab%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_oswa=de_downsp_res_oswa%>% filter(padj_underrepresented<0.05),
                           de_downsp_res_sex=de_downsp_res_sex%>% filter(padj_underrepresented<0.05) )

# saveRDS(sp_go_sig_res_under, file = "yourpath/spleen_go_underrepresented.rds")


