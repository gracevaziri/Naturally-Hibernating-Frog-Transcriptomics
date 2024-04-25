#get normalized counts for plotting----
vs_normalized_counts<- as.data.frame(counts(vs_ddsClean, normalized=TRUE))
sp_normalized_counts<- as.data.frame(counts(sp_ddsClean, normalized=TRUE))
all_normalized_counts<- as.data.frame(counts(all_ddsClean, normalized=TRUE))

# DGE ----

# all----
resultsNames(all_ddsClean)
# vs vs spleen, main effect
vs_v_sp <-results(all_ddsClean, name="tissue_vs_vs_sp", alpha = 0.05)
vs_WvF_shrink <- lfcShrink(all_ddsClean,type="ashr",coef="tissue_vs_vs_sp")

# the tissue vs. period interaction. What is the differnce in how tissues respond to the fall v winter contrast?
summary(results(all_dds, name="tissuevs.period3"), alpha = 0.05)
tiss_int_winter_v_fall<- results(all_dds, name="tissuevs.period3", alpha = 0.05)


# the tissue vs. period interaction. What is the differnce in how tissues respond to the  winter v emergence contrast?
summary(results(all_dds, contrast=list("tissuevs.period4", "tissuevs.period3")))
tiss_int_emergence_v_winter<- results(all_dds, contrast=list("tissuevs.period4", "tissuevs.period3"), alpha = 0.05)

# the tissue vs. period interaction. What is the differnce in how tissues respond to the emergence v spring contrast?
summary(results(all_dds, contrast=list("tissuevs.period6", "tissuevs.period4")))
tiss_int_spring_v_emergence<- results(all_dds, contrast=list("tissuevs.period6", "tissuevs.period4"), alpha = 0.05)

# save results into one big object
all_results <- list(
  
  vs_v_sp_res = cbind(gene_id=rownames(vs_v_sp),vs_v_sp) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  
  tiss_int_winter_v_fall_res = cbind(gene_id=rownames(tiss_int_winter_v_fall),tiss_int_winter_v_fall) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  tiss_int_emergence_v_winter_res = cbind(gene_id=rownames(tiss_int_emergence_v_winter),tiss_int_emergence_v_winter) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  tiss_int_spring_v_emergence_res = cbind(gene_id=rownames(tiss_int_spring_v_emergence),tiss_int_spring_v_emergence) %>% 
    data.frame() %>%
    filter(!is.na(padj))
)

# save the results as R object
# saveRDS(all_results, file = "yourpath/all_de_results.rds")

# ventral skin ----

# winter versus fall, main effect
vs_WvF <-results(vs_ddsClean, name="period_3_vs_1", alpha = 0.05)
vs_WvF_shrink <- lfcShrink(vs_ddsClean,type="ashr",coef="period_3_vs_1")

# post-winter1 versus fall
vs_PH1vF <- results(vs_ddsClean, name="period_4_vs_1", alpha = 0.05)
vs_PH1vF_shrink <- lfcShrink(vs_ddsClean,type="ashr",coef="period_4_vs_1")

# post-winter2 versus fall
vs_PH2vF <- results(vs_ddsClean, name="period_6_vs_1", alpha = 0.05)
vs_PH2vF_shrink <- lfcShrink(vs_ddsClean,type="ashr",coef="period_6_vs_1")

# post-winter1 versus winter
vs_PH1vW <-results(vs_ddsClean, contrast=list(c("period_4_vs_1"), c("period_3_vs_1")), alpha = 0.05)
vs_PH1vW_shrink <- lfcShrink(vs_ddsClean, contrast=list(c("period_4_vs_1"), c("period_3_vs_1")),type="ashr")

# post-winter2 versus winter
vs_PH2vW<-results(vs_ddsClean, contrast=list(c("period_6_vs_1"), c("period_3_vs_1")), alpha = 0.05)
vs_PH2vW_shrink <- lfcShrink(vs_ddsClean,contrast=list(c("period_6_vs_1"), c("period_3_vs_1")), type="ashr")

# post-winter2 versus post-winter1
vs_PH1vPH2 <-results(vs_ddsClean, contrast=list(c("period_6_vs_1"),c("period_4_vs_1")), alpha = 0.05)
vs_PH1vPH2_shrink <- lfcShrink(vs_ddsClean,contrast=list(c("period_6_vs_1"),c("period_4_vs_1")), type="ashr")

# male versus female
vs_MvF <- results(vs_ddsClean, name= "sex_M_vs_F", alpha = 0.05)
vs_MvF_shrink<-lfcShrink(vs_ddsClean,contrast=list(c("sex_M_vs_F")), type ="ashr")

# rhabdias
vs_rhab <- results(vs_ddsClean, name= "rhabdias_log", alpha = 0.05)
vs_rhab_shrink<-lfcShrink(vs_ddsClean,contrast=list(c("rhabdias_log")), type ="ashr")

# oswaldocruzia
vs_oswa <- results(vs_ddsClean, name= "oswaldocruzia_log", alpha = 0.05)
vs_oswa_shrink<-lfcShrink(vs_ddsClean,contrast=list(c("oswaldocruzia_log")), type ="ashr")

# save results into one object
vs_results <- list(
  vs_res_WvF = cbind(gene_id=rownames(vs_WvF),vs_WvF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_PH1vF = cbind(gene_id=rownames(vs_PH1vF),vs_PH1vF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_PH2vF = cbind(gene_id=rownames(vs_PH2vF),vs_PH2vF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_PH1vW = cbind(gene_id=rownames(vs_PH1vW),vs_PH1vW) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_PH2vW = cbind(gene_id=rownames(vs_PH2vW),vs_PH2vW) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_PH1vPH2 = cbind(gene_id=rownames(vs_PH1vPH2),vs_PH1vPH2) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  vs_res_rhab = cbind(gene_id = rownames(vs_rhab), vs_rhab) %>%
    data.frame() %>%
    filter(!is.na(padj)), 
  
  vs_res_oswa = cbind(gene_id = rownames(vs_oswa), vs_oswa) %>%
    data.frame() %>%
    filter(!is.na(padj)), 
  
  vs_res_MvF = cbind(gene_id = rownames(vs_MvF), vs_MvF) %>%
    data.frame() %>%
    filter(!is.na(padj))
  
)

#save the results as R object
# saveRDS(vs_results, file = "yourpath/skin_de_results.rds")

# spleen ----

# winter versus fall
sp_WvF <-results(sp_ddsClean, name="period_3_vs_1", alpha = 0.05)
sp_WvF_shrink <- lfcShrink(sp_dds,type="ashr",coef="period_3_vs_1")

# post-winter1 versus fall
sp_PH1vF <- results(sp_ddsClean, name="period_4_vs_1", alpha = 0.05)
sp_PH1vF_shrink <- lfcShrink(sp_ddsClean,type="ashr",coef="period_4_vs_1")

# post-winter2 versus fall
sp_PH2vF <- results(sp_ddsClean, name="period_6_vs_1", alpha = 0.05)
sp_PH2vF_shrink <- lfcShrink(sp_ddsClean,type="ashr",coef="period_6_vs_1")

# post-winter1 versus winter
sp_PH1vW <-results(sp_ddsClean, contrast=list(c("period_4_vs_1"), c("period_3_vs_1")), alpha = 0.05)
sp_PH1vW_shrink <- lfcShrink(sp_ddsClean, contrast=list(c("period_6_vs_1"), c("period_3_vs_1")),type="ashr")

# post-winter2 versus winter
sp_PH2vW<-results(sp_ddsClean, contrast=list(c("period_6_vs_1"), c("period_3_vs_1")), alpha = 0.05)
sp_PH2vW_shrink <- lfcShrink(sp_ddsClean,contrast=list(c("period_6_vs_1"), c("period_3_vs_1")), type="ashr")

# post-winter2 versus post-winter1
sp_PH1vPH2 <-results(sp_ddsClean, contrast=list(c("period_6_vs_1"),c("period_4_vs_1")), alpha = 0.05)
sp_PH1vPH2_shrink <- lfcShrink(sp_ddsClean,contrast=list(c("period_6_vs_1"),c("period_4_vs_1")), type="ashr")

# male versus female
sp_MvF <- results(sp_ddsClean, name= "sex_M_vs_F", alpha = 0.05)
sp_MvF_shrink <- lfcShrink(sp_ddsClean,contrast=list(c("sex_M_vs_F")), type="ashr")

# rhabdias
sp_rhab <- results(sp_ddsClean, name= "rhabdias_log", alpha = 0.05)
sp_rhab_shrink <- lfcShrink(sp_ddsClean,contrast=list(c("rhabdias_log")), type="ashr")

# oswaldocruzia
sp_oswa <- results(sp_ddsClean, name= "oswaldocruzia_log", alpha = 0.05)
sp_oswa_shrink <- lfcShrink(sp_ddsClean,contrast=list(c("oswaldocruzia_log")), type="ashr")

# save the results 
sp_results <- list(
  sp_res_WvF = cbind(gene_id=rownames(sp_WvF),sp_WvF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_PH1vF = cbind(gene_id=rownames(sp_PH1vF),sp_PH1vF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_PH2vF = cbind(gene_id=rownames(sp_PH2vF),sp_PH2vF) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_PH1vW = cbind(gene_id=rownames(sp_PH1vW),sp_PH1vW) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_PH2vW = cbind(gene_id=rownames(sp_PH2vW),sp_PH2vW) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_PH1vPH2 = cbind(gene_id=rownames(sp_PH1vPH2),sp_PH1vPH2) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_rhab = cbind(gene_id=rownames(sp_rhab),sp_rhab) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_oswa = cbind(gene_id=rownames(sp_oswa),sp_oswa) %>% 
    data.frame() %>%
    filter(!is.na(padj)),
  
  sp_res_sex = cbind(gene_id=rownames(sp_MvF),sp_MvF) %>% 
    data.frame() %>%
    filter(!is.na(padj))
  
)

#save the results as R object
# saveRDS(sp_results, file = "yourpath/spleen_de_results.rds")


