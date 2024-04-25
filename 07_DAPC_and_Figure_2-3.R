#DAPC
#install.packages("adegenet") 
library(adegenet)
library(tidyr)
library(tibble)
library(ggpubr)
#*** Note: the code used for this section is taken and modified from the code published in: 
# Kenkel and Matz (2016).
# Gene expression plasticity as a mechanism of coral adaptation to a variable environment.
# Nature Ecology and Evolution 1(1). 
# 10.1038/s41559-016-0014
# https://www.nature.com/articles/s41559-016-0014



# All ----
########################################################################################################################
#  discriminant function analysis using adegenet package on vst transformed expr data----
########################################################################################################################
dat_all=as.data.frame(assay(all_vsd))
names(dat_all)
all_sampleTable = sampleTable_all 

head(dat_all) 
degs10_all<-rownames(dat_all)

#######################Some genes have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe=(dat_all)


# Initialize an empty vector to store row names
equal_min_max_rows <- c()

for (col in rownames(dframe)) {
  min_val <- min(dframe[col,])
  max_val <- max(dframe[col,])
  
  if (min_val == max_val) {
    print(col)
    equal_min_max_rows <- c(equal_min_max_rows, col)
  }
}

# Print the captured row names
print(equal_min_max_rows)

######################now use appropriate dataset for analysis

pcp_all=prcomp(t(dat_all[degs10_all,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp_all$x
screeplot(pcp_all,bstick=T) 

# adegenet: finding clusters (even though we know what clusters we want) - 
clus_all=adegenet::find.clusters(t(dat_all[degs10_all,]), n.clust = 8, n.pca = 40) #[degs10,]
names(dat_all[degs10_all,])


clus_all$grp=c(rep("Fall spleen",2), "Fall ventral skin", "Fall spleen",
               rep("Fall ventral skin",3),"Fall spleen","Fall ventral skin",
               rep("Winter spleen",2), "Winter ventral skin", "Winter spleen", "Winter ventral skin",
               rep("Emergence spleen",2), "Emergence ventral skin", "Emergence spleen", "Emergence ventral skin",
               "Emergence spleen", "Emergence ventral skin","Emergence spleen", "Emergence ventral skin", 
               "Spring spleen", "Spring ventral skin", rep("Spring spleen", 2),rep("Spring ventral skin", 2),
               "Spring spleen", "Spring ventral skin", "Spring ventral skin",  "Fall spleen", rep("Winter spleen", 2),
               "Spring spleen", "Fall ventral skin", rep("Winter ventral skin", 2), "Emergence ventral skin", "Winter ventral skin")


# now lets build a discriminant function for these 4 groups: I started with using all 40 PCs and 3 DF
dp_all=adegenet::dapc(t(dat_all[degs10_all,]),clus_all$grp) #[degs10,]

#now checking what is the optimal nmber of PCs to retain so I don't overfit
temp <- optim.a.score(dp_all)
names(temp)

#redo DAPC with 6 pcs and 2 dfs
dp_all=adegenet::dapc(t(dat_all[degs10_all,]),clus_all$grp) 

desired_levels <- c("Fall spleen", "Fall ventral skin", "Winter spleen", "Winter ventral skin",
                    "Emergence spleen", "Emergence ventral skin", "Spring spleen", "Spring ventral skin")
dp_all$grp<-factor_variable <- factor(dp_all$grp, levels = desired_levels)
dp_all$means <- dp_all$means[desired_levels,]
dp_all$grp.coord<-dp_all$grp.coord[desired_levels,]
# Create your plot here
par(mar=c(10,10,4,2) + 0.1, mgp=c(5,1,0), oma=c(2, 2, 2, 2))

#mult all coords by -1 so fall shows up on the left of the plot for visualizing
dp_all$ind.coord<-(dp_all$ind.coord*-1)
dp_all$grp.coord<-(dp_all$grp.coord*-1)

scatter(dp_all, bg="white", ratio.pca=0.9,  pch=c(20,18,20,18,20,18,20,18), cell=0, cstar=0,
        col = c("#925A44","#925A44", "#155289","#155289","#489FA7","#489FA7","#335A30","#335A30"),
        solid=1,  cex=4, clab=0, mstree=FALSE, scree.da=FALSE,
        posi.pca="bottomright", leg=FALSE)



# Ventral Skin ----
dat_vs=as.data.frame(assay(vs_vsd))
names(dat_vs)
vs_sampleTable = vs_sampleTable 


########################################################################################################################
#  discriminant function analysis using adegenet package on vst transformed expr data----
########################################################################################################################
degs10_vs<-rownames(dat_vs)

#######################Some genes may have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe_vs=(dat_vs)

# Initialize an empty vector to store row names
equal_min_max_rows <- c()

for (col in rownames(dframe_vs)) {
  min_val <- min(dframe_vs[col,])
  max_val <- max(dframe_vs[col,])
  
  if (min_val == max_val) {
    print(col)
    equal_min_max_rows <- c(equal_min_max_rows, col)
  }
}

# Print the captured row names
print(equal_min_max_rows)

######################now use appropriate dataset for analysis

pcp_vs=prcomp(t(dat_vs[degs10_vs,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores_vs=pcp_vs$x
screeplot(pcp_vs,bstick=T)

# adegenet: finding clusters (even though we know what clusters we want) - 
clus_vs=adegenet::find.clusters(t(dat_vs[degs10_vs,]), n.clust = 4, n.pca = 20) #[degs10,]
names(dat_vs[degs10_vs,])
clus_vs$grp

# define our a priori groups
clus_vs$grp=c(rep("Fall",5),rep("Winter",2), rep ("Emergence", 4),
              rep("Spring",5), "Fall", rep("Winter", 2),  "Emergence", "Winter") 

# now lets build a discriminant function for these 4 groups: I started with using all 22 PCs and 3 DFs
dp_vs=adegenet::dapc(t(dat_vs[degs10_vs,]),clus_vs$grp) #[degs10,]

#now checking what is the optimal nmber of PCs to retain so I don't overfit?
temp <- optim.a.score(dp_vs) #~1 but I'm going to use 3 anyway just so I can visualize
names(temp)

#redo DAPC with 3 pcs and 2 dfs
dp_vs=adegenet::dapc(t(dat_vs[degs10_vs,]),clus_vs$grp) 

desired_levels <- c("Fall", "Winter", "Emergence", "Spring")
dp_vs$grp<-factor_variable <- factor(dp_vs$grp, levels = desired_levels)
dp_vs$means <- dp_vs$means[desired_levels,]
dp_vs$grp.coord<-dp_vs$grp.coord[desired_levels,]

# Recreate the named numeric vector for dp_vs$prior
prior <- c((6/21), rep((5/21),3))
# Set the names as attributes
attr(prior, "names") <- c("Fall", "Winter", "Emergence", "Spring")
# Print the resulting named numeric vector
print(prior)
#reassign dp_vs$prior to be the new vector prior 
dp_vs$prior <- prior
dp_vs$posterior<-dp_vs$posterior[,desired_levels]
dp_vs$assign<-factor(dp_vs$assign, levels= desired_levels)

#mult all coords by -1 so fall shows up on the left of the plot for visualizing
dp_vs$ind.coord<-(dp_vs$ind.coord*-1)
dp_vs$grp.coord<-(dp_vs$grp.coord*-1)

# Create your plot here
par(mar=c(10,10,4,2) + 0.1, mgp=c(5,1,0), oma=c(2, 2, 2, 2))

# Figure 2A----
# 2 axes
scatter(dp_vs, bg="white", ratio.pca=0.5,  pch=20, cell=0, cstar=0,
        col = c("#925A44", "#155289","#489FA7","#335A30"),
        solid=1,  cex=4, clab=0, mstree=FALSE, scree.da=FALSE,
        posi.pca="bottomright", leg=FALSE)

# DAPC 1
scatter(dp_vs,1,1, ratio.pca=0.3, 
                   col= c("#925A44", "#155289","#489FA7","#335A30"), 
                   bg="white", scree.da=FALSE, legend=F, solid=.8)

# DAPC 2
scatter(dp_vs,2,2, ratio.pca=0.3,
                   col= c("#925A44", "#155289","#489FA7","#335A30"),
                   bg="white", scree.da=FALSE, legend=F, solid=.8)





# Spleen ----
dat_sp=as.data.frame(assay(sp_vsd))
names(dat_sp)
sp_sampleTable = sp_sampleTable 

########################################################################################################################
#  discriminant function analysis using adegenet package on vst transformed expr data----
########################################################################################################################
head(dat_sp) 
degs10_sp<-rownames(dat_sp)

#######################Some genes may have insufficient variance for DFA in data subsets!...find those genes using loop below
dframe_sp=(dat_sp)
for(col in rownames(dframe_sp))
{  min=min(dframe_sp[col,])
max=max(dframe_sp[col,])
if(min == max){print(col)}
}

######################now use appropriate dataset for analysis

pcp_sp=prcomp(t(dat_sp[degs10_sp,]), retx=TRUE, center=TRUE, scale.=TRUE) #scale.=TRUE
scores=pcp_sp$x
screeplot(pcp_sp,bstick=T) # only 1st PC is non-random when using diff exp genes only; higher PCs non-random when using all data...

# adegenet: finding clusters (even though we know what clusters we want) - 
clus_sp=adegenet::find.clusters(t(dat_sp[degs10_sp,]), n.clust = 4, n.pca = 20)

#Use clus$grp to rename 
clus_sp$grp=c(rep("Fall",4),rep("Winter",3), rep ("Post-winter 1", 5),
              rep("Post-winter 2",4),"Fall", rep("Winter", 2),  "Post-winter 2") #tell the DF groups you want to cluster; in this case in2in and off2off

# now lets build a discriminant function for these 4 groups: I started with using all 20 PCs and 3 DFs
dp_sp=adegenet::dapc(t(dat_sp[degs10_sp,]),clus_sp$grp) #[degs10,]

#now checking what is the optimal nmber of PCs to retain so I don't overfit?
temp <- optim.a.score(dp_sp) #~3
names(temp)

#redo DAPC with 3 pcs and 2 dfs
dp_sp=adegenet::dapc(t(dat_sp[degs10_sp,]),clus_sp$grp) 

desired_levels <- c("Fall", "Winter", "Post-winter 1", "Post-winter 2")
dp_sp$grp<-factor_variable <- factor(dp_sp$grp, levels = desired_levels)
dp_sp$means <- dp_sp$means[desired_levels,]
dp_sp$grp.coord<-dp_sp$grp.coord[desired_levels,]
# Recreate the named numeric vector for dp_sp$prior
prior <- c(rep(0.25, 4))
# Set the names as attributes
attr(prior, "names") <- c("Fall", "Winter", "Post-winter 1", "Post-winter 2")
# Print the resulting named numeric vector
print(prior)
#reassign dp_sp$prior to be the new vector prior 
dp_sp$prior <- prior
dp_sp$posterior<-dp_sp$posterior[,desired_levels]
dp_sp$assign<-factor(dp_sp$assign, levels= desired_levels)

# Create your plot here
par(mar=c(10,10,4,2) + 0.1, mgp=c(5,1,0), oma=c(2, 2, 2, 2))

# Figure 2B----
#DAPC 2 axes
scatter(dp_sp, bg="white", ratio.pca=0.4,   cell=0, cstar=0,
        col = c("#925A44", "#155289","#489FA7","#335A30"),
        solid=.9,  cex=7, clab=0, mstree=FALSE, scree.da=FALSE,
        posi.pca="bottomright", leg=FALSE,)
# txt.leg=paste("Sample Period",1:4)


#DAPC 1
scatter(dp_sp,1,1, ratio.pca=0.3, 
                   col= c("#925A44", "#155289","#489FA7","#335A30"), 
                   bg="white", scree.da=FALSE, legend=F, solid=.8)

#DAPC 2
scatter(dp_sp,2,2, ratio.pca=0.3,
                   col= c("#925A44", "#155289","#489FA7","#335A30"),
                   bg="white", scree.da=FALSE, legend=F, solid=.8)


#Highload genes Figure 3 ----
#look at which genes contributed most to a particular discriminant function
highload_genes_sp<- data.frame(dp_sp$var.contr) %>% filter(LD1 > 0.001 | LD2 > 0.002)

top_ld1_sp <- data.frame(dp_sp$var.contr) %>%
  arrange(desc(LD1)) %>% 
  slice_head(n = 5)%>%
  mutate(tissue = "spleen", DF = "LD1")

top_ld2_sp <- data.frame(dp_sp$var.contr) %>%
  arrange(desc(LD2)) %>%
  slice_head(n = 5)%>%
  mutate(tissue = "spleen", DF = "LD2")


top_ld1_vs <- data.frame(dp_vs$var.contr) %>%
  arrange(desc(LD1)) %>% 
  slice_head(n = 5)%>%
  mutate(tissue = "ventralskin",DF = "LD1")


top_ld2_vs <- data.frame(dp_vs$var.contr) %>%
  arrange(desc(LD2)) %>%
  slice_head(n = 5)%>%
  mutate(tissue = "ventralskin", DF = "LD2")

highloads <- bind_rows(top_ld1_sp, top_ld2_sp,top_ld1_vs,top_ld2_vs )
highloads <- highloads %>% mutate(gene = m$Predicted.Gene[match(rownames(highloads), m$gene_id)])


# Now result_dataset contains the top 5 observations from LD1 and the top 5 observations from LD2

highloads<-highloads %>%
  mutate(Predicted_gene = m$Predicted.Gene[match(rownames(highloads), m$gene_id)],
         peptide_status = m$peptide_status[match(rownames(highloads), m$gene_id)],
         Description = m$Description[match(rownames(highloads), m$gene_id)],
         GOBP = m$GO.Biological[match(rownames(highloads), m$gene_id)],
         GOMF = m$GO.Molecular[match(rownames(highloads), m$gene_id)],
         GO.CC = m$GO.Cellular[match(rownames(highloads), m$gene_id)],
         imm = case_when(str_detect(GOBP, "imm") ~ 1, TRUE~ 0| case_when(str_detect(GOMF, "imm") ~1, TRUE ~ 0|
                                                                           case_when(str_detect(GO.CC, "imm") ~1, TRUE ~ 0))))

highloads %>% filter(tissue =="spleen", DF == "LD1") %>% dplyr::select(Predicted_gene, peptide_status, Description, GOBP) %>% head(2)
highloads %>% filter(tissue =="spleen", DF == "LD2") %>% dplyr::select(Predicted_gene, peptide_status, Description)
highloads %>% filter(tissue =="ventralskin", DF == "LD1") %>% dplyr::select(Predicted_gene, peptide_status, Description)
highloads %>% filter(tissue =="ventralskin", DF == "LD2") %>% dplyr::select(Predicted_gene, peptide_status, Description)


highloads_sp<-highloads %>% filter(tissue =="spleen")
highloads_vs<-highloads %>% filter(tissue =="ventralskin")


highload_sp_df<-sp_normalized_counts %>%as.data.frame() %>% filter(rownames(sp_normalized_counts)%in%rownames(highloads_sp))
highload_sp_df<-data.frame(t(highload_sp_df))

highload_sp_df=highload_sp_df %>% rownames_to_column() %>% pivot_longer(cols = starts_with("MSTRG"),names_to="gene_id") %>%
  mutate(period = substr(rowname, 4,5)) %>%
  mutate(period = case_when(period == "01" ~ "Fall",
                            period == "03" ~ "Winter",
                            period == "04" ~ "Emergence",
                            period == "06" ~ "Spring"))

# Define the desired order of the levels
desired_order <- c("Fall", "Winter", "Emergence", "Spring")

# Use factor with levels to reorder the factor variable
highload_sp_df$period <- factor(highload_sp_df$period, levels = desired_order)

highload_sp_df=highload_sp_df %>% 
  mutate(predicted_gene = m$Predicted.Gene[match(highload_sp_df$gene_id, m$gene_id)])%>%
  mutate(predicted_gene= case_when(gene_id == "MSTRG.23145"~"MSTRG.23145",
                                   gene_id == "MSTRG.541"~"MSTRG.541", 
                                   TRUE ~ predicted_gene))%>%
  mutate(tissue = "Spleen",
         DF = highloads$DF[match(highload_sp_df$gene_id, rownames(highloads))])

highload_vs_df<-vs_normalized_counts %>%as.data.frame() %>% filter(rownames(vs_normalized_counts)%in%rownames(highloads_vs))
highload_vs_df<-data.frame(t(highload_vs_df))

highload_vs_df=highload_vs_df %>% rownames_to_column() %>% pivot_longer(cols = starts_with("MSTRG"),names_to="gene_id") %>%
  mutate(period = substr(rowname, 4,5)) %>%
  mutate(period = case_when(period == "01" ~ "Fall",
                            period == "03" ~ "Winter",
                            period == "04" ~ "Emergence",
                            period == "06" ~ "Spring"))

highload_vs_df=highload_vs_df %>% 
  mutate(predicted_gene = m$Predicted.Gene[match(highload_vs_df$gene_id, m$gene_id)])%>%
  mutate(predicted_gene= case_when(gene_id == "MSTRG.34173"~"ranaspumin-like \nprotein",
                                   TRUE ~ predicted_gene)) %>%
  mutate(tissue = "Ventral Skin",
         DF = highloads$DF[match(highload_vs_df$gene_id, rownames(highloads))])

highload_vs_df$period <- factor(highload_vs_df$period, levels = desired_order)


(HL_sp<-ggplot(highload_sp_df %>% filter(gene_id == "MSTRG.541" |gene_id == "MSTRG.23145"|
                                           gene_id == "MSTRG.36218"|gene_id == "MSTRG.38809") , #|gene_id =="MSTRG.148"
               aes(x=period, y=log(value), fill = period))+
    scale_fill_manual(values = c("#925A44","#155289", "#489FA7","#335A30"))+
    geom_jitter(alpha = 0.4)+
    geom_boxplot()+
    facet_grid(DF+predicted_gene  ~  ., scales = "free") +
    theme_bw()+
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_text(size = 13),
          axis.title.y  = element_blank(),
          plot.title = element_text(size = 18),
          axis.text.x = element_text(size = 13))+
    labs(x = NULL, y = "log(Normalized count)", fill= "Sample Period", title = "Spleen")+
    guides(fill="none"))

distinct(highload_vs_df, gene_id, predicted_gene)
distinct(highload_sp_df, gene_id, predicted_gene)

vs_gene_order<-c("ranaspumin-like \nprotein","SLPI","IRG1", "S100A12" )
# Use factor with levels to reorder the factor variable
highload_vs_df$predicted_gene <- factor(highload_vs_df$predicted_gene, levels = vs_gene_order)

highload_vs_df<-highload_vs_df %>% mutate(period= factor(period, levels = c("Fall", "Winter", "Emergence", "Spring")))


(HL_vs<-ggplot(highload_vs_df %>% filter(gene_id == "MSTRG.34173"|gene_id == "MSTRG.15881"|
                                           gene_id == "MSTRG.35005"|gene_id == "MSTRG.47854") , #|gene_id =="MSTRG.148"
               aes(x=period, y=log(value), fill = period))+
    scale_fill_manual(values = c("#925A44","#155289", "#489FA7","#335A30"))+
    geom_jitter(alpha = 0.4)+
    geom_boxplot()+
    facet_grid(predicted_gene ~  ., scales = "free") +
    theme_bw()+
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(),
          strip.text = element_text(size = 13),
          axis.title.y  = element_text(size = 13),
          plot.title = element_text(size = 18),
          axis.text.x = element_text(size = 13))+
    labs(x = NULL, y = "log(Normalized count)", fill= "Sample Period",title = "Ventral Skin")+
    guides(fill="none")
)


hl_boxplots<-ggarrange(HL_vs, HL_sp, 
                       ncol = 2, nrow = 1)

hl_boxplots



