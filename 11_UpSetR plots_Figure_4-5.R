# UpSetR----
# ****Note, the plots "vs_main" and "sp_main" were exported to a graphics editing software to be combined into Figure 4 of the paper. 
# ****Note, the plots "vs_main_go_over" and "sp_main_go_over" were exported to a graphics editing software to be combined into Figure 5 of the paper. 
#graphing the DE results

# Set up some vectors of gene annotations related to immune processes
immune_system_process_genes<-tab%>% dplyr::select(gene_id, GO.Biological) %>% filter(str_detect(GO.Biological, "GO:0002376") ) %>% 
  dplyr::select(gene_id)


adap_immune_genes<-tab%>% dplyr::select(gene_id, GO.Biological) %>% filter(str_detect(GO.Biological, "GO:0002250")) %>% 
  dplyr::select(gene_id)

inna_immune_genes<-tab%>% dplyr::select(gene_id, GO.Biological) %>% filter(str_detect(GO.Biological, "GO:0045087")) %>% 
  dplyr::select(gene_id)

`%notin%` <- function(x, table) !(x %in% table)


# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("krassowski/complex-upset")
library(ggplot2)
library(ComplexUpset)

# Figure 4----
#ventral skin ----
#main effects graph
vs_results<-readRDS(file = "youtpath/skin_de_results.rds")
# ---- preparing a matrix for my UpSetR plot
# https://research.libd.org/rstatsclub/2018/07/27/hacking-our-way-through-upsetr/

# Start with the rows being all the expressed genes in the ventral skin, 
# Columns represent each comparison, but I will actually have double the number of comparisons I'm including, because I want to 
# separate out by up regulated and down regulated

vs_results_up <- list(vs_up_WvF = vs_results$vs_res_WvF %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      vs_up_PH1vW =vs_results$vs_res_PH1vW %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      vs_up_PH1vPH2 =vs_results$vs_res_PH1vPH2 %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      vs_up_rhab =vs_results$vs_res_rhab %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      vs_up_oswa =vs_results$vs_res_oswa %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      vs_up_sex =vs_results$vs_res_MvF %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id)
                      )

vs_results_down <- list(vs_down_WvF = vs_results$vs_res_WvF %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        vs_down_PH1vW =vs_results$vs_res_PH1vW %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        vs_down_PH1vPH2 =vs_results$vs_res_PH1vPH2 %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        vs_down_rhab =vs_results$vs_res_rhab %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        vs_down_oswa =vs_results$vs_res_oswa %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        vs_down_sex =vs_results$vs_res_MvF %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id)
                        )

#Start making the UpSetR dataframe
upsetr_df<- data.frame(gene= rownames(vs_ddsClean))



upsetr_df <-upsetr_df %>% mutate(`Up in winter vs. fall` = case_when(gene %in% vs_results_up$vs_up_WvF$gene_id ~ 1,
                                                                     TRUE ~ 0),
                                 `Down in winter vs. fall` = case_when(gene %in% vs_results_down$vs_down_WvF$gene_id ~ 1,
                                                                       TRUE ~ 0),
                                 `Up in emergence vs. winter` = case_when(gene %in% vs_results_up$vs_up_PH1vW$gene_id ~ 1,
                                                                          TRUE ~ 0),
                                 `Down in emergence vs. winter` = case_when(gene %in% vs_results_down$vs_down_PH1vW$gene_id ~ 1,
                                                                            TRUE ~ 0),
                                 `Up in spring vs. emergence` = case_when(gene %in% vs_results_up$vs_up_PH1vPH2$gene_id ~ 1,
                                                                          TRUE ~ 0),
                                 `Down in spring vs. emergence` = case_when(gene %in% vs_results_down$vs_down_PH1vPH2$gene_id ~ 1,
                                                                            TRUE ~ 0),
                                 `Up in Rhabdias sp. infected vs. uninfected` = case_when(gene %in% vs_results_up$vs_up_rhab$gene_id ~ 1,
                                                                                          TRUE ~ 0),
                                 `Down in Rhabdias sp. infected vs. uninfected` = case_when(gene %in% vs_results_down$vs_down_rhab$gene_id ~ 1,
                                                                                            TRUE ~ 0),
                                 `Up in Oswaldocruzia sp. infected vs. uninfected` = case_when(gene %in% vs_results_up$vs_up_oswa$gene_id ~ 1,
                                                                                               TRUE ~ 0),
                                 `Down in Oswaldocruzia sp. infected vs. uninfected` = case_when(gene %in% vs_results_down$vs_down_owsa$gene_id ~ 1,
                                                                                                 TRUE ~ 0),
                                 `Up in male vs. female` = case_when(gene %in% vs_results_up$vs_up_sex$gene_id ~ 1,
                                                                     TRUE ~ 0),
                                 `Down in male vs. female` = case_when(gene %in% vs_results_down$vs_down_sex$gene_id ~ 1,
                                                                       TRUE ~ 0),
)



upsetr_df<-upsetr_df %>% mutate(`GO Term` = case_when(gene %in% adap_immune_genes$gene_id & gene %in% immune_system_process_genes$gene_id ~ "Adaptive",
                                                      
                                                      gene %in% inna_immune_genes$gene_id & gene %in% immune_system_process_genes$gene_id ~ "Innate",
                                                      
                                                      gene %in% immune_system_process_genes$gene_id &  
                                                        gene %notin%adap_immune_genes$gene_id & 
                                                        gene%notin%inna_immune_genes$gene_id~ "Other immune",
                                                      
                                                      TRUE ~ "Not immune related"))

upsetr_df$`GO Term` <- factor(upsetr_df$`GO Term`, levels = c("Adaptive", "Innate", "Other immune", "Not immune related"))


upsetr_df_up <- upsetr_df %>% dplyr::select(starts_with("Up"))
names(upsetr_df_up)<-substr(names(upsetr_df_up), 6, 70)

upsetr_df_down <- upsetr_df %>% dplyr::select(starts_with("Down"))
names(upsetr_df_down)<-substr(names(upsetr_df_down), 8, 80)


(vs_main<-ComplexUpset::upset(upsetr_df, name = "Comparison",intersect = names(upsetr_df)[c(2:13)],
                              width_ratio=0.3,height_ratio = 0.8, min_size=10, max_size = 5000, wrap = TRUE,
                              themes=upset_modify_themes(
                                list(
                                  'intersections_matrix'=theme(text=element_text(size=12)))),
                              queries=list(
                                
                                upset_query(set='Up in Rhabdias sp. infected vs. uninfected', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in Rhabdias sp. infected vs. uninfected', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in spring vs. emergence', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in spring vs. emergence', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in emergence vs. winter', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in emergence vs. winter', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in winter vs. fall', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in winter vs. fall', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set="Down in Oswaldocruzia sp. infected vs. uninfected", fill = "darkblue", only_components = c('intersections_matrix')),
                                upset_query(set="Up in Oswaldocruzia sp. infected vs. uninfected", fill = "red", only_components = c('intersections_matrix')),
                                upset_query(set='Up in male vs. female', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in male vs. female', fill="darkblue", only_components=c('intersections_matrix'))
                              ),
                              base_annotations=list(
                                'Intersection size'=(intersection_size(text=list(size=4),bar_number_threshold=1,  # show all numbers on top of bars
                                                                       width=0.6,   # reduce width of the bars
                                                                       mapping=aes(fill=`GO Term`)) 
                                                     + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                                     
                                                     # add some space on the top of the bars
                                                     + scale_y_continuous(expand=expansion(mult=c(0, 0.2)))
                                                     # +scale_y_break(breaks = c(150,1500), scales = "fixed")
                                                     
                                                     + ggtitle('Ventral Skin Differential Expression')
                                                     
                                                     + guides(fill = "none")
                                                     + theme(
                                                       # hide grid lines
                                                       panel.grid.major=element_blank(),
                                                       panel.grid.minor=element_blank(),
                                                       # show axis lines
                                                       axis.line=element_line(colour='black'),
                                                       axis.text = element_text(size = 12),
                                                       axis.title.y = element_text(size = 14),
                                                       plot.title = element_text(size = 16),
                                                       legend.title = element_text(size = 14),
                                                       legend.text = element_text(size = 14)
                                                       
                                                     )
                                )
                                # 'Intersection ratio'=intersection_ratio()
                                
                              ),
                              stripes=upset_stripes(
                                geom=geom_segment(size=7),  # make the stripes larger
                                colors=c('gray92', 'white')
                              ),
                              # to prevent connectors from getting the colorured
                              # use `fill` instead of `color`, together with `shape='circle filled'`
                              matrix=intersection_matrix(
                                geom=geom_point(
                                  shape='circle filled',
                                  size=3.5,
                                  stroke=0.45
                                )
                              ),
                              
                              set_sizes=(
                                upset_set_size(geom=geom_bar(aes(fill = `GO Term`),width=0.4),position = "right" )
                                + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                + geom_text(aes(label=after_stat(count)), hjust=0, stat='count')
                                + expand_limits(y=2200)
                                
                                + theme(
                                  axis.line.x=element_line(colour='black'),
                                  axis.ticks.x=element_line(),
                                  axis.text.x=element_text(angle=90, size =12),
                                  # axis.text = element_text(size = 14)
                                  
                                )
                                
                              ),
                              guides = 'over',
                              
                              sort_sets=FALSE,
                              sort_intersections='descending'
)
)

vs_main


# spleen ----
# separate out by up regulated and down regulated

sp_results<-readRDS(file = "../../../Final bioinformatics/DEgene results/spleen_de_results.rds")

sp_results_up <- list(sp_up_WvF = sp_results$sp_res_WvF %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      sp_up_PH1vW =sp_results$sp_res_PH1vW %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      sp_up_PH1vPH2 =sp_results$sp_res_PH1vPH2 %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      sp_up_rhab =sp_results$sp_res_rhab %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      sp_up_oswa =sp_results$sp_res_oswa %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id),
                      sp_up_sex =sp_results$sp_res_sex %>% as.data.frame() %>% 
                        filter(padj<0.05 & log2FoldChange >2) %>% dplyr::select(gene_id)
)

sp_results_down <- list(sp_down_WvF = sp_results$sp_res_WvF %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        sp_down_PH1vW =sp_results$sp_res_PH1vW %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                       sp_down_PH1vPH2 =sp_results$sp_res_PH1vPH2 %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        sp_down_rhab =sp_results$sp_res_rhab %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        sp_down_oswa =sp_results$sp_res_oswa %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id),
                        sp_down_sex =sp_results$sp_res_sex %>% as.data.frame() %>% 
                          filter(padj<0.05 & log2FoldChange <(-2)) %>% dplyr::select(gene_id)
                       )


upsetr_df_sp<- data.frame(gene= rownames(sp_ddsClean))


upsetr_df_sp <-upsetr_df_sp %>% mutate(`Up in winter vs. fall` = case_when(gene %in% sp_results_up$sp_up_WvF$gene_id ~ 1,
                                                                           TRUE ~ 0),
                                       `Down in winter vs. fall` = case_when(gene %in% sp_results_down$sp_down_WvF$gene_id ~ 1,
                                                                             TRUE ~ 0),
                                       `Up in emergence vs. winter` = case_when(gene %in% sp_results_up$sp_up_PH1vW$gene_id ~ 1,
                                                                                TRUE ~ 0),
                                       `Down in emergence vs. winter` = case_when(gene %in% sp_results_down$sp_down_PH1vW$gene_id ~ 1,
                                                                                  TRUE ~ 0),
                                       `Up in spring vs. emergence` = case_when(gene %in% sp_results_up$sp_up_PH1vPH2$gene_id ~ 1,
                                                                                TRUE ~ 0),
                                       `Down in spring vs. emergence` = case_when(gene %in% sp_results_down$sp_down_PH1vPH2$gene_id ~ 1,
                                                                                  TRUE ~ 0),
                                       `Up in Rhabdias sp. infected vs. uninfected` = case_when(gene %in% sp_results_up$sp_up_rhab$gene_id ~ 1,
                                                                                                TRUE ~ 0),
                                       `Down in Rhabdias sp. infected vs. uninfected` = case_when(gene %in% sp_results_down$sp_down_rhab$gene_id ~ 1,
                                                                                                  TRUE ~ 0),
                                       `Up in Oswaldocruzia sp. infected vs. uninfected` = case_when(gene %in% sp_results_up$sp_up_oswa$gene_id ~ 1,
                                                                                                     TRUE ~ 0),
                                       `Down in Oswaldocruzia sp. infected vs. uninfected` = case_when(gene %in% sp_results_down$sp_down_owsa$gene_id ~ 1,
                                                                                                       TRUE ~ 0),
                                       `Up in male vs. female` = case_when(gene %in% sp_results_up$sp_up_sex$gene_id ~ 1,
                                                                           TRUE ~ 0),
                                       `Down in male vs. female` = case_when(gene %in% sp_results_down$sp_down_sex$gene_id ~ 1,
                                                                             TRUE ~ 0),
)

names(upsetr_df_sp)


upsetr_df_sp<-upsetr_df_sp %>% mutate(`GO Term` = case_when(gene %in% adap_immune_genes$gene_id & gene %in% immune_system_process_genes$gene_id ~ "Adaptive",
                                                            
                                                            gene %in% inna_immune_genes$gene_id & gene %in% immune_system_process_genes$gene_id ~ "Innate",
                                                            
                                                            gene %in% immune_system_process_genes$gene_id &  
                                                              gene %notin%adap_immune_genes$gene_id & 
                                                              gene%notin%inna_immune_genes$gene_id~ "Other immune",
                                                            
                                                            TRUE ~ "Not immune related"))

upsetr_df_sp$`GO Term` <- factor(upsetr_df_sp$`GO Term`, levels = c("Adaptive", "Innate", "Other immune", "Not immune related"))



(sp_main<-ComplexUpset::upset(upsetr_df_sp, name = "Comparison",intersect = names(upsetr_df)[c(2:13)],
                              width_ratio=0.3,height_ratio = 0.8, min_size=10, max_size = 5000, wrap = TRUE,
                              themes=upset_modify_themes(
                                list(
                                  'intersections_matrix'=theme(text=element_text(size=12)))),
                              queries=list(
                                upset_query(set='Up in Rhabdias sp. infected vs. uninfected', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in Rhabdias sp. infected vs. uninfected', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in spring vs. emergence', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in spring vs. emergence', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in emergence vs. winter', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in emergence vs. winter', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set='Up in winter vs. fall', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in winter vs. fall', fill="darkblue", only_components=c('intersections_matrix')),
                                upset_query(set="Up in Oswaldocruzia sp. infected vs. uninfected", fill = "red", only_components = c('intersections_matrix')),
                                upset_query(set="Down in Oswaldocruzia sp. infected vs. uninfected", fill = "darkblue", only_components = c('intersections_matrix')),
                                upset_query(set='Up in male vs. female', fill="red", only_components=c('intersections_matrix')),
                                upset_query(set='Down in male vs. female', fill="darkblue", only_components=c('intersections_matrix'))
                              ),
                              base_annotations=list(
                                'Intersection size'=(intersection_size(text=list(size=4),bar_number_threshold=1,  # show all numbers on top of bars
                                                                       width=0.6,   # reduce width of the bars
                                                                       mapping=aes(fill=`GO Term`)) 
                                                     + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                                     
                                                     # add some space on the top of the bars
                                                     + scale_y_continuous(expand=expansion(mult=c(0, 0.2)))
                                                     # +scale_y_break(breaks = c(150,1500), scales = "fixed")
                                                     
                                                     + ggtitle('Spleen Differential Expression')
                                                     
                                                     + guides(fill = "none")
                                                     + theme(
                                                       # hide grid lines
                                                       panel.grid.major=element_blank(),
                                                       panel.grid.minor=element_blank(),
                                                       # show axis lines
                                                       axis.line=element_line(colour='black'),
                                                       axis.text = element_text(size = 12),
                                                       axis.title.y = element_text(size = 14),
                                                       plot.title = element_text(size = 16),
                                                       legend.title = element_text(size = 14),
                                                       legend.text = element_text(size = 14)
                                                       
                                                     )
                                )
                                # 'Intersection ratio'=intersection_ratio()
                                
                              ),
                              stripes=upset_stripes(
                                geom=geom_segment(size=6),  # make the stripes larger
                                colors=c('gray92', 'white')
                              ),
                              # to prevent connectors from getting the colorured
                              # use `fill` instead of `color`, together with `shape='circle filled'`
                              matrix=intersection_matrix(
                                geom=geom_point(
                                  shape='circle filled',
                                  size=3.5,
                                  stroke=0.45
                                )
                              ),
                              
                              set_sizes=(
                                upset_set_size(geom=geom_bar(aes(fill = `GO Term`),width=0.4),position = "right" )
                                + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                + geom_text(aes(label=after_stat(count)), hjust=0, stat='count')
                                + expand_limits(y=1500)
                                + guides(fill = "none")
                                + theme(
                                  axis.line.x=element_line(colour='black'),
                                  axis.ticks.x=element_line(),
                                  axis.text.x=element_text(angle=90, size =12),
                                  # axis.text = element_text(size = 14)
                                  
                                )
                                
                              ),
                              guides = 'over',
                              
                              sort_sets=FALSE,
                              sort_intersections='descending'
)
)




###Figure 5----
sp_go_sig_res_over<- readRDS("yourpath/spleen_go_underrepresented.rds")
vs_go_sig_res_over<- readRDS("youtpath/ventralskin_go_overrepresented.rds")

upsetr_df_go<- data.frame(GO.term=unique(go_map$GO.term))
names(upsetr_df_go)
length(unique(go_map$GO.term))

adap_offspring_bp = GOBPOFFSPRING[["GO:0002250"]]
inna_offspring_bp = GOBPOFFSPRING[["GO:0045087"]]
immun_offspring_bp = GOBPOFFSPRING[["GO:0002376"]]

adap<-which(immun_offspring_bp %in% adap_offspring_bp)
inna<-which(immun_offspring_bp %in% inna_offspring_bp)
adap_inna<- c(adap,inna)

immun_other_bp = immun_offspring_bp[-adap_inna]


# ventral skin enriched ----


upsetr_df_go_vs_over <-upsetr_df_go %>%
  mutate(`Up in winter vs. fall`  = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_WvF$category ~ 1,
                                              TRUE ~ 0),
         `Down in winter vs. fall` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_WvF$category ~ 1,
                                               TRUE ~ 0),
         `Up in emergence vs. winter` = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_PH1vW$category ~ 1,
                                                  TRUE ~ 0),
         `Down in emergence vs. winter` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_PH1vW$category ~ 1,
                                                    TRUE ~ 0),
         `Up in spring vs. emergence` = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_PH1vPH2$category ~ 1,
                                                  TRUE ~ 0),
         `Down in spring vs. emergence` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_PH1vPH2$category ~ 1,
                                                    TRUE ~ 0),
         `Up in Rhabdias sp. infected vs. uninfected` = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_rhab$category ~ 1,
                                                                  TRUE ~ 0),
         `Down in Rhabdias sp. infected vs. uninfected` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_rhab$category ~ 1,
                                                                    TRUE ~ 0),
         `Up in Oswaldocruzia sp. infected vs. uninfected` = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_oswa$category ~ 1,
                                                                       TRUE ~ 0),
         `Down in Oswaldocruzia sp. infected vs. uninfected` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_oswa$category ~ 1,
                                                                         TRUE ~ 0),
         `Up in male vs. female` = case_when(GO.term %in% vs_go_sig_res_over$de_upvs_res_sex$category ~ 1,
                                             TRUE ~ 0),
         `Down in male vs. female` = case_when(GO.term %in% vs_go_sig_res_over$de_downvs_res_sex$category ~ 1,
                                               TRUE ~ 0))




upsetr_df_go_vs_over<-upsetr_df_go_vs_over %>% 
  mutate(`Immune annotation` = case_when(GO.term %in% adap_offspring_bp ~ "Adaptive",
                                         GO.term %in% inna_offspring_bp ~ "Innate",
                                         GO.term %in% immun_other_bp ~ "Other immune",
                                         TRUE ~ "Not immune related"))

upsetr_df_go_vs_over$`Immune annotation` <- factor(upsetr_df_go_vs_over$`Immune annotation`, levels = c("Adaptive", "Innate", "Other immune", "Not immune related"))


unique(upsetr_df_go_vs_over$`Immune annotation`)
(vs_main_go_over<-ComplexUpset::upset(upsetr_df_go_vs_over, name = "Comparison",intersect = names(upsetr_df_go_vs_over)[c(2:13)],
                                      width_ratio=0.3,height_ratio = 0.8, min_size=10, max_size = 5000, wrap = TRUE,
                                      themes=upset_modify_themes(
                                        list(
                                          'intersections_matrix'=theme(text=element_text(size=12)))),
                                      queries=list(
                                        upset_query(set=names(upsetr_df_go_vs_over[8]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[9]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[2]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[3]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[4]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[5]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[6]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[7]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[10]), fill = "darkblue", only_components = c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[11]), fill = "red", only_components = c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[12]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_vs_over[13]), fill="darkblue", only_components=c('intersections_matrix'))
                                      ),
                                      base_annotations=list(
                                        'Intersection size'=(intersection_size(text=list(size=3),bar_number_threshold=1,  # show all numbers on top of bars
                                                                               width=0.6,   # reduce width of the bars
                                                                               mapping=aes(fill=`Immune annotation`)) 
                                                             + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                                             
                                                             # add some space on the top of the bars
                                                             + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
                                                             # +scale_y_break(breaks = c(150,1500), scales = "fixed")
                                                             
                                                             + ggtitle('Enriched ventral skin GO terms')
                                                             
                                                             + guides(fill = "none")
                                                             + theme(
                                                               # hide grid lines
                                                               panel.grid.major=element_blank(),
                                                               panel.grid.minor=element_blank(),
                                                               # show axis lines
                                                               axis.line=element_line(colour='black'),
                                                               axis.text = element_text(size = 12),
                                                               axis.title.y = element_text(size = 14),
                                                               plot.title = element_text(size = 16),
                                                               legend.title = element_text(size = 14),
                                                               legend.text = element_text(size = 14)
                                                             )
                                        )
                                        # 'Intersection ratio'=intersection_ratio()
                                        
                                      ),
                                      stripes=upset_stripes(
                                        geom=geom_segment(size=8),  # make the stripes larger
                                        colors=c('gray92', 'white')
                                      ),
                                      # to prevent connectors from getting the colorured
                                      # use `fill` instead of `color`, together with `shape='circle filled'`
                                      matrix=intersection_matrix(
                                        geom=geom_point(
                                          shape='circle filled',
                                          size=3.5,
                                          stroke=0.45
                                        )
                                      ),
                                      
                                      set_sizes=(
                                        upset_set_size(geom=geom_bar(aes(fill = `Immune annotation`),width=0.4),position = "right" )
                                        + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                        + geom_text(aes(label=after_stat(count)), hjust=0, stat='count')
                                        + expand_limits(y=300)
                                        + theme(
                                          axis.line.x=element_line(colour='black'),
                                          axis.ticks.x=element_line(),
                                          axis.text.x=element_text(angle=90)
                                        )
                                        
                                      ),
                                      guides = 'over',
                                      
                                      sort_sets=FALSE,
                                      sort_intersections='descending',
)
)


# spleen enriched ----
upsetr_df_go_sp_over <-upsetr_df_go %>%
  mutate(`Up in winter vs. fall`  = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_WvF$category ~ 1,
                                              TRUE ~ 0),
         `Down in winter vs. fall` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_WvF$category ~ 1,
                                               TRUE ~ 0),
         `Up  in emergence vs. winter` = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_PH1vW$category ~ 1,
                                                   TRUE ~ 0),
         `Down in emergence vs. winter` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_PH1vW$category ~ 1,
                                                    TRUE ~ 0),
         `Up in spring vs. emergence` = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_PH1vPH2$category ~ 1,
                                                  TRUE ~ 0),
         `Down in spring vs. emergence` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_PH1vPH2$category ~ 1,
                                                    TRUE ~ 0),
         `Up in Rhabdias sp. infected vs. uninfected` = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_rhab$category ~ 1,
                                                                  TRUE ~ 0),
         `Down in Rhabdias sp. infected vs. uninfected` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_rhab$category ~ 1,
                                                                    TRUE ~ 0),
         `Up in Oswaldocruzia sp. infected vs. uninfected` = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_oswa$category ~ 1,
                                                                       TRUE ~ 0),
         `Down in Oswaldocruzia sp. infected vs. uninfected` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_oswa$category ~ 1,
                                                                         TRUE ~ 0),
         `Up in male vs. female` = case_when(GO.term %in% sp_go_sig_res_over$de_upsp_res_sex$category ~ 1,
                                             TRUE ~ 0),
         `Down in male vs. female` = case_when(GO.term %in% sp_go_sig_res_over$de_downsp_res_sex$category ~ 1,
                                               TRUE ~ 0))

names(upsetr_df_go_sp_over)
head(Gomap)






upsetr_df_go_sp_over<-upsetr_df_go_sp_over %>% 
  mutate(`Immune annotation` = case_when(GO.term %in% adap_offspring_bp ~ "Adaptive",
                                         GO.term %in% inna_offspring_bp ~ "Innate",
                                         GO.term %in% immun_other_bp ~ "Other immune",
                                         TRUE ~ "Not immune related"))
unique(upsetr_df_go_sp_over$`Immune annotation`)

head(upsetr_df_go_sp_over)

upsetr_df_go_sp_over$`Immune annotation` <- factor(upsetr_df_go_sp_over$`Immune annotation`, levels = c("Adaptive", "Innate", "Other immune", "Not immune related"))


(sp_main_go_over<-ComplexUpset::upset(upsetr_df_go_sp_over, name = "Comparison",intersect = names(upsetr_df_go_sp_over)[c(2:13)],
                                      width_ratio=0.3,height_ratio = 0.8, min_size=10, max_size = 5000, wrap = TRUE,
                                      themes=upset_modify_themes(
                                        list(
                                          'intersections_matrix'=theme(text=element_text(size=12)))),
                                      queries=list(
                                        upset_query(set=names(upsetr_df_go_sp_over[8]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[9]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[2]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[3]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[4]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[5]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[6]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[7]), fill="darkblue", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[10]), fill = "darkblue", only_components = c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[11]), fill = "red", only_components = c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[12]), fill="red", only_components=c('intersections_matrix')),
                                        upset_query(set=names(upsetr_df_go_sp_over[13]), fill="darkblue", only_components=c('intersections_matrix'))
                                      ),
                                      base_annotations=list(
                                        'Intersection size'=(intersection_size(text=list(size=3),bar_number_threshold=1,  # show all numbers on top of bars
                                                                               width=0.6,   # reduce width of the bars
                                                                               mapping=aes(fill=`Immune annotation`)) 
                                                             + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                                             
                                                             # add some space on the top of the bars
                                                             + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
                                                             # +scale_y_break(breaks = c(150,1500), scales = "fixed")
                                                             
                                                             + ggtitle('Enriched spleen GO terms')
                                                             
                                                             + guides(fill = "none")
                                                             + theme(
                                                               # hide grid lines
                                                               panel.grid.major=element_blank(),
                                                               panel.grid.minor=element_blank(),
                                                               # show axis lines
                                                               axis.line=element_line(colour='black'),
                                                               axis.text = element_text(size = 12),
                                                               axis.title.y = element_text(size = 14),
                                                               plot.title = element_text(size = 16),
                                                               legend.title = element_text(size = 14),
                                                               legend.text = element_text(size = 14)
                                                             )
                                        )
                                        # 'Intersection ratio'=intersection_ratio()
                                        
                                      ),
                                      stripes=upset_stripes(
                                        geom=geom_segment(size=8),  # make the stripes larger
                                        colors=c('gray92', 'white')
                                      ),
                                      # to prevent connectors from getting the colorured
                                      # use `fill` instead of `color`, together with `shape='circle filled'`
                                      matrix=intersection_matrix(
                                        geom=geom_point(
                                          shape='circle filled',
                                          size=3.5,
                                          stroke=0.45
                                        )
                                      ),
                                      
                                      set_sizes=(
                                        upset_set_size(geom=geom_bar(aes(fill = `Immune annotation`),width=0.4),position = "right" )
                                        + scale_fill_manual(values = c("Not immune related"= "darkgray", "Adaptive" = "#1b9e77","Innate"="#d95f02","Other immune" ="#7570b3"))
                                        + geom_text(aes(label=after_stat(count)), hjust=0, stat='count')
                                        + expand_limits(y=300)
                                        + theme(
                                          axis.line.x=element_line(colour='black'),
                                          axis.ticks.x=element_line(),
                                          axis.text.x=element_text(angle=90)
                                        )
                                        
                                      ),
                                      guides = 'over',
                                      
                                      sort_sets=FALSE,
                                      sort_intersections='descending',
)
)

