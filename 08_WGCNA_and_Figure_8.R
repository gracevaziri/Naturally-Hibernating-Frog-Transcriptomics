library(forcats)

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ComplexHeatmap")
}
library(magrittr)
library(clue)
library(WGCNA)
####################################################################################

#WGCNA----

####################################################################################


# setting up a custom plotting theme for text elements ----
my_theme<-theme(
  text = element_text(size = 14),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 12),
  legend.text = element_text(size = 12)
  # Add more customization as needed
)

# All----

# Ventral Skin----
datExpr_vs =t(dat_vs)

gsg_vs = goodSamplesGenes(datExpr_vs, verbose = 3);
gsg_vs$allOK

if (!gsg_vs$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_vs$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr_vs)[!gsg_vs$goodGenes], collapse = ", ")));
  if (sum(!gsg_vs$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr_vs)[!gsg_vs$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr_vs = datExpr_vs[gsg_vs$goodSamples, gsg_vs$goodGenes]
}

sampleTree_vs = hclust(dist(datExpr_vs), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_vs, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Define numbers of genes and samples
nGenes_vs = ncol(datExpr_vs)
nSamples_vs = nrow(datExpr_vs)


traitData_vs = sampleTable[sampleTable$sampleName %in% rownames(datExpr_vs),]

traitData_vs = traitData_vs  %>% mutate(Fall = case_when(period == "1" ~ 1,
                                                         TRUE ~ 0),
                                        Winter = case_when(period == "3" ~ 1,
                                                           TRUE ~ 0),
                                        Emergence = case_when(period=="4" ~ 1,
                                                              TRUE ~0),
                                        Spring = case_when(period=="6" ~ 1,
                                                           TRUE ~0),
                                        SexBin = case_when(sex == "F" ~ 0,
                                                           sex == "M" ~1),
                                        parasitizedbin = case_when( parasitized== "Y" ~ 1,
                                                                    parasitized == "N" ~ 0))
# remove columns that hold information we do not need.
allTraits_vs = traitData_vs[, -c(1:5,8:11,14:16)];

collectGarbage();

# Re-cluster samples
sampleTree2_vs = hclust(dist(datExpr_vs), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_vs = numbers2colors(allTraits_vs, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_vs, traitColors_vs,
                    groupLabels = names(allTraits_vs), 
                    main = "Sample dendrogram and trait heatmap")

enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft_vs = pickSoftThreshold(datExpr_vs, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_vs$fitIndices[,1], -sign(sft_vs$fitIndices[,3])*sft_vs$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_vs$fitIndices[,1], -sign(sft_vs$fitIndices[,3])*sft_vs$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_vs$fitIndices[,1], sft_vs$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_vs$fitIndices[,1], sft_vs$fitIndices[,5], labels=powers, cex=cex1,col="red")

cor <- WGCNA::cor
set.seed(123213)
bwnet_vs = WGCNA::blockwiseModules(datExpr_vs, maxBlockSize = 15500,
                                   power = 6, TOMType = "signed",networkType = "signed",
                                   minModuleSize = 100,
                                   reassignThreshold = 0, mergeCutHeight = 0.25,
                                   numericLabels = TRUE,
                                   saveTOMs = TRUE,
                                   # randomSeed = 54321,
                                   saveTOMFileBase = "VentralSkinTOM",
                                   verbose = 3)


bwLabels_vs =bwnet_vs$colors
bwColors_vs = labels2colors(bwLabels_vs)

table(bwnet_vs$colors)

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors_vs = labels2colors(bwnet_vs$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet_vs$dendrograms[[1]], mergedColors_vs[bwnet_vs$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels_vs = bwnet_vs$colors
moduleColors_vs = labels2colors(bwnet_vs$colors)

MEs_vs = bwnet_vs$MEs;
colnames(MEs_vs) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs_vs),"ME",""))))

MEs_vs = orderMEs(MEs_vs)

# The correlation map of each module obtained by clustering according to the expression of genes
# marDendro/marHeatmap sets the bottom, left, top, and right margins
plotEigengeneNetworks(MEs_vs, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
## If there is phenotypic data, it can also be put together with ME data and plotted together
MEs_colpheno_vs = orderMEs(cbind(MEs_vs, allTraits_vs))
plotEigengeneNetworks(MEs_colpheno_vs, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)

geneTree_vs = bwnet_vs$dendrograms[[1]];

# Recalculate MEs with color labels
MEs0_vs = moduleEigengenes(datExpr_vs, moduleColors_vs)$eigengenes
MEs_vs = orderMEs(MEs0_vs)
moduleTraitCor_vs = cor(MEs_vs, allTraits_vs, use = "p");
moduleTraitPvalue_vs = corPvalueStudent(moduleTraitCor_vs, nSamples_vs);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix_vs =  paste(signif(moduleTraitCor_vs, 2), "\n(",
                       signif(moduleTraitPvalue_vs, 1), ")", sep = "");
dim(textMatrix_vs) = dim(moduleTraitCor_vs)


yLabels_vs <- names(MEs_vs)
yLabels_vs <- paste0("\\textcolor{black}{", yLabels_vs, "}")
rownames(textMatrix_vs) <- yLabels_vs


# Display the correlation values within a heatmap plot
par(mar = c(6, 8.5, 3, 3));labeledHeatmap.multiPage(Matrix = moduleTraitCor_vs,
                                                    xLabels = names(allTraits_vs),
                                                    yLabels = names(MEs_vs),
                                                    ySymbols = names(MEs_vs),
                                                    colorLabels = TRUE,
                                                    colors = blueWhiteRed(50),
                                                    textMatrix = textMatrix_vs,
                                                    setStdMargins = FALSE,
                                                    cex.text = 0.6,
                                                    rowsPerPage = NULL, maxRowsPerPage = 15, 
                                                    colsPerPage = NULL, maxColsPerPage = 10, 
                                                    addPageNumberToMain = TRUE, 
                                                    zlim = c(-1,1),
                                                    main = paste("Module-trait relationships"))

#identify modules that might be interesting based on high correlation with experimental variable
int_modules_vs<- moduleTraitCor_vs %>% as.data.frame() %>% filter(abs(svl) >=0.7 | abs(mass) >=0.7 | abs(rhabdias_log) >=0.7|
                                                                    abs(oswaldocruzia_log) >=0.7 | abs(Winter) >=0.7 |
                                                                    abs(Fall) >=0.7 | abs(Emergence) >=0.7 | 
                                                                    abs(Spring) >=0.7 | abs(SexBin) >= 0.7 ) %>% rownames() 

int_modules_vs

modules_vs<- substr(rownames(moduleTraitCor_vs), 3,49)


# Read in the probe annotation
m
# Match probes in the data set to the probe IDs in the annotation file
probes_vs = names(data.frame(datExpr_vs))
probes2annot_vs = match(probes_vs, m$gene_id)
# The following is the number or probes without annotation:
sum(is.na(probes2annot_vs))

# Get the corresponding Locus Link IDs
allLLIDs_vs = m$entrezID[probes2annot_vs];
str(allLLIDs_vs)
allLLIDs_vs<-as.integer(allLLIDs_vs)
sum(is.na(allLLIDs_vs))
sum(!is.na(allLLIDs_vs))

# $ Choose interesting modules

# Red was also interesting, downreg with both osw and rhab infections (lots to do with lipid metabolism - upreg in postwinter2 and fall, downreg in winter adn postwinter1, )
int_modules_vs
# [1] "MElightsteelblue1" "MEturquoise"       "MEpink"            "MEdarkolivegreen"  "MEblue"            "MEgreen"          

intModules_vs = c("lightsteelblue", "turquoise", "pink", "darkolivegreen", "blue", "green")



# make lists named after the modules that show all the genes present in each module
ventralskin_module_genelists<- list()
for (i in intModules_vs)
{
  # Select module probes
  modGenes = (moduleColors_vs==i)
  # Get their entrez ID codes
  modProbes = probes_vs[modGenes];
  #remove NAs
  modProbes = modProbes[which(!is.na(modProbes))]
  # Write them into a file
  fileName = paste("ProbeIDs-", i, ".txt", sep="");
  write.table(as.data.frame(modProbes), file = fileName,
              row.names = FALSE, col.names = FALSE)
  ventralskin_module_genelists[[i]] <- modProbes
  
}


# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("VS_ProbeIDs-all.txt", sep="");
write.table(as.data.frame(probes_vs), file = fileName,
            row.names = FALSE, col.names = FALSE)

getwd()
list.files()
colorh_vs = labels2colors(bwnet_vs$colors)

#find the hub gene for each module
vs_hub_genes = chooseTopHubInEachModule(datExpr_vs, colorh_vs  )
intModules_vs

#turquoise
tab %>% filter(gene_id=="MSTRG.38521") %>% dplyr::select(Predicted.Gene) #NAA40
#  Acts as a negative regulator of apoptosis (PubMed:26666750). May play a role in hepatic lipid metabolism (By similarity). 

# pink
tab %>% filter(gene_id=="MSTRG.31141") %>% dplyr::select(Predicted.Gene) #XPO6

# blue
tab %>% filter(gene_id=="MSTRG.39761") %>% dplyr::select(Predicted.Gene) #GTF2H1

# lightsteelblue
tab %>% filter(gene_id=="MSTRG.13944") %>% dplyr::select(Predicted.Gene) #MANEA

# darkolivegreen
tab %>% filter(gene_id=="MSTRG.24103") %>% dplyr::select(Predicted.Gene) #ZBED6

# green
tab %>% filter(gene_id=="MSTRG.25330") %>% dplyr::select(Predicted.Gene) #ZNF317

hub_vs_df<-vs_normalized_counts %>%as.data.frame() %>% filter(rownames(vs_normalized_counts) %in% data.frame(vs_hub_genes)$vs_hub_genes)
hub_vs_df<-data.frame(t(hub_vs_df))

hub_vs_df=hub_vs_df %>% rownames_to_column() %>% pivot_longer(cols = starts_with("MSTRG"),names_to="gene_id") %>%
  mutate(period = substr(rowname, 4,5)) %>%
  mutate(period = case_when(period == "01" ~ "Fall",
                            period == "03" ~ "Winter",
                            period == "04" ~ "Emergence",
                            period == "06" ~ "Spring"))

hub_vs_df=hub_vs_df %>% 
  mutate(predicted_gene = m$Predicted.Gene[match(hub_vs_df$gene_id, m$gene_id)])%>%
    mutate(tissue = "Ventral Skin",
         module =rownames(data.frame(vs_hub_genes))[match(hub_vs_df$gene_id, data.frame(vs_hub_genes)$vs_hub_genes)])



###---- complex heatmap for WGCNA
module_eigengenes_vs <- bwnet_vs$MEs

# Print out a preview
head(module_eigengenes_vs)

#see if our eigengenes relate to our metadata labels. 
# First we double check that our samples are still in order.
all.equal(traitData_vs$sampleName, rownames(MEs_vs))

# Create the design matrix from the `period` variable
des_mat_vs <- model.matrix(~ traitData_vs$period +traitData_vs$sex+ traitData_vs$rhabdias_log + traitData_vs$oswaldocruzia_log)

# Run linear model on each module. Limma wants our tests to be per row, 
# so we also need to transpose so the eigengenes are rows

# lmFit() needs a transposed version of the matrix
fit_vs <- limma::lmFit(t(module_eigengenes_vs), design = des_mat_vs)

# Apply empirical Bayes to smooth standard errors
fit_vs <- limma::eBayes(fit_vs)

# Apply multiple testing correction and obtain stats
stats_df_vs <- limma::topTable(fit_vs, number = ncol(module_eigengenes_vs)) %>%
  tibble::rownames_to_column("module")

head(stats_df_vs,49)
vs_col_mod_size = data.frame(table(moduleColors_vs))
vs_col_mod_size =vs_col_mod_size %>% arrange(Freq)
vs_mod_size = data.frame(table(moduleLabels_vs))
vs_mod_size =vs_mod_size %>% arrange(Freq)
vs_mods<-cbind(vs_col_mod_size, vs_mod_size)

#add module size to stats_df_sp
stats_df_vs=stats_df_vs %>% mutate(module_num = substr(module, 3,5)) 
stats_df_vs=stats_df_vs %>% mutate( module_size = vs_mod_size$Freq[match(stats_df_vs$module_num, vs_mod_size$moduleLabels_vs)])
stats_df_vs=stats_df_vs %>% mutate(mod_col = vs_mods$moduleColors_vs[match(
  stats_df_vs$module_size, vs_mods$Freq)])


# # let’s use ggplot to see what module turquoise eigengene 
# # looks like among periods
# Module_vs_df <- module_eigengenes_vs %>%
#   tibble::rownames_to_column("sampleID") %>%
#   # Here we are performing an inner join with a subset of metadata
#   dplyr::inner_join(traitData_vs %>%
#                       dplyr::select(sampleName, period),
#                     by = c("sampleID" = "sampleName")
#   )
# 
# 
# Module_vs_df_long<-Module_vs_df %>% pivot_longer(cols=!c(sampleID ,period),names_to="eigengene",values_to="Average log2-expression")


gene_module_key_vs <- tibble::enframe(bwnet_vs$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))


make_module_heatmap_vs <- function(module_name,
                                   expression_mat = vs_normalized_counts,
                                   metadata_df = traitData_vs %>% 
                                     mutate(period =fct_recode(period,"Fall"= "1",
                                                               "Winter"="3",
                                                               "Emergence"="4", 
                                                               "Spring"="6")),
                                   gene_module_key_df = gene_module_key_vs,
                                   module_eigengenes_df = module_eigengenes_vs) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sampleName")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sampleName, period) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sampleName") %>%
    # Arrange by patient and time point
    dplyr::arrange(period, sampleName) %>%
    # Store sample
    tibble::column_to_rownames("sampleName")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    `Sample Period` = col_annot_df$period,
    # Add annotation barplot
    `Module Eigengene` = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(`Sample Period` = c("Fall" = "#925A44", "Winter" = "#155289", "Emergence" = "#489FA7", "Spring"="#335A30"))
  )
  
  # c("#925A44", "#155289","#489FA7","#335A30")
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    # t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = "Normalized gene expression",
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
                                     # show_heatmap_legend	= FALSE
                                     
  )
  
  # Return heatmap
  return(heatmap)
}

# Supplemental Information Figure 1----
make_module_heatmap_vs(module_name = "ME1")


# Spleen----
datExpr_sp =t(dat_sp)

gsg_sp = goodSamplesGenes(datExpr_sp, verbose = 3);
gsg_sp$allOK

if (!gsg_sp$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_sp$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr_sp)[!gsg_sp$goodGenes], collapse = ", ")));
  if (sum(!gsg_sp$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr_sp)[!gsg_sp$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr_sp = datExpr_sp[gsg_sp$goodSamples, gsg_sp$goodGenes]
}

sampleTree_sp = hclust(dist(datExpr_sp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)


par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_sp, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Define numbers of genes and samples

nGenes_sp = ncol(datExpr_sp)
nSamples_sp = nrow(datExpr_sp)

# traitData_sp = metadata[metadata$sampleNames %in% rownames(datExpr_sp),]
traitData_sp = sampleTable[sampleTable$sampleName %in% rownames(datExpr_sp),]

traitData_sp = traitData_sp %>% mutate(Fall = case_when(period == "1" ~ 1,
                                                        TRUE ~ 0),
                                       Winter = case_when(period == "3" ~ 1,
                                                          TRUE ~ 0),
                                       PostWinter1 = case_when(period=="4" ~ 1,
                                                               TRUE ~0),
                                       PostWinter2 = case_when(period=="6" ~ 1,
                                                               TRUE ~0),
                                       SexBin = case_when(sex == "F" ~ 0,
                                                          sex == "M" ~1),
                                       parasitizedbin = case_when( parasitized== "Y" ~ 1,
                                                                   parasitized == "N" ~ 0))


# remove columns that hold information we do not need.
allTraits_sp = traitData_sp[, -c(1:5,8:11,14:16)];

collectGarbage();


# Re-cluster samples
sampleTree2_sp = hclust(dist(datExpr_sp), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_sp = numbers2colors(allTraits_sp, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_sp, traitColors_sp,
                    groupLabels = names(allTraits_sp), 
                    main = "Sample dendrogram and trait heatmap")



enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft_sp = pickSoftThreshold(datExpr_sp, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_sp$fitIndices[,1], -sign(sft_sp$fitIndices[,3])*sft_sp$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_sp$fitIndices[,1], -sign(sft_sp$fitIndices[,3])*sft_sp$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_sp$fitIndices[,1], sft_sp$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_sp$fitIndices[,1], sft_sp$fitIndices[,5], labels=powers, cex=cex1,col="red")

set.seed(12345)
bwnet_sp = WGCNA::blockwiseModules(datExpr_sp, maxBlockSize = 14500,
                                   power = 10, TOMType = "signed",networkType = "signed",
                                   minModuleSize = 100,
                                   reassignThreshold = 0, mergeCutHeight = 0.25,
                                   randomSeed = 54321,
                                   numericLabels = TRUE,
                                   saveTOMs = TRUE,
                                   saveTOMFileBase = "SpleenTOM",
                                   verbose = 3)

bwLabels_sp =bwnet_sp$colors
bwColors_sp = labels2colors(bwLabels_sp)

table(bwnet_sp$colors)
 

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors_sp = labels2colors(bwnet_sp$colors)
table(mergedColors_sp)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet_sp$dendrograms[[1]], mergedColors_sp[bwnet_sp$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



moduleLabels_sp = bwnet_sp$colors
sp_mod_size = data.frame(table(moduleLabels_sp))

moduleColors_sp = labels2colors(bwnet_sp$colors)
sp_col_mod_size = data.frame(table(moduleColors_sp))

MEs_sp = bwnet_sp$MEs;
geneTree_sp = bwnet_sp$dendrograms[[1]];
# save(MEs_sp, moduleLabels_sp, moduleColors_sp, geneTree_sp,
#      file = "Spleen-networkConstruction-auto.RData")

# Recalculate MEs with color labels
MEs0_sp = moduleEigengenes(datExpr_sp, moduleColors_sp)$eigengenes
MEs_sp = orderMEs(MEs0_sp)
moduleTraitCor_sp = cor(MEs_sp, allTraits_sp, use = "p");
moduleTraitPvalue_sp = corPvalueStudent(moduleTraitCor_sp, nSamples_sp);

sizeGrWindow(10,6)
# Will display correlations and their p-values

textMatrix_sp =  paste(signif(moduleTraitCor_sp, 2), "\n(",
                       signif(moduleTraitPvalue_sp, 1), ")", sep = "");
dim(textMatrix_sp) = dim(moduleTraitCor_sp)
par(mar = c(6, 10.5, 3, 3));

# Display the correlation values within a heatmap plot
dev.off()
labeledHeatmap.multiPage(Matrix = moduleTraitCor_sp,
                         xLabels = names(allTraits_sp),
                         yLabels = names(MEs_sp),
                         ySymbols = names(MEs_sp),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
                         textMatrix = textMatrix_sp,
                         setStdMargins = FALSE,
                         cex.text = 0.6,
                         rowsPerPage = NULL, maxRowsPerPage = 20, 
                         colsPerPage = NULL, maxColsPerPage = 10, 
                         addPageNumberToMain = TRUE, 
                         zlim = c(-1,1),
                         main = paste("Module-trait relationships"))

# identify modules that might be interesting based on high correlation with experimental variable
int_modules_sp<-moduleTraitCor_sp %>% as.data.frame() %>% filter(abs(svl) >=0.7 | abs(mass) >=0.7 | abs(rhabdias_log) >=0.7|
                                                                   abs(oswaldocruzia_log) >=0.7 | abs(Winter) >=0.7 |
                                                                   abs(Fall) >=0.7 | abs(SexBin) >= 0.7 ) %>%  rownames() 

# Read in the probe annotation
m
# Match probes in the data set to the probe IDs in the annotation file
probes_sp = colnames(datExpr_sp)
probes2annot_sp = match(probes_sp, m$gene_id)
# Get the corresponding Locus Link IDs
allLLIDs_sp = m$entrezID[probes2annot_sp];
str(allLLIDs_sp)
allLLIDs_sp<-as.integer(allLLIDs_sp)
sum(is.na(allLLIDs_sp))
sum(!is.na(allLLIDs_sp))

# $ Choose interesting modules


int_modules_sp
int_modules_sp = c("yellow","turquoise","darkturquoise", "pink")




setwd("../../MakeYourOwnWGCNAOutputDirectory/Module Gene lists/Spleen")
getwd()

# make lists named after the modules that show all the genes present in each module
spleen_module_genelists<- list()
for (i in int_modules_sp)
{
  # Select module probes
  modGenes = (moduleColors_sp==i)
  # Get their entrez ID codes
  modProbes = probes_sp[modGenes];
  #remove NAs
  modProbes = modProbes[which(!is.na(modProbes))]
  # Write them into a file
  fileName = paste("ProbeIDs-", i, ".txt", sep="");
  write.table(as.data.frame(modProbes), file = fileName,
              row.names = FALSE, col.names = FALSE)
   spleen_module_genelists[[i]] <- modProbes
  
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("SP_ProbeIDs-all.txt", sep="");
write.table(as.data.frame(probes_sp), file = fileName,
            row.names = FALSE, col.names = FALSE)

getwd()
list.files()

colorh_sp = labels2colors(bwnet_sp$colors)

#find the hub gene for each module
sp_hub_genes = chooseTopHubInEachModule(datExpr_sp, colorh_sp  )

# yellow
tab %>% filter(gene_id=="MSTRG.39691") %>% dplyr::select(Predicted.Gene) #ARNTL

# darkturquoise
tab %>% filter(gene_id=="MSTRG.5055") %>% dplyr::select(peptide_status) #no peptide

# turquoise
tab %>% filter(gene_id=="MSTRG.8874") %>% dplyr::select(Predicted.Gene) #HNRNPR

# pink
tab %>% filter(gene_id=="MSTRG.49458") %>% dplyr::select(Predicted.Gene) #MEGF8

hub_sp_df<-sp_normalized_counts %>%as.data.frame() %>% filter(rownames(sp_normalized_counts) %in% data.frame(sp_hub_genes)$sp_hub_genes)
hub_sp_df<-data.frame(t(hub_sp_df))

hub_sp_df=hub_sp_df %>% rownames_to_column() %>% pivot_longer(cols = starts_with("MSTRG"),names_to="gene_id") %>%
  mutate(period = substr(rowname, 4,5)) %>%
  mutate(period = case_when(period == "1" ~ "Fall",
                            period == "3" ~ "Winter",
                            period == "4" ~ "Post-winter 1",
                            period == "6" ~ "Post-winter 2"))

hub_sp_df=hub_sp_df %>% 
  mutate(predicted_gene = m$Predicted.Gene[match(hub_sp_df$gene_id, m$gene_id)])%>%
  mutate(tissue = "Spleen",
         module =rownames(data.frame(sp_hub_genes))[match(hub_sp_df$gene_id, data.frame(sp_hub_genes)$sp_hub_genes)])


###---- complex heatmap for WGCNA
module_eigengenes <- bwnet_sp$MEs

# Print out a preview
head(module_eigengenes)

#see if our eigengenes relate to our metadata labels. 
# First we double check that our samples are still in order.
all.equal(traitData_sp$sampleName, rownames(MEs_sp))

# Create the design matrix from the `period` variable

des_mat_sp <- model.matrix(~ traitData_sp$period +traitData_sp$sex+ traitData_sp$rhabdias_log + traitData_sp$oswaldocruzia_log)

# Run linear model on each module. Limma wants our tests to be per row, 
# so we also need to transpose so the eigengenes are rows

# lmFit() needs a transposed version of the matrix
fit_sp <- limma::lmFit(t(module_eigengenes), design = des_mat_sp)

# Apply empirical Bayes to smooth standard errors
fit_sp <- limma::eBayes(fit_sp)

# Apply multiple testing correction and obtain stats
stats_df_sp <- limma::topTable(fit_sp, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df_sp,36)

#add module size to stats_df_sp
stats_df_sp=stats_df_sp %>% mutate(module_num = substr(module, 3,5)) 
stats_df_sp=stats_df_sp %>% mutate( module_size = sp_mod_size$Freq[match(stats_df_sp$module_num, sp_mod_size$moduleLabels_sp)])
stats_df_sp=stats_df_sp %>% mutate(mod_col = sp_col_mod_size$moduleColors_sp[match(
  stats_df_sp$module_size, sp_col_mod_size$Freq)])



# let’s use ggplot to see what module turquoise eigengene 
# looks like among periods
Module_sp_df <- module_eigengenes %>%
  tibble::rownames_to_column("sampleID") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(traitData_sp %>%
                      dplyr::select(sampleName, period),
                    by = c("sampleID" = "sampleName")
  )

Module_sp_df_long<-Module_sp_df %>% pivot_longer(cols=!c(sampleID ,period),names_to="eigengene",values_to="Average log2-expression")

# What genes are a part of Turquoise module?

gene_module_key_sp <- tibble::enframe(bwnet_sp$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key_sp %>%
  dplyr::filter(module == "ME1")


#module
make_module_heatmap_sp <- function(module_name,
                                   expression_mat = sp_normalized_counts,
                                   metadata_df = traitData_sp %>% 
                                     mutate(period =fct_recode(period,"Fall"= "1",
                                                               "Winter"="3",
                                                               "Emergence"="4", 
                                                               "Spring"="6")),
                                   gene_module_key_df = gene_module_key_sp,
                                   module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sampleName")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sampleName, period) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sampleName") %>%
    # Arrange by patient and time point
    dplyr::arrange(period, sampleName) %>%
    # Store sample
    tibble::column_to_rownames("sampleName")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    `Sample Period` = col_annot_df$period,
    # Add annotation barplot
    `Module Eigengene` = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(`Sample Period` = c("Fall" = "#925A44", "Winter" = "#155289", "Emergence" = "#489FA7", "Spring"="#335A30"))
  )
  
  # c("#925A44", "#155289","#489FA7","#335A30")
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    # t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = "Spleen \nTurquoise",
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
                                     # show_heatmap_legend	= FALSE
                                     
  )
  
  # Return heatmap
  return(heatmap)
}
# Figure 8---- 
make_module_heatmap_sp(module_name = "ME1")








