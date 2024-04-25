# All----
# normalized, variance-stabilized transformed counts for visualization
all_vsd <- vst(all_ddsClean, blind=FALSE)

#  plot using ggplot
all_dat <- plotPCA(all_vsd,returnData=TRUE,intgroup=c("period", "tissue"))

(all_p <- ggplot(all_dat,aes(x=PC1,y=PC2,col=paste(period), shape=paste(tissue))))
all_p <- all_p + geom_point(size =5) + 
  xlab(paste("PC1: ", round(attr(all_dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(all_dat,"percentVar")[2],2)*100, "% variation explained", sep="")) 

all_p


# Ventral Skin----
# normalized, variance-stabilized transformed counts for visualization
vs_vsd <- vst(vs_ddsClean, blind=FALSE)

#  plot using ggplot
vs_dat <- plotPCA(vs_vsd,returnData=TRUE,intgroup=c("period", "sex"))


(vs_p <- ggplot(vs_dat,aes(x=PC1,y=PC2,col=paste(period), shape=paste(sex))))
vs_p <- vs_p + geom_point(size =5) + 
  xlab(paste("PC1: ", round(attr(vs_dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(vs_dat,"percentVar")[2],2)*100, "% variation explained", sep="")) 
vs_p


# Spleen----
# normalized, variance-stabilized transformed counts for visualization

sp_vsd <- vst(sp_ddsClean, blind=FALSE)

#  plot using ggplot
sp_dat <- plotPCA(sp_vsd,returnData=TRUE,intgroup=c("period", "sex"))

(sp_p <- ggplot(sp_dat,aes(x=PC1,y=PC2,col=paste(period), shape=paste(sex))))
sp_p <- sp_p + geom_point(size = 5) + 
  xlab(paste("PC1: ", round(attr(sp_dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(sp_dat,"percentVar")[2],2)*100, "% variation explained", sep="")) 
sp_p


