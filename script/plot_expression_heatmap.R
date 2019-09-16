library(pheatmap)
setwd("/Users/slbai/Documents/git_local/HET_identity/小立碗藓HET/小立碗藓自身/表达量")

tree_order = c("Pp3c27_8690", "Pp3c27_8695", "Pp3c27_8710","Pp3c13_970", "Pp3c9_560",
               "Pp3c8_25331", "Pp3c8_25330", "Pp3c12_12000", "Pp3c13_190", "Pp3c27_8560",
               "Pp3c8_22200", "Pp3c8_22202", "Pp3c15_6400", "Pp3c27_8820", "Pp3c27_8810", 
               "Pp3c14_10", "Pp3c27_8580", "Pp3c5_11620", "Pp3c27_8740", "Pp3c13_14910", 
               "Pp3c27_8610", "Pp3c4_30400", "Pp3c13_930", "Pp3c13_650", "Pp3c13_670", 
               "Pp3c20_19770","Pp3c22_17530" )
#HET FPKM
all_fpkm<-read.table('het_exp_t.txt',header = T,sep = '\t' ,
                     row.names = 1, stringsAsFactors = F)
#remove FPKM lower than 0.1 gene
remain_data <- all_fpkm[apply(all_fpkm, MARGIN = 1, FUN = function(x) max(x) > 0.1),]
remove_data <- all_fpkm[apply(all_fpkm, MARGIN = 1, FUN = function(x) max(x) <= 0.1),]
remove_data <- 0 * remove_data
# myGene_fpkm_log=log2(all_fpkm + 1)

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_data <- scale_rows(remain_data)

plot_data <- rbind(scale_data, remove_data)[tree_order, ]

bk <- seq(-max(abs(plot_data))+0.1, max(abs(plot_data))+0.1, by=0.01)

pheatmap(plot_data, cluster_cols = F, cluster_rows = F, scale = "none",
         # color = c(colorRampPalette(c("blue","white"))(100),colorRampPalette(c("white","red"))(100)),
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-4,4,2),
         breaks=bk,
         cellwidth = 8, cellheight = 8,
         # border_color = NA,
        filename = 'het_FPKM_row_scale.pdf')

#HET RMA
all_fpkm<-read.table('het_RMA_t.txt',header = T,sep = '\t' ,
                     row.names = 1, stringsAsFactors = F)
#remove FPKM lower than 0.1 gene
remain_data <- all_fpkm[apply(all_fpkm, MARGIN = 1, FUN = function(x) max(x) > 0.1),]
remove_data <- all_fpkm[apply(all_fpkm, MARGIN = 1, FUN = function(x) max(x) <= 0.1),]
remove_data <- 0 * remove_data
# myGene_fpkm_log=log2(all_fpkm + 1)

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_data <- scale_rows(remain_data)

plot_data <- rbind(scale_data, remove_data)[tree_order, ]

bk <- seq(-max(abs(plot_data))+0.1, max(abs(plot_data))+0.1, by=0.01)

pheatmap(plot_data, cluster_cols = F, cluster_rows = F, scale = "none",
         # color = c(colorRampPalette(c("blue","white"))(100),colorRampPalette(c("white","red"))(100)),
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         cellwidth = 8, cellheight = 8,
         # border_color = NA,
         filename = 'het_RMA_row_scale.pdf')
