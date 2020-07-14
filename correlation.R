library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(patchwork)

# set workspace
setwd('/Users/orson/Documents/bioinfo/Covid19-CCLE')
rm(list=ls())

# import cell line info
cell.info = read.csv('CCLE_data/cell_line_info/CCLE_info.csv', stringsAsFactors = F)
head(cell.info)

# dissect out the info
colnames(cell.info)
df.info = cell.info[,c(1,2,5,6,9)]
head(df.info)

# import datasets
gene = c('ace2', 'tmprss2', 'tmprss4')
quant = c('Affy', 'RNAseq')

# collect dataframe of Affy and RNAseq 
df.ac = data.frame()
for (g in gene){
  for (q in quant){
    df = read.csv(paste0('CCLE_data/dataset/', q, '_', g,'.csv'), header = T)
    assay = rep(q, nrow(df))
    gn = rep(g, nrow(df))
    df = cbind(df, assay, gn)
    colnames(df) = c('CCLE.name','expr','assay','gene')
    df.ac = rbind(df.ac, df)
  }
}
head(df.ac)

# dcast dataframe
df.ac$expr = df.ac$expr %>% as.numeric()
class(df.ac$exp)
df.dc = dcast(df.ac, formula = CCLE.name + assay ~ gene, value.var = 'expr')
head(df.dc)

# dissect dfaffy and dfrna
df.affy = df.dc[which(df.dc$assay=='Affy'),c(1,3:5)] %>% na.omit()
rownames(df.affy) = NULL
head(df.affy)
df.rnaseq = df.dc[which(df.dc$assay=='RNAseq'),c(1,3:5)] %>% na.omit()
rownames(df.rnaseq) = NULL
head(df.rnaseq)
## start with df.rna
# merge cell info into df.rna
df.rnaseq = merge(df.affy, df.info)
head(df.rnaseq)

# plot correlation of ace2, tmprss2 and tmprss4
# pairs(df.affy[3:5])
# pairs(df.rnaseq[3:5])
# test = log2(df.affy[3:5])
# pairs(test)

# # plot histogram of tmprss2
# hist(df.rnaseq$tmprss2)

# 2^ normalization
# df.2 = 2^df.rnaseq[,2:4] 
df.2 = df.rnaseq[,2:4] 

rownames(df.2) = make.names(df.rnaseq$CCLE.name, unique = TRUE)
head(df.2)

 # test out the score with assumed a
 b = mean(df.2$tmprss2)/10
 a = 0.5
 FUN = function(x){x[1] * (x[2] + a*x[3] + b)}

df.sc = apply(df.2, 1, FUN)
df = df.sc
df = data.frame(score = df, cell.line = df.rnaseq$Cell.line.primary.name, tissue = df.rnaseq$Site.Primary, source = df.rnaseq$Source) %>% na.omit()
index.ord = order(df$score, decreasing = TRUE)
df.ord = df[index.ord,]
label.cell = rep('', nrow(df))
label.tissue = rep(NA, nrow(df))
label.source = rep(NA, nrow(df))
n = 30
label.cell[1:n] = df.ord$cell.line[1:n] %>% as.character()
label.tissue[1:n] = df.ord$tissue[1:n] %>% as.character()
label.source[1:n] = df.ord$source[1:n] %>% as.character()
# ##
 # label.cell = rep('', nrow(df))
 # label.tissue = rep(NA, nrow(df))
 # calu3 = which(df.ord$cell.line == 'CALU3')
 # caco2 = which(df.ord$cell.line == 'CACO2')
 # index.ca = c(calu3,caco2)
 # label.cell[index.ca] = df.ord$cell.line[index.ca]
 # label.tissue[index.ca] = df.ord$tissue[index.ca]
 # #

 jpeg(filename = paste0('Affy_score_source_log',a,'.jpg'), width = 24, height = 8, units = 'in', res = 300)
 ggplot(data = df.ord, aes(x = cell.line, y = log(score),
                           color = label.tissue, 
                           label = label.cell)) +
   facet_grid(cols = vars(source)) +
   geom_point(size = 1, alpha = 0.7) +
   labs(title = paste0('alpha = ',a),
        x = 'cell lines', y = 'Log2(score)', color = 'tissue') +
   geom_text_repel(show.legend = FALSE,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines")) +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
  dev.off()



alist = (1:100)/100
b = mean(df.2$tmprss2)/10
a = 0.5

df.sc.ac = data.frame(cell.line = df.rnaseq$Cell.line.primary.name, 
  tissue = df.rnaseq$Site.Primary)
head(df.sc.ac)
for (a in alist){
  print(a)
  FUN = function(x){x[1] * (x[2] + a*x[3] + b)}
  arr = apply(df.2, 1, FUN)
  df.sc.ac = cbind(df.sc.ac, arr)
}
rownames(df.sc.ac) = NULL
dim(df.sc.ac)
head(df.sc.ac)
df.sc.ac = df.sc.ac %>% na.omit()
df.ar = t(df.sc.ac[,-c(1:2)])
colnames(df.ar) = df.sc.ac$cell.line
rownames(df.ar) = c(1:100)/100
head(df.ar)

stat.mean = apply(df.ar, 2, mean)
index.stat = order(stat.mean, decreasing = TRUE)
df.sc.ac.od = df.sc.ac[index.stat,1:2]
df.od.final = cbind(df.sc.ac.od, 
                    source = df.rnaseq$Source[index.stat],
                    avg.score  = stat.mean[index.stat])
write.csv(df.od.final, file = 'Affy_Score_in_order.csv' ,row.names = F)



#
# index.stat.top = order(stat.mean, decreasing = TRUE) %>% head(20)
# cell.top = colnames(df.ar)[index.stat.top]
# df.ar.melt = melt(df.ar)
# colnames(df.ar.melt) = c('alpha','cell.line','score')
# col.label = rep(NA, nrow(df.ar.melt))
# index.col.label = which(df.ar.melt$cell.line %in% cell.top)
# col.label[index.col.label] = as.character(df.ar.melt$cell.line[index.col.label])
# col.label = col.label%>% as.factor()
#
# jpeg(filename = 'scores/score_simulation.jpg', width = 10, height = 6, units = 'in', res = 300)
# ggplot(data = df.ar.melt, aes(x = alpha, y = score, color = col.label)) +
# geom_point() +
# theme_minimal() +
# labs(title = 'Simulate Scores to Alpha', color = 'cell lines')
# dev.off()
# 
# 
# 
# 
# # visualize back to dot plot
# cell.index = which(df.rnaseq$cell.line %in% cell.top)
# cell.label = rep(NA, nrow(df.rnaseq))
# cell.label[cell.index] = df.rnaseq$cell.line[cell.index]
# color.label = rep(NA, nrow(df.rnaseq))
# color.label[cell.index] = df.rnaseq$tissue[cell.index]
# 
# 
# p1 = ggplot(data = df.rnaseq, aes(x = ace2, y = tmprss2,
#                            color = color.label, label = cell.label)) +
# geom_point() +
# geom_text_repel(show_guide=FALSE) +
# theme_minimal() +
# theme(legend.position = "none")
# 
# p2 = ggplot(data = df.rnaseq, aes(x = ace2, y = tmprss4,
#                                 color = color.label, label = cell.label)) +
# geom_point() +
# geom_text_repel(show_guide=FALSE) +
# theme_minimal() +
# theme(legend.position = "none")
# 
# p3 = ggplot(data = df.rnaseq, aes(x = tmprss2, y = tmprss4,
#                                 color = color.label, label = cell.label)) +
# geom_point() +
# geom_text_repel(show_guide=FALSE) +
# theme_minimal() +
# labs(color = 'tissue')
# 
# 
# jpeg(filename = 'scores/mapping_paired_expression.jpg', width = 15, height = 10, units = 'in', res = 300)
# p1 + p2 + p3 + guide_area() +
# plot_layout(guides = 'collect') +
# plot_annotation(tag_levels = 'A',
#                 title = 'ACE2/TMPRSS2/TMPRSS4 paired expression',
#                 theme = theme(plot.title = element_text(hjust=0.5)))
# dev.off()



# basically do the same thing on Microarray








