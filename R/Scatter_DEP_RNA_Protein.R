library(ggplot2)
library(hrbrthemes)
#install.packages("hrbrthemes")
library(ggpubr)
# library(ggpmisc)
library(ggrepel)

matrix <- read.csv('/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/FGFR2_DEP_pro_phos_scater.csv', header = T)
colnames(matrix)
matrix$Significant <- ifelse(matrix$display == 'display', "display", "not_display")
# matrix$size <- ifelse(matrix$display == 'display', 2, 1)
ggplot(data = matrix, aes(x = logFC_Pro, y = logFC_phos,color=pathway,size = size)) +
  geom_point() + ## aes(color = Significant)
  scale_size_continuous(range = c(1,2))+
  scale_color_brewer(palette="Paired") +
  # scale_color_manual(values = c("#515151", "red","#7AA9CE")) + # ("gray", "#7AA9CE") c( "#EA686B","#7AA9CE","gray"))
  theme_bw(base_size = 12) + theme(legend.position = "right") +
  geom_text_repel(
    data = subset(matrix, Significant == "display"),
    aes(label = phosphosite),
    size = 5,
    max.overlaps = Inf,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines",
    )
  )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("FGFR2 alteration on RNA and Protein")+
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 16,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12)) +
        # panel.grid.major = element_blank(),   #不显示网格线
        # panel.grid.minor = element_blank())+  #不显示网格线
  geom_hline(yintercept = -0.35,color = 'gray', linetype="dashed") +
  geom_hline(yintercept = 0.35,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = -0.25,color = 'gray', linetype="dashed") +
  geom_vline(xintercept = 0.25,color = 'gray', linetype="dashed") +
  # scale_x_continuous(limits = c(-1.8, 1.8))+
  # scale_y_continuous(limits = c(-1.8, 1.8))
  scale_x_continuous(breaks=seq(-2, 2, 0.5)) +
  scale_y_continuous(breaks=seq(-2.5, 2.5, 1)) +
  ylab("Log2FC (FGFR2/others)-RNA")+xlab("Log2FC (FGFR2/others)-protein")  #设置x轴和y轴的标题



matrix$Significant <- ifelse(matrix$display == 'display', "display", "not_display")
p1= ggplot(data=matrix, aes(x=logFC_pro, y=logFC_RNA))+    #aes(x=Protein, y=Phosphosite)分别代表log2FC
  geom_point(aes(colour = Significant))+
  stat_smooth(method="lm",se=TRUE)+
  xlim(-1.8,1.8)+        #删除则默认适合的坐标轴
  ylim(-1.8,1.8)+        #删除则默认适合的坐标轴
  theme_bw() +         #加边框
  scale_color_manual(values = c('blue',"gray",'#C42D43',"#F9813A","#1B96BC"))+
  geom_text_repel(
    data = subset(matrix, Significant == "display"),
    aes(label = gene),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines") )+
  stat_cor(data=matrix, method = "spearman")   #spearman比pearson相关性好
#stat_cor(data=dat, method = "pearson")

#以下三选一即可：出图方式一
p1 + geom_text_repel(
  data = subset(matrix, matrix$display == 'display'),
  aes(label = gene),   #label = phosphosite代表基因标签
  size = 6,
  max.overlaps = Inf,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines")
) +
  theme(axis.text.x=element_text(angle=0,hjust = 0.6,colour="black",family="Arial",size=16), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.title.x=element_text(angle=0,hjust = 0.5,colour="black",family="Arial",size=16),
        axis.text.y=element_text(family="Arial",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 16,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.8), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=12),
        legend.title=element_text(face="plain", family="Arial", colour="black", #设置图例的总标题的字体属性
                                  size=12))+
        # panel.grid.major = element_blank(),   #不显示网格线
        # panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Log2FC (FGFR2/others)-RNA")+xlab("Log2FC (FGFR2/others)-protein")  #设置x轴和y轴的标题


#以下三选一即可：出图方式二
p1 + geom_label_repel(
  data = subset(matrix, matrix$Significant == 'display'),
  aes(label = Gene,
      fill = factor(color)),
  color = 'red', size = 3.5) +
  theme(legend.position = "right")  

#以下三选一即可：出图方式三
#scale_fill_manual(values = c('#AFAFAF','#3BA236','#CC3333'))
p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),axis.line = element_line(colour = "black"))
#geom_text(aes(label = matrix$X))



matrix <- read.csv('/Users/ranpeng/Desktop/CCA/Data/Fig3/RNA/RNA_protein_scatter.csv', header = T)
colnames(matrix)
matrix$color <- ifelse(matrix$RNA > 0.5 & matrix$Protein > 0.5 ,'both_up',
                       ifelse(matrix$RNA < -0.5 & matrix$Protein < -0.5 ,'both_down',
                              ifelse(matrix$RNA > 0.5 & matrix$Protein < -0.5 ,'RNA_Up',
                                     ifelse(matrix$RNA < -0.5 & matrix$Protein > 0.5 ,'Protein_Up', "not_sig"))))

# matrix$Significant <- ifelse(matrix$RNA > 5 & matrix$Protein > 5 ,'both_up',
#                        ifelse(matrix$RNA < -5 & matrix$Protein < -5 ,'both_down',
#                               ifelse(matrix$RNA > 5 & matrix$Protein < -5 ,'RNA_Up',
#                                      ifelse(matrix$RNA < -5 & matrix$Protein > 5 ,'Protein_Up', "not_sig"))))

matrix$Significant <- ifelse(matrix$RNA > 2 & matrix$Protein > 2|matrix$RNA < -2 & matrix$Protein < -1 |matrix$RNA > 1 & matrix$Protein < -1|matrix$RNA < -1 & matrix$Protein > 1
                             ,'display',"not_display")


write.csv(matrix,"/Users/ranpeng/Desktop/CCA/Data/Fig3/RNA/RNA_protein_scatter_1.csv")                             



