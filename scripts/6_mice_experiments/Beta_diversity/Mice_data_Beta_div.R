set.seed(12345)
library(phyloseq)
library("ggplot2")
library(ggpubr)
library(gridExtra)


###### Read the phyloseq object ######

physeq<-readRDS("~/Downloads/physeq_sv_decontam_Qfemto.rds")
physeq_genus <- tax_glom(physeq, taxrank = 'Genus')


###### Bay Curtis distance ######

dist<-phyloseq::distance(physeq_genus,"bray")

sampledf <- data.frame(sample_data(physeq_genus))
Adonis<-adonis2(dist ~ time_enterotype, data = sampledf)
Adonis

###### PCoA Plot ######

ord <- ordinate(physeq_genus, method = "PCoA", distance = "bray")
p <- plot_ordination(physeq_genus, ord, color = "time_enterotype")
p <- p + theme_bw() + ggtitle(paste("Bray distance (Adonis p=0.001)"))
p<-p + scale_colour_manual(values = c("darkgreen","brown4","darkolivegreen4","red3","springgreen3" ,"tomato1" ))
p <- p + stat_ellipse() + theme_bw()
p <- p +theme(legend.title=element_blank(),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15),axis.text.y = element_text(size=13),axis.text.x = element_text(size=13),legend.text = element_text(size=15))
p<-p+theme(plot.title = element_text(size=18))

p<-p+theme(legend.spacing.y = unit(0.1, 'mm'))


pdf("Beta_time_enterotype.pdf")
print(p)
dev.off()