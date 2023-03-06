# Venn diagrams (Kallisto, DeSeq2)
```R
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done')
library(ggvenn)

# load contig names
MF_ccdc_DeSeq2 = read.table("MF_Kallisto_ccdc_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_ccdc_DeSeq2)[1] <- "gene"
MF_dmrt1L_DeSeq2 = read.table("MF_Kallisto_dmrt1L_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_dmrt1L_DeSeq2)[1] <- "gene"
MF_dmrt1S_DeSeq2 = read.table("MF_Kallisto_dmrt1S_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_dmrt1S_DeSeq2)[1] <- "gene"
wtko_dmw_DeSeq2 = read.table("wtko_Kallisto_dmw_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtko_dmw_DeSeq2)[1] <- "gene"
wtko_scan_DeSeq2 = read.table("wtko_Kallisto_scan_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtko_scan_DeSeq2)[1] <- "gene"
wtko_ccdc_DeSeq2 = read.table("wtko_Kallisto_ccdc_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtko_ccdc_DeSeq2)[1] <- "gene"


# Make list for comparison
x <- list(MF_1=MF_ccdc_DeSeq2$gene, 
          MF_2=MF_dmrt1L_DeSeq2$gene, 
          MF_3=MF_dmrt1S_DeSeq2$gene,
          dmw = wtko_dmw_DeSeq2$gene)

# plot
pdf("dmw_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot1 <- ggvenn(x, 
       show_elements = F,
       label_sep = "\n", 
       text_size = 3,
       fill_color = c("black","grey70", "grey80", "grey90"),
       set_name_size = 3,
       show_percentage = F)
dev.off()


# Make list for comparison
x <- list(MF_1=MF_ccdc_DeSeq2$gene, 
          MF_2=MF_dmrt1L_DeSeq2$gene, 
          MF_3=MF_dmrt1S_DeSeq2$gene,
          scan = wtko_scan_DeSeq2$gene)

# plot

pdf("scan_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot2 <- ggvenn(x, 
                show_elements = F,
                label_sep = "\n", 
                text_size = 3,
                fill_color = c("black","grey70", "grey80", "grey90"),
                set_name_size = 3,
                show_percentage = F)
dev.off()

# Make list for comparison
x <- list(MF_1=MF_ccdc_DeSeq2$gene, 
          MF_2=MF_dmrt1L_DeSeq2$gene, 
          MF_3=MF_dmrt1S_DeSeq2$gene,
          ccdc = wtko_ccdc_DeSeq2$gene)

# plot
pdf("ccdc_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot3 <- ggvenn(x, 
       show_elements = F,
       label_sep = "\n", 
       text_size = 3,
       fill_color = c("black","grey70", "grey80", "grey90"),
       set_name_size = 3,
       show_percentage = F)
dev.off()






# Make list for comparison
x <- list(dmw=wtko_dmw_DeSeq2$gene,
          scan = wtko_scan_DeSeq2$gene,
          ccdc=wtko_ccdc_DeSeq2$gene)

# plot
#pdf("dmw_scan_ccdc_wtko_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot4 <- ggvenn(x, 
       show_elements = F,
       label_sep = "\n", 
       text_size = 3,
       fill_color = c("black","grey70", "grey80", "grey90"),
       set_name_size = 3,
       show_percentage = F)
#dev.off()
library(gridExtra)
#merge all three plots within one grid (and visualize this)
t <- textGrob("")
g <- grid.arrange(plot1, plot2, plot3,t, plot4, 
                  nrow=2, ncol=3,
                  top = "Kallisto and DeSeq2") #arranges plots within grid

#save
ggsave(file="Kallisto_DeSeq2_Venn.pdf", w=8, h=4, g)
```
