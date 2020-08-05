######################################################################################################
# R programs for creating PCoA plots in Figures 3a, S3, & S4
# 1) Purposes: Biogeography of microbial communities 
#              Associations of Mediterranean diet with overall structure of microbial communities
# 2) Microbial variables: beta-diversity matrix of microbial communities
# 3) Dietary exposures: Mediterranean diet index
######################################################################################################
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)

# read in taxonomy data
tax_rpk_name <-   read.table(    file = '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/taxonomy/bugs_dna_929_unFilt.tsv',
                                 sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)

# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]

# read in metadata
meta_med<-read.csv(file="/udd/nhdow/microbiome/meta_ffq_8612.csv")
# id list for those with dietary data
med_id<-meta_med %>% select(mid) 
ttax_rpk <- as.data.frame(t(tax_rpk_species))
ttax_rpk$mid<-rownames(ttax_rpk)
# keep participants with both dietary and microbiome data
ttax_rpk<-inner_join(med_id, ttax_rpk, by="mid")
rownames(ttax_rpk)<-ttax_rpk$mid
ttax_rpk <- subset(ttax_rpk, select = -mid)
# calculate relative abundance
ttax_rpk_rel <-
  sweep(ttax_rpk,
        1,
        STATS = rowSums(ttax_rpk),
        FUN = "/")
# calculate bray-crutis
ttax_rpk_rel_bc <- vegdist(ttax_rpk_rel, "bray")
# PCOA using cmdscale
mod <- cmdscale(ttax_rpk_rel_bc, eig = T)
pcoap <- data.frame(mod$points)
pcoap$mid<-rownames(pcoap)
pcoap<-pcoap[order(pcoap$mid),]
# percentages of variation explained by PCO1 & 2
mod$eig[1]/sum(mod$eig)
mod$eig[2]/sum(mod$eig)
# create dataset for PCo scores
med_score<-subset(meta_med, select = c(mid, emed122ch))
pcoap<-inner_join(pcoap, med_score, by="mid")
write.csv(pcoap, file="/udd/nhdow/microbiome/review/data_generated/pcoa_species.csv")
# PCOA figure colored by meddiet score
png(
  file = "/udd/nhdow/microbiome/review/fig3a_pcoa_bray_curtis_sample.png",
  width = 2400,
  height = 2175,
  pointsize = 50
)
ggplot(data = pcoap, aes(X1, X2, colour = pcoap$emed122ch)) +
  geom_point(size=22, alpha=0.8) + 
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(
    legend.title = element_text(size = 80),
    legend.text = element_text(size = 80),
    axis.title.x = element_text(size = 80),
    axis.title.y = element_text(size = 80),
    axis.text.y = element_text(size = 80),
    axis.text.x = element_text(size = 80),
    axis.line = element_line(color = "Black", size=3)
  ) +
  scale_color_viridis(name = "MedDiet \nscore", trans = "reverse")
dev.off()

######################################################################################################
# only keep phylum-level feature
tax_rpk_phylum <- subset(tax_rpk_name,!is.na(phylum) & is.na(class) & is.na(order) & is.na(family) & is.na(genus) & is.na(species) & is.na(strain))
rownames(tax_rpk_phylum)<-tax_rpk_phylum$phylum
tax_rpk_phylum<-tax_rpk_phylum[,-1:-9]
tax_rpk_phylum <-
  sweep(tax_rpk_phylum,
        2,
        STATS = colSums(tax_rpk_phylum),
        FUN = "/")

# only keep Bacteroidetes & Firmicutes 
tax_rpk_phylum_fnb <- tax_rpk_phylum[c("p__Bacteroidetes", "p__Firmicutes"),] 
ttax_rpk_phylum_fnb <- as.data.frame(t(tax_rpk_phylum_fnb))
ttax_rpk_phylum_fnb$mid<-rownames(ttax_rpk_phylum_fnb)
ttax_rpk_phylum_fnb <-
  inner_join(ttax_rpk_phylum_fnb, med_id, by = "mid")
# merge with species-level PCo scores
pcoa_fnb<-inner_join(pcoap, ttax_rpk_phylum_fnb, by="mid")

# create PCOA colored by Firmicutes
p1<-ggplot(data = pcoa_fnb, aes(X1, X2)) +
  geom_point(aes(colour =pcoa_fnb$p__Firmicutes), size=15, alpha=0.7) + 
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +theme_classic()+
  ggtitle("Firmicutes")+
  scale_color_distiller(palette='PuRd', name = "Relative \nabundance", trans = "reverse")+
  theme(legend.title = element_text(size = 80),
        legend.text = element_text(size = 60),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3),
        title =element_text(size=100, face='bold')
  ) 
# create PCOA colored by Bacteroidetes
p2<-ggplot(data = pcoa_fnb, aes(X1, X2)) +
  geom_point(aes(colour =p__Bacteroidetes), size=15, alpha=0.7) + 
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +theme_classic()+
  ggtitle("Bacteroidetes")+
  scale_color_distiller(palette='PuRd', name = "Relative \nabundance", trans = "reverse")+
  theme(legend.title = element_text(size = 80),
        legend.text = element_text(size = 60),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3),
        title =element_text(size=100, face='bold')
  ) 

png(
  file = "/udd/nhdow/microbiome/review/sfig1_pcoa_phylum.png",
  width = 3200,
  height = 1450,
  pointsize = 50
)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")
dev.off()

#######################################################################################################################
# create PCOA figures colored by most abundant species
species_4_plot<-ttax_rpk_rel[,order(-colMeans(ttax_rpk_rel))]
species_4_plot$mid<-rownames(species_4_plot)

pcoap_sp <- function(num) {
  sp<-species_4_plot[, c(num, ncol(species_4_plot))]
  pcoap_sp<-inner_join(sp, pcoap, by="mid")
  a <- colnames(sp)[1]
  name <- gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", substring(a, 4))
  ggplot(data = pcoap_sp, aes(X1, X2)) +
    geom_point(aes(colour = pcoap_sp[,1]), size = 6, alpha=0.7) +  
    ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +  theme_classic() +
    ggtitle(paste0(name)) +
    scale_color_viridis(name = "Relative \nabundance", trans = "reverse") +
    theme(
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.line = element_line(color = "Black", size = 3),
      title = element_text(size = 25, face = 'bold'),
      plot.title = element_text(face = "bold.italic")
    )
}
p1<-pcoap_sp(num=1)
p2<-pcoap_sp(num=2)
p3<-pcoap_sp(num=3)
p4<-pcoap_sp(num=4)
p5<-pcoap_sp(num=5)
p6<-pcoap_sp(num=6)
p7<-pcoap_sp(num=7)
p8<-pcoap_sp(num=8)
p9<-pcoap_sp(num=9)
p10<-pcoap_sp(num=10)
p11<-pcoap_sp(num=11)
p12<-pcoap_sp(num=12)

pp1<-ggarrange(p1,
               p2,
               p3,
               p4,
               ncol=4, nrow=1, common.legend = TRUE, legend="right")

pp2<-ggarrange(p5,
               p6,
               p7,
               p8, ncol=4, nrow=1, common.legend = TRUE, legend="right")

pp3<-ggarrange(p9,
               p10,
               p11,
               p12, ncol=4, nrow=1, common.legend = TRUE, legend="right")

png(file="/udd/nhdow/microbiome/review/sfig2_pcoa_12species_combined.png",  width = 2400,
    height = 1800, pointsize=50)
ggarrange(pp1, pp2, pp3, ncol=1, nrow=3)
dev.off()