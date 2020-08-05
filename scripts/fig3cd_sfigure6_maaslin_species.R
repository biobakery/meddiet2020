##############################################################################################################################################
# R programs for creating Figures 3c, 3d & S6
# 1) Purpose: Associations of Mediterranean diet and its components with microbial species
# 2) Microbial variables: Microbial species
# 3) Dietary variables: Main exposure: Mediterranean diet index
#                       Secondary exposures: fruit, vegetables, whole grains, monounsaturated fat/saturated fat ratio, 
#                                            red meat, alcohol, legume, nuts, & fish
# 4) Covariates: age, total energy intake, physical activity, smoking status, antibiotic use, 
#                probiotic use, stool type, PPI use, metformin use
#############################################################################################################################################
library(Maaslin2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(cowplot)
# read in taxonomy data
tax_rpk_name <-   read.table(file = '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/taxonomy/bugs_dna_929_unFilt.tsv',
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
rownames(meta_med)<-meta_med$mid
meta_med$age_fecal<-as.numeric(as.character(meta_med$age_fecal))
meta_med$calor122cn<-as.numeric(as.character(meta_med$calor122cn))
meta_med$totMETs_paq<-as.numeric(as.character(meta_med$totMETs_paq))
# standardize dietary intake data
transform_diet <- function(var1, name) {
  var1 <- as.numeric(as.character(var1))
  i_var <-
    qnorm((rank(var1, na.last = "keep") - 0.5) / sum(!is.na(var1)))
  ires_var <-
    as.data.frame(lm(i_var ~ meta_med$calor122cn)$residuals)
  colnames(ires_var) <- name
  new.df <- cbind(meta_med[!is.na(var1),], ires_var)
  return(new.df)
}
meta_med<-transform_diet(var1=meta_med$emed122ch, name="ires_emed122ch")
meta_med<-transform_diet(var1=meta_med$ms122c, name="ires_ms122c")
meta_med<-transform_diet(var1=meta_med$fish122c, name="ires_fish122c")
meta_med<-transform_diet(var1=meta_med$meat122c, name="ires_meat122c")
meta_med<-transform_diet(var1=meta_med$wgrain122c, name="ires_wgrain122c")
meta_med<-transform_diet(var1=meta_med$legu122c, name="ires_legu122c")
meta_med<-transform_diet(var1=meta_med$fruit122c, name="ires_fruit122c")
meta_med<-transform_diet(var1=meta_med$veg122c, name="ires_veg122c")
meta_med<-transform_diet(var1=meta_med$nut122c, name="ires_nut122c")
meta_med<-transform_diet(var1=meta_med$alco122cn, name="ires_alco122cn")
# maaslin regressions using maaslin default normalization and transformation (NORMALIZATION = TSS, TRANSFORM = LOG)
maaslin_diet<-function(dir, exposure){
  Maaslin2(tax_rpk_species, 
           meta_med, 
           dir,
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="SubjectID", 
           fixed_effects = c(exposure,"calor122cn", "age_fecal", "probio_2mo_qu", "totMETs_paq", "stool_type", "smk", "ant_12mo_qu", "pril12", "metfo12"))
}
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_emed122ch_taxonomy",   exposure = "ires_emed122ch")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_ms122c_taxonomy",     exposure = "ires_ms122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_fish122c_taxonomy",   exposure = "ires_fish122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_meat122c_taxonomy",   exposure = "ires_meat122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_wgrain122c_taxonomy", exposure = "ires_wgrain122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_legu122c_taxonomy",   exposure = "ires_legu122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_fruit122c_taxonomy",  exposure = "ires_fruit122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_veg122c_taxonomy",    exposure = "ires_veg122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_nut122c_taxonomy",    exposure = "ires_nut122c")
maaslin_diet(dir="/udd/nhdow/microbiome/maaslin/ires_alco122cn_taxonomy",  exposure = "ires_alco122cn")
# combine results from maasline regressions across different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}
sigamed <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_emed122ch_taxonomy/significant_results.tsv",   exposure = "ires_emed122ch")
sigmsra <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_ms122c_taxonomy/significant_results.tsv",     exposure = "ires_ms122c")
sigfish <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_fish122c_taxonomy/significant_results.tsv",   exposure = "ires_fish122c")
sigmeat <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_meat122c_taxonomy/significant_results.tsv",   exposure = "ires_meat122c")
sigwhlg <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_wgrain122c_taxonomy/significant_results.tsv", exposure = "ires_wgrain122c")
siglegu <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_legu122c_taxonomy/significant_results.tsv",   exposure = "ires_legu122c")
sigfrui <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_fruit122c_taxonomy/significant_results.tsv",  exposure = "ires_fruit122c")
sigvegg <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_veg122c_taxonomy/significant_results.tsv",    exposure = "ires_veg122c")
signuts <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_nut122c_taxonomy/significant_results.tsv",    exposure = "ires_nut122c")
sigalco <- sig_results(dir="/udd/nhdow/microbiome/maaslin/ires_alco122cn_taxonomy/significant_results.tsv",  exposure = "ires_alco122cn")
all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}
sigamedall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_emed122ch_taxonomy/all_results.tsv",   exposure = "ires_emed122ch", label="MedDiet")
sigmsraall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_ms122c_taxonomy/all_results.tsv",     exposure = "ires_ms122c", label="M/S ratio")
sigfishall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_fish122c_taxonomy/all_results.tsv",   exposure = "ires_fish122c", label="Fish")
sigmeatall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_meat122c_taxonomy/all_results.tsv",   exposure = "ires_meat122c", label="R/P meat")
sigwhlgall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_wgrain122c_taxonomy/all_results.tsv", exposure = "ires_wgrain122c", label="Whole grains")
sigleguall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_legu122c_taxonomy/all_results.tsv",   exposure = "ires_legu122c", label="Legume")
sigfruiall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_fruit122c_taxonomy/all_results.tsv",  exposure = "ires_fruit122c", label="Fruits")
sigveggall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_veg122c_taxonomy/all_results.tsv",    exposure = "ires_veg122c", label="Vegetables")
signutsall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_nut122c_taxonomy/all_results.tsv",    exposure = "ires_nut122c", label="Nuts")
sigalcoall <- all_results(dir="/udd/nhdow/microbiome/maaslin/ires_alco122cn_taxonomy/all_results.tsv",  exposure = "ires_alco122cn", label="Alcohol")
join1<-full_join(sigamed, sigwhlg, by="feature")
join2<-full_join(join1, sigfrui, by="feature")
join3<-full_join(join2, sigvegg, by="feature")
join4<-full_join(join3, siglegu, by="feature")
join5<-full_join(join4, signuts, by="feature")
join6<-full_join(join5, sigmeat, by="feature")
join7<-full_join(join6, sigfish, by="feature")
join8<-full_join(join7, sigmsra, by="feature")
join9<-full_join(join8, sigalco, by="feature")
sigfeature<-subset(join9, select=feature)
join1<-left_join(sigfeature, sigamedall, by="feature")
join2<-left_join(sigfeature, sigwhlgall, by="feature")
join3<-left_join(sigfeature, sigveggall, by="feature")
join4<-left_join(sigfeature, sigleguall, by="feature")
join5<-left_join(sigfeature, signutsall, by="feature")
join6<-left_join(sigfeature, sigmeatall, by="feature")
join7<-left_join(sigfeature, sigfishall, by="feature")
join8<-left_join(sigfeature, sigmsraall, by="feature")
join9<-left_join(sigfeature, sigalcoall, by="feature")
join10<-left_join(sigfeature, sigfruiall, by="feature")
bind9<-rbind(join1, join2, join3, join4, join5, join6, join7, join8, join9, join10)
# read in table of species matched to phyla
phylum<-read.csv(file = '/udd/nhdow/microbiome/phylum_species.csv', 
                 check.names=FALSE)
phylum<-phylum[phylum$species!="" & phylum$strain=="",]
phylum<-subset(phylum, select=-strain)
colnames(phylum)[which(colnames(phylum) == 'species')] <- 'feature'
# format specie & phylum names
bind9<-left_join(bind9, phylum, by="feature")
bind9$feature <- substring(bind9$feature, 4)
bind9$phylum <- substring(bind9$phylum, 4)
bind9<-bind9[order(bind9$meta, bind9$phylum, bind9$feature),]
bind9$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind9$feature)
bug_names<-bind9[1:(nrow(bind9)/10), (ncol(bind9)-1):ncol(bind9)]
# create heatmap for diet-taxonomy associations
level_x_order <- factor(bind9$meta, level = c('MedDiet', 'Whole grains', 'Vegetables', 'Fruits', 
                                              'Nuts', 'Legume', 'Fish', 'R/P meat', 'M/S ratio', 'Alcohol'))
level_y_order <- factor(bind9$name, levels = bug_names$name)
bind9$stars <- cut(bind9$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
p1<-ggplot(bind9, aes(level_x_order, level_y_order)) +
           geom_tile(aes(fill = coef), color = "white") +
           scale_fill_gradient2(low = "red", high = "blue", space = "Lab" ) +
           geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
           xlab("MedDiet Score & Components") +
           theme(legend.title = element_text(size = 50),
                 legend.text = element_text(size = 50),
                 legend.position = "left",
                 plot.title = element_text(size=50),
                 axis.title=element_text(size=50,face="bold"),
                 axis.title.y= element_blank(),
                 axis.text.y=element_text(size=50, face="bold.italic"),
                 axis.text.x = element_text(angle = 90, hjust = 1, size=50, face="bold")) +
  labs(fill = "Relative Abundance")
# create left bar for categorizing species to phyla
bind9$category<-as.factor(bind9$phylum)
p2<-ggplot(bind9, aes(level_x_order, level_y_order))+
           scale_y_discrete(position = "right")+
           geom_tile(aes(fill = category)) +
           xlab("MedDiet Score & Components") +
           theme(legend.title = element_text(size = 50),
                 legend.text = element_text(size = 50),
                 axis.title.y= element_blank(),
                 axis.title=element_text(size=50,face="bold"),
                 axis.text.x = element_text(angle = 90, hjust = 1, size=50)) +
  labs(fill = "Phyla")
png(file="/udd/nhdow/microbiome/review/figS6_maaslin_taxonomy.png",width=4500,height=3000, pointsize=50)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()
# create scatter plots for figure 3d
tax_rpk<-tax_rpk_species+1
tax_rpk_rel <-
  sweep(tax_rpk,
        2,
        STATS = colSums(tax_rpk),
        FUN = "/")
tax_rpk_rel<-log10(tax_rpk_rel)
ttax_rel <- as.data.frame(t(tax_rpk_rel))
ttax_rel<- ttax_rel[,order(-colMeans(ttax_rel))]
ttax_rel$mid<-rownames(ttax_rel)
ttax_rel<-inner_join(meta_med, ttax_rel, by="mid")
rownames(ttax_rel)<-ttax_rel$mid

scatter_plot<-function(diet, bug, color){
  pp<-ggplot(ttax_rel, aes(x=diet, y=bug)) + 
    geom_point(col="#4490ba", alpha=0.5, size=15)+
    geom_smooth(method=lm, size=5, color=color)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ theme_classic()+
    theme(axis.title= element_blank(),
          axis.text.y=element_text(size=50,),
          axis.text.x = element_text(size=50))
  return(pp)
}
p1<-scatter_plot(ttax_rel$ires_emed122ch, 
                 ttax_rel$s__Faecalibacterium_prausnitzii, "#E69F00")
p2<-scatter_plot(ttax_rel$ires_emed122ch, 
                 ttax_rel$s__Eubacterium_eligens, "#E69F00")
p3<-scatter_plot(ttax_rel$ires_emed122ch, 
                 ttax_rel$s__Bacteroides_cellulosilyticus, "#E69F00")
p4<-scatter_plot(ttax_rel$ires_emed122ch, 
                 ttax_rel$s__Collinsella_aerofaciens, "#E69F00")
p5<-scatter_plot(ttax_rel$ires_meat122c, 
                 ttax_rel$s__Faecalibacterium_prausnitzii, "#E69F00")
p6<-scatter_plot(ttax_rel$ires_meat122c, 
                 ttax_rel$s__Eubacterium_eligens, "#E69F00")
p7<-scatter_plot(ttax_rel$ires_meat122c, 
                 ttax_rel$s__Bacteroides_cellulosilyticus, "#E69F00")
p8<-scatter_plot(ttax_rel$ires_meat122c, 
                 ttax_rel$s__Collinsella_aerofaciens, "#E69F00")
png(file="/udd/nhdow/microbiome/revision/fig3d_scatter1.png",width=2400,height=1125, pointsize=50)
plot_grid(p1, p5, align = "h", ncol =  2, nrow=1)
dev.off()

png(file="/udd/nhdow/microbiome/revision/fig3d_scatter2.png",width=2400,height=1125, pointsize=50)
plot_grid(p2, p6, align = "h", ncol =  2, nrow=1)
dev.off()

png(file="/udd/nhdow/microbiome/revision/fig3d_scatter3.png",width=2400,height=1125, pointsize=50)
plot_grid(p3, p7, align = "h", ncol =  2, nrow=1)
dev.off()

png(file="/udd/nhdow/microbiome/revision/fig3d_scatter4.png",width=2400,height=1125, pointsize=50)
plot_grid(p4, p8, align = "h", ncol =  2, nrow=1)
dev.off()

