##############################################################################################################################################
# R programs for creating Figure 4a
# 1) Purpose: Associations of Mediterranean diet and its components with DNA pathways
# 2) Microbial variables: DNA pathways
# 3) Dietary variables: Main exposure: Mediterranean diet index
#                       Secondary exposures: fruit, vegetables, whole grains, monounsaturated fat/saturated fat ratio, 
#                                            red meat, alcohol, legume, nuts, & fish
# 4) Covariates: age, total energy intake, physical activity, smoking status, antibiotic use, 
#                probiotic use, stool type, PPI use, metformin use
#############################################################################################################################################

library(tidyverse)
library(Maaslin2)
library(grid)
library(ggpubr)
library(haven)
# read in DNA pathway data
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s4.sas7bdat")
# link cohort ID to aliasid ID then to sample ID
id_alias<-read.csv(file="/udd/nhdow/microbiome/id_aliasID.csv")
id_alias<-as.data.frame(unique(id_alias[, c(1,3)]))
alias_subject<-read.csv(file="/udd/nhdow/microbiome/subjectID_aliasid_key.csv")
id_subject<-inner_join(id_alias, alias_subject, by="aliasid")
# combine datasets of different times in to one "wide format" database
transform_data<-function(data1, label1, label2){
  data1<-inner_join(data1, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data2<-as.data.frame(t(data1[,-(ncol(data1)-3):-(ncol(data1))]))
  colnames(data2)<-data1$sampleid
  rownames(data2)<-gsub(label2, "", rownames(data2))
  data2<-data2[order(rownames(data2)),]
  return(data2)
}
t_s1<-transform_data(s1, "_SF05", "_s1")
t_s2<-transform_data(s2, "_SF06", "_s2")
t_s3<-transform_data(s3, "_SF07", "_s3")
t_s4<-transform_data(s4, "_SF08", "_s4")
ec_table<-cbind(t_s1, t_s2, t_s3, t_s4)
t_relab1<-as.data.frame(t(ec_table))
t_relab1<-t_relab1[, colSums(t_relab1)!=0]
t_relab1<-t_relab1[,order(-colSums(t_relab1))]
# QC DNA pathway data: orthoganal filtering
dna.pwy.qc <- t_relab1[,c(),drop=F]

for (i in seq_len(ncol(t_relab1)))
{
  good <- T
  for (j in seq_len(ncol(dna.pwy.qc)))
  {
    if (abs(cor.test(t_relab1 [, i], dna.pwy.qc[, j], na.rm=TRUE)$estimate) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    dna.pwy.qc <- cbind(dna.pwy.qc, t_relab1 [, i])
    colnames(dna.pwy.qc)[ncol(dna.pwy.qc)] <-
      colnames(t_relab1)[i]
  }
}

# read in meatdata
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
meta_med<-transform_diet(var1=meta_med$dairy122c, name="ires_dairy122c")
# maaslin regressions using maaslin default normalization and transformation (NORMALIZATION = TSS, TRANSFORM = LOG)
setwd("/udd/nhdow/microbiome/maaslin/")
maaslin_diet<-function(dir, exposure){
  Maaslin2(dna.pwy.qc, 
           meta_med, 
           dir,
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="SubjectID", 
           fixed_effects = c(exposure,"calor122cn", "age_fecal", "probio_2mo_qu", "totMETs_paq", "stool_type", "smk", "ant_12mo_qu", "pril12", "metfo12"))
}
maaslin_diet(dir="./ires_emed122ch_pwy_dna",   exposure = "ires_emed122ch")
maaslin_diet(dir="./ires_ms122c_pwy_dna",     exposure = "ires_ms122c")
maaslin_diet(dir="./ires_fish122c_pwy_dna",   exposure = "ires_fish122c")
maaslin_diet(dir="./ires_meat122c_pwy_dna",   exposure = "ires_meat122c")
maaslin_diet(dir="./ires_wgrain122c_pwy_dna", exposure = "ires_wgrain122c")
maaslin_diet(dir="./ires_legu122c_pwy_dna",   exposure = "ires_legu122c")
maaslin_diet(dir="./ires_fruit122c_pwy_dna",  exposure = "ires_fruit122c")
maaslin_diet(dir="./ires_veg122c_pwy_dna",    exposure = "ires_veg122c")
maaslin_diet(dir="./ires_nut122c_pwy_dna",    exposure = "ires_nut122c")
maaslin_diet(dir="./ires_alco122cn_pwy_dna",  exposure = "ires_alco122cn")
maaslin_diet(dir="./ires_dairy122c_pwy_dna",  exposure = "ires_dairy122c")
# combine resutls from maaslin regressions of different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}
sigamed <- sig_results(dir="./ires_emed122ch_pwy_dna/significant_results.tsv",   exposure = "ires_emed122ch")
sigmsra <- sig_results(dir="./ires_ms122c_pwy_dna/significant_results.tsv",     exposure = "ires_ms122c")
sigfish <- sig_results(dir="./ires_fish122c_pwy_dna/significant_results.tsv",   exposure = "ires_fish122c")
sigmeat <- sig_results(dir="./ires_meat122c_pwy_dna/significant_results.tsv",   exposure = "ires_meat122c")
sigwhlg <- sig_results(dir="./ires_wgrain122c_pwy_dna/significant_results.tsv", exposure = "ires_wgrain122c")
siglegu <- sig_results(dir="./ires_legu122c_pwy_dna/significant_results.tsv",   exposure = "ires_legu122c")
sigfrui <- sig_results(dir="./ires_fruit122c_pwy_dna/significant_results.tsv",  exposure = "ires_fruit122c")
sigvegg <- sig_results(dir="./ires_veg122c_pwy_dna/significant_results.tsv",    exposure = "ires_veg122c")
signuts <- sig_results(dir="./ires_nut122c_pwy_dna/significant_results.tsv",    exposure = "ires_nut122c")
sigalco <- sig_results(dir="./ires_alco122cn_pwy_dna/significant_results.tsv",  exposure = "ires_alco122cn")
all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}
sigamedall <- all_results(dir="./ires_emed122ch_pwy_dna/all_results.tsv",   exposure = "ires_emed122ch", label="MedDiet")
sigmsraall <- all_results(dir="./ires_ms122c_pwy_dna/all_results.tsv",     exposure = "ires_ms122c", label="M/S ratio")
sigfishall <- all_results(dir="./ires_fish122c_pwy_dna/all_results.tsv",   exposure = "ires_fish122c", label="Fish")
sigmeatall <- all_results(dir="./ires_meat122c_pwy_dna/all_results.tsv",   exposure = "ires_meat122c", label="R/P meat")
sigwhlgall <- all_results(dir="./ires_wgrain122c_pwy_dna/all_results.tsv", exposure = "ires_wgrain122c", label="Whole grains")
sigleguall <- all_results(dir="./ires_legu122c_pwy_dna/all_results.tsv",   exposure = "ires_legu122c", label="Legume")
sigfruiall <- all_results(dir="./ires_fruit122c_pwy_dna/all_results.tsv",  exposure = "ires_fruit122c", label="Fruits")
sigveggall <- all_results(dir="./ires_veg122c_pwy_dna/all_results.tsv",    exposure = "ires_veg122c", label="Vegetables")
signutsall <- all_results(dir="./ires_nut122c_pwy_dna/all_results.tsv",    exposure = "ires_nut122c", label="Nuts")
sigalcoall <- all_results(dir="./ires_alco122cn_pwy_dna/all_results.tsv",  exposure = "ires_alco122cn", label="Alcohol")
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
bind9<-rbind(join1, join2, join3, join4, join5, join6,
             join7, join8, join9, join10)
colnames(bind9)[which(colnames(bind9)=="feature")]<-"pwy_channing"
# add pathway names 
pwy_label<-read.csv(file="/udd/nhdow/microbiome/dna_pwy_channing_label.csv")
bind9_label<-left_join(bind9, pwy_label, by="pwy_channing")
write.csv(bind9_label, file="/udd/nhdow/microbiome/revision/data_generated/fig4a_maaslin.csv")
# only keep pathways used in plotting
bind9display<-subset(bind9_label, pwy_name=="d-fructuronate degradation"|pwy_name=="mannan degradation" |
                                  pwy_name=="glutaryl-coa degradation"|pwy_name=="lactose and galactose degradation I" | 
                                  pwy_name=="pyruvate fermentation to acetate and lactate II" )
level_x_order <- factor(bind9display$meta, level = c('MedDiet', 'Whole grains', 'Vegetables', 'Fruits', 
                                                     'Nuts', 'Legume', 'Fish', 'R/P meat', 'M/S ratio', 'Alcohol'))
level_y_order <- factor(bind9display$pwy_name, levels =c('lactose and galactose degradation I',
                                                         'glutaryl-coa degradation',
                                                         'pyruvate fermentation to acetate and lactate II',
                                                         'mannan degradation',
                                                         'd-fructuronate degradation'))
bind9display$starss <- cut(bind9display$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
# create heatmap for diet-pathway associations
png(file="/udd/nhdow/microbiome/review/fig4a_dna_pwy.png",width=2800,height=800, pointsize=50)
ggplot(bind9display, aes(level_x_order, level_y_order)) +
       geom_tile(aes(fill = coef), color = "white") +
       scale_fill_gradient2(low = "red", high = "blue", space = "Lab" ) +
       geom_text(aes(label=starss), color="black", size=16, show.legend = TRUE) +
       theme(legend.title = element_text(size = 50),
             legend.text = element_text(size = 50),
             legend.position = "right",
             plot.title = element_text(size=50),
             axis.title=element_blank(),
             axis.title.y= element_blank(),
             axis.text.y=element_text(size=50),
             axis.text.x = element_text(angle = 90, hjust = 1, size=50)) +
  labs(fill = "Beta coefficient")
dev.off()
