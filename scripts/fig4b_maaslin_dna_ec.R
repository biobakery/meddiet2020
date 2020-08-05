##############################################################################################################################################
# R programs for creating Figure 4a
# 1) Purpose: Associations of Mediterranean diet and its components with DNA enzymes
# 2) Microbial variables: DNA enzymes
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
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(haven)
library(reshape2)
library(Hmisc)

# read in DNA enzyme data
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s4.sas7bdat")
# link cohort ID to aliasid ID then to sample ID
id_alias<-read.csv(file="/udd/nhdow/microbiome/review/data_generated/id_aliasID.csv")
id_alias<-as.data.frame(unique(id_alias[, c(1,3)]))
alias_subject<-read.csv(file="/udd/nhdow/microbiome/review/data_generated/subjectID_aliasid_key.csv")
id_subject<-inner_join(id_alias, alias_subject, by="aliasid")
# combine datasets of different times in to one "wide format" database
transform_data<-function(data1, label1, label2){
  data1<-inner_join(data1, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data2<-as.data.frame(t(data1[,-1912:-1915]))
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
ec_table<-ec_table[-1910:-1911,]
t_relab1<-as.data.frame(t(ec_table))
t_relab1<-t_relab1[, colSums(t_relab1)!=0]
t_relab1<-t_relab1[,order(-colSums(t_relab1))]
# QC DNA enzyme data: orthoganal filtering
dna.ec.qc <- t_relab1[,c(),drop=F]
for (i in seq_len(ncol(t_relab1)))
{
  good <- T
  for (j in seq_len(ncol(dna.ec.qc)))
  {
    if (abs(cor.test(t_relab1 [, i], dna.ec.qc[, j], na.rm=TRUE)$estimate) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    dna.ec.qc <- cbind(dna.ec.qc, t_relab1 [, i])
    colnames(dna.ec.qc)[ncol(dna.ec.qc)] <-
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
  Maaslin2(dna.ec.qc, 
           meta_med, 
           dir,
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="SubjectID", 
           fixed_effects = c(exposure,"calor122cn", "age_fecal", "probio_2mo_qu", "totMETs_paq", "stool_type", "smk", "ant_12mo_qu", "pril12", "metfo12"))
}
maaslin_diet(dir="./ires_emed122ch_ec_dna",   exposure = "ires_emed122ch")
maaslin_diet(dir="./ires_ms122c_ec_dna",     exposure = "ires_ms122c")
maaslin_diet(dir="./ires_fish122c_ec_dna",   exposure = "ires_fish122c")
maaslin_diet(dir="./ires_meat122c_ec_dna",   exposure = "ires_meat122c")
maaslin_diet(dir="./ires_wgrain122c_ec_dna", exposure = "ires_wgrain122c")
maaslin_diet(dir="./ires_legu122c_ec_dna",   exposure = "ires_legu122c")
maaslin_diet(dir="./ires_fruit122c_ec_dna",  exposure = "ires_fruit122c")
maaslin_diet(dir="./ires_veg122c_ec_dna",    exposure = "ires_veg122c")
maaslin_diet(dir="./ires_nut122c_ec_dna",    exposure = "ires_nut122c")
maaslin_diet(dir="./ires_alco122cn_ec_dna",  exposure = "ires_alco122cn")
maaslin_diet(dir="./ires_dairy122c_ec_dna",  exposure = "ires_dairy122c")
# combine resutls from maaslin regression of different dietary variables
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}
sigamed <- sig_results(dir="./ires_emed122ch_ec_dna/significant_results.tsv",   exposure = "ires_emed122ch")
sigmsra <- sig_results(dir="./ires_ms122c_ec_dna/significant_results.tsv",     exposure = "ires_ms122c")
sigfish <- sig_results(dir="./ires_fish122c_ec_dna/significant_results.tsv",   exposure = "ires_fish122c")
sigmeat <- sig_results(dir="./ires_meat122c_ec_dna/significant_results.tsv",   exposure = "ires_meat122c")
sigwhlg <- sig_results(dir="./ires_wgrain122c_ec_dna/significant_results.tsv", exposure = "ires_wgrain122c")
siglegu <- sig_results(dir="./ires_legu122c_ec_dna/significant_results.tsv",   exposure = "ires_legu122c")
sigfrui <- sig_results(dir="./ires_fruit122c_ec_dna/significant_results.tsv",  exposure = "ires_fruit122c")
sigvegg <- sig_results(dir="./ires_veg122c_ec_dna/significant_results.tsv",    exposure = "ires_veg122c")
signuts <- sig_results(dir="./ires_nut122c_ec_dna/significant_results.tsv",    exposure = "ires_nut122c")
sigalco <- sig_results(dir="./ires_alco122cn_ec_dna/significant_results.tsv",  exposure = "ires_alco122cn")
all_results <- function(dir, exposure, label){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}
sigamedall <- all_results(dir="./ires_emed122ch_ec_dna/all_results.tsv",   exposure = "ires_emed122ch", label="MedDiet")
sigmsraall <- all_results(dir="./ires_ms122c_ec_dna/all_results.tsv",     exposure = "ires_ms122c", label="M/S ratio")
sigfishall <- all_results(dir="./ires_fish122c_ec_dna/all_results.tsv",   exposure = "ires_fish122c", label="Fish")
sigmeatall <- all_results(dir="./ires_meat122c_ec_dna/all_results.tsv",   exposure = "ires_meat122c", label="R/P meat")
sigwhlgall <- all_results(dir="./ires_wgrain122c_ec_dna/all_results.tsv", exposure = "ires_wgrain122c", label="Whole grains")
sigleguall <- all_results(dir="./ires_legu122c_ec_dna/all_results.tsv",   exposure = "ires_legu122c", label="Legume")
sigfruiall <- all_results(dir="./ires_fruit122c_ec_dna/all_results.tsv",  exposure = "ires_fruit122c", label="Fruits")
sigveggall <- all_results(dir="./ires_veg122c_ec_dna/all_results.tsv",    exposure = "ires_veg122c", label="Vegetables")
signutsall <- all_results(dir="./ires_nut122c_ec_dna/all_results.tsv",    exposure = "ires_nut122c", label="Nuts")
sigalcoall <- all_results(dir="./ires_alco122cn_ec_dna/all_results.tsv",  exposure = "ires_alco122cn", label="Alcohol")
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
colnames(bind9)[which(colnames(bind9)=="feature")]<-"ec_channing"
# add enzyme names 
ec_label <- read.csv(file = "/udd/nhdow/microbiome/dna_ec_channing_label.csv")
bind9_label<-left_join(bind9, ec_label, by="ec_channing")
write.csv(bind9_label, file="/udd/nhdow/microbiome/review/fig4b_maaslin.csv")
# read in stratified DNA enzyme data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s4.sas7bdat")
# read in enzyme names for stratified DNA enzyme data
str_ec_label<-read.csv(file="/udd/nhdow/microbiome/review/data_generated/stratified_dna_ec_channing_label.csv")
meddiet<-subset(meta_med, select = c(mid, emed122ch))
transform_data<-function(data, ec, label1, label2){
  data_sub<-data[ , grepl( ec , names(data) ) ]
  data_sub$id<-data$id
  data1<-inner_join(data_sub, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data11<-data1[,-c((ncol(data1)-3):ncol(data1))]
  rownames(data11)<-data1$sampleid
  data2<-as.data.frame(t(data11))
  rownames(data2)<-gsub(label2, "", rownames(data2))
  data2<-data2[order(rownames(data2)),]
  return(data2)
}
# creat composition figures for DNA enzymes
stack_area <- function(ec_name, name) {
  all_dna_ec1$ec_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_ec_label, by="ec_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$ec_label
  numex<-nchar(ec_name)
  all_dna_ec1$jud<-substr(all_dna_ec1$ec_label, (numex+1), (numex+1))
  all_dna_ec1<-all_dna_ec1[all_dna_ec1$jud==":",]
  all_dna_ec1<-subset(all_dna_ec1, select = -c(jud, ec_channing, ec_label))
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1, na.rm = T),
                     FUN = "/")
  all_dna_ec1_noun<-all_dna_ec1[!grepl( "unclassified" , rownames(all_dna_ec1) ), ]
  mean_rel<-as.data.frame(rowSums(all_dna_ec1_noun, na.rm = T)/rowSums(all_dna_ec1_noun==0, na.rm = T))
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm=T),
                  FUN = "/")
  colnames(mean_rel)<-"mean_rel"
  dim_count<-mean_rel[mean_rel>0.1]
  num<-length(dim_count)
  all_dna_ec1_noun<-cbind(all_dna_ec1_noun, mean_rel)
  all_dna_ec1_noun<-all_dna_ec1_noun[order(-all_dna_ec1_noun$mean_rel),]
  
  stratum1 <- all_dna_ec1_noun[1:num, -ncol(all_dna_ec1_noun)]
  stratum2 <-
    all_dna_ec1_noun[(num + 1):nrow(all_dna_ec1_noun), -ncol(all_dna_ec1_noun)]
  stratum3<- all_dna_ec1[grepl( "unclassified" , rownames(all_dna_ec1) ), ]
  
  rownames(stratum3)<-"s__Unclassified"
  other <- as.data.frame(t(colSums(stratum2)))
  rownames(other) <- "s__Other_species"
  stratum_cleared <- rbind(stratum1, other, stratum3)
  t_stratum <-
    as.data.frame(t(stratum_cleared))
  
  t_stratum$mid <- rownames(t_stratum)
  t_stratum1 <-inner_join(t_stratum, meddiet, by="mid")
  t_stratum1<-t_stratum1[order(t_stratum1$emed122ch),]
  t_stratum1<-subset(t_stratum1, select=-emed122ch)
  p.stratum <- melt(t_stratum1, id = 'mid')
  p.stratum$variable <- as.factor(p.stratum$variable)
  p.stratum$value <- as.numeric(as.character(p.stratum$value))
  
  p.stratum <- p.stratum %>% separate(variable, c(NA, "bugs"), "s__")
  p.stratum$bugs <- gsub("_", " ", p.stratum$bugs)
  p.stratum$cat<-name
  return(p.stratum)
}

t_s1_sub<-transform_data(s1, "ecd2_7_1_45", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd2_7_1_45", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd2_7_1_45", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd2_7_1_45", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p1<-stack_area(ec_name = "2.7.1.45", name="EC 2.7.1.45: 2-dehydro-3-deoxygluconokinase")

t_s1_sub<-transform_data(s1, "ecd5_4_2_8", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd5_4_2_8", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd5_4_2_8", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd5_4_2_8", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p2<-stack_area(ec_name = "5.4.2.8", name="EC 5.4.2.8: phosphomannomutase")

t_s1_sub<-transform_data(s1, "ecd1_2_7_1", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd1_2_7_1", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd1_2_7_1", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd1_2_7_1", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p3<-stack_area(ec_name = "1.2.7.1", name="EC 1.2.7.1: pyruvate:ferredoxin oxidoreductase")

t_s1_sub<-transform_data(s1, "ecd3_2_1_85", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd3_2_1_85", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd3_2_1_85", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd3_2_1_85", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p4<-stack_area(ec_name = "3.2.1.85", name="EC 3.2.1.85: 6-phospho-beta-galactosidase")

t_s1_sub<-transform_data(s1, "ecd5_5_1_1", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd5_5_1_1", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd5_5_1_1", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd5_5_1_1", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p5<-stack_area(ec_name = "5.5.1.1", name="EC 5.5.1.1: muconate cycloisomerase")
# 1.17.99.5 renamed to 1.17.98.1
t_s1_sub<-transform_data(s1, "ecd1_17_99_5", "_SF05", "_s1")
t_s2_sub<-transform_data(s2, "ecd1_17_99_5", "_SF06", "_s2")
t_s3_sub<-transform_data(s3, "ecd1_17_99_5", "_SF07", "_s3")
t_s4_sub<-transform_data(s4, "ecd1_17_99_5", "_SF08", "_s4")
all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)

p6<-stack_area(ec_name = "1.17.99.5", name="EC 1.17.98.1: bile-acid 7alpha-dehydroxylase")

all_ec_composition<-rbind(p1, p2, p3, p4, p5, p6)
all_ec_composition$bugs<-capitalize(all_ec_composition$bugs)
all_ec_composition$value[all_ec_composition$value=="NaN"]<-0
all_ec_composition$cat_f<-factor(all_ec_composition$cat, levels = c("EC 2.7.1.45: 2-dehydro-3-deoxygluconokinase",
                                        "EC 5.4.2.8: phosphomannomutase",
                                        "EC 1.2.7.1: pyruvate:ferredoxin oxidoreductase",
                                        "EC 3.2.1.85: 6-phospho-beta-galactosidase",
                                        "EC 5.5.1.1: muconate cycloisomerase",
                                        "EC 1.17.98.1: bile-acid 7alpha-dehydroxylase"))
bugs<-unique(all_ec_composition$bugs)
color<-c('Faecalibacterium prausnitzii'='#48A462',
         'Other species'='#ffffbf',
         'Unclassified'='#d3d3d3',
         'Lachnospiraceae bacterium 5 1 63faa'='#F88A89',
         'Roseburia intestinalis'='#F06C45',
         'Anaerostipes hadrus'='#FE870D',
         'Methanobrevibacter smithii'='#B294C7',
         'Collinsella aerofaciens'='#b20058',
         'Bifidobacterium adolescentis'='#FF7F00')

png(file="/udd/nhdow/microbiome/review/fig4b_composition.png",width=1000,height=3000, pointsize=50)
ggplot(all_ec_composition,
       aes(
         x = mid,
         y = value,
         fill = bugs,
         group = bugs
       )) +
  geom_area(position="fill") + theme_bw() + 
  scale_fill_manual(values = color, na.value = "white")+
  facet_wrap(. ~ cat_f, ncol=1) +
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(color="black", size=45),
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
dev.off()

# read in unstratified DNA enzyme data for creating scatter plots
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_unstratified_raw_s4.sas7bdat")
# link cohort ID to aliasid ID then to sample ID
id_alias<-read.csv(file="/udd/nhdow/microbiome/id_aliasID.csv")
id_alias<-as.data.frame(unique(id_alias[, c(1,3)]))
alias_subject<-read.csv(file="/udd/nhdow/microbiome/subjectID_aliasid_key.csv")
id_subject<-inner_join(id_alias, alias_subject, by="aliasid")
# combine datasets of different times in to one "wide format" database
transform_data<-function(data1, label1, label2){
  data1<-inner_join(data1, id_subject, by="id")
  data1$sampleid<-paste0(data1$SubjectID, label1)
  data2<-as.data.frame(t(data1[,-1912:-1915]))
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
# add pseudo count of 1 for plotting purpose only
ec_table<-ec_table[-1910:-1911,]+1
ec_table <-
  sweep(ec_table,
        2,
        STATS = colSums(ec_table, na.rm = T),
        FUN = "/")
t_relab1<-as.data.frame(t(ec_table))
t_relab1<-t_relab1[, colSums(t_relab1)!=0]
t_relab1$mid<-rownames(t_relab1)
tax_meta_combined<-inner_join(t_relab1, meta_med, by="mid")

scatter_plot<-function(diet, bug, title, color, xtitle){
  pp<-ggplot(tax_meta_combined, aes(x=diet, y=log10(bug))) + 
    geom_point(col=color, alpha=0.5, size=5)+
    geom_smooth(method=lm, size=3)+
    ggtitle(title)+
    xlab(xtitle)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ theme_classic()+
    theme(plot.title = element_text(size=20, face = "bold"),
          axis.title.x = element_text(size=20),
          axis.title.y= element_blank(),
          axis.text.y=element_text(size=20,),
          axis.text.x = element_text(size=20))
  return(pp)
}

# D-fructuronate degradation
p1<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd2_7_1_45, 
                 "FDR P=0.03", "olivedrab4", "MedDiet")
# mannan degradation
p2<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd5_4_2_8, 
                 "FDR P=0.13", "olivedrab4", "MedDiet")
# pyruvate fermentation to acetate and lactate II
p3<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd1_2_7_1, 
                 "FDR P=0.13", "olivedrab4", "MedDiet")
# lactose and galactose degradation I
p4<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd3_2_1_85, 
                 "FDR P=0.05", "olivedrab4", "MedDiet")
# lignin degradation
p5<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd5_5_1_1, 
                 "FDR P=0.0008", "olivedrab4", "MedDiet")
# bile acid metabolism
p6<-scatter_plot(tax_meta_combined$ires_emed122ch, 
                 tax_meta_combined$ecd1_17_99_5, 
                 "FDR P=0.0008", "olivedrab4", "MedDiet")
png(file="/udd/nhdow/microbiome/revision/fig4b_scatter_plot1.png",width=500,height=1500, pointsize=50)
plot_grid(p1, p2, p3, align = "h", ncol =  1, nrow=3)
dev.off()
png(file="/udd/nhdow/microbiome/revision/fig4b_scatter_plot2.png",width=500,height=1500, pointsize=50)
plot_grid(p4, p5, p6, align = "h", ncol =  1, nrow=3)
dev.off()
