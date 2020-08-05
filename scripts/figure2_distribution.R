##############################################################################################################################################
# R programs for creating Figure 2
# 1) Purposes: Descriptive analysis for distributions of diet, microbial species, 
#             functional potential (DNA enzymes and pathways) and activity (RNA/DNA ratio of enzymes and pathways)
# 2) Variables: Mediterranean diet and its components, top 10 most abundant microbial species, 
#               top 5 most abundant DNA enzymes & pathways and top 5 most abundant species that contribute to each functional feature
#               top 5 most transcribed RNA enzymes & pathways and top 3 most abundant species that contribute to each functional feature
#############################################################################################################################################
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(diverse)
library(grid)
library(ggpubr)
library(viridis)
library(gridExtra)
library(cowplot)
library(haven)

################################ Distribution of MedDiet and its components ######################################
#read in dietary and covariate data
meta_med<-read.csv(file="/udd/nhdow/microbiome/meta_ffq_8612.csv")
rownames(meta_med)<-meta_med$mid
#keep only MedDiet and its components
med_distribution<-meta_med[order(meta_med$emed122ch), c("emed122ch", "ms122c", "fish122c", "meat122c", "wgrain122c", "legu122c", "fruit122c", "veg122c", "nut122c","alco122cn")]
colnames(med_distribution)<-c('Mediterranean diet index', 'Monounsaturated fat/saturated fat ratio', 'Fish', 'Red/processed meat', 
                              'Whole grains', 'Legume', 'Fruits', 'Vegetables', 
                              'Nuts', 'Alcohol')
med_distribution$mid<-rownames(med_distribution)
# standardize dietary data
for (i in 1:(ncol(med_distribution)-1)) {
  df<-med_distribution[,c(i,ncol(med_distribution))]
  colnames(df)<-c("intake", "mid")
  df$intake <- qnorm((rank(df$intake, na.last = "keep") - 0.5) / sum(!is.na(df$intake)))
  df$cat<-colnames(med_distribution)[i]
  assign(paste0("component",i), df )
}

meddiet<-rbind(component1, component2, component3, component4, component5, component6, component7, component8, component9, component10)
level_x_order <- factor(meddiet$mid, levels = component1$mid)
level_y_order <- factor(meddiet$cat, levels = c('Alcohol', 'Monounsaturated fat/saturated fat ratio', 'Red/processed meat', 
                                                'Fish', 'Legume', 'Nuts', 'Fruits', 
                                                'Vegetables', 'Whole grains', 'Mediterranean diet index'))
#create plot for distributions of MedDiet and its components
png(file="/udd/nhdow/microbiome/review/fig2a_diet.png",width=4000,height=1600, pointsize=50)
ggplot(meddiet, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = intake)) +
  scale_fill_viridis_c(name ="Intake Level") +
  theme(legend.title = element_text(size = 80),
        legend.text = element_text(size = 80),
        legend.position = "left",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=100),
        axis.text.x = element_blank())
dev.off()

################################ distributions of species #######################################################################################
#read in taxonomy data & keep only species-level features
tax_rpk_name <-   read.table(    file = '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/taxonomy/bugs_dna_929_unFilt.tsv',
                                 sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))

rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]
# ID list of 925 participants who have both diet and gut microbiome data
med_id<-meta_med %>% select(mid, emed122ch) 
# add pseudo-count of 1 only for plotting purpose
tax_rpk_species<-tax_rpk_species+1
# calculate relative abundance
tax_rpk_rel <-
  sweep(tax_rpk_species,
        2,
        STATS = colSums(tax_rpk_species),
        FUN = "/")
ttax_rel <- as.data.frame(t(tax_rpk_rel))
ttax_rel<- ttax_rel[,order(-colMeans(ttax_rel))]
ttax_rel$mid<-rownames(ttax_rel)
ttax_rel<-inner_join(med_id, ttax_rel, by="mid")
rownames(ttax_rel)<-ttax_rel$mid
# sort samples by MedDiet index
ttax_rel<-ttax_rel[order(ttax_rel$emed122ch),]
# keep top 10 most abundant species
for (i in 3:12) {
  df<-ttax_rel[,c(1,i)]
  colnames(df)<-c("mid", "species")
  df$cat<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", substring(colnames(ttax_rel)[i], 4))
  assign(paste0("sp",i), df )
}
spall<-rbind(sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10, sp11, sp12)

tax_y_order <- factor(spall$cat, levels = c( "Bacteroides stercoris",  "Ruminococcus bromii", "Bacteroides vulgatus", 
                                             "Eubacterium siraeum", "Alistipes putredinis", "Prevotella copri",
                                             "Bacteroides uniformis", "Subdoligranulum unclassified",  "Faecalibacterium prausnitzii", "Eubacterium rectale"))
tax_x_order <-
  factor(spall$mid, levels = ttax_rel$mid)
spall$lgsp<-log10(spall$species)

png(file="/udd/nhdow/microbiome/review/fig2b_species.png",width=4000,height=1600, pointsize=50)
ggplot(spall, aes(tax_x_order, tax_y_order)) +
       geom_tile(aes(fill = lgsp)) +
       scale_fill_gradientn(colours = c("black", "#fc4e2a", "#fed976", "#ffffcc"), name ="Relative abundance\nlog10 scale", na.value = "white") +
       theme(legend.title = element_text(size = 80),
             legend.text = element_text(size = 80),
             legend.position = "left",
             plot.title = element_blank(),
             axis.title=element_blank(),
             axis.text.y=element_text(size=100, face = "italic"),
             axis.text.x = element_blank())
dev.off()

#################################### Distribution of DNA ECs ##########################################################
# read in unstratified DNA EC data
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
ec_table<-ec_table[-1910:-1911,]
ec_table$ec_channing<-rownames(ec_table)
# add name of EC to database
ec_label <- read.csv(file = "/udd/nhdow/microbiome/dna_ec_channing_label.csv")
ec_table<-left_join(ec_table, ec_label, by="ec_channing")
rownames(ec_table)<-ec_table$ec_label

t_dna_ec1<-as.data.frame(t(ec_table[,-c(926:930)]))
# add pseudo count of 1 only for plotting purpose
t_dna_ec1<-t_dna_ec1+1
t_dna_ec1 <-
  sweep(t_dna_ec1,
        1,
        STATS = rowSums(t_dna_ec1),
        FUN = "/")
t_dna_ec1 <- t_dna_ec1[, order(-colMeans(t_dna_ec1))]
dna_ec2 <- t_dna_ec1[,c(),drop=F]
# orthogonal filtering
for (i in seq_len(ncol(t_dna_ec1)))
{
  good <- T
  for (j in seq_len(ncol(dna_ec2)))
  {
    if (abs(cor(t_dna_ec1 [, i], dna_ec2[, j])) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    dna_ec2 <- cbind(dna_ec2, t_dna_ec1 [, i])
    colnames(dna_ec2)[ncol(dna_ec2)] <-
      colnames(t_dna_ec1)[i]
  }
}

dna_ec_rel<-dna_ec2[,order(-colMeans(dna_ec2))]
dna_ec_rel$mid<-rownames(dna_ec_rel)
dna_ec_rel<-inner_join(med_id, dna_ec_rel, by="mid")
rownames(dna_ec_rel)<-dna_ec_rel$mid
# sort samples by MedDiet index
dna_ec_rel<-dna_ec_rel[order(dna_ec_rel$emed122ch),]
# keep top 5 most abundant DNA ECs
for (i in 3:12) {
  df<-dna_ec_rel[,c(1,i)]
  colnames(df)<-c("mid", "level")
  df$cat<-colnames(dna_ec_rel)[i]
  assign(paste0("dnaec",i), df )
}
dna_ec_all<-rbind(dnaec3, dnaec4, dnaec5, dnaec6, dnaec7)
ec_label<-unique(dna_ec_all[,3])
dna_ec_y_order <- factor(dna_ec_all$cat, levels = rev(ec_label))
dna_ec_x_order <- factor(dna_ec_all$mid, levels = dna_ec_rel$mid)
dna_ec_all$lglevel<-log10(dna_ec_all$level)

# read in stratified DNA EC data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna/ecd_stratified_raw_s4.sas7bdat")

str_ec_label<-read.csv(file="/udd/nhdow/microbiome/stratified_dna_ec_channing_label.csv")
ec_labe_st<-as.data.frame(ec_label)
ec_label <- read.csv(file = "/udd/nhdow/microbiome/dna_ec_channing_label.csv")
ec_labe_st<-left_join(ec_labe_st, ec_label, by="ec_label")

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
# create data for stack plots for top 5 most abundant DNA EC
stack_area <- function(ec_name, num) {
  all_dna_ec1$ec_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_ec_label, by="ec_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$ec_label
  all_dna_ec1<-all_dna_ec1[, -c(926,927)]
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1),
                     FUN = "/")
  mean_rel<-as.data.frame(rowSums(all_dna_ec1)/rowSums(all_dna_ec1!=0))
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm=T),
                  FUN = "/")
  colnames(mean_rel)<-"mean_rel"
  mean_rel$bug_name<-rownames(mean_rel)
  stratum <- mean_rel[order(-mean_rel$mean_rel), ]
  stratum <- stratum %>% separate(bug_name, c(NA, "bugs"), "s__")
  stratum$bugs[is.na(stratum$bugs)]<-"Unclassified"
  stratum_cleared <- stratum[1:num, ]
  stratum_cleared$bugs <- gsub("_", " ", stratum_cleared$bugs)
  stratum_cleared$bugs<-paste(toupper(substr(stratum_cleared$bugs, 1, 1)), substr(stratum_cleared$bugs, 2, nchar(stratum_cleared$bugs)), sep="")
  stratum_cleared$ec<-ec_name
  return(stratum_cleared)
}

for (i in 1:5){
  if(i==1){
    ec_name<-as.character(ec_labe_st[i,1])
    ec<-as.character(ec_labe_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 5)
    bugall_ec<-bug1
  } else if (i>1){  
    ec_name<-as.character(ec_labe_st[i,1])
    ec<-as.character(ec_labe_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 5)
    bugall_ec<-rbind(bugall_ec, bug1)
  }
}
########################  Distribution of DNA pathways  #####################################################
# read in unstratified DNA pathway data
s1<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s1.sas7bdat")
s2<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s2.sas7bdat")
s3<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s3.sas7bdat")
s4<-read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_unstratified_raw_s4.sas7bdat")
# combind databases of different time points to one database with "wide format"
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
pwy_table<-cbind(t_s1, t_s2, t_s3, t_s4)
# add pseudo-count of 1 for plotting purpose only
pwy_table<-pwy_table+1
pwy_table<-  sweep(pwy_table,
                   2,
                   STATS = colSums(pwy_table),
                   FUN = "/")
pwy_label<-read.csv(file="/udd/nhdow/microbiome/dna_pwy_channing_label.csv")
pwy_table$pwy_channing<-rownames(pwy_table)
pwy_table<-left_join(pwy_table, pwy_label, by="pwy_channing")
rownames(pwy_table)<-paste(toupper(substr(pwy_table$pwy_name, 1, 1)), substr(pwy_table$pwy_name, 2, nchar(as.character(pwy_table$pwy_name))), sep="")

t_relab1<-as.data.frame(t(pwy_table[, -c(926,927)]))
t_relab1<-t_relab1[,order(-colSums(t_relab1))]
t_relab1<-t_relab1[, colMeans(t_relab1)!=0]
# orthogonal filtering
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
dna_pwy_rel<-dna.pwy.qc[,order(-colMeans(dna.pwy.qc))]
dna_pwy_rel$mid<-rownames(dna_pwy_rel)
dna_pwy_rel<-inner_join(med_id, dna_pwy_rel, by="mid")
rownames(dna_pwy_rel)<-dna_pwy_rel$mid
# sort by MedDiet index
dna_pwy_rel<-dna_pwy_rel[order(dna_pwy_rel$emed122ch),]
# keep top 5 most abundant DNA pathways
for (i in 3:12) {
  df<-dna_pwy_rel[,c(1,i)]
  colnames(df)<-c("mid", "level")
  df$cat<-colnames(dna_pwy_rel)[i]
  assign(paste0("dnapwy",i), df )
}
dna_pwy_all<-rbind(dnapwy3, dnapwy4, dnapwy5, dnapwy6, dnapwy7)
dna_pwy_label<-as.data.frame(unique(dna_pwy_all$cat))
colnames(dna_pwy_label)<-"pwy_name"
dna_pwy_y_order <- factor(dna_pwy_all$cat, levels = rev(dna_pwy_label$pwy_name))
dna_pwy_x_order <- factor(dna_pwy_all$mid, levels = dna_pwy_rel$mid)
dna_pwy_all$lglevel<-log10(dna_pwy_all$level)

pwy_label$pwy_name<-paste(toupper(substr(pwy_label$pwy_name, 1, 1)), substr(pwy_label$pwy_name, 2, nchar(as.character(pwy_label$pwy_name))), sep="")
dna_pwy_label_st<-left_join(dna_pwy_label, pwy_label, by="pwy_name")

# read in stratified DNA pathway data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpd_stratified_raw_s4.sas7bdat")

str_ec_label<-read.csv(file="/udd/nhdow/microbiome/stratified_dna_pwy_channing_label.csv")

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
# creat data for stack plots for top 5 most abundant DNA pathways
stack_area <- function(ec_name, num) {
  all_dna_ec1$pwy_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_ec_label, by="pwy_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$pwy_name
  all_dna_ec1<-all_dna_ec1[, -c(926,927)]
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1),
                     FUN = "/")
  mean_rel<-as.data.frame(rowSums(all_dna_ec1, na.rm = T)/rowSums(all_dna_ec1!=0, na.rm = T))
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm = T),
                  FUN = "/")
  colnames(mean_rel)<-"mean_rel"
  mean_rel$bug_name<-rownames(mean_rel)
  stratum <- mean_rel[order(-mean_rel$mean_rel), ]
  stratum <- stratum %>% separate(bug_name, c(NA, "bugs"), "s__")
  stratum$bugs[is.na(stratum$bugs)]<-"Unclassified"
  stratum_cleared <- stratum[1:num, ]
  stratum_cleared$bugs <- gsub("_", " ", stratum_cleared$bugs)
  stratum_cleared$bugs<-paste(toupper(substr(stratum_cleared$bugs, 1, 1)), substr(stratum_cleared$bugs, 2, nchar(stratum_cleared$bugs)), sep="")
  stratum_cleared$ec<-ec_name
  return(stratum_cleared)
}

for (i in 1:5){
  if(i==1){
    ec_name<-as.character(dna_pwy_label_st[i,1])
    ec<-as.character(dna_pwy_label_st[i,2]) 
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 5)
    bugall_pwy<-bug1
  } else if (i>1){  
    ec_name<-as.character(dna_pwy_label_st[i,1])
    ec<-as.character(dna_pwy_label_st[i,2]) 
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 5)
    bugall_pwy<-rbind(bugall_pwy, bug1)
  }
}
# bug color code
bug_color_value<-c('Bacteroides uniformis'='#A6CEE3',
                   'Bacteroides vulgatus'='#5FA0CA',
                   'Bacteroides cellulosilyticus'='#257CB2',
                   'Unclassified'='#D3D3D3',
                   'Bacteroides plebeius'='#A5D981',
                   'Faecalibacterium prausnitzii'='#63B84F',
                   'Methanobrevibacter smithii'='#4F9F3B',
                   'Roseburia intestinalis'='#B89B74',
                   'Alistipes onderdonkii'='#F68181',
                   'Alistipes shahii'='#E93E3F',
                   'Eubacterium siraeum'='#E9412F',
                   'Acidaminococcus sp bv3l6'='#F6975B',
                   'Synergistes sp 3 1 syn1'='#FDAC4F',
                   'Akkermansia muciniphila'='#FE8B15',
                   'Bacteroides coprocola'='#ED8F47',
                   'Escherichia coli'='#D1AAB7',
                   'Prevotella copri'='#A585BF',
                   'Eubacterium rectale'='#73489F',
                   'Prevotella buccae'='#A99099',
                   'Prevotella denticola'='#F7F599',
                   'Anaerococcus obesiensis'='#D9AF63',
                   'Streptococcus lutetiensis'='#B15928')
# barplot of DNA EC contributions by species
level_ec <- factor(bugall_ec$ec, levels = rev(ec_labe_st$ec_label))
dna_ec_bar<-ggplot(bugall_ec, aes(x = level_ec, y = mean_rel, fill = bugs)) +
                   geom_bar(stat = "identity")+theme_bw()+ylim(0,1)+
                   scale_fill_manual(values = bug_color_value)+ coord_flip()+
                   theme(legend.position = "none",
                         plot.title = element_blank(),
                         axis.title=element_blank(),
                         axis.text.y=element_blank(),
                         axis.text.x = element_text(size = 50))
# barplot of DNA pathway contributions by species
level_pwy <- factor(bugall_pwy$ec, levels = rev(dna_pwy_label$pwy_name))
dna_pwy_bar<-ggplot(bugall_pwy, aes(x = level_pwy, y = mean_rel, fill = bugs)) +
                    geom_bar(stat = "identity")+theme_bw()+ylim(0,1)+
                    scale_fill_manual(values = bug_color_value)+ coord_flip()+
                    theme(legend.position = "none",
                         plot.title = element_blank(),
                         axis.title=element_blank(),
                         axis.text.y=element_blank(),
                         axis.text.x = element_text(size = 50))
# plot for DNA EC distribution
dna_ec_plot<-ggplot(dna_ec_all, aes(dna_ec_x_order, dna_ec_y_order)) +
                    geom_tile(aes(fill =lglevel)) +
                    scale_fill_gradientn(colours = c("black", "#b2182b", "#d6604d", 
                                                     "#f4a582", "#fddbc7", "#92c5de", "#2166ac"), 
                    name ="Relative abundance\nlog10 scale", na.value = "white")+
                    theme(legend.position = "left",
                          legend.text = element_text(size=30),
                          legend.title = element_text(size=30),
                          plot.title = element_blank(),
                          axis.title=element_blank(),
                          axis.text.y=element_text(size=100),
                          axis.text.x = element_blank())
# plot for DNA pathway distribution
dna_pwy_plot<-ggplot(dna_pwy_all, aes(dna_pwy_x_order, dna_pwy_y_order)) +
                     geom_tile(aes(fill =lglevel)) +
                     scale_fill_gradientn(colours = c("black", "#b2182b", "#d6604d", 
                                                      "#f4a582", "#fddbc7", "#92c5de", "#2166ac"), 
                     name ="Relative abundance\nlog10 scale", na.value = "white")+
                     theme(legend.position = "left",
                           legend.text = element_text(size=30),
                           legend.title = element_text(size=30),
                           plot.title = element_blank(),
                           axis.title=element_blank(),
                           axis.text.y=element_text(size=100),
                           axis.text.x = element_blank())

png(file="/udd/nhdow/microbiome/review/fig2d_dna_pwy.png",width=7000,height=800, pointsize=50)
plot_grid(dna_pwy_plot, dna_pwy_bar, align = "h", ncol =  2, rel_widths = c(6/7, 1/7))
dev.off()

png(file="/udd/nhdow/microbiome/review/fig2d_dna_ec.png",width=7000,height=800, pointsize=50)
plot_grid(dna_ec_plot, dna_ec_bar, align = "h", ncol =  2, rel_widths = c(6/7, 1/7))
dev.off()
# create a common legend for barplots
png(file="/udd/nhdow/microbiome/review/fig2d_bug_dna_legend.png",width=3500,height=2000, pointsize=50)
ggplot(bug_all, aes(x = ec, y = mean_rel, fill = bugs)) +
       geom_bar(stat = "identity")+theme_bw()+ ylim(0,0.8)+
       scale_fill_manual(values = bug_color_value)+ coord_flip()+
       theme(legend.position = "right",
             legend.text = element_text(size=50, face = "italic"),
             legend.title = element_blank(),
             plot.title = element_blank(),
             axis.title=element_blank(),
             axis.text.y=element_blank(),
             axis.text.x = element_text(size = 30))
dev.off()
########################  Distribution of RNA/DNA ratio for ECs ###########################################################################################
# read in RNA/DNA ratio for ECs
ec_rnadna<- read.csv(file = "/udd/nhdow/microbiome/rna-dna-relativeExpression_ec_noBUGS.csv", 
                     header = TRUE, check.names=FALSE, row.names = 1)

# orthogonal filtering
ec_rnadna<-ec_rnadna[,order(-colMeans(ec_rnadna))]
ratio_ec <- ec_rnadna[,c(),drop=F]

for (i in seq_len(ncol(ec_rnadna)))
{
  good <- T
  for (j in seq_len(ncol(ratio_ec)))
  {
    if (abs(cor.test(ec_rnadna [, i], ratio_ec[, j])$estimate) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    ratio_ec <- cbind(ratio_ec, ec_rnadna [, i])
    colnames(ratio_ec)[ncol(ratio_ec)] <-
      colnames(ec_rnadna)[i]
  }
}
ratio_ec<-ratio_ec[,order(-colMeans(ratio_ec))]
ratio_ec$mid<-rownames(ratio_ec)
ratio_ec<-left_join(med_id, ratio_ec, by="mid")
rownames(ratio_ec)<-ratio_ec$mid
ratio_ec<-ratio_ec[order(ratio_ec$emed122ch),]
# keep top 5 most transcribed ECs
for (i in 3:7) {
  df<-ratio_ec[,c(1,i)]
  colnames(df)<-c("mid", "level")
  df$cat<-colnames(ratio_ec)[i]
  assign(paste0("ratioec",i), df )
}
ratioec6$cat<-"1.3.8.1: short-chain acyl-CoA dehydrogenase"
ratio_ec_all<-rbind(ratioec3, ratioec4, ratioec5, ratioec6, ratioec7)

ratio_ec_label<-as.data.frame(unique(ratio_ec_all[,3]))
colnames(ratio_ec_label)<-"ec_label"

ratio_ec_y_order <- factor(ratio_ec_all$cat, levels = rev(ratio_ec_label$ec_label))
ratio_ec_x_order <- factor(ratio_ec_all$mid, levels = ratio_ec$mid)
ratio_ec_all$lglevel<-log2(ratio_ec_all$level)
# read in stratified RNA data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s2.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s3.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecrna/ecr_stratified_raw_s4.sas7bdat")

str_rna_ec_label<-read.csv(file="/udd/nhdow/microbiome/stratified_rna_ec_channing_label.csv")
rna_ec_label <- read.csv(file = "/udd/nhdow/microbiome/rna_ec_channing_label.csv")
rna_ec_label_st<-left_join(ratio_ec_label, rna_ec_label, by="ec_label")

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
# create data for stack plots for RNA EC contributions by species
stack_area <- function(ec_name, num) {
  all_dna_ec1$ec_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_rna_ec_label, by="ec_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$ec_label
  all_dna_ec1<-all_dna_ec1[, -c(373,374)]
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1),
                     FUN = "/")
  mean_rel<-as.data.frame(rowSums(all_dna_ec1, na.rm = T)/rowSums(all_dna_ec1!=0, na.rm = T))
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm = T),
                  FUN = "/")
  
  colnames(mean_rel)<-"mean_rel"
  mean_rel$bug_name<-rownames(mean_rel)
  stratum <- mean_rel[order(-mean_rel$mean_rel), ]
  stratum <- stratum %>% separate(bug_name, c(NA, "bugs"), "s__")
  stratum$bugs[is.na(stratum$bugs)]<-"Unclassified"
  stratum_cleared <- stratum[1:num, ]
  stratum_cleared$bugs <- gsub("_", " ", stratum_cleared$bugs)
  stratum_cleared$bugs<-paste(toupper(substr(stratum_cleared$bugs, 1, 1)), substr(stratum_cleared$bugs, 2, nchar(stratum_cleared$bugs)), sep="")
  stratum_cleared$ec<-ec_name
  return(stratum_cleared)
}

for (i in 1:5){
  if(i==1){
    ec_name<-as.character(rna_ec_label_st[i,1])
    ec<-as.character(rna_ec_label_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 3)
    bugall_ratioec<-bug1
  } else if (i>1){  
    ec_name<-as.character(rna_ec_label_st[i,1])
    ec<-as.character(rna_ec_label_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 3)
    bugall_ratioec<-rbind(bugall_ratioec, bug1)
  }
}
bugall_ratioec$ec[bugall_ratioec$ec=="1.3.99.2: transferred entry: 1.3.8.1"]<-"1.3.8.1: short-chain acyl-CoA dehydrogenase"

#################################  Distribution of RNA/DNA ratio of pathways ###################################
# read in unstratified 
pwy_rnadna <-read.csv(file = "/udd/nhdow/microbiome/rna-dna-relativeExpression_pwy_noBUGS.csv", 
                     header = TRUE, check.names=FALSE, row.names = 1)

# orthogonal filtering
pwy_rnadna<-pwy_rnadna[,order(-colMeans(pwy_rnadna))]
ratio_pwy <- pwy_rnadna[,c(),drop=F]

for (i in seq_len(ncol(pwy_rnadna)))
{
  good <- T
  for (j in seq_len(ncol(ratio_pwy)))
  {
    if (abs(cor.test(pwy_rnadna [, i], ratio_pwy[, j])$estimate) > 0.9)
    {
      good <- F
    }
  }
  if (good)
  {
    ratio_pwy <- cbind(ratio_pwy, pwy_rnadna [, i])
    colnames(ratio_pwy)[ncol(ratio_pwy)] <-
      colnames(pwy_rnadna)[i]
  }
}

ratio_pwy<-ratio_pwy[, order(-colMeans(ratio_pwy))]
ratio_pwy$mid<-rownames(ratio_pwy)
ratio_pwy<-left_join(med_id, ratio_pwy, by="mid")
rownames(ratio_pwy)<-ratio_pwy$mid
ratio_pwy<-ratio_pwy[order(ratio_pwy$emed122ch),]
# keep top 5 most transcribed pathways
for (i in 3:7) {
  df<-ratio_pwy[,c(1,i)]
  colnames(df)<-c("mid", "level")
  df$cat<-colnames(ratio_pwy)[i]
  assign(paste0("ratiopwy",i), df )
}

ratio_pwy_all<-rbind(ratiopwy3, ratiopwy4, ratiopwy5, ratiopwy6, ratiopwy7)
ratio_pwy_label<-as.data.frame(unique(ratio_pwy_all$cat))
colnames(ratio_pwy_label)<-"pwy_name"
ratio_pwy_y_order <- factor(ratio_pwy_all$cat, levels = rev(ratio_pwy_label$pwy_name))
ratio_pwy_x_order <- factor(ratio_pwy_all$mid, levels = ratio_pwy$mid)
ratio_pwy_all$lglevel<-log2(ratio_pwy_all$level)
# read in stratified RNA pathway data
s1<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpr_stratified_raw_s1.sas7bdat")
s2<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpr_stratified_raw_s1.sas7bdat")
s3<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpr_stratified_raw_s1.sas7bdat")
s4<- read_sas("/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/gapfill_pathway/gpr_stratified_raw_s1.sas7bdat")

str_rna_pwy_label<-read.csv(file="/udd/nhdow/microbiome/stratified_rna_pwy_channing_label.csv")

rna_pwy_label <- read.csv(file = "/udd/nhdow/microbiome/rna_pwy_channing_label.csv")
rna_pwy_label$pwy_name<-paste(toupper(substr(rna_pwy_label$pwy_name, 1, 1)), substr(rna_pwy_label$pwy_name, 2, nchar(as.character(rna_pwy_label$pwy_name))), sep="")
ratio_pwy_label_st<-left_join(ratio_pwy_label, rna_pwy_label, by="pwy_name")

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
# create data for stack plots for RNA/DNA ratio pathway contributions by species
stack_area <- function(ec_name, num) {
  all_dna_ec1$pwy_channing<-rownames(all_dna_ec1)
  all_dna_ec1<-left_join(all_dna_ec1, str_rna_pwy_label, by="pwy_channing")
  rownames(all_dna_ec1)<-all_dna_ec1$pwy_name
  all_dna_ec1<-all_dna_ec1[, -c(381,382)]
  all_dna_ec1<-sweep(all_dna_ec1,
                     2,
                     STATS = colSums(all_dna_ec1),
                     FUN = "/")
  mean_rel<-as.data.frame(rowSums(all_dna_ec1, na.rm = T)/rowSums(all_dna_ec1!=0, na.rm = T))
  
  mean_rel<-sweep(mean_rel,
                  2,
                  STATS = colSums(mean_rel, na.rm = T),
                  FUN = "/")
  colnames(mean_rel)<-"mean_rel"
  mean_rel$bug_name<-rownames(mean_rel)
  stratum <- mean_rel[order(-mean_rel$mean_rel), ]
  stratum <- stratum %>% separate(bug_name, c(NA, "bugs"), "s__")
  stratum$bugs[is.na(stratum$bugs)]<-"Unclassified"
  stratum_cleared <- stratum[1:num, ]
  stratum_cleared$bugs <- gsub("_", " ", stratum_cleared$bugs)
  stratum_cleared$bugs<-paste(toupper(substr(stratum_cleared$bugs, 1, 1)), substr(stratum_cleared$bugs, 2, nchar(stratum_cleared$bugs)), sep="")
  stratum_cleared$ec<-ec_name
  return(stratum_cleared)
}

for (i in 1:5){
  if(i==1){
    ec_name<-as.character(ratio_pwy_label_st[i,1])
    ec<-as.character(ratio_pwy_label_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 3)
    bugall_ratiopwy<-bug1
  } else if (i>1){  
    ec_name<-as.character(ratio_pwy_label_st[i,1])
    ec<-as.character(ratio_pwy_label_st[i,2])
    t_s1_sub<-transform_data(s1, ec, "_SF05", "_s1")
    t_s2_sub<-transform_data(s2, ec, "_SF06", "_s2")
    t_s3_sub<-transform_data(s3, ec, "_SF07", "_s3")
    t_s4_sub<-transform_data(s4, ec, "_SF08", "_s4")
    all_dna_ec1<-cbind(t_s1_sub, t_s2_sub, t_s3_sub, t_s4_sub)
    bug1<-stack_area(ec_name = ec_name,num = 3)
    bugall_ratiopwy<-rbind(bugall_ratiopwy, bug1)
  }
}

bug_ratioall<-rbind(bugall_ratioec, bugall_ratiopwy, bugall_ec, bugall_pwy)
length(unique(bug_ratioall$bugs))
bug_ratio_color<-data.frame(bugs=unique(bug_ratioall$bugs), colors=colorRampPalette(brewer.pal(12, "Paired"))(22))

bug_ratio_color_value<-c('Bacteroides uniformis'='#A6CEE3',
                         'Bacteroides vulgatus'='#5FA0CA',
                         'Bacteroides cellulosilyticus'='#257CB2',
                         'Unclassified'='#D3D3D3',
                         'Bacteroides plebeius'='#A5D981',
                         'Faecalibacterium prausnitzii'='#63B84F',
                         'Methanobrevibacter smithii'='#4F9F3B',
                         'Roseburia intestinalis'='#B89B74',
                         'Alistipes onderdonkii'='#F68181',
                         'Alistipes shahii'='#E93E3F',
                         'Eubacterium siraeum'='#E9412F',
                         'Acidaminococcus sp bv3l6'='#F6975B',
                         'Synergistes sp 3 1 syn1'='#FDAC4F',
                         'Akkermansia muciniphila'='#FE8B15',
                         'Bacteroides coprocola'='#ED8F47',
                         'Escherichia coli'='#D1AAB7',
                         'Prevotella copri'='#A585BF',
                         'Eubacterium rectale'='#73489F',
                         'Prevotella buccae'='#A99099',
                         'Prevotella denticola'='#F7F599',
                         'Anaerococcus obesiensis'='#D9AF63',
                         'Streptococcus lutetiensis'='#B15928')
# barplot of RNA/DNA ratio EC contributions by species
level_ratioec <- factor(bugall_ratioec$ec, levels = rev(ratio_ec_label$ec_label))
ratio_ec_bar<-ggplot(bugall_ratioec, aes(x = level_ratioec, y = mean_rel, fill = bugs)) +
                     geom_bar(stat = "identity")+theme_bw()+ylim(0,1)+
                     scale_fill_manual(values = bug_ratio_color_value)+ coord_flip()+
                     theme(legend.position = "none",
                           plot.title = element_blank(),
                           axis.title=element_blank(),
                           axis.text.y=element_blank(),
                           axis.text.x = element_text(size = 50))
# barplot of RNA/DNA ratio pathway contributions by species
level_ratiopwy <- factor(bugall_ratiopwy$ec, levels = rev(ratio_pwy_label$pwy_name))
ratio_pwy_bar<-ggplot(bugall_ratiopwy, aes(x = level_ratiopwy, y = mean_rel, fill = bugs)) +
                      geom_bar(stat = "identity")+theme_bw()+ylim(0,1)+
                      scale_fill_manual(values = bug_ratio_color_value)+ coord_flip()+
                      theme(legend.position = "none",
                            plot.title = element_blank(),
                            axis.title=element_blank(),
                            axis.text.y=element_blank(),
                            axis.text.x = element_text(size = 50))
# plot for RNA/DNA ratio EC distribution
ratio_ec_plot<-ggplot(ratio_ec_all, aes(ratio_ec_x_order, ratio_ec_y_order)) +
                      geom_tile(aes(fill =lglevel)) +
                      scale_fill_gradientn(colours = c("black", "#d73027", "#fc8d59", "#6eb7c8", 
                                                       "#56889f", "#28647c"), 
                      name ="RNA/DNA ratio\nlog2 scale", na.value = "white")+
                      theme(legend.position = "left",
                            legend.text = element_text(size=80),
                            legend.title = element_text(size=80),
                            plot.title = element_blank(),
                            axis.title=element_blank(),
                            axis.text.y=element_text(size=100),
                            axis.text.x = element_blank())
# plot for RNA/DNA ratio pathway distribution
ratio_pwy_plot<-ggplot(ratio_pwy_all, aes(ratio_pwy_x_order, ratio_pwy_y_order)) +
                       geom_tile(aes(fill =lglevel)) +
                       scale_fill_gradientn(colours = c("black", "#d73027", "#fc8d59", 
                                                        "#6eb7c8", "#56889f", "#28647c"), 
                       name ="RNA/DNA ratio\nlog2 scale", na.value = "white")+
                       theme(legend.position = "left",
                             legend.text = element_text(size=80),
                             legend.title = element_text(size=80),
                             plot.title = element_blank(),
                             axis.title=element_blank(),
                             axis.text.y=element_text(size=100),
                             axis.text.x = element_blank())

png(file="/udd/nhdow/microbiome/review/fig2e_ratio_ec.png",width=7000,height=800, pointsize=50)
plot_grid(ratio_ec_plot, ratio_ec_bar, align = "h", ncol =  2, rel_widths = c(6/7, 1/7))
dev.off()

png(file="/udd/nhdow/microbiome/review/fig2e_ratio_pwy.png",width=7000,height=800, pointsize=50)
plot_grid(ratio_pwy_plot, ratio_pwy_bar, align = "h", ncol =  2, rel_widths = c(6/7, 1/7))
dev.off()
# create a common legend for contribution barplots
png(file="/udd/nhdow/microbiome/review/fig2e_bug_ratio_legend.png",width=3500,height=2000, pointsize=50)
ggplot(bug_ratioall, aes(x = ec, y = mean_rel, fill = bugs)) +
       geom_bar(stat = "identity")+theme_bw()+ylim(0,1)+
       scale_fill_manual(values = bug_ratio_color_value)+ coord_flip()+
       theme(legend.position = "right",
             legend.text = element_text(size=50, face="italic"),
             legend.title = element_blank(),
             plot.title = element_blank(),
             axis.title=element_blank(),
             axis.text.y=element_blank(),
             axis.text.x = element_blank())
dev.off()