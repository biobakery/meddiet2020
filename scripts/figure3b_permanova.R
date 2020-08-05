##############################################################################################################################################
# R programs for creating Figure 3b
# 1) Purpose: PERMANOVA Test: Associations (% of total variation explained) of Mediterranean diet and its components,   
#                            covariables & biomarkers with overall structure of microbial communities
# 2) Microbial variable: species-level beta-diversity of microbial communities
# 3) Dietary exposures: Mediterranean diet index, fruit, vegetables, whole grains, monounsaturated fat/saturated fat ratio, 
#                       red meat, alcohol, legume, nuts, & fish
# 4) Covariates: age, total energy intake, physical activity, smoking status, 
#                antibiotic use, probiotic use, stool type, PPI use, metformin use
# 5) Biomarkers: HbA1c, TC, TG, HDL, & CRP
#############################################################################################################################################
library(scales)
library(tidyverse)
library(vegan)

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
med_id<-meta_med %>% select(mid) 
ttax_rpk <- as.data.frame(t(tax_rpk_species))
ttax_rpk$mid<-rownames(ttax_rpk)
ttax_rpk<-inner_join(med_id, ttax_rpk, by="mid")
rownames(ttax_rpk)<-ttax_rpk$mid
ttax_rpk <- subset(ttax_rpk, select = -mid)
# calculate relative abundance
ttax_rpk_rel <-
  sweep(ttax_rpk,
        1,
        STATS = rowSums(ttax_rpk),
        FUN = "/")
# calculate bray-curtis 
ttax_rpk_rel_bc <- vegdist(ttax_rpk_rel, "bray")
# PERMANOVA using adonis
x <- adonis2(formula = ttax_rpk_rel_bc ~ emed122ch, strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
emed122ch <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ ms122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
ms122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ fish122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
fish122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ meat122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
meat122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ wgrain122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
wgrain122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ legu122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
legu122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ fruit122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
fruit122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ veg122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
veg122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ nut122c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
nut122c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ alco122cn,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
alco122cn <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ lgcrp,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
lgcrp <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ lgtg,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
lgtg<- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ lgtc,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
lgtc<- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ lghdl,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
lghdl<- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ lghba1c,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
lghba1c <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ calor122cn,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
calor122cn <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ totMETs_paq,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
totMETs_paq <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ ant_12mo_qu,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
ant_12mo_qu <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ probio_2mo_qu,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
probio_2mo_qu <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ stool_type,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
stool_type <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ age_fecal,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
age_fecal <- x[!is.na(x$F), ]
x <- adonis2(formula = ttax_rpk_rel_bc ~ smk,  strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
smk <- x[!is.na(x$F), ]
x<-adonis2(formula = ttax_rpk_rel_bc ~ metfo12, strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999) 
metfo12 <- x[!is.na(x$F), ]
x<-adonis2(formula = ttax_rpk_rel_bc ~pril12, strata=as.factor(meta_med$SubjectID),  data = meta_med, permutations = 999)													
pril12<- x[!is.na(x$F), ]
# merge PERMANOVA results
all_taxonomy <-rbind(emed122ch, ms122c, fish122c, meat122c, wgrain122c, legu122c, fruit122c, veg122c, nut122c, alco122cn, 
                     lgcrp, lghdl, lgtc, lgtg, lghba1c, calor122cn, totMETs_paq, probio_2mo_qu, stool_type, age_fecal, 
                     smk, pril12, metfo12, ant_12mo_qu)
# label for variables
component<- c('Mediterranean diet index',
              'Monounsaturated fat/saturated fat ratio',
              'Fish',
              'Red/processed meat',
              'Whole grains',
              'Legume',
              'Fruits',
              'Vegetables',
              'Nuts',
              'Alcohol',
              'C-reactive protein',
              'High-density lipoprotein',
              'Total cholesterol',
              'Triglyceride',
              'Hemoglobin A1c',
              'Total energy intake',
              'Physical activity',
              'Probiotics use',
              'Bristol stool scale',
              'Age',
              'Smoking',
              'Proton pump inhibitors',
              'Metformin',
              'Antibiotics')  
all_taxonomy<-cbind(all_taxonomy, component)
all_taxonomy_diet<-all_taxonomy[1:10,]
all_taxonomy_diet$cat<-"Diet"
all_taxonomy_bio<-all_taxonomy[11:15,]
all_taxonomy_bio$cat<-"Plasma biomarkers"
all_taxonomy_cov<-all_taxonomy[16:24,]
all_taxonomy_cov$cat<-"Covariables"
all_taxonomy_tax<-rbind(all_taxonomy_diet, all_taxonomy_bio, all_taxonomy_cov)
all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)
# FDR adjustment for p-values
all_taxonomy_tax$fdr <-
  p.adjust(
    all_taxonomy_tax$`Pr(>F)`,
    method = "BH",
    n = length(all_taxonomy_tax$`Pr(>F)`)
  )
# create figure to show PERMANOVA results
all_taxonomy_tax$stars <- cut(all_taxonomy_tax$fdr, breaks=c(-Inf, 0.005, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))
color_code<-c("Diet"="#56B4E9",
              "Plasma biomarkers"="#999999",
              "Covariables"="#E69F00")
png(file="/udd/nhdow/microbiome/review/fig3b_barplot_permanova_change_font.png",width=2200, height=2000, pointsize=50)
ggplot(all_taxonomy_tax, aes(level_y_order, R2, fill=cat)) +
       geom_bar(stat="identity")+
       theme_light()+
       scale_y_continuous(labels=scales::percent, limits = c(0, 0.01))+
       geom_text(aes(label=stars), color="black", size=30) +
       scale_fill_manual(values = color_code)+ coord_flip()+
       theme(legend.position = "bottom",
             legend.text = element_text(size = 50),
             legend.title = element_blank(),
             plot.title = element_blank(),
             axis.title.x=element_blank(),
             axis.title.y= element_blank(),
             axis.text.y=element_text(size=60),
             axis.text.x = element_text(size=60))  
dev.off()