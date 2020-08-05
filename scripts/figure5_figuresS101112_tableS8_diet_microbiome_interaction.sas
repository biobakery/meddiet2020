/*************************************************************************************************************************
  SAS programs for creating Figures 5a, S10, S11, S12 and Table S8 using linear mixed models
  1) Purpose: a) p-values for interaction between PCo1 /P. copri and meddiet index in relation
                 to cadiometabolic disease risk score and biomarker levels
              b) p-values for interaction between other abundant species and meddiet index in relation
                 to cadiometabolic disease risk score and biomarker levels
              c) stratified analysis in strata defined by PCo1 and P. Copri
              d) independent association of meddiet index with cadiometabolic disease risk score and biomarker levels 
  2) Microbial variables: PCo1 score, species relative abundances
  3) Dietary variables: Mediterranean diet index
  4) Biomarkers: cardiometabolic disease risk score, HDL, TG, TC, HbA1c, CRP
  5) Covariates: age, total energy intake, physical activity, smoking status, antibiotic use, 
                 probiotic use, stool type, PPI use, metformin use
*************************************************************************************************************************/

libname meddiet "/udd/nhdow/microbiome/";
* read in PCo 1 & 2 Scores;
proc sort data=meddiet.pcoap_2t;by id_2v;run;
* read in species-level feature relative abundance;
proc sort data=meddiet.species_2t;by id_2v;RUN;

data species_2t;
merge meddiet.species_2t (in=a) meddiet.pcoap_2t (in=b);
by id_2v;
if a and b;
subjectid=substr(id_2v, 1, 8)*1;
cvisit=substr(id_2v, 11, 1)*1;
proc sort data=species_2t;by subjectid cvisit;
run;
* read in meta data;
proc sort data=meddiet.meta_ffq_8612;by subjectid cvisit;run;
* read in plasma biomarker data;
proc sort data=meddiet.biomarker_mlvs;by subjectid cvisit;run;
* merge all data together;
data meta_biomarker_spp;
merge meddiet.biomarker_mlvs (in=a) meddiet.meta_ffq_8612 species_2t (in=b);
by subjectid cvisit;
if a and b;
run;
* calculate quinitles of biomarkers & P copri;
proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=5;
var tc_plasma hdlc_plasma tg_plasma crp_plasma hba1cp s__Prevotella_copri;
ranks tc_plasmar hdlc_plasmar tg_plasmar crp_plasmar hba1cpr pcopriq;
run;
* dichotomize most abundant species and PCo1 score by median;
proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=2;
var s__Eubacterium_rectale s__Faecalibacterium_prausnitzii s__Subdoligranulum_unclassified 
    s__Bacteroides_uniformis s__Alistipes_putredinis s__Eubacterium_siraeum s__Bacteroides_stercoris s__Ruminococcus_bromii 
    s__Bacteroides_vulgatus s__Akkermansia_muciniphila s__Eubacterium_eligens x1;
ranks erectaleq fprauq sunclassq bunifq aputrq esiraq bstercq rbromq bvulgq amucinq eeligenq x1median;
run;

data meta_biomarker_spp;
set meta_biomarker_spp;
* calculate cardiometabolic disease risk score;
mscore=(tc_plasmar+1)+(6-hdlc_plasmar)+(crp_plasmar+1)+(hba1cpr+1)+(tg_plasmar+1);
* log transform biomarkers;
lgtg=log(tg_plasma);
lgtc=log(tc_plasma);
lgcrp=log(crp_plasma);
lghdl=log(hdlc_plasma);
lghba1c=log(hba1cp);
* dichotomize P. copri by 20th percentile;
if pcopriq=4 then pcopriq=1;
else pcopriq=0;
run;
* calculate energy residuals of meddiet score;
proc glm data=meta_biomarker_spp;
model emed122ch=calor122cn;
output out=meta_biomarker_spp r=emed122chr;
run;
quit;
* macro for calculate p-value for interaction using LRT;
%macro lrt(biomarker, var);
ods output FitStatistics=model_inter;
proc mixed data=meta_biomarker_spp method=ml;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq;
Model &biomarker = &var &var*emed122chr emed122chr age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION ; 
random intercept/subject=subjectid  TYPE=cs;
run;
ods output close;

ods output FitStatistics=model_nointer;
proc mixed data=meta_biomarker_spp method=ml;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq;
Model &biomarker = &var emed122chr age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION ; 
random intercept/subject=subjectid  TYPE=cs;
run;
ods output close;
data model_inter;
set model_inter;
if descr="-2 Log Likelihood";
rename value=lgld_inter;
run;

data model_nointer;
set model_nointer;
if descr="-2 Log Likelihood";
rename value=lgld_nointer;
run;

data &biomarker.&var.;
merge model_inter model_nointer;
by descr;
lrt_p=1-CDF ('CHISQUARE', abs(lgld_nointer-lgld_inter), 1);
biomarker="&biomarker.&var.";
run;
%mend;
* p-values for interaction between PCo1 and meddiet index in relation
  to cadiometabolic disease risk score and biomarker levels in Figure 5a & Figure S12;
%lrt(mscore, x1);
%lrt(lgcrp, x1);
%lrt(lgtc, x1);
%lrt(lgtg, x1);
%lrt(lghba1c, x1);
%lrt(lghdl, x1);

data pvalue_lrt;
set 
mscorex1
lgcrpx1
lgtcx1
lgtgx1
lghba1cx1
lghdlx1;
run;
* macro for calculating multivariate-adjusted values from mixied linear models;
%macro inter(biomarker, var1, var2);
proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq;
Model &biomarker = &var1 &var1*emed122chr emed122chr age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION outpm=meta_biomarker_spp; 
random intercept/subject=subjectid  TYPE=cs;
run; 

data meta_biomarker_spp;
set meta_biomarker_spp;
&var2=Pred;
drop Pred StdErrPred DF Alpha Lower Upper Resid;
run;
%mend;
* multivariate-adjusted values for cardiometabolic disease risk score and biomarkers;
%inter(mscore, x1, mscore_p);
%inter(lgcrp, x1, lgcrp_p);
%inter(lgtc, x1, lgtc_p);
%inter(lgtg, x1, lgtg_p);
%inter(lghba1c, x1, lghba1c_p);
%inter(lghdl, x1, lghdl_p);

data exp;
set meta_biomarker_spp;
crp=exp(lgcrp_p);
hba1c=exp(lghba1c_p);
tg=exp(lgtg_p);
tc=exp(lgtc_p);
hdl=exp(lghdl_p);
keep x1 crp mscore_p hba1c tg tc hdl emed122ch subjectid cvisit;
run;
* output values for plotting in Figure 5a & Figure S12;
proc export data=exp
     outfile="/udd/nhdow/microbiome/pco1_diet_cardiometabolic.csv"
     dbms=csv 
     replace;
run;
* p-values for interaction between microbial species and meddiet index in relation
  to cadiometabolic disease risk score in Figure S11;
%lrt(mscore, pcopriq);
%lrt(mscore, erectaleq);
%lrt(mscore, fprauq);
%lrt(mscore, sunclassq);
%lrt(mscore, bunifq);
%lrt(mscore, aputrq);
%lrt(mscore, esiraq);
%lrt(mscore, bstercq);
%lrt(mscore, rbromq);
%lrt(mscore, bvulgq);
%lrt(mscore, amucinq);
%lrt(mscore, eeligenq);

data pvalue_lrt_spp;
set 
mscorepcopriq
mscoreerectaleq
mscorefprauq
mscoresunclassq
mscorebunifq
mscoreaputrq
mscoreesiraq
mscorebstercq
mscorerbromq
mscorebvulgq
mscoreamucinq
mscoreeeligenq;
run;
* p-values for interaction between P. copri and meddiet index in relation
  to cadiometabolic disease risk score and biomarker levels in Table S8;
%lrt(mscore, pcopriq);
%lrt(lgcrp, pcopriq);
%lrt(lgtc, pcopriq);
%lrt(lgtg, pcopriq);
%lrt(lghba1c, pcopriq);
%lrt(lghdl, pcopriq);

data lrt_pvalue_pcopri;
set 
mscorepcopriq
lgcrppcopriq
lgtcpcopriq
lgtgpcopriq
lghba1cpcopriq
lghdlpcopriq;
run;
* output multivariate-adjusted values for plotting Figure S11;
%inter(mscore, pcopriq, mscore_pcopriq);
%inter(mscore, erectaleq, mscore_erectaleq);
%inter(mscore, fprauq, mscore_fprauq);
%inter(mscore, sunclassq, mscore_sunclassq);
%inter(mscore, bunifq, mscore_bunifq);
%inter(mscore, aputrq, mscore_aputrq);
%inter(mscore, esiraq, mscore_esiraq);
%inter(mscore, bstercq, mscore_bstercq);
%inter(mscore, rbromq, mscore_rbromq);
%inter(mscore, bvulgq, mscore_bvulgq);
%inter(mscore, amucinq, mscore_amucinq);
%inter(mscore, eeligenq, mscore_eeligenq);

data exp_spp;
set meta_biomarker_spp;
keep x1 emed122ch subjectid cvisit  mscore_pcopriq pcopriq
     mscore_erectaleq erectaleq mscore_fprauq fprauq  mscore_sunclassq sunclassq
     mscore_bunifq bunifq mscore_aputrq aputrq mscore_esiraq esiraq mscore_bstercq bstercq
     mscore_rbromq rbromq mscore_bvulgq bvulgq mscore_amucinq amucinq mscore_eeligenq eeligenq;
run;

proc export data=exp_spp
     outfile="/udd/nhdow/microbiome/species_diet_cardiometabolic.csv"
     dbms=csv 
     replace;
run;
* calculate quartiles of meddiet index for Table S8;
proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=4;
var emed122chr;ranks emed122chrq;run;
* macro for calculating multivariate-adjusted means across quintiles
  of meddiet index in strata defined by PCo1 score / P. copri relative abundance;
%macro strata(biomarker, stratum);
proc sort data=meta_biomarker_spp;by &stratum;run;
ods output LSMeans=adj_mean;
proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq emed122chrq (ref=first);
Model &biomarker = emed122chrq age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION; 
random intercept/subject=subjectid  TYPE=cs;
by &stratum;
lsmeans emed122chrq/cl;
run; 
ods output close;

data &biomarker.&stratum.;
set adj_mean;
rename &stratum=strata;
group="&stratum.";
biomarker="&biomarker.";
run;
%mend;
%strata(mscore, pcopriq);
%strata(lgcrp, pcopriq);
%strata(lgtc, pcopriq);
%strata(lgtg, pcopriq);
%strata(lghba1c, pcopriq);
%strata(lghdl, pcopriq);
%strata(mscore, x1median);
%strata(lgcrp, x1median);
%strata(lgtc, x1median);
%strata(lgtg, x1median);
%strata(lghba1c, x1median);
%strata(lghdl, x1median);
* output adjusted mean for creating Table S8;
data adj_mean_stratified;
set mscorepcopriq
lgcrppcopriq
lgtcpcopriq
lgtgpcopriq
lghba1cpcopriq
lghdlpcopriq
mscorex1median
lgcrpx1median
lgtcx1median
lgtgx1median
lghba1cx1median
lghdlx1median;
exp_mean=exp(Estimate);
exp_lower=exp(Lower);
exp_upper=exp(Upper);
if biomarker="mscore" then do;
exp_mean=Estimate;
exp_lower=Lower;
exp_upper=Upper;
end;
mean_ci= put(round(exp_mean,0.1),4.2)||' ('||put(round(exp_lower,0.1), 4.2)||', '||put(round(exp_upper,0.1), 4.2)||')';
run;
proc export data=adj_mean_stratified
     outfile="/udd/nhdow/microbiome/species_diet_cardiometabolic_adj_mean_stratified.csv"
     dbms=csv 
     replace;
run;
* create score variable (quartile-specific median) for quartiles of meddiet index used in trend test;
data meta_biomarker_spp;
set meta_biomarker_spp;
if emed122chrq=0 then emed122chm= 2.6666667 ;
if emed122chrq=1 then emed122chm= 3.8888889 ;
if emed122chrq=2 then emed122chm= 4.8888889 ;
if emed122chrq=3 then emed122chm= 6.4444444 ;
run;
* macro for trend test in Table S8;
%macro ptrend(biomarker, stratum);
proc sort data=meta_biomarker_spp;by &stratum;
ods output SolutionF=solution;
proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq emed122chrq;
Model &biomarker = emed122chm age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION; 
random intercept/subject=subjectid  TYPE=cs;
by &stratum;
run; 
ods output close;

data &biomarker.&stratum.;
set solution;
if effect="emed122chm";
rename &stratum=strata;
group="&stratum.";
biomarker="&biomarker.";
run;
%mend;
%ptrend(mscore, pcopriq);
%ptrend(lgcrp, pcopriq);
%ptrend(lgtc, pcopriq);
%ptrend(lgtg, pcopriq);
%ptrend(lghba1c, pcopriq);
%ptrend(lghdl, pcopriq);
%ptrend(mscore, x1median);
%ptrend(lgcrp, x1median);
%ptrend(lgtc, x1median);
%ptrend(lgtg, x1median);
%ptrend(lghba1c, x1median);
%ptrend(lghdl, x1median);

data p_trend_stratified;
set mscorepcopriq
lgcrppcopriq
lgtcpcopriq
lgtgpcopriq
lghba1cpcopriq
lghdlpcopriq
mscorex1median
lgcrpx1median
lgtcx1median
lgtgx1median
lghba1cx1median
lghdlx1median;
run;

* Meddiet-biomarker independent associations in Figure S10;
%macro ind_effect(biomarker, var2);
proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq;
Model &biomarker = emed122chr age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 x1/ SOLUTION outpm=meta_biomarker_spp; 
random intercept/subject=subjectid  TYPE=cs;
run; 

data meta_biomarker_spp;
set meta_biomarker_spp;
&var2=Pred;
drop Pred StdErrPred DF Alpha Lower Upper Resid;
run;
%mend;
%ind_effect(mscore, mscore_p);
%ind_effect(lgcrp, lgcrp_p);
%ind_effect(lgtc, lgtc_p);
%ind_effect(lgtg, lgtg_p);
%ind_effect(lghba1c, lghba1c_p);
%ind_effect(lghdl, lghdl_p);
* output multivariate-adjusted values for plotting Figure S10;
data exp_ind;
set meta_biomarker_spp;
crp=exp(lgcrp_p);
hba1c=exp(lghba1c_p);
tg=exp(lgtg_p);
tc=exp(lgtc_p);
hdl=exp(lghdl_p);
keep x1 crp mscore_p hba1c tg tc hdl emed122ch subjectid cvisit;
run;

proc export data=exp_ind
     outfile="/udd/nhdow/microbiome/meddiet_biomarker_independent_association.csv"
     dbms=csv 
     replace;
run;
endsas;
