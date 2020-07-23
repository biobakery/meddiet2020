proc import datafile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\meta_biomarker_468.csv"
     out=meta_biomarker
     dbms=csv
     replace;
     getnames=yes;
run;

proc import datafile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\species_relative_abundance_2t.csv"
     out=species_2t
     dbms=csv
     replace;
     getnames=yes;
run;

proc import datafile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\pcoap_2t.csv"
     out=pcoap_2t
     dbms=csv
     replace;
     getnames=yes;
run;

proc sort data=pcoap_2t;by VAR1;run;
proc sort data=species_2t;by VAR1;RUN;

data species_2t;
merge species_2t (in=a) pcoap_2t (in=b);
by VAR1;
if a and b;
run;

data species_2t;
set species_2t;
subjectid=substr(VAR1, 1, 8)*1;
cvisit=substr(VAR1, 11, 1)*1;
proc sort data=species_2t;by subjectid cvisit;
run;

proc sort data=meta_biomarker;by subjectid cvisit;run;
data meta_biomarker_spp;
merge meta_biomarker (in=a) species_2t (in=b);
by subjectid cvisit;
if a and b;
run;

proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=5;
var ldlc_plasma tc_plasma hdlc_plasma tg_plasma crp_plasma hba1cp age_fecal;
ranks ldlc_plasmar tc_plasmar hdlc_plasmar tg_plasmar crp_plasmar hba1cpr age_fecalr;
run;

data meta_biomarker_spp;
set meta_biomarker_spp;
mscore=(tc_plasmar+1)+(5-hdlc_plasmar)+(crp_plasmar+1)+(hba1cpr+1)+(tg_plasmar+1);
lgtg=log(tg_plasma);
lgtc=log(tc_plasma);
lgcrp=log(crp_plasma);
lgldl=log(ldlc_plasma);
lghdl=log(hdlc_plasma);
lghba1c=log(hba1cp);
thratio=tc_plasma/hdlc_plasma;
run;

proc rank data=meta_biomarker_spp out=meta_biomarker_spp normal=blom;
var crp_plasma tg_plasma tc_plasma hdlc_plasma hba1cp ldlc_plasma mscore;
ranks crp_int tg_int tc_int hdl_int hba1c_int ldl_int mscore_int;
run;

proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=10;var totmets_paq calor122cn;ranks totmets_paqq calor122cnr;run;

proc glm data=meta_biomarker_spp;
model emed122ch=calor122cn;
output out=meta_biomarker_spp r=emed122chr;run;

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
%lrt(mscore_int, x1);
%lrt(crp_int, x1);
%lrt(tc_int, x1);
%lrt(tg_int, x1);
%lrt(hba1c_int, x1);
%lrt(hdl_int, x1);

%lrt(mscore, x1);
%lrt(lgcrp, x1);
%lrt(lgtc, x1);
%lrt(lgtg, x1);
%lrt(lghba1c, x1);
%lrt(lghdl, x1);

data pvalue_lrt;
set mscore_intx1
crp_intx1
tc_intx1
tg_intx1
hba1c_intx1
hdl_intx1
mscorex1
lgcrpx1
lgtcx1
lgtgx1
lghba1cx1
lghdlx1;
run;
proc export data=pvalue_lrt
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\pco1_diet_cardiometabolic_lrt_pvalues.csv"
     dbms=csv 
     replace;
run;
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

proc export data=exp
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\pco1_diet_cardiometabolic_predicted_values.csv"
     dbms=csv 
     replace;
run;
proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=2;
var s__Eubacterium_rectale s__Faecalibacterium_prausnitzii s__Subdoligranulum_unclassified 
    s__Bacteroides_uniformis s__Alistipes_putredinis s__Eubacterium_siraeum s__Bacteroides_stercoris s__Ruminococcus_bromii 
    s__Bacteroides_vulgatus s__Akkermansia_muciniphila s__Eubacterium_eligens;
ranks erectaleq fprauq sunclassq bunifq aputrq esiraq bstercq rbromq bvulgq amucinq eeligenq;
run;

proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=5;var s__Prevotella_copri;ranks pcopriq;run;
data meta_biomarker_spp; set meta_biomarker_spp;if pcopriq=4 then pcopriq=1;else pcopriq=0;run;
/*proc contents data=species_2t short order=varnum;run;
s__Eubacterium_rectale s__Faecalibacterium_prausnitzii s__Subdoligranulum_unclassified 
s__Bacteroides_uniformis s__Alistipes_putredinis s__Prevotella_copri
s__Eubacterium_siraeum s__Bacteroides_stercoris s__Ruminococcus_bromii 
s__Bacteroides_vulgatus s__Akkermansia_muciniphila s__Eubacterium_eligens       */

%lrt(mscore_int, pcopriq);
%lrt(mscore_int, erectaleq);
%lrt(mscore_int, fprauq);
%lrt(mscore_int, sunclassq);
%lrt(mscore_int, bunifq);
%lrt(mscore_int, aputrq);
%lrt(mscore_int, esiraq);
%lrt(mscore_int, bstercq);
%lrt(mscore_int, rbromq);
%lrt(mscore_int, bvulgq);
%lrt(mscore_int, amucinq);
%lrt(mscore_int, eeligenq);

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
set mscore_intpcopriq
mscore_interectaleq
mscore_intfprauq
mscore_intsunclassq
mscore_intbunifq
mscore_intaputrq
mscore_intesiraq
mscore_intbstercq
mscore_intrbromq
mscore_intbvulgq
mscore_intamucinq
mscore_inteeligenq

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

proc export data=pvalue_lrt_spp
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\species_diet_cardiometabolic_lrt_pvalues.csv"
     dbms=csv 
     replace;
run;

%lrt(mscore_int, pcopriq);
%lrt(crp_int, pcopriq);
%lrt(tc_int, pcopriq);
%lrt(tg_int, pcopriq);
%lrt(hba1c_int, pcopriq);
%lrt(hdl_int, pcopriq);

%lrt(mscore, pcopriq);
%lrt(lgcrp, pcopriq);
%lrt(lgtc, pcopriq);
%lrt(lgtg, pcopriq);
%lrt(lghba1c, pcopriq);
%lrt(lghdl, pcopriq);

data lrt_pvalue_pcopri;
set mscore_intpcopriq
crp_intpcopriq
tc_intpcopriq
tg_intpcopriq
hba1c_intpcopriq
hdl_intpcopriq

mscorepcopriq
lgcrppcopriq
lgtcpcopriq
lgtgpcopriq
lghba1cpcopriq
lghdlpcopriq;
run;

proc export data=lrt_pvalue_pcopri
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\pcopri_diet_cardiometabolic_lrt_pvalues.csv"
     dbms=csv 
     replace;
run;

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
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\species_diet_cardiometabolic_predicted_values.csv"
     dbms=csv 
     replace;
run;

proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=4;
var emed122chr;ranks emed122chrq;run;
proc rank data=meta_biomarker_spp out=meta_biomarker_spp groups=2;
var x1;ranks x1median;run;

/*proc means data=meta_biomarker_spp mean median;class emed122chrq;var emed122ch;run;

data meta_biomarker_spp; set meta_biomarker_spp; emed122ch4=emed122ch/4;run;

proc sort data=meta_biomarker_spp;by pcopriq;run;

proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq emed122chrq;
Model mscore = emed122ch4 age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION; 
random intercept/subject=subjectid  TYPE=cs;
by pcopriq;
run; */


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
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\species_diet_cardiometabolic_adj_mean_stratified.csv"
     dbms=csv 
     replace;
run;
proc means data=meta_biomarker_spp median;class emed122chrq;var emed122ch;run;

data meta_biomarker_spp;
set meta_biomarker_spp;
if emed122chrq=0 then emed122chm= 2.6666667 ;
if emed122chrq=1 then emed122chm= 3.8888889 ;
if emed122chrq=2 then emed122chm= 4.8888889 ;
if emed122chrq=3 then emed122chm= 6.4444444 ;
run;

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

proc export data=p_trend_stratified
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\species_diet_cardiometabolic_p_trend_stratified.csv"
     dbms=csv 
     replace;
run;
proc means data=meta_biomarker_spp std;var emed122chr;run;
data meta_biomarker_spp;
set meta_biomarker_spp;
emed122chrsd=emed122chr/1.4324244 ;
run;
proc sort data=meta_biomarker_spp;by pcopriq;run;

proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq emed122chrq (ref=first);
Model mscore = emed122chrq age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION; 
random intercept/subject=subjectid  TYPE=cs;
by pcopriq;
run; 

proc sort data=meta_biomarker_spp;by x1median;run;

proc mixed data=meta_biomarker_spp;
Class subjectid cvisit calor122cnr age_fecalr totmets_paqq emed122chrq (ref=first);
Model mscore = emed122chrq age_fecalr calor122cnr totmets_paqq
               smk probio_2mo_qu stool_type ant_12mo_qu pril12 metfo12 / SOLUTION; 
random intercept/subject=subjectid  TYPE=cs;
by x1median;
run; 

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
%ind_effect(mscore, mscore_p);*0.0306;
%ind_effect(lgcrp, lgcrp_p);*0.0866;
%ind_effect(lgtc, lgtc_p);*0.1165;
%ind_effect(lgtg, lgtg_p);*0.0421;
%ind_effect(lghba1c, lghba1c_p);*0.1503;
%ind_effect(lghdl, lghdl_p);*0.0076;

%ind_effect(mscore_int, mscore_int_p);*0.0321;
%ind_effect(crp_int, lgcrp_int_p);*0.0343;
%ind_effect(tc_int, lgtc_int_p);*0.1231;
%ind_effect(tg_int, lgtg_int_p);*0.0672;
%ind_effect(hba1c_int, lghba1c_int_p);*0.2706;
%ind_effect(hdl_int, lghdl_int_p);*0.0083;

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
     outfile="\\Mac\Dropbox\hpfs_meddiet\revision_1st\data\meddiet_biomarker_independent_association.csv"
     dbms=csv 
     replace;
run;
