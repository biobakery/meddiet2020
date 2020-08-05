/*************************************************************************************************************************
  SAS programs for creating Figures 5c and Table S7 using data from a case-control study on myocardial infarction
  1) Purpose: calculate multivariate-adjusted risk ratio of MI associated with cardiometabolic disease risk score
              and biomarkers 
  2) exposure variables: cardiometabolic disease risk score, HDL, TG, TC, HbA1c, CRP
  3) outcome variables: incident MI
  4) Covariates: age, smoking status, month of blood sampling, parental history of MI before the age of 60 years
                 alcohol intake, level of physical activity, and body mass index 
*************************************************************************************************************************/
libname mi "/udd/nhdow/microbiome/mi_caco/";

data caco_mi;
set mi.caco_final;
/* quintiles for biomarkers using cutoffs calculated from controls only*/
/* HbA1c */
if .<hba1c<=5.180 then hba1cq=1;
else if 5.180<hba1c<=5.414 then hba1cq=2;
else if 5.414<hba1c<=5.587 then hba1cq=3;
else if 5.587<hba1c<=5.769 then hba1cq=4;
else if hba1c>5.769 then hba1cq=5;
/* CRP */
if .<crp<=0.41 then crpq=1;
else if 0.41<crp<=0.74 then crpq=2;
else if 0.74<crp<=1.31 then crpq=3;
else if 1.31<crp<=2.45 then crpq=4;
else if crp>2.45 then crpq=5;
/* HDL-C */
if .<hdl<=36.1 then hdlq=1;
else if 36.1<hdl<=41.5 then hdlq=2;
else if 41.5<hdl<=47.8 then hdlq=3;
else if 47.8<hdl<=55.7 then hdlq=4;
else if hdl>55.7 then hdlq=5;
/* TC */
if .<chol<=173 then cholq=1;
else if 173<chol<=193 then cholq=2;
else if 193<chol<=209 then cholq=3;
else if 209<chol<=230 then cholq=4;
else if chol>230 then cholq=5;
/* TG */
if .<trig<=78 then trigq=1;
else if 78<trig<=102 then trigq=2;
else if 102<trig<=132 then trigq=3;
else if 132<trig<=188 then trigq=4;
else if trig>188 then trigq=5;
* create score variable (quartile-specific median) for quintiles 
of biomarkers used in trend test;
if crpq=1 then crpmed=0.26; 
        else if crpq=2 then crpmed=0.56; 
        else if crpq=3 then crpmed=0.97; 
        else if crpq=4 then crpmed=1.74; 
	else if crpq=5 then crpmed=4.15; 

if cholq=1 then cholmed=162; 
else if cholq=2 then cholmed=184; 
        else if cholq=3 then cholmed=201; 
        else if cholq=4 then cholmed=220; 
        else if cholq=5 then cholmed=249; 

if trigq=1 then trigmed=63; 
else if trigq=2 then trigmed=90; 
        else if trigq=3 then trigmed=116; 
        else if trigq=4 then trigmed=156; 
        else if trigq=5 then trigmed=251.5; 

if hdlq=1 then hdlmed=31.9; 
else if hdlq=2 then hdlmed=38.8; 
        else if hdlq=3 then hdlmed=44.8; 
        else if hdlq=4 then hdlmed=51.1; 
        else if hdlq=5 then hdlmed=62.5; 

if hba1cq=1 then hba1cmed=5.017; 
else if hba1cq=2 then hba1cmed=5.313; 
        else if hba1cq=3 then hba1cmed=5.502; 
        else if hba1cq=4 then hba1cmed=5.6785; 
        else if hba1cq=5 then hba1cmed=5.9025;
/* calculate cardiometabolic disease risk score */
mscore=hba1cq+crpq+cholq+trigq+(6-hdlq);
/* quintiles for cardiometabolic disease risk score using 
cutoffs calculated from controls only*/
if .<mscore=<12 then mscoreq=1;
else if 12<mscore=<14 then mscoreq=2;
else if 14<mscore=<17 then mscoreq=3;
else if 17<mscore=<19 then mscoreq=4;
else if mscore>19 then mscoreq=5;
* create score variable (quartile-specific median) for quintiles of 
cardiometabolic disease risk score used in trend test;
if mscoreq=1 then mscoremed=10;
else if mscoreq=2 then mscoremed=13;
else if mscoreq=3 then mscoremed=15;
else if mscoreq=4 then mscoremed=17;
else if mscoreq=5 then mscoremed=20;
* create per-unit and per-SD variables for 
cardiometabolic disease risk score and biomarkers;
mscore2=mscore/2;
mscoresd=mscore/3.9760251;
hba1csd=hba1cmed/0.3477154;
trigsd=trigmed/101.4668594;
hdlsd=hdlmed/12.7048773;
cholsd=cholmed/36.1669291;
crpsd=crpmed/5.4667112;
run;
/* logistic regressions for calculating RRs */
/* cardiometabolic disease risk score */
/* model 1 */
proc logistic data=caco_mi;
class mscoreq (ref=first);
model cc= mscoreq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= mscoremed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= mscore2 &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= mscoresd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class mscoreq (ref=first);
model cc= mscoreq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= mscoremed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= mscore2 &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= mscoresd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;
/* CRP */
/* model 1 */
proc logistic data=caco_mi;
class crpq (ref=first);
model cc= crpq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= crpmed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= crpsd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class crpq (ref=first);
model cc= crpq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= crpmed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= crpsd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;
/* TG */
/* model 1 */
proc logistic data=caco_mi;
class trigq (ref=first);
model cc= trigq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= trigmed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= trigsd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class trigq (ref=first);
model cc= trigq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= trigmed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= trigsd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;
/* TC */
/* model 1 */
proc logistic data=caco_mi;
class cholq (ref=first);
model cc= cholq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= cholmed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= cholsd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class cholq (ref=first);
model cc= cholq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= cholmed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= cholsd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;
/* HbA1c */
/* model 1 */
proc logistic data=caco_mi;
class hba1cq (ref=first);
model cc= hba1cq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= hba1cmed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= hba1csd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class hba1cq (ref=first);
model cc= hba1cq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= hba1cmed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= hba1csd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;
/* HDL-C */
/* model 1 */
proc logistic data=caco_mi;
class hdlq (ref=first);
model cc= hdlq &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= hdlmed &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;

proc logistic data=caco_mi;
model cc= hdlsd &agecat_ &smoke_ &bldq_  &hbatch_ ;
run;
/* model 2 */
proc logistic data=caco_mi;
class hdlq (ref=first);
model cc= hdlq &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= hdlmed &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

proc logistic data=caco_mi;
model cc= hdlsd &agecat_ &smoke_ &bldq_  &hbatch_ 
&fammi_ &actq_ &nsaid_  &alcog_ &bmicat_;
run;

/* Calculate case numbers in each quintile of cardiometabolic disease risk score and biomarkers */
proc freq data=total;tables cc*mscoreq;run;
proc freq data=total;tables cc*crpq;run;
proc freq data=total;tables cc*trigq;run;
proc freq data=total;tables cc*cholq;run;
proc freq data=total;tables cc*hba1cq;run;
proc freq data=total;tables cc*hdlq;run;
endsas;
