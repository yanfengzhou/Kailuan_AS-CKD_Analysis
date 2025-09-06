/*import data*/
data  PWV_CKD;set 'C:\Users\Documents\Kailuan\PWV_CKD.sas7bdat';run;

data PWV_CKD; set PWV_CKD;
map=dbp+ (sbp - dbp)/3;
if sbp >= 140 or dbp >= 90 or self_hypertension=1 or anti_hy=1 then hypertension = 1;
    else hypertension = 0;
if fbg >=7.0 or self_diabetes=1 or anti_dia=1 then diabetes=1;
    else diabetes=0;
if sex=1 and hemoglobin<12.0 then anemia=1;
if sex=2 and hemoglobin<13.0 then anemia=1;
    else anemia=0;
pwv_measdate=measdatesy;
start=pwv_measdate;
stop=ckd_measdate;
format pwv_measdate start stop date9.;
t_ckd=(stop-start)/365.25; 
age=(pwv_measdate-birthdate)/365.25;
run;

proc rank data=PWV_CKD out=PWV_CKD group=5;
var pwv;
ranks pwv5; 
run;

/*table 1 baseline characteristics of participants*/
proc sort data=PWV_CKD;
by ID pwv_measdate;
run;

/*Remove duplicates by ID and generate the baseline database*/
proc sort data=PWV_CKD nodupkey out=PWV_baseline;
by ID;run;

/*Import ggBaseline macro for baseline table analysis*/
%include "C:\Users\Documents\Kailuan\ggBaseline2.sas";
%ggBaseline2(data=PWV_baseline,
var=pwv|anova|pwv\ckd|chisq|ckd\
sex|chisq|sex\age|anova|age\bmi|anova|bmi\
edu|chisq|edu\marri|chisq|marri\
physical|chisq|physical\smoking|chisq|smoking\alcohol|chisq|alcohol\
salt|chisq|salt\hypertension|chisq|hypertension\diabetes|chisq|diabetes\anti_hy|chisq|anti_hy\anti_dia|chisq|anti_dia\statin|chisq|statin\
sbp|anova|sbp\dbp|anova|dbp\map|anova|map\
fbg|anova|fbg\tg|krswls|tg\
ldl|anova|ldl\
hdl|anova|hdl\
egfr|anova|egfr\
crp|krswls|crp\anemia|chisq|anemia\,

grp=pwv5,
grplabel=0|1|2|3|4,
stdiff=y,
totcol=y,
pctype=col,
exmisspct=y,
showp=y,
filetype=rtf,
file=C:\Users\Documents\Kailuan\Baseline characteristics\table1,
title=baseline characteristics of participants,
footnote=,
fnspace=,
page=portrait,
deids=N
);
ods pdf close;
ods pdf file="C:\Users\Documents\Kailuan\Baseline characteristics\table1"; 


/*Table 2 baPWV at baseline*/
proc means data=PWV_baseline n nmiss mean std min max sum;
var pwv;
class pwv5;
run;
proc freq data=PWV_baseline;
table pwv5*ckd;
run;

/*model 1*/
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class pwv5(ref='4')  ;
model (start,stop)*ckd(0)=pwv5 
/rl;id ID;
output out=res  ressch=sch;
run;
/*model 2*/
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class pwv5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start,stop)*ckd(0)=pwv5 age map bmi sex marri edu
smoke alcohol physical salt tg ldl hdl fbg crp egfr
/rl;id ID;
output out=res  ressch=sch;
run;
/*model 3*/
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class pwv5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1');
model (start,stop)*ckd(0)=pwv5 age map bmi sex marri edu
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes anti_hy anti_dia statin anemia
/rl;id ID;
output out=res  ressch=sch;
run; 

/*P trend*/
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
model (start,stop)*ckd(0)=pwv5 
/rl;id ID;
output out=res  ressch=sch;
run;
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start,stop)*ckd(0)=pwv5 age map bmi  sex marri  edu
smoke alcohol physical salt tg ldl hdl fbg crp egfr
/rl;id ID;
output out=res  ressch=sch;
run;
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv5 age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia  marri statin anemia
/rl;id ID;
output out=res  ressch=sch;
run;
/*per SD increment*/
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
model (start,stop)*ckd(0)=pwvsd 
/rl;id ID;
output out=res  ressch=sch;
run;
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start,stop)*ckd(0)=pwvsd age map bmi  sex marri  edu
smoke alcohol physical salt tg ldl hdl fbg crp egfr
/rl;id ID;
output out=res  ressch=sch;
run;
proc phreg data=PWV_baseline covout noprint outest=res covs(aggregate); by _imputation_;
class  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwvsd age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia  marri statin anemia
/rl;id ID;
output out=res  ressch=sch;
run;


/*Table 3. Time-dependent association*/
proc sort data=PWV_CKD;
by ID decending pwv_measdate;
run;

data PWV_CKD1;set PWV_CKD;
ckd_measdate2=lag(ckd_measdate);
pwv2=lag(pwv);
pwv_measdate2=ckd_measdate;
format ckd_measdate2 pwv_measdate2 date9.;
run;

data PWV_CKD2; set PWV_CKD1;
if first.ID then pwv_measdate2=.;
by ID;
run;

data PWV_CKD3; set PWV_CKD2;
if pwv_measdate2=. then delete; 
run;
/*Generate progression of arterial stiffness*/
data pwv_pro; set PWV_CKD3;
pwv_progression=(pwv2-pwv)/((pwv_measdate2-pwv_measdate)/365.25);
run;

data pwv_pro; set pwv_pro;
start=pwv_measdate2;
stop=ckd_measdate2;
t_ckd=(stop-start)/365.25; 
format start stop date9.;
run;

/*quintiles*/
proc rank data=pwv_pro out=pwv_pro group=5;
var pwv_progression;
ranks pwv_progression5; 
run;
proc means data=pwv_pro n nmiss mean std min max;
var  pwv_progression;
class pwv_progression5;
run;
proc freq data=pwv_pro;
table pwv_progression5*ckd;
run;

/*Change in baPWV (category, cm/s)*/
data pwv_pro;set pwv_pro;
if pwv<1400 and pwv2<1400 then pwv_category=1;
if pwv<1400 and pwv2>=1400 then pwv_category=2;
if pwv>=1400 and pwv2<1400 then pwv_category=3;
if pwv>=1400 and pwv2>=1400 then pwv_category=4;
pwv_diff=pwv2-pwv; 
run;

proc means data=pwv_pro n nmiss mean std min max;
var pwv_diff;
class pwv_category;
run;

proc freq data=pwv_pro;
table pwv_category*ckd;
run;

/*multiple-imputation*/
proc mi data=pwv_pro seed=20221111 nimpute=5 out=pwv_mi;
  class ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia ;
  fcs nbiter=50;
      discrim(ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia / details)
      reg(pwv age map bmi tg ldl hdl fbg crp egfr/details);
  var ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia pwv age map bmi tg ldl hdl fbg crp egfr;
mcmc plots=trace;
run;

/*Changes in baPWV (progression of baPWV, cm/s/year)*/
/*model 1*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')   ;
model (start, stop)*ckd(0)=pwv_progression5 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3;
run;

/*model 2*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
 marri edu smoke alcohol physical salt tg ldl hdl fbg crp egfr
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi  sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr;
run;

/*model 3*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') hypertension(ref='1') diabetes(ref='1') 
anti_hy(ref='1') anti_dia(ref='1') statin(ref='1') anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi  sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*P trend*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
model (start, stop)*ckd(0)=pwv_progression5 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5;
run;

proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
  edu marri 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5 pwv age map bmi  sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr;
run;

proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5 pwv age map bmi  sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*Per SD*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
model (start, stop)*ckd(0)=pwv_progression_sd 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd;
run;

proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression_sd  pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
  edu marri 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr;
run;

proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression_sd pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*Change in baPWV (category, cm/s)*/
/*model 1*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_category(ref='4')   ;
model (start, stop)*ckd(0)=pwv_category 
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_category1 pwv_category2 pwv_category3;
run;

/*model 2*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_category(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') ;
model (start, stop)*ckd(0)=pwv_category pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
  edu   marri
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_category1 pwv_category2 pwv_category3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr;
run;

/*model 3*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_category(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_category pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia statin  marri anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_category1 pwv_category2 pwv_category3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*Table 4 subgroup analyses*/
/*baseline pwv*/
data g1;set pwv_mi;if pwv<1400 then yz1=1;if pwv>=1400 then yz1=2;run;
proc freq data=g1;table yz1*ckd;run;
/*P for interaction*/
proc phreg data=g1 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*yz1 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*yz1 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1a;set g1;if yz1=1;run;
proc phreg data=g1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1');
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1b;set g1;if yz1=2;run;
proc phreg data=g1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*bmi*/
data g1;set pwv_mi;if bmi<23 then yz2=1;if bmi>=23 then yz2=2;run;
proc freq data=g1;table yz2*ckd;run;
/*P for interaction*/
proc phreg data=g1 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*yz2 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*yz2 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1a;set g1;if yz2=1;run;
proc phreg data=g1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1b;set g1;if yz2=2;run;
proc phreg data=g1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*sex*/
proc freq data=pwv_pro;table sex*ckd;run;
/*P for interaction*/
proc phreg data=pwv_pro  covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*sex pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*sex pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1a;set pwv_mi;if sex=1;run;
proc phreg data=g1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')   marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1b;set pwv_mi;if sex=2;run;
proc phreg data=g1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*age*/
data g1;set pwv_mi;if age<60 then yz4=1;if age>=60 then yz4=2;run;
proc freq data=g1;table yz4*ckd;run;
/*P for interaction*/
proc phreg data=g1 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*yz4 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*yz4 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1a;set g1;if yz4=1;run;
proc phreg data=g1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1b;set g1;if yz4=2;run;
proc phreg data=g1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*pre_dia*/
data h1;set pwv_mi;
if diabetes=0 and fbg>=6.1 and fbg<7 then pre_dia=2;else pre_dia=0;
if diabetes=1 then pre_dia=1;
run;
proc freq data=h1;table pre_dia*ckd;run;
/*P for interaction*/
proc phreg data=h1 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*pre_dia pwv age map
bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension   edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*pre_dia pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data h1a;set h1;if pre_dia=0;run;
proc phreg data=h1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension   edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data h1b;set h1;if pre_dia=2;run;
proc phreg data=h1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension   edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data h1c;set h1;if pre_dia=1;run;
proc phreg data=h1c covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension   edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*hypertension*/
proc freq data=pwv_pro;table hypertension*ckd;run;
/*P for interaction*/
proc phreg data=pwv_pro covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1')  diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv_progression5*hypertension pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv_progression5*hypertension pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1a;set pwv_mi;if hypertension=0;run;
proc phreg data=g1a covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

data g1b;set pwv_mi;if hypertension=1;run;
proc phreg data=g1b covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start, stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*Sensitivity analysis by adjusting SBP rather than MAP*/
/*model 3*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1')  anemia(ref='1')
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start,stop)*ckd(0)=pwv_progression5 pwv age sbp bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age sbp bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P trend*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv_progression5 pwv age sbp bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5 pwv age sbp bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P SD*/
proc phreg data=pwv_mi covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv_progression_sd  pwv age sbp bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd pwv age sbp bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;


/*Sensitivity analysis by excluding participants with less than one year of follow-up*/
data pwv_pro1;set pwv_pro;
if t_ckd<1 then delete;
run;

proc rank data=pwv_pro1 out=pwv_pro1 group=5;
var pwv_progression;
ranks pwv_progression5; 
run;

proc means data=pwv_pro1 n nmiss std mean min max;
var pwv_progression;
run;

proc freq data=pwv_pro1;
table pwv_progression5*ckd;
run;

/*multiple-imputation*/
proc mi data=pwv_pro1 seed=20221111 nimpute=5 out=pwv_mi2;
  class ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia ;
  fcs nbiter=50;
      discrim(ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia / details)
      reg(pwv age map bmi tg ldl hdl fbg crp egfr/details);
  var ckd pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia pwv age map bmi tg ldl hdl fbg crp egfr;
mcmc plots=trace;
run;

/*model 3*/
proc phreg data=pwv_mi2 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1') anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia statin marri anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P trend*/
proc phreg data=pwv_mi2 covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P SD*/
proc phreg data=pwv_mi2 covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd(0)=pwv_progression_sd  pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*Sensitivity analysis by using a definition of CKD that 
requires at least two confirmed values of eGFR <60 ml/min/1.73 m2 or albuminuria*/

/*import data, the CKD_data database includes the definition of CKD at each follow-up visit*/
data CKD_data;set 'C:\Users\Documents\Kailuan\CKD_data.sas7bdat';
run;

proc means data=CKD_data sum;
var ckd_definiton;
class ID;
output out=CKD_data2 sum=ckd_definitonl; 
run;

proc sql;
create table pwv_pro2 as select * from pwv_pro a left join CKD_data2 b on a.ID=b.ID;
quit;

data pwv_pro2; set pwv_pro2;
if ckd=1 and ckd_definitonl>=2 then ckd1=1;
else ckd1=0;
run;

proc rank data=pwv_pro2 out=pwv_pro2 group=5;
var pwv_progression;
ranks pwv_progression5; 
run;

proc means data=pwv_pro2 n nmiss mean std min max;
var  pwv_progression;
run;

proc freq data=pwv_pro2;
table pwv_progression5*ckd1 ;
run;

/*multiple-imputation*/
proc mi data=pwv_pro2 seed=20221111 nimpute=5 out=pwv_mi3;
  class ckd1 pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia ;
  fcs nbiter=50;
      discrim(ckd1 pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia / details)
      reg(pwv age map bmi tg ldl hdl fbg crp egfr/details);
  var ckd1 pwv_progression5 sex marri edu salt smoke alcohol physical 
      hypertension diabetes anti_hy anti_dia statin anemia pwv age map bmi tg ldl hdl fbg crp egfr;
mcmc plots=trace;
run;

/*model 3*/
proc phreg data=pwv_mi3 covout noprint outest=res covs(aggregate); by _imputation_;
class pwv_progression5(ref='4')  sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') anemia(ref='1')
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  ;
model (start,stop)*ckd1(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression0 pwv_progression1 pwv_progression2 pwv_progression3 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P trend*/
proc phreg data=pwv_mi3 covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd1(0)=pwv_progression5 pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression5 pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;

/*P SD*/
proc phreg data=pwv_pro2 covout noprint outest=res covs(aggregate); by _imputation_;
class   sex(ref='1') marri(ref='1') edu(ref='1') 
salt(ref='1') smoke(ref='1') alcohol(ref='1') physical(ref='1') 
hypertension(ref='1') diabetes(ref='1') anti_hy(ref='1') anti_dia(ref='1') statin(ref='1')  anemia(ref='1') ;
model (start,stop)*ckd1(0)=pwv_progression_sd  pwv age map bmi  sex 
smoke alcohol physical salt tg ldl hdl fbg crp egfr
hypertension diabetes  edu anti_hy anti_dia marri statin anemia
/rl ties=exact risklimits;id ID;
run;
/*combined data*/
proc mianalyze data=res;
modeleffects pwv_progression_sd pwv age map bmi sex2 
 marri2 edu2 smoke0 alcohol0 physical0 salt2 salt3 tg ldl hdl fbg crp egfr
hypertension0 diabetes0 anti_hy0 anti_dia0 statin0 anemia0;
run;
