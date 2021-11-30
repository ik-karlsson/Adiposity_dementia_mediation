////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: BMI_dem_mediation_Cox_20210224_mid
Study: FORTE postdoc aim 2: Identify potential mediator of the associatin between BMI and
dementia
Purpose: Survival analysis of BMI on dementia, and cholesterol and inflammation as mediators
NOTE! This code is for midlife BMI - separate code for late-life

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20201207 by Ida Karlsson (IK)
	Updated: 	20210224 by IK - New data, include PRS for potential mediators
				20210506 by IK - Drop PRS analyses, add medaition
				20210707 by IK - Fix error in fasting variable
STATA v15.1		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
log using session, replace
set more off

use "P:\Dementia_IK\FORTE_postdoc\Aim_2\Data\Derived_data\BMI_dem_mediation_20210707_mid.dta", clear
rename *, lower

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										// DATA MANAGEMENT //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

summarize bmi
summarize whr

// Setting BMI below 15 or above 55 as missing
replace bmi=. if bmi < 15
replace bmi=. if bmi > 55
// 2 changes

//graph box whr
list twinnr bmi whr if(whr!=. & whr>1.35)
// These 3 are outliers, and the BMI does not correspond. Set to missing.
replace whr=. if whr > 1.35

drop if bmi==. & whr==.
// n=0 

// drop those with no mediator measures
drop if crp==. & tc==. & hdl==. & ldl==. & tg==.
// 18 dropped

/* Look at distribution of mediators 
graph box crp
graph box tc 
graph box hdl
graph box ldl
graph box tg

histogram crp
histogram tc
histogram hdl
histogram ldl
histogram tg
*/
// TC, HDL, and LDL look ok, but not TG and CRP. 
list twinnr tc hdl ldl tg fast if(tg>tc & tg!=.)
// n=24 - will recode all lipids to missing..
list twinnr tc hdl ldl tg fast if(hdl>tc & hdl!=.)
list twinnr tc hdl ldl tg fast if(ldl>tc & ldl!=.)
// n=0
gen tc_dup=tc
replace tc=. if tg>tc_dup & tg!=.
replace hdl=. if tg>tc_dup & tg!=.
replace ldl=. if tg>tc_dup & tg!=.
replace tg=. if tg>tc_dup & tg!=.
drop tc_dup

sum tc
sum hdl 
sum ldl 
sum tg
// Took care of a lot of TG outliers!

// Setting CRP >100 to missing, as it indicates bacterial infection
replace crp=. if crp > 100
//graph box crp
sum crp
// Gen logged value
gen crp_l = log(crp)
//histogram crp_l

// And logged TG
gen tg_l = log(tg)
//histogram tg_l

/// Standardizing variables for analyses
foreach y of varlist bmi tc hdl ldl crp_l tg_l   {
	egen z_`y' = std(`y')
	}

/// WHR is standardized by sex
egen whr_sexm= mean(whr), by(sex)
egen whr_sexsd= sd(whr), by(sex)
gen z_whr= (whr-whr_sexm)/whr_sexsd
	
// study as numeric
// NOTE! Only TG and SATSA, so generating binary
gen study_n=0
replace study_n=1 if study_main=="TwinGene"

// sex as binary
gen bsex=sex-1

// recode fasting to binary (and assign missing as fasting..)
recode fasting 2=0
recode fasting .=0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								///// SETTING UP THE SURVIVAL DATA /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate age at end of follow-up
gen age_exit=min(lastage, dem_onset)
list lastage dem_onset age_exit in 1/20
// Drop individuals who were already demented or has no registry info at baseline
drop if age >= age_exit
// 13 deleted

// Check remaining nr of cases
tab dementia
// 110 dementia cases remaining, 5889 ctrls

stset age_exit, failure(dementia==1) id(twinnr) enter(age)

list twinnr age dem_onset age_exit _st _t0 _t _d dementia in 1/20, noobs
stdes

drop if _st!=1
// n=0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									/////  DESCRIPTIVE STATISTICS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Convert education
gen educ_r=1 if educ==0
replace educ_r=0 if educ==1

cap postclose descriptive
	postfile descriptive str20 (exp_var all all_ps ctrl ctrl_ps dem_case dem_case_ps dem_p_value) ///
	using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210707_mid_descriptive.dta, replace

	foreach exp_var in sex educ_r smoke {
			set more off
			// All
			estpost tab `exp_var'
				mat b = e(b)
				mat pct = e(pct)
				local all = string(b[1,2],"%9.0f")  
				local all_ps = "("+string(pct[1,2],"%9.2f")+")" 
			// Dementia
			estpost tab `exp_var' dementia, chi2 nototal
				mat b = e(b)
				mat colpct = e(colpct)
				local ctrl = string(b[1,2],"%9.0f")  
				local ctrl_ps = "("+string(colpct[1,2],"%9.2f")+")"    
				local dem_case = string(b[1,4],"%9.0f")  
				local dem_case_ps = "("+string(colpct[1,4],"%9.2f")+")"    
				local dem_p_value = string(e(p))  
				
         post descriptive 	("`exp_var'") ("`all'") ("`all_ps'") ///
							("`ctrl'") ("`ctrl_ps'") ("`dem_case'") ("`dem_case_ps'") ("`dem_p_value'") 
        }
	
		foreach exp_var in age lastage age_death bmi whr crp tc hdl ldl tg {
			set more off
			// All
			mean `exp_var'
				mat b = e(b)
				mat sd  = e(sd)
				local all = string(b[1,1],"%9.2f")  
				local all_ps = "("+string(sd[1,1],"%9.2f")+")"    
			// Dementia
			ttest `exp_var', by(dementia)
				local ctrl = string(r(mu_1),"%9.2f")  
				local ctrl_ps = "("+string(r(sd_1),"%9.2f")+")"    
				local dem_case = string(r(mu_2),"%9.2f")  
				local dem_case_ps = "("+string(r(sd_2),"%9.2f")+")"    
				local dem_p_value = string(r(p))  

		post descriptive ("`exp_var'") ("`all'") ("`all_ps'") ("`ctrl'") ("`ctrl_ps'") ///
						 ("`dem_case'") ("`dem_case_ps'") ("`dem_p_value'") 
		} 
postclose descriptive

// Separate look at dementia onset
summarize dem_onset
/*  Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
   dem_onset |        110    74.89091    7.211542         60         88 */
   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										/////  SURVIVAL MODELS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// test ph assumption
// K-M curves and log-log plots
// By BMI
xtile tert_bmi = bmi, nquantiles(3)
sts graph, by(tert_bmi)
stphplot, by(tert_bmi)
// Bouncy but ok
// WHR
xtile tert_whr = whr, nquantiles(3)
sts graph, by(tert_whr)
stphplot, by(tert_whr)
// Also bouncy but ok

// Check with ph-test
stcox z_bmi z_whr educ smoke i.sex, vce(cluster pairid)
estat phtest, detail
// Looks good
stcox z_bmi z_whr z_crp_l z_tc z_hdl z_ldl z_tg_l educ smoke i.sex, vce(cluster pairid)
estat phtest, detail
// all look good½

// --> Happy with the PH assumption


// Base model of BMI and WHR, full sample
foreach ad_measure in z_bmi z_whr {
	stcox `ad_measure' i.sex educ, vce(cluster pairid) strata(study_main)
// N of individuals with adiposity + biomarker measure
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		sum `exp_var' if(`ad_measure'!=.)
	}
}
/*No. of subjects      =        5,997             Number of obs    =       5,997
No. of failures      =          110
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.128218   .0934562     1.46   0.145     .9591432    1.327096
------------------------------------------------------------------------------

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      5,396   -.0021856    .9989616  -1.973062    3.80082
        z_tc |      5,972    .0011623    1.000352  -2.957563   4.996911
       z_hdl |      5,931   -.0005055    .9988354  -2.434451   5.137273
       z_ldl |      5,527    .0015814    .9998895  -3.146436   5.318139
      z_tg_l |      5,972    .0000252    .9999838  -5.515275   3.852978
-------------+---------------------------------------------------------
	  
No. of subjects      =        5,803             Number of obs    =       5,803
No. of failures      =           92
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_whr |   1.248261   .1271342     2.18   0.029     1.022377    1.524051
       2.sex |    .860053   .1960027    -0.66   0.508     .5502243    1.344345
        educ |   1.074567   .2377284     0.33   0.745     .6964996    1.657855
------------------------------------------------------------------------------
                                                      
    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      5,311    .0022404     .998102  -1.973062    3.80082
        z_tc |      5,778   -.0107576    .9903179  -2.957563   4.996911
       z_hdl |      5,741   -.0011473     .999305  -2.434451   5.137273
       z_ldl |      5,439    .0035976    .9989248  -3.146436   5.318139
      z_tg_l |      5,779   -.0026878    1.001183  -5.515275   3.852978
-------------+---------------------------------------------------------
*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr i.sex educ, vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.101635   .1078665     0.99   0.323     .9092689    1.334698
       z_whr |   1.204796   .1316217     1.71   0.088     .9725712     1.49247
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom inter_bmi inter_biom inter_bmixbiom tertbiom_c1 tertbiom_c2 tertbiom_c3 tertdiff_pval) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210707_mid_Cox_`ad_measure'.dta, replace
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		set more off
		
		// BMI/WHR on the biomarker
		regress `exp_var' `ad_measure' i.sex educ fasting i.study_n, vce(cluster pairid)
		lincom `ad_measure' 
			local m_lin = r(estimate)
			local se_lin = r(se)
			local Lci_lin = `m_lin'-1.96*`se_lin'
			local Uci_lin = `m_lin'+1.96*`se_lin'
			local base_lin= string(`m_lin',"%9.2f")+" ("+string(`Lci_lin',"%9.2f")+"-"+string(`Uci_lin',"%9.2f")+")"  
			
		// Base model
		stcox `ad_measure' i.sex educ if(`exp_var'!=.), vce(cluster pairid) strata(study_main)
		lincom `ad_measure' 
			local m_m0 = r(estimate)
			local se_m0 = r(se)
			local Lci_m0 = exp(`m_m0'-1.96*`se_m0')
			local Uci_m0 = exp(`m_m0'+1.96*`se_m0')
			local base_bmi= string(exp(`m_m0'),"%9.2f")+" ("+string(`Lci_m0',"%9.2f")+"-"+string(`Uci_m0',"%9.2f")+")"  
			
		// Base model, biomarker
		stcox `exp_var' i.sex educ fasting  if(`ad_measure'!=.), vce(cluster pairid) strata(study_main)
		lincom `exp_var' 
			local m_m1 = r(estimate)
			local se_m1 = r(se)
			local Lci_m1 = exp(`m_m1'-1.96*`se_m1')
			local Uci_m1 = exp(`m_m1'+1.96*`se_m1')
			local base_biom= string(exp(`m_m1'),"%9.2f")+" ("+string(`Lci_m1',"%9.2f")+"-"+string(`Uci_m1',"%9.2f")+")"  
		
		// Adjusted for the biomarker
		stcox `ad_measure' `exp_var' i.sex educ fasting, vce(cluster pairid) strata(study_main)
		lincom `ad_measure' 
			local m_m2 = r(estimate)
			local se_m2 = r(se)
			local Lci_m2 = exp(`m_m2'-1.96*`se_m2')
			local Uci_m2 = exp(`m_m2'+1.96*`se_m2')
			local adj_bmi= string(exp(`m_m2'),"%9.2f")+" ("+string(`Lci_m2',"%9.2f")+"-"+string(`Uci_m2',"%9.2f")+")"  
		lincom `exp_var'
			local m_m3 = r(estimate)
			local se_m3 = r(se)
			local Lci_m3 = exp(`m_m3'-1.96*`se_m3')
			local Uci_m3 = exp(`m_m3'+1.96*`se_m3')
			local adj_biom= string(exp(`m_m3'),"%9.2f")+" ("+string(`Lci_m3',"%9.2f")+"-"+string(`Uci_m3',"%9.2f")+")"  

		// Interaction model
		stcox c.`ad_measure'##c.`exp_var' i.sex educ fasting, vce(cluster pairid) strata(study_main)
		lincom `ad_measure' 
			local m_m4 = r(estimate)
			local se_m4 = r(se)
			local Lci_m4 = exp(`m_m4'-1.96*`se_m4')
			local Uci_m4 = exp(`m_m4'+1.96*`se_m4')
			local inter_bmi= string(exp(`m_m4'),"%9.2f")+" ("+string(`Lci_m4',"%9.2f")+"-"+string(`Uci_m4',"%9.2f")+")"  
		lincom `exp_var'
			local m_m5 = r(estimate)
			local se_m5 = r(se)
			local Lci_m5 = exp(`m_m5'-1.96*`se_m5')
			local Uci_m5 = exp(`m_m5'+1.96*`se_m5')
			local inter_biom= string(exp(`m_m5'),"%9.2f")+" ("+string(`Lci_m5',"%9.2f")+"-"+string(`Uci_m5',"%9.2f")+")"  
		lincom c.`ad_measure'#c.`exp_var'
			local m_m6 = r(estimate)
			local se_m6 = r(se)
			local Lci_m6 = exp(`m_m6'-1.96*`se_m6')
			local Uci_m6 = exp(`m_m6'+1.96*`se_m6')
			local inter_bmixbiom= string(exp(`m_m6'),"%9.2f")+" ("+string(`Lci_m6',"%9.2f")+"-"+string(`Uci_m6',"%9.2f")+")"  
		
		/// By tertiles of the biomarker
		xtile tert_biom = `exp_var', nquantiles(3)
		stcox i.tert_biom c.`ad_measure'#i.tert_biom educ i.sex fasting, vce(cluster pairid) strata(study_main) 
		/// 1st quartile
		lincom 1.tert_biom#c.`ad_measure'
			local m_m7 = r(estimate)
			local se_m7 = r(se)
			local Lci_m7 = exp(`m_m7'-1.96*`se_m7')
			local Uci_m7 = exp(`m_m7'+1.96*`se_m7')
			local tertbiom_c1= string(exp(`m_m7'),"%9.2f")+" ("+string(`Lci_m7',"%9.2f")+"-"+string(`Uci_m7',"%9.2f")+")"  
		/// 2nd quartile
		lincom 2.tert_biom#c.`ad_measure'
			local m_m8 = r(estimate)
			local se_m8 = r(se)
			local Lci_m8 = exp(`m_m8'-1.96*`se_m8')
			local Uci_m8 = exp(`m_m8'+1.96*`se_m8')
			local tertbiom_c2= string(exp(`m_m8'),"%9.2f")+" ("+string(`Lci_m8',"%9.2f")+"-"+string(`Uci_m8',"%9.2f")+")"  
		/// 3rd quartile
		lincom 3.tert_biom#c.`ad_measure'
			local m_m9 = r(estimate)
			local se_m9 = r(se)
			local Lci_m9 = exp(`m_m9'-1.96*`se_m9')
			local Uci_m9 = exp(`m_m9'+1.96*`se_m9')
			local tertbiom_c3= string(exp(`m_m9'),"%9.2f")+" ("+string(`Lci_m9',"%9.2f")+"-"+string(`Uci_m9',"%9.2f")+")"  
		/// P-value for group difference
		stcox c.`ad_measure'##c.tert_biom educ i.sex fasting, vce(cluster pairid) strata(study_main) 
		lincom c.`ad_measure'#c.tert_biom
		local tertdiff_pval = string(r(p),"%9.2f")
		drop tert_biom

		post coxreg_dem ("`exp_var'") ("`base_lin'") ("`base_bmi'") ("`base_biom'") ("`adj_bmi'") ("`adj_biom'") ("`inter_bmi'") ("`inter_biom'") ///
						("`inter_bmixbiom'") ("`tertbiom_c1'") ("`tertbiom_c2'") ("`tertbiom_c3'") ("`tertdiff_pval'")  
	}
	postclose coxreg_dem
}         

/// Adding CRP without fasting adjustment..
regress z_crp_l z_bmi i.sex educ i.study_n, vce(cluster pairid)
regress z_crp_l z_whr i.sex educ i.study_n, vce(cluster pairid)
/*      z_bmi |   .3774341   .3519805    .4028878
        z_whr |   .2986336   .2715189    .3257484*/
	   
stcox z_bmi i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
stcox z_whr i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.161156   .9479229    1.422355
       z_whr |   1.270244   .9631204    1.675306*/
	   
stcox z_crp_l i.sex educ if(z_bmi!=.), vce(cluster pairid) strata(study_main)
stcox z_crp_l i.sex educ if(z_whr!=.), vce(cluster pairid) strata(study_main)
/*      z_crp_l |   .8699591   .6265946    1.207844
		z_crp_l |   .9161087   .6563016    1.278764*/

stcox z_bmi z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.268072   .9924599    1.620224
     z_crp_l |   .7920433   .5451674    1.150716
       z_whr |   1.328741   1.003598    1.759223
     z_crp_l |    .845485   .5932842    1.204895*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// MEDIATION ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

med4way z_whr z_crp_l bsex educ study_n, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0)
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0113233   .0087658    -1.29   0.196    -.0285038    .0058573
   ereri_cde |  -.0134992   .0081512    -1.66   0.098    -.0294752    .0024768
ereri_intref |  -.0001619   .0021256    -0.08   0.939    -.0043281    .0040042
ereri_intmed |  -.0001806   .0023646    -0.08   0.939    -.0048152     .004454
   ereri_pie |   .0025184   .0025686     0.98   0.327     -.002516    .0075528
------------------------------------------------------------------------------*/

med4way z_whr z_tc bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0173618   .0091572    -1.90   0.058    -.0353095    .0005859
   ereri_cde |   -.019265   .0094669    -2.03   0.042    -.0378198   -.0007103
ereri_intref |   .0018632   .0014255     1.31   0.191    -.0009307    .0046571
ereri_intmed |  -.0000459   .0001546    -0.30   0.767    -.0003489    .0002572
   ereri_pie |   .0000859    .000287     0.30   0.765    -.0004766    .0006484
------------------------------------------------------------------------------*/

med4way z_bmi z_hdl bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0120775   .0099641    -1.21   0.225    -.0316068    .0074518
   ereri_cde |  -.0136277   .0092416    -1.47   0.140     -.031741    .0044856
ereri_intref |  -.0005721   .0019214    -0.30   0.766    -.0043379    .0031936
ereri_intmed |  -.0005785   .0019363    -0.30   0.765    -.0043737    .0032167
   ereri_pie |   .0027008   .0023833     1.13   0.257    -.0019705     .007372
------------------------------------------------------------------------------*/

med4way z_whr z_ldl bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   -.007259   .0088869    -0.82   0.414     -.024677    .0101591
   ereri_cde |   -.010183    .007663    -1.33   0.184    -.0252022    .0048361
ereri_intref |   .0029544   .0034701     0.85   0.395    -.0038468    .0097557
ereri_intmed |   -.000019   .0000977    -0.19   0.846    -.0002106    .0001725
   ereri_pie |  -.0000113   .0000606    -0.19   0.852    -.0001301    .0001075
------------------------------------------------------------------------------*/

med4way z_whr z_tg_l bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0147951   .0091693    -1.61   0.107    -.0327665    .0031763
   ereri_cde |   -.014105   .0092524    -1.52   0.127    -.0322394    .0040294
ereri_intref |   .0005951   .0021873     0.27   0.786     -.003692    .0048822
ereri_intmed |   .0007453   .0027356     0.27   0.785    -.0046164     .006107
   ereri_pie |  -.0020305    .003144    -0.65   0.518    -.0081926    .0041316
------------------------------------------------------------------------------*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20210707_mid.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

