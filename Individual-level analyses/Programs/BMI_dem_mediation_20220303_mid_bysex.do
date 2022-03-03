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
				20210824 by IK - Run sex-stratified analyses
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

/// NOTE! Removing any non-essential code, as everything is set up as in main analyses

// Setting BMI below 15 or above 55 as missing
replace bmi=. if bmi < 15
replace bmi=. if bmi > 55
replace whr=. if whr > 1.35
drop if bmi==. & whr==.

// drop those with no mediator measures
drop if crp==. & tc==. & hdl==. & ldl==. & tg==.

gen tc_dup=tc
replace tc=. if tg>tc_dup & tg!=.
replace hdl=. if tg>tc_dup & tg!=.
replace ldl=. if tg>tc_dup & tg!=.
replace tg=. if tg>tc_dup & tg!=.
drop tc_dup

// Setting CRP >100 to missing, as it indicates bacterial infection
replace crp=. if crp > 100
gen crp_l = log(crp)

// And logged TG
gen tg_l = log(tg)
replace tg_l=. if tg_l<-2

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
										/////  SURVIVAL MODELS, MEN /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Base model of BMI and WHR, full sample
stcox z_bmi educ if(sex==1), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(sex==1), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.217668   .1927734     1.24   0.214     .8928367    1.660679
       z_whr |   1.496664   .2521531     2.39   0.017     1.075763    2.082246
------------------------------------------------------------------------------*/


// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr educ if(sex==1), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.036808   .1597869     0.23   0.815     .7665071    1.402429
       z_whr |   1.474113   .2519215     2.27   0.023     1.054539    2.060623
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210824_mid_men_`ad_measure'.dta, replace
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		set more off
		
		// BMI/WHR on the biomarker
		regress `exp_var' `ad_measure' educ fasting i.study_n if(sex==1), vce(cluster pairid)
		lincom `ad_measure' 
			local m_lin = r(estimate)
			local se_lin = r(se)
			local Lci_lin = `m_lin'-1.96*`se_lin'
			local Uci_lin = `m_lin'+1.96*`se_lin'
			local base_lin= string(`m_lin',"%9.2f")+" ("+string(`Lci_lin',"%9.2f")+"-"+string(`Uci_lin',"%9.2f")+")"  
			
		// Base model
		stcox `ad_measure'  educ if(`exp_var'!=. & sex==1), vce(cluster pairid) strata(study_main)
		lincom `ad_measure' 
			local m_m0 = r(estimate)
			local se_m0 = r(se)
			local Lci_m0 = exp(`m_m0'-1.96*`se_m0')
			local Uci_m0 = exp(`m_m0'+1.96*`se_m0')
			local base_bmi= string(exp(`m_m0'),"%9.2f")+" ("+string(`Lci_m0',"%9.2f")+"-"+string(`Uci_m0',"%9.2f")+")"  
			
		// Base model, biomarker
		stcox `exp_var'  educ fasting  if(`ad_measure'!=. & sex==1), vce(cluster pairid) strata(study_main)
		lincom `exp_var' 
			local m_m1 = r(estimate)
			local se_m1 = r(se)
			local Lci_m1 = exp(`m_m1'-1.96*`se_m1')
			local Uci_m1 = exp(`m_m1'+1.96*`se_m1')
			local base_biom= string(exp(`m_m1'),"%9.2f")+" ("+string(`Lci_m1',"%9.2f")+"-"+string(`Uci_m1',"%9.2f")+")"  
		
		// Adjusted for the biomarker
		stcox `ad_measure' `exp_var'  educ fasting if(sex==1), vce(cluster pairid) strata(study_main)
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

		post coxreg_dem ("`exp_var'") ("`base_lin'") ("`base_bmi'") ("`base_biom'") ("`adj_bmi'") ("`adj_biom'") 
	}
	postclose coxreg_dem
}         

/// Adding CRP without fasting adjustment..
regress z_crp_l z_bmi educ i.study_n if(sex==1), vce(cluster pairid)
regress z_crp_l z_whr educ i.study_n if(sex==1), vce(cluster pairid)
/*     z_bmi |   .3257541     .2835498    .3679584
       z_whr |   .2742516     .2329577    .3155455*/
	   
stcox z_bmi educ if(z_crp_l!=. & sex==1), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(z_crp_l!=. & sex==1), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.238713      .884262    1.735244
       z_whr |   1.508397     1.029996    2.209001*/
	   
stcox z_crp_l educ if(z_bmi!=. & sex==1), vce(cluster pairid) strata(study_main)
stcox z_crp_l educ if(z_whr!=. & sex==1), vce(cluster pairid) strata(study_main)
/*      z_crp_l |   1.035075     .5925954    1.807946
        z_crp_l |   1.123358     .6440828    1.959271*/

stcox z_bmi z_crp_l educ if(sex==1), vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l educ if(sex==1), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.254233      .798062     1.97115
     z_crp_l |   .9670577     .4995634    1.872036
       z_whr |   1.503372      1.00534    2.248121
     z_crp_l |   1.015949     .5471336    1.886471*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										/////  SURVIVAL MODELS, WOMEN /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Base model of BMI and WHR, full sample
stcox z_bmi educ if(sex==2), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(sex==2), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |    1.10091   .1066903     0.99   0.321     .9104606    1.331198
       z_whr |   1.131922   .1447541     0.97   0.333      .880972    1.454358
------------------------------------------------------------------------------*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr educ if(sex==2), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.141172   .1391912     1.08   0.279     .8985212    1.449353
       z_whr |   1.083331   .1525142     0.57   0.570     .8221045    1.427564
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210824_mid_women_`ad_measure'.dta, replace
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		set more off
		
		// BMI/WHR on the biomarker
		regress `exp_var' `ad_measure' educ fasting i.study_n if(sex==2), vce(cluster pairid)
		lincom `ad_measure' 
			local m_lin = r(estimate)
			local se_lin = r(se)
			local Lci_lin = `m_lin'-1.96*`se_lin'
			local Uci_lin = `m_lin'+1.96*`se_lin'
			local base_lin= string(`m_lin',"%9.2f")+" ("+string(`Lci_lin',"%9.2f")+"-"+string(`Uci_lin',"%9.2f")+")"  
			
		// Base model
		stcox `ad_measure'  educ if(`exp_var'!=. & sex==2), vce(cluster pairid) strata(study_main)
		lincom `ad_measure' 
			local m_m0 = r(estimate)
			local se_m0 = r(se)
			local Lci_m0 = exp(`m_m0'-1.96*`se_m0')
			local Uci_m0 = exp(`m_m0'+1.96*`se_m0')
			local base_bmi= string(exp(`m_m0'),"%9.2f")+" ("+string(`Lci_m0',"%9.2f")+"-"+string(`Uci_m0',"%9.2f")+")"  
			
		// Base model, biomarker
		stcox `exp_var'  educ fasting  if(`ad_measure'!=. & sex==2), vce(cluster pairid) strata(study_main)
		lincom `exp_var' 
			local m_m1 = r(estimate)
			local se_m1 = r(se)
			local Lci_m1 = exp(`m_m1'-1.96*`se_m1')
			local Uci_m1 = exp(`m_m1'+1.96*`se_m1')
			local base_biom= string(exp(`m_m1'),"%9.2f")+" ("+string(`Lci_m1',"%9.2f")+"-"+string(`Uci_m1',"%9.2f")+")"  
		
		// Adjusted for the biomarker
		stcox `ad_measure' `exp_var'  educ fasting if(sex==2), vce(cluster pairid) strata(study_main)
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

		post coxreg_dem ("`exp_var'") ("`base_lin'") ("`base_bmi'") ("`base_biom'") ("`adj_bmi'") ("`adj_biom'") 
	}
	postclose coxreg_dem
}         

/// Adding CRP without fasting adjustment..
regress z_crp_l z_bmi educ i.study_n if(sex==2), vce(cluster pairid)
regress z_crp_l z_whr educ i.study_n if(sex==2), vce(cluster pairid)
/*     z_bmi |   .4078047     .3765329    .4390765
       z_whr |   .3171314     .2814394    .3528234*/
	   
stcox z_bmi educ if(z_crp_l!=. & sex==2), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(z_crp_l!=. & sex==2), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.151245     .9050642    1.464388
       z_whr |   1.152808     .8042227    1.652486*/
	   
stcox z_crp_l educ if(z_bmi!=. & sex==2), vce(cluster pairid) strata(study_main)
stcox z_crp_l educ if(z_whr!=. & sex==2), vce(cluster pairid) strata(study_main)
/*   z_crp_l |   .7986668      .547838    1.164338
     z_crp_l |   .8278943     .5634444    1.216462*/

stcox z_bmi z_crp_l educ if(sex==2), vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l educ if(sex==2), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.318108     1.001405    1.734972
     z_crp_l |   .7048029     .4655697    1.066966
       z_whr |     1.2424     .8560617    1.803091
     z_crp_l |   .7748874     .5145235    1.167003
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20210824_mid_bysex.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

