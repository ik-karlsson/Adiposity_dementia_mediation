////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: BMI_dem_mediation_Cox_20210224_mid
Study: FORTE postdoc aim 2: Identify potential mediator of the associatin between BMI and
dementia
Purpose: Survival analysis of BMI on dementia, and cholesterol and inflammation as mediators
NOTE! This code is for late-life BMI - separate code for midlife

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20201207 by Ida Karlsson (IK)
	Updated: 	20210224 by IK - New data, include PRS for potential mediators
				20210707 by IK - Fix error in fasting variable
				20210824 by IK - Run sex-stratified analyses
				20220303 by IK - Refine outlier removal
				
STATA v15.1		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
log using session, replace
set more off

use "P:\Dementia_IK\FORTE_postdoc\Aim_2\Data\Derived_data\BMI_dem_mediation_20210707_late.dta", clear
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

// Setting CRP >100 to missing, as it indicates bacterial infection
replace crp=. if crp > 100
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
encode study_main, gen(study_n)

// study as dummy
gen satsa=0
gen octo=0
gen gender=0
gen harmony=0
gen twingene=0
replace satsa=1 if study_main=="SATSA"
replace octo=1 if study_main=="OCTOtw"
replace gender=1 if study_main=="GENDER"
replace twingene=1 if study_main=="TwinGene"

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
// 226 deleted

// Check remaining nr of cases
tab dementia
// 1000 dementia cases remaining, 6257 ctrls

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
       z_bmi |   .9721786   .0630324    -0.44   0.663     .8561648    1.103913
       z_whr |   1.055938    .050592     1.14   0.256     .9612925    1.159901
------------------------------------------------------------------------------*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr educ if(sex==1), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .9338655   .0679731    -0.94   0.347     .8097074    1.077062
       z_whr |   1.084585   .0573274     1.54   0.124     .9778496    1.202971
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20220303_late_men_`ad_measure'.dta, replace
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
regress z_crp_l z_bmi educ i.study_n fasting if(sex==1), vce(cluster pairid)
regress z_crp_l z_whr educ i.study_n if(sex==1), vce(cluster pairid)
/*----------------------------------------------------------------------------
             |               Robust
     z_crp_l | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .1839802   .0216514     8.50   0.000     .1415228    .2264376
       z_whr |   .1940886   .0184592    10.51   0.000     .1578908    .2302865
------------------------------------------------------------------------------*/

stcox z_bmi educ if(z_crp_l!=. & sex==1), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(z_crp_l!=. & sex==1), vce(cluster pairid) strata(study_main)
/*------------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .9921833   .0678155    -0.11   0.909     .8677857    1.134413
       z_whr |   1.075768   .0545789     1.44   0.150     .9739421     1.18824
------------------------------------------------------------------------------*/
	   
stcox z_crp_l educ if(z_bmi!=. & sex==1), vce(cluster pairid) strata(study_main)
stcox z_crp_l educ if(z_whr!=. & sex==1), vce(cluster pairid) strata(study_main)
/*------------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
     z_crp_l |   .8835286    .049438    -2.21   0.027     .7917563    .9859383
     z_crp_l |   .8825776   .0496541    -2.22   0.026     .7904311    .9854662
------------------------------------------------------------------------------*/

stcox z_bmi z_crp_l educ if(sex==1), vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l educ if(sex==1), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.013491   .0699407     0.19   0.846     .8852765    1.160276
     z_crp_l |    .881993   .0499077    -2.22   0.026     .7894049    .9854406
       z_whr |   1.102197   .0565922     1.90   0.058     .9966767    1.218889
     z_crp_l |   .8659607   .0496653    -2.51   0.012     .7738902     .968985
------------------------------------------------------------------------------*/

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
       z_bmi |   .9002047   .0368104    -2.57   0.010      .830873    .9753217
       z_whr |   .9720087    .045328    -0.61   0.543     .8871065    1.065037
 ------------------------------------------------------------------------------*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr educ if(sex==2), vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .8880573   .0387228    -2.72   0.006     .8153146    .9672901
       z_whr |   1.017454   .0505319     0.35   0.728     .9230813    1.121476
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20220303_late_women_`ad_measure'.dta, replace
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
/*------------------------------------------------------------------------------
             |               Robust
     z_crp_l | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .3066194   .0150354    20.39   0.000     .2771366    .3361021
       z_whr |   .2478852   .0171799    14.43   0.000     .2141972    .2815733
------------------------------------------------------------------------------*/
	   
stcox z_bmi educ if(z_crp_l!=. & sex==2), vce(cluster pairid) strata(study_main)
stcox z_whr educ if(z_crp_l!=. & sex==2), vce(cluster pairid) strata(study_main)
/*------------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .8944408     .04159    -2.40   0.016       .81653    .9797855
       z_whr |   .9849069   .0522988    -0.29   0.775     .8875569    1.092935
------------------------------------------------------------------------------*/
	   
stcox z_crp_l educ if(z_bmi!=. & sex==2), vce(cluster pairid) strata(study_main)
stcox z_crp_l educ if(z_whr!=. & sex==2), vce(cluster pairid) strata(study_main)
/*------------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
     z_crp_l |   .9207776   .0463469    -1.64   0.101     .8342762    1.016248
     z_crp_l |   .9241309   .0465996    -1.56   0.118     .8371656     1.02013
------------------------------------------------------------------------------*/

stcox z_bmi z_crp_l educ if(sex==2), vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l educ if(sex==2), vce(cluster pairid) strata(study_main)
/*------------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .9080419    .043952    -1.99   0.046     .8258575    .9984049
     z_crp_l |   .9529305   .0498592    -0.92   0.357     .8600519    1.055839
       z_whr |    1.00035    .054153     0.01   0.995     .8996487    1.112323
     z_crp_l |   .9240684   .0474736    -1.54   0.124     .8355531    1.021961
------------------------------------------------------------------------------*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20220303_late_bysex.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
