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
				20220302 by IK - Sensitivity analysis dropping controls with lastage <70
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

// Setting BMI below 15 or above 55 as missing (as default)
replace bmi=. if bmi < 15
replace bmi=. if bmi > 55
// 2 changes
// These 3 are outliers, and the BMI does not correspond. Set to missing.
replace whr=. if whr > 1.35

drop if bmi==. & whr==.
// n=0 
// drop those with no mediator measures
drop if crp==. & tc==. & hdl==. & ldl==. & tg==.
// 18 dropped

// Setting CRP >100 to missing, as it indicates bacterial infection
replace crp=. if crp > 100
// Gen logged value
gen crp_l = log(crp)

// And logged TG
gen tg_l = log(tg)

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

/// Dropping controls with last age < 70
drop if dementia==0 & lastage < 70
// n=3,232

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
// 110 dementia cases remaining, 2,667  ctrls

// Compare age at follow-up and death between cases and controls
bysort dementia: summarize age_exit
/* -> dementia = 0
    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
    age_exit |      2,667    72.74766    3.849629         70         95

-----------------------------------------------------------------------
-> dementia = 1
    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
    age_exit |        110    74.89091    7.211542         60         88	*/

bysort dementia: summarize age_death
/*-> dementia = 0
    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
   age_death |        345    77.95072    6.997542         70         97

-----------------------------------------------------------------------
-> dementia = 1
    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
   age_death |         80     81.3875    6.809327         65         92	*/
   
stset age_exit, failure(dementia==1) id(twinnr) enter(age)

list twinnr age dem_onset age_exit _st _t0 _t _d dementia in 1/20, noobs
stdes

drop if _st!=1
// n=0

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										/////  SURVIVAL MODELS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Base model of BMI and WHR, full sample
foreach ad_measure in z_bmi z_whr {
	stcox `ad_measure' i.sex educ, vce(cluster pairid) strata(study_main)
// N of individuals with adiposity + biomarker measure
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		sum `exp_var' if(`ad_measure'!=.)
	}
}
/*No. of subjects =       2,777                           Number of obs =  2,777
No. of failures =         110
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_bmi |    1.12473   .0937562     1.41   0.159      .955197    1.324352
------------------------------------------------------------------------------

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      2,338    .0572838    1.002135  -1.973062    3.80082
        z_tc |      2,776    .0554001    1.040199   -2.94952    4.98959
       z_hdl |      2,740   -.0184406    .9658394  -2.423786   4.895744
       z_ldl |      2,393     .018569    1.020166  -2.828485   4.582056
      z_tg_l |      2,776    .0533257    .9597477  -2.484768    4.77826
	  
No. of subjects =       2,642                           Number of obs =  2,642
No. of failures =          92
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_whr |   1.237645   .1248862     2.11   0.035     1.015557      1.5083
------------------------------------------------------------------------------
                                                      
    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      2,306    .0662496    .9995689  -1.973062    3.80082
        z_tc |      2,641    .0346488    1.023546   -2.94952    4.98959
       z_hdl |      2,609   -.0216509    .9651019  -2.423786   4.895744
       z_ldl |      2,359    .0262869    1.019283  -2.828485   4.582056
      z_tg_l |      2,642    .0527873    .9611413  -2.484768    4.77826			*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr i.sex educ, vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   1.100807   .1071348     0.99   0.324     .9096388    1.332151
       z_whr |   1.196209   .1279029     1.68   0.094     .9700489    1.475097
------------------------------------------------------------------------------*/

// Loop over models for BMI, WHR, and biomarkers
foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom inter_bmi inter_biom inter_bmixbiom tertbiom_c1 tertbiom_c2 tertbiom_c3 tertdiff_pval) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20220303_mid_sens70_Cox_`ad_measure'.dta, replace
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
/*      z_bmi |   .3600998   .3220258    .3981737
        z_whr |   .2670211   .2242283    .3098139 */
	   
stcox z_bmi i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
stcox z_whr i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.149201   .9358427    1.411202
       z_whr |    1.25742   .9579568    1.650498	*/
	   
stcox z_crp_l i.sex educ if(z_bmi!=.), vce(cluster pairid) strata(study_main)
stcox z_crp_l i.sex educ if(z_whr!=.), vce(cluster pairid) strata(study_main)
/*      z_crp_l |    .856662   .6238614    1.176335
		z_crp_l |   .8997015   .6520883    1.241339	*/

stcox z_bmi z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
/*     z_bmi |   1.262149   .9886192    1.611359
     z_crp_l |   .7831341   .5459782    1.123303
       z_whr |   1.321117   1.007319    1.732669
     z_crp_l |   .8309857   .5900974    1.170209	*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// MEDIATION ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

med4way z_whr z_crp_l bsex educ study_n, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0)
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0206239   .0170598    -1.21   0.227    -.0540604    .0128126
   ereri_cde |  -.0239898    .016233    -1.48   0.139    -.0558058    .0078263
ereri_intref |  -.0005334   .0039462    -0.14   0.892    -.0082677    .0072009
ereri_intmed |  -.0005221   .0038372    -0.14   0.892    -.0080429    .0069986
   ereri_pie |   .0044214    .004422     1.00   0.317    -.0042456    .0130884
------------------------------------------------------------------------------*/

med4way z_whr z_tc bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0201017   .0108425    -1.85   0.064    -.0413526    .0011492
   ereri_cde |  -.0226592    .011255    -2.01   0.044    -.0447186   -.0005998
ereri_intref |   .0023262   .0018832     1.24   0.217    -.0013647    .0060172
ereri_intmed |   -.000249   .0003086    -0.81   0.420    -.0008538    .0003558
   ereri_pie |   .0004802   .0005459     0.88   0.379    -.0005897    .0015502
------------------------------------------------------------------------------*/

med4way z_bmi z_hdl bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |    -.00328   .0100096    -0.33   0.743    -.0228984    .0163384
   ereri_cde |  -.0121068   .0090141    -1.34   0.179     -.029774    .0055604
ereri_intref |   .0015895   .0018238     0.87   0.383    -.0019851    .0051641
ereri_intmed |   .0021863   .0023711     0.92   0.356     -.002461    .0068337
   ereri_pie |    .005051   .0027318     1.85   0.064    -.0003033    .0104053
------------------------------------------------------------------------------*/

med4way z_whr z_ldl bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |  -.0122152   .0159002    -0.77   0.442     -.043379    .0189485
   ereri_cde |  -.0168778   .0142668    -1.18   0.237    -.0448403    .0110847
ereri_intref |   .0042206   .0057217     0.74   0.461    -.0069937    .0154349
ereri_intmed |   .0002755   .0004045     0.68   0.496    -.0005173    .0010683
   ereri_pie |   .0001665   .0003474     0.48   0.632    -.0005143    .0008473
------------------------------------------------------------------------------*/

med4way z_whr z_tg_l bsex educ study_n fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   -.017648   .0107676    -1.64   0.101    -.0387521     .003456
   ereri_cde |  -.0168949   .0109726    -1.54   0.124    -.0384008    .0046109
ereri_intref |    .000432   .0022185     0.19   0.846    -.0039161    .0047802
ereri_intmed |   .0005663   .0029049     0.19   0.845    -.0051271    .0062598
   ereri_pie |  -.0017515   .0033436    -0.52   0.600    -.0083048    .0048019
------------------------------------------------------------------------------*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20220303_mid_sens70.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

