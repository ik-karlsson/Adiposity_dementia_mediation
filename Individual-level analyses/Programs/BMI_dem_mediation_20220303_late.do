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

summarize bmi
summarize whr

graph box bmi, name(raw_bmi, replace)
graph box whr, name(raw_whr, replace)
graph combine raw_bmi raw_whr
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BOX_raw_bmi_whr_late.tif", replace

// Setting BMI below 15 or above 55 as missing
replace bmi=. if bmi < 15
replace bmi=. if bmi > 55
// 2 changes

list twinnr bmi whr if(whr!=. & whr>1.35)
// 1 outlier, with missing BMI (dropped in QC?). Set to missing.
replace whr=. if whr > 1.35

drop if bmi==. & whr==.
// n=1

graph box bmi, name(clean_bmi, replace)
graph box whr, name(clean_whr, replace)
graph combine clean_bmi clean_whr
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BOX_clean_bmi_whr_late.tif", replace

// drop those with no mediator measures
drop if crp==. & tc==. & hdl==. & ldl==. & tg==.
// 198 dropped

/* Look at distribution of mediators */
graph box crp, name(box_raw_crp, replace)
graph box tc, name(box_raw_tc, replace) 
graph box hdl, name(box_raw_hdl, replace)
graph box ldl, name(box_raw_ldl, replace)
graph box tg, name(box_raw_tg, replace)

histogram crp, name(hist_raw_crp, replace)
histogram tc, name(hist_raw_tc, replace)
histogram hdl, name(hist_raw_hdl, replace)
histogram ldl, name(hist_raw_ldl, replace)
histogram tg, name(hist_raw_tg, replace)

graph combine box_raw_crp box_raw_tc box_raw_hdl box_raw_ldl box_raw_tg ///
			  hist_raw_crp hist_raw_tc hist_raw_hdl hist_raw_ldl hist_raw_tg, rows(2) xsize(20) ysize(10)
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BOX_raw_biomarkers_late.tif", replace

// TC, HDL, and LDL look ok, but not TG and CRP. 
list twinnr tc hdl ldl tg fast if(tg>tc & tg!=.)
// n=16 
list twinnr tc hdl tg fast if(hdl>tc & hdl!=.)
list twinnr tc hdl tg fast if(ldl>tc & ldl!=.)
// n=0


// Setting CRP >100 to missing, as it indicates bacterial infection
replace crp=. if crp > 100
//n=4
graph box crp
sum crp
// Gen logged value
gen crp_l = log(crp)
histogram crp_l

// And logged TG
gen tg_l = log(tg)
histogram tg_l
graph box tg_l
// One questionable, but not a clear outlier so will leave..

// Plot clean values
graph box crp_l, name(box_clean_crp, replace)
graph box tc, name(box_clean_tc, replace) 
graph box hdl, name(box_clean_hdl, replace)
graph box ldl, name(box_clean_ldl, replace)
graph box tg_l, name(box_clean_tg, replace)

histogram crp_l, name(hist_clean_crp, replace)
histogram tc, name(hist_clean_tc, replace)
histogram hdl, name(hist_clean_hdl, replace)
histogram ldl, name(hist_clean_ldl, replace)
histogram tg_l, name(hist_clean_tg, replace)

graph combine box_clean_crp box_clean_tc box_clean_hdl box_clean_ldl box_clean_tg ///
			  hist_clean_crp hist_clean_tc hist_clean_hdl hist_clean_ldl hist_clean_tg, rows(2) xsize(20) ysize(10)
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BOX_clean_biomarkers_late.tif", replace

// Boxplot by study 
graph box bmi, by(study_main) name(bmi_by_s, replace) 
graph box whr, by(study_main) name(whr_by_s, replace) 

graph box crp_l, by(study_main) name(crp_by_s, replace) 
graph box tc, by(study_main) name(tc_by_s, replace) 
graph box hdl, by(study_main) name(hdl_by_s, replace) 
graph box ldl, by(study_main) name(ldl_by_s, replace) 
graph box tg_l, by(study_main) name(tg_by_s, replace) 

graph combine bmi_by_s whr_by_s crp_by_s tc_by_s hdl_by_s ldl_by_s tg_by_s , rows(2) xsize(20) ysize(10)
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\Measures_by_study_late.tif", replace

// Look ok

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

// Look at the cumulative hazard function
sts graph, cumhaz
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\Cumulative_haz_late.tif", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									/////  DESCRIPTIVE STATISTICS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
// Convert education
gen educ_r=1 if educ==0
replace educ_r=0 if educ==1

cap postclose descriptive
	postfile descriptive str20 (exp_var all all_ps ctrl ctrl_ps dem_case dem_case_ps dem_p_value) ///
	using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20220303_late_descriptive.dta, replace

	foreach exp_var in sex educ_r smoke vitalstatus {
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
				mat V  = e(V)
				local all = string(b[1,1],"%9.2f")  
				local all_ps = "("+string(b[1,1],"%9.2f")+")"    
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
   dem_onset |      1,000      82.354    6.236793         67        101 */
   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										/////  SURVIVAL MODELS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
// test ph assumption
// K-M curves and log-log plots
// By BMI
xtile tert_bmi = bmi, nquantiles(3)
sts graph, by(tert_bmi) name(stsgraph_bmi, replace)
stphplot, by(tert_bmi) name(phplot_bmi, replace)
// Bouncy but ok
// WHR
xtile tert_whr = whr, nquantiles(3)
sts graph, by(tert_whr) name(stsgraph_whr, replace)
stphplot, by(tert_whr) name(phplot_whr, replace)
// Also bouncy but ok
graph combine stsgraph_bmi phplot_bmi stsgraph_whr phplot_whr, rows(2) xsize(30) ysize(30)
graph export "P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\PHassumption_test_late.tif", replace

// Check with ph-test
stcox z_bmi z_whr educ i.sex, vce(cluster pairid)
estat phtest, detail
// whr close at 0.08
stcox z_bmi z_whr z_crp_l z_tc z_hdl z_ldl z_tg_l educ i.sex, vce(cluster pairid)
estat phtest, detail
// all look fine

// --> Happy with the PH assumption

// Base model of BMI and WHR, full sample
foreach ad_measure in z_bmi z_whr {
	stcox `ad_measure' i.sex educ, vce(cluster pairid) strata(study_main)
// N of individuals with adiposity + biomarker measure
	foreach exp_var in z_crp_l z_tc z_hdl z_ldl z_tg_l {
		sum `exp_var' if(`ad_measure'!=.)
	}
}
/*
No. of subjects      =        7,256             Number of obs    =       7,256
No. of failures      =        1,000
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_bmi |   .9217965   .0318319    -2.36   0.018     .8614717    .9863456
       2.sex |   .9659459   .0638825    -0.52   0.600     .8485137     1.09963
        educ |   .8880413   .0575582    -1.83   0.067     .7821007    1.008332
------------------------------------------------------------------------------

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      6,470    -.007533    .9945532  -2.132165   3.603468
        z_tc |      7,252    .0004261    .9954632  -2.954753    5.02368
       z_hdl |      6,810    .0036899    1.000585  -2.423104   6.158111
       z_ldl |      6,061      .00515    .9981037  -2.900693   5.792454
      z_tg_l |      6,891   -.0018964    .9995702  -5.838809   5.126534
-------------+---------------------------------------------------------

No. of subjects      =        7,125             Number of obs    =       7,125
No. of failures      =          971
------------------------------------------------------------------------------
             |               Robust
          _t | Haz. Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
       z_whr |   1.006681    .033668     0.20   0.842     .9428087    1.074879
       2.sex |    .967148   .0647731    -0.50   0.618     .8481743     1.10281
        educ |   .9109322    .059935    -1.42   0.156     .8007207    1.036313
------------------------------------------------------------------------------

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
     z_crp_l |      6,416   -.0094702    .9930969  -2.132165   3.603468
        z_tc |      7,121   -.0063122    .9899812  -2.954753    5.02368
       z_hdl |      6,687     .005802     1.00045  -2.423104   6.158111
       z_ldl |      6,005    .0063554    .9984901  -2.900693   5.792454
      z_tg_l |      6,760   -.0070755    1.000334  -5.838809   5.126534
-------------+---------------------------------------------------------*/

// BMI, adjusted for WHR, and vice versa 
stcox z_bmi z_whr i.sex educ, vce(cluster pairid) strata(study_main)
/*----------------------------------------------------------------------------
             |               Robust
          _t | Haz. ratio   std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
       z_bmi |   .9021741   .0336227    -2.76   0.006      .838624      .97054
       z_whr |   1.047553   .0374396     1.30   0.194      .976684    1.123565
------------------------------------------------------------------------------*/


foreach ad_measure in z_bmi z_whr {

	cap postclose coxreg_dem
		postfile coxreg_dem str20  (exp_var base_lin base_bmi base_biom adj_bmi adj_biom inter_bmi inter_biom inter_bmixbiom tertbiom_c1 tertbiom_c2 tertbiom_c3 tertdiff_pval) ///
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20220303_late_Cox_`ad_measure'.dta, replace
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
		stcox `exp_var' i.sex educ fasting if(`ad_measure'!=.), vce(cluster pairid) strata(study_main)
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
/*     z_bmi |   .2649925   .2406988    .2892861
       z_whr |   .2224152   .1977302    .2471002*/
	   
stcox z_bmi i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
stcox z_whr i.sex educ if(z_crp_l!=.), vce(cluster pairid) strata(study_main)
/*     z_bmi |      .9253   .8582151    .9976287
       z_whr |   1.023286   .9524094    1.099437*/
	   
stcox z_crp_l i.sex educ if(z_bmi!=.), vce(cluster pairid) strata(study_main)
stcox z_crp_l i.sex educ if(z_whr!=.), vce(cluster pairid) strata(study_main)
/*   z_crp_l |   .9026683   .8391072     .971044
     z_crp_l |    .904104   .840272     .972785*/

stcox z_bmi z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
stcox z_whr z_crp_l i.sex educ, vce(cluster pairid) strata(study_main)
/*     z_bmi |   .9473872   .8772773      1.0231
     z_crp_l |   .9153312   .8495876    .9861623
       z_whr |   1.043641   .9704347    1.122369
     z_crp_l |   .8967574   .8325555    .9659101*/

	 
regress z_hdl z_bmi i.sex educ fasting i.study_n, vce(cluster pairid)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									///// MEDIATION ANALYSES /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0) fulloutput

/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0055768   .0034423     1.62   0.105      -.00117    .0123236
   ereri_cde |   .0036778   .0034099     1.08   0.281    -.0030055    .0103611
ereri_intref |  -.0000271   .0001751    -0.15   0.877    -.0003703    .0003161
ereri_intmed |  -.0001319   .0008587    -0.15   0.878    -.0018149    .0015511
   ereri_pie |   .0020579    .000889     2.31   0.021     .0003155    .0038004
      terira |   1.005577   .0034423   292.12   0.000       .99883    1.012324
       p_cde |    .659486   .2975246     2.22   0.027     .0763485    1.242624
    p_intref |  -.0048576   .0324133    -0.15   0.881    -.0683865    .0586713
    p_intmed |  -.0236488   .1589018    -0.15   0.882    -.3350906    .2877929
       p_pie |   .3690204   .2709886     1.36   0.173    -.1621075    .9001484
        op_m |   .3453716   .2844234     1.21   0.225    -.2120881    .9028313
      op_ati |  -.0285064   .1912563    -0.15   0.882    -.4033619     .346349
        op_e |    .340514   .2975246     1.14   0.252    -.2426234    .9236514
------------------------------------------------------------------------------*/

med4way z_bmi z_tc bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0067094    .003627     1.85   0.064    -.0003994    .0138182
   ereri_cde |   .0075136   .0032385     2.32   0.020     .0011662    .0138609
ereri_intref |  -.0005309   .0012172    -0.44   0.663    -.0029165    .0018547
ereri_intmed |  -.0000872   .0002007    -0.43   0.664    -.0004806    .0003063
   ereri_pie |  -.0001861   .0001979    -0.94   0.347    -.0005739    .0002017
------------------------------------------------------------------------------*/

med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0) fulloutput
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0032707   .0040637     0.80   0.421    -.0046941    .0112355
   ereri_cde |   .0104221    .003416     3.05   0.002     .0037268    .0171174
ereri_intref |  -.0029965   .0014766    -2.03   0.042    -.0058905   -.0001025
ereri_intmed |  -.0017097   .0008457    -2.02   0.043    -.0033673    -.000052
   ereri_pie |  -.0024452   .0009645    -2.54   0.011    -.0043356   -.0005548
      terira |   1.003271   .0040637   246.88   0.000     .9953059    1.011235
       p_cde |   3.186498   3.263562     0.98   0.329    -3.209967    9.582962
    p_intref |  -.9161714   1.454639    -0.63   0.529    -3.767212    1.934869
    p_intmed |  -.5227171   .8308865    -0.63   0.529    -2.151225    1.105791
       p_pie |   -.747609   1.015814    -0.74   0.462    -2.738568     1.24335
        op_m |  -1.270326   1.822693    -0.70   0.486     -4.84274    2.302087
      op_ati |  -1.438889   2.285104    -0.63   0.529     -5.91761    3.039833
        op_e |  -2.186498   3.263562    -0.67   0.503    -8.582962    4.209967
------------------------------------------------------------------------------*/

med4way z_bmi z_ldl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0073932   .0035221     2.10   0.036       .00049    .0142964
   ereri_cde |    .007567   .0034627     2.19   0.029     .0007801    .0143538
ereri_intref |  -.0002201   .0006127    -0.36   0.719    -.0014211    .0009808
ereri_intmed |  -.0000732   .0002048    -0.36   0.721    -.0004745    .0003282
   ereri_pie |   .0001195   .0001894     0.63   0.528    -.0002517    .0004907
------------------------------------------------------------------------------*/

med4way z_bmi z_tg_l bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0059837   .0032444     1.84   0.065    -.0003753    .0123427
   ereri_cde |    .008845   .0034722     2.55   0.011     .0020396    .0156503
ereri_intref |   .0000447   .0000969     0.46   0.644    -.0001452    .0002346
ereri_intmed |  -.0011349   .0009814    -1.16   0.247    -.0030584    .0007885
   ereri_pie |   -.001771    .000968    -1.83   0.067    -.0036682    .0001262
------------------------------------------------------------------------------*/

/* tereri=total excess relative risk; ereri_cde=excess relative risk due to controlled direct effect; ereri_intref=excess
   relative risk due to reference interaction; ereri_intmed=excess relative risk due to mediated interaction;
   ereri_pie=excess relative risk due to pure indirect effect
*/


/// CDE of +1 SD BMI at different levels of CRP
med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(-2) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0)
med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(-1) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0)
med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0)
med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(1) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0)
med4way z_bmi z_crp_l bsex educ satsa octo gender twingene, a0(0) a1(1) m(2) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0)
/* m= -2 to +2
   ereri_cde |   .0046552   .0073897     0.63   0.529    -.0098283    .0191387
   ereri_cde |   .0041703   .0048082     0.87   0.386    -.0052536    .0135941
   ereri_cde |   .0036778   .0034099     1.08   0.281    -.0030055    .0103611
   ereri_cde |   .0031777   .0045843     0.69   0.488    -.0058074    .0121628
   ereri_cde |   .0026698   .0072243     0.37   0.712    -.0114895    .0168292
*/   

   
/// CDE of +1 SD BMI at different levels of HDL
med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(-2) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(-1) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(1) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(2) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)

/* m= -2 to +2
   ereri_cde |  -.0020323   .0068585    -0.30   0.767    -.0154747      .01141
   ereri_cde |   .0041198   .0045029     0.91   0.360    -.0047057    .0129454
   ereri_cde |   .0104221    .003416     3.05   0.002     .0037268    .0171174
   ereri_cde |   .0168773   .0048203     3.50   0.000     .0074295     .026325
   ereri_cde |   .0234882   .0075595     3.11   0.002     .0086719    .0383044
*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20220303_late.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
