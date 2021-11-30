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

// Setting BMI below 15 or above 55 as missing
replace bmi=. if bmi < 15
replace bmi=. if bmi > 55
// 2 changes

//graph box whr
list twinnr bmi whr if(whr!=. & whr>1.35)
// 1 outlier, with missing BMI (dropped in QC?). Set to missing.
replace whr=. if whr > 1.35

drop if bmi==. & whr==.
// n=1

// drop those with no mediator measures
drop if crp==. & tc==. & hdl==. & ldl==. & tg==.
// 198 dropped

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
// n=16 - will recode all lipids to missing..
list twinnr tc hdl tg fast if(hdl>tc & hdl!=.)
list twinnr tc hdl tg fast if(ldl>tc & ldl!=.)
// n=0
gen tc_dup=tc
replace tc=. if tg>tc_dup & tg!=.
replace hdl=. if tg>tc_dup & tg!=.
replace ldl=. if tg>tc_dup & tg!=.
replace tg=. if tg>tc_dup & tg!=.
drop tc_dup

// One LDL outlier
replace ldl=. if ldl>=9
//graph box ldl

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
// One outlier..
replace tg_l=. if tg_l<-2

/* Boxplot by study
graph box bmi, by(study_main)
graph box whr, by(study_main)

graph box crp_l, by(study_main)
graph box tc, by(study_main)
graph box hdl, by(study_main)
graph box ldl, by(study_main)
graph box tg_l, by(study_main)
*/
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									/////  DESCRIPTIVE STATISTICS /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
// Convert education
gen educ_r=1 if educ==0
replace educ_r=0 if educ==1

cap postclose descriptive
	postfile descriptive str20 (exp_var all all_ps ctrl ctrl_ps dem_case dem_case_ps dem_p_value) ///
	using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210707_late_descriptive.dta, replace

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
sts graph, by(tert_bmi)
stphplot, by(tert_bmi)
// Looks ok
// WHR
xtile tert_whr = whr, nquantiles(3)
sts graph, by(tert_whr)
stphplot, by(tert_whr)
// Looks ok

// Check with ph-test
stcox z_bmi z_whr educ smoke i.sex, vce(cluster pairid)
estat phtest, detail
// whr close at 0.08
stcox z_bmi z_whr z_crp_l z_tc z_hdl z_ldl z_tg_l educ smoke i.sex, vce(cluster pairid)
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
        z_tc |      7,231    .0001444    .9956915  -2.957685   5.026659
       z_hdl |      6,791    .0036212     1.00058  -2.429217   6.161885
       z_ldl |      6,057    .0047144    .9986047  -2.915636   4.429988
      z_tg_l |      6,869   -.0018905    .9992465  -3.261021   3.770288
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
        z_tc |      7,100   -.0066104    .9901801  -2.957685   5.026659
       z_hdl |      6,668    .0058017    1.000417  -2.429217   6.161885
       z_ldl |      6,001    .0059302    .9989539  -2.915636   4.429988
      z_tg_l |      6,738   -.0073913    .9995723  -3.261021   3.770288
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
		using P:\Dementia_IK\FORTE_postdoc\Aim_2\Output\BMI_dem_mediation_20210707_late_Cox_`ad_measure'.dta, replace
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
      tereri |   .0068133   .0036341     1.87   0.061    -.0003093     .013936
   ereri_cde |   .0075915   .0032414     2.34   0.019     .0012384    .0139447
ereri_intref |  -.0005094   .0012223    -0.42   0.677     -.002905    .0018863
ereri_intmed |  -.0000835   .0002013    -0.42   0.678     -.000478    .0003109
   ereri_pie |  -.0001853   .0001985    -0.93   0.350    -.0005743    .0002037
------------------------------------------------------------------------------*/

med4way z_bmi z_hdl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0) fulloutput
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0035831   .0040662     0.88   0.378    -.0043865    .0115526
   ereri_cde |   .0105743   .0034185     3.09   0.002     .0038741    .0172745
ereri_intref |  -.0029212   .0014758    -1.98   0.048    -.0058138   -.0000287
ereri_intmed |  -.0016631   .0008435    -1.97   0.049    -.0033163   -9.96e-06
   ereri_pie |  -.0024069   .0009607    -2.51   0.012    -.0042899    -.000524
      terira |   1.003583   .0040662   246.81   0.000     .9956135    1.011553
       p_cde |     2.9512   2.719178     1.09   0.278     -2.37829     8.28069
    p_intref |  -.8152854   1.217473    -0.67   0.503    -3.201489    1.570918
    p_intmed |  -.4641612   .6939976    -0.67   0.504    -1.824371    .8960491
       p_pie |  -.6717532   .8452624    -0.79   0.427    -2.328437    .9849307
        op_m |  -1.135914   1.515363    -0.75   0.453     -4.10597    1.834142
      op_ati |  -1.279447    1.91107    -0.67   0.503    -5.025076    2.466182
        op_e |    -1.9512   2.719177    -0.72   0.473    -7.280688    3.378289
------------------------------------------------------------------------------*/

med4way z_bmi z_ldl bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0) 
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |   .0073723   .0035187     2.10   0.036     .0004759    .0142687
   ereri_cde |   .0075553   .0034633     2.18   0.029     .0007673    .0143432
ereri_intref |  -.0002236   .0006018    -0.37   0.710    -.0014031    .0009559
ereri_intmed |  -.0000741   .0002008    -0.37   0.712    -.0004678    .0003195
   ereri_pie |   .0001148    .000186     0.62   0.537    -.0002499    .0004794
------------------------------------------------------------------------------*/

med4way z_bmi z_tg_l bsex educ satsa octo gender twingene fasting, a0(0) a1(1) m(0) yreg(aft, weibull) yregoptions(tratio) mreg(linear) c(0 0 0 0 0 0 0)
/*----------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      tereri |    .006093    .003247     1.88   0.061    -.0002711    .0124571
   ereri_cde |   .0089592   .0034839     2.57   0.010     .0021309    .0157875
ereri_intref |   .0000532   .0000998     0.53   0.594    -.0001424    .0002487
ereri_intmed |  -.0011291   .0009839    -1.15   0.251    -.0030575    .0007993
   ereri_pie |  -.0017903   .0009712    -1.84   0.065    -.0036939    .0001132
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
   ereri_cde |   -.001576   .0068635    -0.23   0.818    -.0150281    .0118762
   ereri_cde |   .0044272   .0045063     0.98   0.326    -.0044051    .0132594
   ereri_cde |   .0105743   .0034185     3.09   0.002     .0038741    .0172745
   ereri_cde |   .0168681   .0048192     3.50   0.000     .0074226    .0263137
   ereri_cde |   .0233113   .0075531     3.09   0.002     .0085075    .0381151
*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										///// SAVING THE LOG/OUTPUT /////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

log close

translate session.smcl "P:\Dementia_IK\FORTE_postdoc\Aim_2\Logs\BMI_dem_mediation_20210707_late.log", replace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// END OF FILE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
