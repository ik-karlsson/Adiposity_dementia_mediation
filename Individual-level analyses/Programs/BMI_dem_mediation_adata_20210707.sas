/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************
Name: BMI_dem_mediation_20210707
Purpose: Combine data for analysis of the association between BMI and dementia, and the 
mediating effect of CRP and cholesterol fractions

Study: FORTE postdoc aim 2: Identify potential mediator of the associatin between BMI and
dementia

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20201126 by Ida Karlsson (IK)
	Updated: 	20210218 by IK - Expand the study population to GHOST, add PRS for mediators
				20210506 by IK - Drop HARMONY (such a selected sample..), drop the PRSs
				20210707 by IK - Correct error in fasting variable

	Steps:		1. Libnames 
				2. Define the study population
				3. Dementia data
				4. BMI data (long format)
				5. PRS data
				6. Biomarker data
				7. Covariates and additional variables
				8. Save the data

  	SAS v9.4

************************************************************************************************/
/**************************************** 1. Libnames ******************************************/
/***********************************************************************************************/

libname aim2 'P:\Dementia_IK\FORTE_postdoc\Aim_2\Data\Derived_data';

libname tgfast 'P:\twins\twins_Research\internal_base\TwinGene';

libname bmi 'P:\Dementia_IK\BMI\BMI_data\Data\Derived_data';

libname s2 'P:\Dementia_IK\Archived\PhD project\Study II - CHD_genes_dem\Data\Derived_data';

libname aim1 'P:\Dementia_IK\FORTE_postdoc\Aim_1\Data\Derived_data\Aim_1a\GHOST';

/***********************************************************************************************/
/******************************** 2. Define the study population *******************************/
/***********************************************************************************************/

data studypop;
	merge 	satsa.resp0517 (in=a)
			octo.OCTO_INCIDENT_080930 (in=b keep=twinnr)
			gender.INCIDENT_COMPLETE_080812 (in=c keep=twinnr)
			twingene.v_age (in=d);
	by twinnr;

	if a and (iage1 ne . or iage2 ne . or iage3 ne . or iage5 ne . or iage6 ne .) then study_main='SATSA   ';
	else if b then study_main='OCTOtw';
	else if c then study_main='GENDER';
	else if d then study_main='TwinGene';

	if study_main=' ' then delete;
	keep twinnr study_main;
run;
*n=14580;	

proc freq data=studypop;
	table study_main/missprint;
run;

/***********************************************************************************************/
/************************************** 3. Dementia data ***************************************/
/***********************************************************************************************/
data master;
	set "P:\Dementia_IK\Various data across and outside studies\Master_dementia_file\Data\DerivedData\master_dementia_20191018.sas7bdat";
run;
*n=78315;

/* NOTE! Remove twins missing lastage or <60 at last age */
data demdata;
	merge 	studypop (in=a)
			master (in=b keep=twinnr dementia AD dem_onset dem_info lastage mixed flag_early_onset flag_2017);
	by twinnr;
	if a and b and dementia in(0,1) and lastage ge 60; 
	/* Recode mixed as AD */
	if mixed=1 then AD=1;
	/* Recode early onset cases as missing */
	if flag_early_onset=1 then do;
		dementia=.;
		ad=.;
		dem_onset=.;
	end;
	/* Recode cases identified during 2017 as controls (will be time-to-event analysis, so end of follow-up 2016) */
	if flag_2017=1 then do; 
		dementia=0;
		AD=0;
		dem_onset=.;
	end;
	if dementia=1 then dem_onset=int(dem_onset);
	lastage=int(lastage);
	drop mixed flag_early_onset flag_2017;
	if dementia=. then delete;
run;
*n=14033;
*n=12 not in dementia file, n=144 dementia=9, n=6 flag early onset, n=385 lastage <60; 

proc freq data=demdata;
	table dementia AD/missprint;
run;
*n=1360 dementia, 800 AD, 12673 controls;

/***********************************************************************************************/
/**************************** 4. Covariates and additional variables ***************************/
/***********************************************************************************************/

/* Birthyear, age at death, sex, and zygosity from STR */
data demdata_cov1;
	merge 	demdata (in=a)
			stradmin.person_info (in=b keep=twinnr pairid tvab sex birthdate deathdate vitalstatus)
			stradmin.zygosity (in=c keep=twinnr zygosity);
	by twinnr;
	if a and b and c;

	/***** Calculate age at death *****/
	/* Convert birthdate to date format */
	if birthdate='1000000' then delete;
	bdate=input(put(birthdate,8.),yymmdd8.);
	/* 4 twins are dead, but without a date. Assigning date of death from CDR. */
	if twinnr=xxx then deathdate='19920400';
	if twinnr=xxx then deathdate='19930700';
	if twinnr=xxx then deathdate='20060825';
	if twinnr=xxx then deathdate='20020000';
	/* Deal with faulty deathdates */
	if deathdate ne ' ' then do;
		deathyr=substr(deathdate, 1, 4);
		deathmon=substr(deathdate, 5, 2);
		deathday=substr(deathdate, 7,2);

		if deathmon in(' ', '00') then deathmon='12';
		if deathday in (' ', '00') then do;
			if deathmon in ('01', '03', '05', '07', '08', '10', '12') then deathday='31';
			else if deathmon in ('04', '06', '09', '11') then deathday='30';
			else if deathmon='02' then deathday='28';
		end;
		/* Create deathdate in date format */
		ddate=MDY(deathmon, deathday, deathyr);
		/* Calculate age at death */
		age_death=int(YRDIF(bdate, ddate, 'actual'));
	end;
	drop bdate ddate deathyr deathmon deathday;

	if vitalstatus=1 then do;
		if age_death gt lastage then vitalstatus=0;
	end;

run;
*n=14021;
*NOTE! Twinnrs are excluded for ethical reasons;

/* Education */
data demdata_cov2;
	merge 	demdata_cov1 (in=a)
			s2.mr_education_20140626 (in=b keep=twinnr educ)
			;
	by twinnr;
	if a;
	if educ=. then delete;
run;
*n=14002 (31 removed);

/* Smoking */
data demdata_cov3;
	merge	demdata_cov2 (in=a) 
			aim1.smoking_GOSH_20190603 (in=b);
	by twinnr;
	if a;
	*if STR_smoke=. then delete;
	rename STR_smoke=smoke;
	drop smoke_Q61--smoke_gender;
run;
proc freq data=demdata_cov3;
	table smoke/missprint norow nocol nopercent;
run;
*n=24 missing;

/* Indicator for complete twin pairs */
data tvab1 tvab2;
	set demdata_cov3;
	if tvab=1 then output tvab1;
	else if tvab=2 then output tvab2;
run;

data pair;
	merge	tvab1 (in=a keep=pairid zygosity)
			tvab2 (in=b keep=pairid);
	by pairid;
	if a and b;
run;
*n=5700;
proc freq data=pair;
	tables zygosity/missprint norow nocol nopercent;
run;

data demdata_cov4;
	merge 	demdata_cov3 (in=a)
			pair (in=b drop=zygosity);
	by pairid;
	if b then complete_pair=1;
	else complete_pair=0;
run;

/***********************************************************************************************/
/********************************** 5. BMI and biomarker data **********************************/
/***********************************************************************************************/
/* Select out BMI data for those with no biomarker data */
/* BMI measures closest to 50 (between 35 and 65) and closest to 70 (between 65 and 90)*/
data bmi_only;
	merge 	bmi.Bmi_whr_data_long_20201105 (in=a)
			aim2.Biomarkers_gost_20210218 (in=b);
	by twinnr;
	if a and not b;
	keep twinnr bmi age;
run;
	
/* Split in two files: mid- and late-life.
   Select measurements with maximum nr of markers (plus BMI) */
data biom1;
	set aim2.Biomarkers_gost_20210218 bmi_only;
	if BMI ne .;
	/* Combine GOS and TG info on fasting */
	if fasting=. then fasting=fast;
	/* Number of biomarker values (NOTE! Crp gets one and lipids one) */
	crp_m=N(of crp);
	lip_m=(N(of tc,hdl,ldl,tg)/4);
	cnt_m=crp_m+lip_m;
	/* age from 50 and 70 */
	n50=abs(age-50);
	n70=abs(age-70);
	keep twinnr age study bmi whr tc hdl ldl tg crp cholmed fasting cnt_m n50 n70; 
	rename study=study_biom;
run;

data biom_e biom_l;
	set biom1;
	if age ge 40 and age lt 65 then output biom_e;
	else if age ge 65 and age lt 90 then output biom_l;
run;

/* Midlife data */
/* Closest to 50, maximizing biomarker measures */
proc sort data=biom_e;
	by twinnr descending cnt_m n50;
run;

data biom_e1;
	set biom_e;
	by twinnr;
	if first.twinnr;
	drop cnt_m n50 n70;
run;

/* Add to the study population */
data demdata_mid;
	merge 	demdata_cov4 (in=a)
			biom_e1 (in=b);
	by twinnr;
	if a and b;
run;
*n=6030;

/* Late-life */
/* Closest to 70, maximizing biomarker measures */
proc sort data=biom_l;
	by twinnr descending cnt_m n70;
run;

data biom_l1;
	set biom_l;
	by twinnr;
	if first.twinnr;
	drop cnt_m n50 n70;
run;

/* Add to the study population */
data demdata_late;
	merge 	demdata_cov4 (in=a)
			biom_l1 (in=b);
	by twinnr;
	if a and b;
run;
*n=7682;

/* Check N for MS */
data checkmid;
	set	demdata_mid;	
	if (bmi ne . or whr ne .) and (crp ne . or tc ne . or hdl ne . or ldl ne . or tg ne .);
run;
*n=6012;
data checkmid2;
	set	checkmid;	
	if (age lt lastage) and (dem_onset=. or age lt dem_onset);
run;
*n=5999 - good;

data checklate;
	set	demdata_late;	
	if (bmi ne . or whr ne .) and (crp ne . or tc ne . or hdl ne . or ldl ne . or tg ne .);
run;
*n=7483;
data checklate2;
	set	checklate;	
	if (age lt lastage) and (dem_onset=. or age lt dem_onset);
run;
*n=7257 - good;

data checkmidlate;
	merge	checkmid (in=a) 
			checklate (in=b);
	by twinnr;
run;
*13113;
data checkmidlate2;
	merge	checkmid2 (in=a) 
			checklate2 (in=b);
	by twinnr;
	if a and b;
run;
*n=12893;
*n=363 in both;


/***********************************************************************************************/
/************************************* 6. Save the data ****************************************/
/***********************************************************************************************/

/* Midlife */
data aim2.BMI_dem_mediation_20210707_mid;
	retain	twinnr
			pairid
			tvab
			sex 
			zygosity
			complete_pair
			birthdate
			deathdate
			age_death
			vitalstatus
			lastage
			educ
			smoke
			dementia
			AD
			dem_onset
			dem_info
			age
			BMI
			WHR
			CRP
			TC
			HDL
			LDL
			TG
			cholmed
			fasting
			study_biom
			;
	set demdata_mid;
run;

/* STATA format */
PROC EXPORT DATA= aim2.BMI_dem_mediation_20210707_mid
	OUTFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\Data\Derived_data\BMI_dem_mediation_20210707_mid.dta" 
    DBMS=STATA REPLACE;
RUN;


/* Late life */
data aim2.BMI_dem_mediation_20210707_late;
	retain	twinnr
			pairid
			tvab
			sex 
			zygosity
			complete_pair
			birthdate
			deathdate
			age_death
			vitalstatus
			lastage
			educ
			smoke
			dementia
			AD
			dem_onset
			dem_info
			age
			BMI
			WHR
			CRP
			TC
			HDL
			LDL
			TG
			cholmed
			fasting
			study_biom
			;
	set demdata_late;
run;

/* STATA format */
PROC EXPORT DATA= aim2.BMI_dem_mediation_20210707_late
	OUTFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\Data\Derived_data\BMI_dem_mediation_20210707_late.dta" 
    DBMS=STATA REPLACE;
RUN;

/***********************************************************************************************/
/**************************************** END OF FILE ******************************************/
/***********************************************************************************************/
