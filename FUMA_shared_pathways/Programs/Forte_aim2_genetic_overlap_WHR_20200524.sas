/***********************************************************************************************
Name: Forte_aim2_genetic_overlap_20200428
Purpose: Work through FUMA output and generate lists of overlapping genes for pathway analysis
This code is for WHRadjBMI and AD

Study: FORTE postdoc aim 2: Study mediating factors in the effect of BMI on dementia and
	   cognitive abilities

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20200524 by Ida Karlsson (IK)
	Updated: 	

	Steps:	1. Results from gene mapping (position, eQTL, and cromatin interactions)
			2. Results from MAGMA

  	SAS v9.4

************************************************************************************************/
/********************************** 1. Gene mapping results ************************************/
/***********************************************************************************************/

/***** Import the gene mapping results *****/
/* WHR */
PROC IMPORT OUT= WORK.WHR_map
            DATAFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\FUMA\Extracted_results\genes_WHRadjBMI.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;

/* AD */
PROC IMPORT OUT= WORK.AD_map
            DATAFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\FUMA\Extracted_results\genes_AD.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=5000;
RUN;

/***** Look at numbers for venn diagrams for positional, eQTL, and chromatin interaction mapping *****/
/* WHR */
data WHR_venn;
	set WHR_map;
	if posMapSNPs ge 1 then m1='Position';
	if eqtlMapSNPs ge 1 then m2='eQTL';
	if ciMap='Yes' then m3='Chromatin interaction';

	WHR_GeneMap = catx(', ',m1,m2,m3);
	nval=cmiss(of m1,m2,m3);
run;
*n=5589;

proc sort data=WHR_venn;
	by symbol nval;
data WHR_venn_nodup;
	set WHR_venn;
	by symbol;
	if first.symbol;
run;
*n=5568;
	
proc freq data=WHR_venn_nodup;
	table WHR_GeneMap;
run;

/*                                                                   Cumulative    Cumulative
  BMI_GeneMap                              Frequency     Percent     Frequency      Percent
  ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒ
  Chromatin interaction                        2884       51.80          2884        51.80
  Position                                       80        1.44          2964        53.23
  Position, Chromatin interaction               292        5.24          3256        58.48
  Position, eQTL                                252        4.53          3508        63.00
  Position, eQTL, Chromatin interaction        1058       19.00          4566        82.00
  eQTL                                          444        7.97          5010        89.98
  eQTL, Chromatin interaction                   558       10.02          5568       100.00
*/

/* AD */
data AD_venn;
	set AD_map;
	if posMapSNPs ge 1 then m1='Position';
	if eqtlMapSNPs ge 1 then m2='eQTL';
	if ciMap='Yes' then m3='Chromatin interaction';

	AD_GeneMap = catx(', ',m1,m2,m3);
	nval=cmiss(of m1,m2,m3);
run;
*n=295;

proc sort data=AD_venn;
	by symbol nval;
data AD_venn_nodup;
	set AD_venn;
	by symbol;
	if first.symbol;
run;
*n=293;
	
proc freq data=AD_venn_nodup;
	table AD_GeneMap;
run;
/*                                                                   Cumulative    Cumulative
  AD_GeneMap                               Frequency     Percent     Frequency      Percent
  ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒ
  Chromatin interaction                         158       53.92           158        53.92
  Position                                        5        1.71           163        55.63
  Position, Chromatin interaction                15        5.12           178        60.75
  Position, eQTL                                 17        5.80           195        66.55
  Position, eQTL, Chromatin interaction          50       17.06           245        83.62
  eQTL                                           23        7.85           268        91.47
  eQTL, Chromatin interaction                    25        8.53           293       100.00
*/

/***** Overlapping genes *****/
proc sort data=WHR_venn_nodup;
	by symbol;
proc sort data=AD_venn_nodup;
	by symbol;
run;

/* WHR and AD */
data WHR_AD_map;
	merge WHR_venn_nodup (in=a) AD_venn_nodup (in=b keep=symbol AD_GeneMap);
	by symbol;
	if a and b;
	keep symbol chr start end WHR_GeneMap AD_GeneMap;
run;
*n=65;

/***** Output lists *****/
PROC EXPORT DATA= WHR_AD_map
	OUTFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\FUMA\Gene_lists\WHR_AD_mapped_genes.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
/**************************************** END OF FILE ******************************************/
/***********************************************************************************************/
