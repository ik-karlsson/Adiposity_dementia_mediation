/***********************************************************************************************
Name: Forte_aim2_genetic_overlap_20200428
Purpose: Work through FUMA output and generate lists of overlapping genes for pathway analysis

Study: FORTE postdoc aim 2: Study mediating factors in the effect of BMI on dementia and
	   cognitive abilities

Created by Ida Karlsson (IK)
Institutet För Gerontologi, Hälsohögskolan Jönköping
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20200428 by Ida Karlsson (IK)
	Updated: 	

	Steps:	1. Results from gene mapping (position, eQTL, and cromatin interactions)
			2. Results from MAGMA

  	SAS v9.4

************************************************************************************************/
/********************************** 1. Gene mapping results ************************************/
/***********************************************************************************************/

/***** Import the gene mapping results *****/
/* BMI */
PROC IMPORT OUT= WORK.BMI_map
            DATAFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\FUMA\Extracted_results\genes_BMI.txt" 
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
/* BMI */
data BMI_venn;
	set BMI_map;
	if posMapSNPs ge 1 then m1='Position';
	if eqtlMapSNPs ge 1 then m2='eQTL';
	if ciMap='Yes' then m3='Chromatin interaction';

	BMI_GeneMap = catx(', ',m1,m2,m3);
	nval=cmiss(of m1,m2,m3);
run;
*n=8663;

proc sort data=BMI_venn;
	by symbol nval;
data BMI_venn_nodup;
	set BMI_venn;
	by symbol;
	if first.symbol;
run;
*n=8626;
	
proc freq data=BMI_venn_nodup;
	table BMI_GeneMap;
run;

/*                                                                   Cumulative    Cumulative
  BMI_GeneMap                              Frequency     Percent     Frequency      Percent
  ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒ
  Chromatin interaction                        4462       51.73          4462        51.73
  Position                                      111        1.29          4573        53.01
  Position, Chromatin interaction               464        5.38          5037        58.39
  Position, eQTL                                411        4.76          5448        63.16
  Position, eQTL, Chromatin interaction        1658       19.22          7106        82.38
  eQTL                                          553        6.41          7659        88.79
  eQTL, Chromatin interaction                   967       11.21          8626       100.00
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
proc sort data=BMI_venn_nodup;
	by symbol;
proc sort data=AD_venn_nodup;
	by symbol;
run;

/* BMI and AD */
data BMI_AD_map;
	merge BMI_venn_nodup (in=a) AD_venn_nodup (in=b keep=symbol AD_GeneMap);
	by symbol;
	if a and b;
	keep symbol chr start end BMI_GeneMap AD_GeneMap;
run;
*n=175;

/***** Output lists *****/
PROC EXPORT DATA= BMI_AD_map
	OUTFILE= "P:\Dementia_IK\FORTE_postdoc\Aim_2\FUMA\Gene_lists\BMI_AD_mapped_genes.txt" 
    DBMS=TAB REPLACE;
RUN;

/***********************************************************************************************/
/**************************************** END OF FILE ******************************************/
/***********************************************************************************************/
