////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Name: BMI_DirectEffect_by_HDLlev_20210708
Study: FORTE postdoc aim 2: Identify potential mediator of the associatin between BMI and
dementia
Purpose: Create plot of direct effects of BMI at different levels of HDL-c (Figure 2 of paper)

NOTE! Estimates of the direct effects of BMI comes from mediation analyses in BMI_dem_mediation_20210707_late.do

Created by Ida Karlsson (IK)
Department of Medical Epidemiology and Biostatistics (MEB), Karolinska Institutet, Stockholm

   	Created: 	20210708 by Ida Karlsson (IK)
	Updated: 	
	
STATA v16.0		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

clear
input row estimate min95 max95 
1 -.001576 -.0150281 .0118762
2 .0044272 -.0044051 .0132594
3 .0105743 .0038741 .0172745
4 .0168681 .0074226 .0263137
5 .0233113 .0085075 .0381151
end

set scheme s1mono // black and white

twoway ///
(rcap min95 max95 row, vert) /// 
(scatter estimate row, mcolor(black)) /// 
, legend(off) /// 
xlabel(1 "-2" 2 "-1" 3 "0" 4 "1" 5 "2", angle(0) noticks labsize(med)) ///
xscale(r(0.5 5.5)) ///
ylabel(,labsize(med)) ///
title("Direct effects of BMI, by levels of HDL-c", size(vlarge)) ///
xtitle("Standardized HDL-c level", size(large)) ///
ytitle("Direct effect of BMI", size(large)) ///
yline(0, lc(black) lw(thin) lpattern(line))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

