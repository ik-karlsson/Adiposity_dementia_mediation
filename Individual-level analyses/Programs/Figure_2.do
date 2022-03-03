
clear
input row estimate min95 max95 
1 -.0020323 -.0154747 .01141
2 .0041198 -.0047057 .0129454
3 .0104221 .0037268 .0171174
4 .0168773 .0074295 .026325
5 .0234882 .0086719 .0383044
   
end

set scheme s1mono // black and white

twoway ///
(rcap min95 max95 row, vert) /// code for 95% CI
(scatter estimate row, mcolor(black)) /// dot for group 
, legend(off) /// legend at 6 o’clock position
xlabel(1 "-2" 2 "-1" 3 "0" 4 "1" 5 "2", angle(0) noticks labsize(med)) ///
xscale(r(0.5 5.5)) ///
ylabel(,labsize(med)) ///
title("Direct effects of BMI, by levels of HDL-c", size(vlarge)) ///
xtitle("Standardized HDL-c level", size(large)) ///
ytitle("Direct effect of BMI", size(large)) ///
yline(0, lc(black) lw(thin) lpattern(line))

 ///
/// aspect (next line) is how tall or wide the figure is
///aspect(.5)
