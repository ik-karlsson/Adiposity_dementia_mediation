
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
