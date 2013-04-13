set term postscript eps enhanced color "Helvetica" 14
set title "Statistics plot for method FT\nIntensity normalized to maximum of 128 allplanes FT, and a=1"
set size ratio 0.618
set xlabel "Number of random ni out of ni_{max}"
set ylabel "Correlation Coefficient R^2\nIntensity proportionality a, f(x)=a*x+b"
set ytics nomirror
set y2tics border
set y2label "Residual stdev {/Symbol s}"
#set y2label "Covariance {/Symbol s}_{xy}"
set nokey
#set key autotitle columnheader
set key bottom center outside font "Helvetica,7"
set xrange [-128:]
set yrange[-1:1.1]
set y2range[0:0.1]
#set xtics border in scale 1,0.5 mirror rotate by -90 font "Helvetica,5"
#set ytics border in scale 1,0.5 mirror font "Helvetica,5"
set output "allplanes_FT.eps"

#set fit errorvariables
#f(x)=a*x + b
#fit f(x) "allplanes_FT.ser" using ($6*$8 /5376722):($6*$73 /5376722) via a, b
#f_FIT_NDF=FIT_NDF
#f_FIT_STDFIT=FIT_STDFIT
#f_FIT_WSSR=FIT_WSSR
#f_WSSR_NDF=FIT_WSSR/FIT_NDF

#set label sprintf("f(x)=a*x") at graph 0.8,0.80 font "Helvetica,9"
#set label sprintf("a=%1.2f +- %1.4f",a,a_err) at graph 0.8,0.70 font "Helvetica,9"
#set label sprintf("b=%1.2f +- %1.4f",b,b_err) at graph 0.8,0.60 font "Helvetica,9"
#set label sprintf("NDF=%1.2f",f_FIT_NDF) at graph 0.8,0.50 font "Helvetica,9"
#set label sprintf("STDFIT=%1.2f",f_FIT_STDFIT) at graph 0.8,0.40 font "Helvetica,9"
#set label sprintf("WSSR=%1.2f",f_FIT_WSSR) at graph 0.8,0.30 font "Helvetica,9"
#set label sprintf("WSSR/NDF=%1.2f",f_WSSR_NDF) at graph 0.8,0.20 font "Helvetica,9"

plot "allplanes_FT.stats" using (-1*$1):($71):($21):($25) title "Residual average" with yerrorbars,\
"allplanes_FT.stats" using (-1*$1):($47) title "Intensity proportionality a, f(x)=a*x+b",\
"allplanes_FT.stats" using (-1*$1):($45) title "Correlation Coefficient R^2",\
"allplanes_FT.stats" using (-1*$1):($39) title "Y2: Residual stdev {/Symbol s}" axis x1y2,\
"allplanes_FT.stats" using (-1*$1):($41) title "Y2: Covariance {/Symbol s}_{xy}" axis x1y2
