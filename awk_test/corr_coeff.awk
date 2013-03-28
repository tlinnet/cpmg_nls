#http://people.sc.fsu.edu/~sshanbhag/awk-shell/linreg.awk
#http://www.cyberciti.biz/faq/bash-scripting-using-awk/
#http://en.wikipedia.org/wiki/Simple_linear_regression
#awk -v max=1 '{print $1, $2*max}' testdataP225.txt | awk -f corr_coeff.awk
BEGIN{}{ x[NR] = $1; y[NR] = $2;  sx += x[NR]; sy += y[NR]; sxx += x[NR]*x[NR]; sxy += x[NR]*y[NR]; syy += y[NR]*y[NR]; sresi += (x[NR]-y[NR])}
END{ det = NR*sxx - sx*sx; a = (NR*sxy - sx*sy)/det; b = (-sx*sxy+sxx*sy)/det; r = (NR*sxy-sx*sy)/sqrt((NR*sxx-sx*sx)*(NR*syy-sy*sy));
beta=(NR*sxy-sx*sy)/(NR*sxx-sx*sx); alpha=sy/NR-beta/NR*sx; se2=1/(NR*(NR-2))*(NR*syy-sy*sy-beta*beta*(NR*sxx-sx*sx)); sb2=(NR*se2)/(NR*sxx-sx*sx); sa2=sb2/NR*sxx;
for(xi=1;xi<=NR;xi++){xsumsq+=((x[xi]-sx/NR)**2)} ; xvar=xsumsq/(NR); xstdev=sqrt(xvar);
for(yi=1;yi<=NR;yi++){ysumsq+=((y[yi]-sy/NR)**2)} ; yvar=ysumsq/(NR); ystdev=sqrt(yvar);
for(yi=1;yi<=NR;yi++){resisumsq+=((x[yi]-y[yi]-sresi/NR)**2)} ; resivar=resisumsq/(NR); resistdev=sqrt(resivar);
for(xi=1;xi<=NR;xi++){xysumsq+=((x[xi]-sx/NR)*(y[xi]-sy/NR))} ; xycovar=xysumsq/(NR);
print "NR= "NR, "xvar= "xvar, "xstdev= "xstdev, "yvar= "yvar, "ystdev= "ystdev, "xycovar= "xycovar, "r= "r, "r2= "r*r, "a= " a, "b= "b, "resivar= "resivar, "resistdev= "resistdev
#print "sx= "sx, "sy= "sy, "sxx= "sxx, "sxy= "sxy, "syy= "syy, "beta= "beta, "alpha= "alpha;
#print "se2= "se2, "se= "sqrt(se2), "sb2= "sb2, "sb= "sqrt(sb2), "sa2= "sa2, "sa= "sqrt(sa2)
}
