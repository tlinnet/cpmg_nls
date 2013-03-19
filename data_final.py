from pylab import *
import scipy.optimize
import os
#
whichtask = 1

#
co = 0.005
#http://www.swvgs.us/rfisher/mathematics/exponentials_and_logarithms/section_1/exponential_translations_1b.htm
#http://www.scipy.org/Cookbook/FittingData
fitfunc = lambda p, x: p[0]**(x-p[1])+p[2]*x # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

#Combine data
D = []
### http://www.scipy.org/NumPy_for_Matlab_Users
D.append([genfromtxt('../kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'b','bbl-19/112','co-mdd',100,500])
D.append([genfromtxt('../kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'r','acbp-17/82','co-mdd',100,500])
D.append([genfromtxt('../kte/080716_cpmgDisp_HEWLpH65_normal_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'g','hewl-75/33','co-mdd',100,500])

P = []
for i in range(len(D)):
    D[i][0][:,1]=1.0-D[i][0][:,0]/max(D[i][0][:,0])
    pb = [1e8, 1.19, 1.] 
    pa, suc = scipy.optimize.leastsq(errfunc, pb[:], args=(D[i][0][:,1], D[i][0][:,38]))
    Y = fitfunc(pa,D[i][0][:,1])
    Yi = next(i for i,v in enumerate(Y) if v > co)
    Yv = Y[Yi]
    Yx = D[i][0][:,1][Yi]
    print Yi, Yv
    P.append([pa,suc,Yx,Yv])

###
figsize = 16
titfont = 26
labfont = 12
figfont = 12
###
if whichtask == 1:
    figA = figure(figsize=(figsize, figsize/1.618))
    figA.suptitle('Tada',fontsize=titfont)
    #figA.gca().invert_xaxis()
    #figA.gca().axes.get_xaxis().set_visible(False)
    #figA.gca().axes.get_yaxis().set_visible(False)
    #Fig 1
    figA1 = figA.add_subplot(111)
    for i in range(len(D)):
        figA1.plot(D[i][0][:,1], D[i][0][:,38],'.',color=D[i][1],linewidth=0.5,label='%s %s %s %s'%(D[i][2],D[i][3],D[i][4],D[i][5]))
        figA1.plot(D[i][0][:,1], fitfunc(P[i][0],D[i][0][:,1]),'-',color=D[i][1],linewidth=0.5,label='fit')
        figA1.annotate('Cut-off:%s,\nx=%3.2f'%(co,P[i][2]), xy=(P[i][2],P[i][3]),xytext=(P[i][2]-i*0.1,20*P[i][3]),arrowprops=dict(facecolor=D[i][1], shrink=0.05, width=2))
    #figA1.set_title("Overlay of all experiments",fontsize=labfont)
    #figA1.set_xlim([5,2.1])
    figA1.legend(loc="best")
    figA1.set_ylabel('Residual stdev',fontsize=labfont)
    figA1.set_xlabel("Sparseness",fontsize=labfont)

    show()
