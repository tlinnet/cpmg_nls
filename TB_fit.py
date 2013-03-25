from pylab import *
import scipy.optimize
import os


###Cut-off
co = 0.005
### figure thing
figsize = 16
titfont = 26
labfont = 12
figfont = 12

#Data FT
hewl_FT = []
hewl_FT.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'r','hewl-75/33','FT',10,50,'y'])
hewl_FT.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'g','hewl-75/33','FT',30,500,'y'])
hewl_FT.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'b','hewl-75/33','FT',100,500,'y'])
#hewl_FT.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'c','hewl-75/33','FT',10,50,'n'])
#hewl_FT.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'y','hewl-75/33','FT',30,500,'n'])

acbp_FT = []
acbp_FT.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'r','acbp-17/82','FT',10,50,'y'])
acbp_FT.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'g','acbp-17/82','FT',30,500,'y'])
acbp_FT.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'b','acbp-17/82','FT',100,500,'y'])

bbl_FT = []
bbl_FT.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'r','bbl-19/112','FT',10,50,'y'])
bbl_FT.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS30_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'g','bbl-19/112','FT',30,500,'y'])
bbl_FT.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/FT/allplanes_FT.stats'),'b','bbl-19/112','FT',100,500,'y'])

#Data MDD
hewl_MDD = []
hewl_MDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'r','hewl-75/33','MDD',10,50,'y'])
hewl_MDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'g','hewl-75/33','MDD',30,500,'y'])
hewl_MDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'b','hewl-75/33','MDD',100,500,'y'])
#hewl_MDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'c','hewl-75/33','MDD',10,50,'n'])
#hewl_MDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'y','hewl-75/33','MDD',30,500,'n'])

acbp_MDD = []
acbp_MDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'r','acbp-17/82','MDD',10,50,'y'])
acbp_MDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'g','acbp-17/82','MDD',30,500,'y'])
acbp_MDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'b','acbp-17/82','MDD',100,500,'y'])

bbl_MDD = []
bbl_MDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'r','bbl-19/112','MDD',10,50,'y'])
bbl_MDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS30_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'g','bbl-19/112','MDD',30,500,'y'])
bbl_MDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/MDD/allplanes_MDD.stats'),'b','bbl-19/112','MDD',100,500,'y'])

#Data coMDD
hewl_coMDD = []
hewl_coMDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'r','hewl-75/33','coMDD',10,50,'y'])
hewl_coMDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'g','hewl-75/33','coMDD',30,500,'y'])
hewl_coMDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'b','hewl-75/33','coMDD',100,500,'y'])
#hewl_coMDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'c','hewl-75/33','coMDD',10,50,'n'])
#hewl_coMDD.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'y','hewl-75/33','coMDD',30,500,'n'])

acbp_coMDD = []
acbp_coMDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'r','acbp-17/82','coMDD',10,50,'y'])
acbp_coMDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'g','acbp-17/82','coMDD',30,500,'y'])
acbp_coMDD.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'b','acbp-17/82','coMDD',100,500,'y'])

bbl_coMDD = []
bbl_coMDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'r','bbl-19/112','coMDD',10,50,'y'])
bbl_coMDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS30_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'g','bbl-19/112','coMDD',30,500,'y'])
bbl_coMDD.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/coMDD/allplanes_coMDD.stats'),'b','bbl-19/112','coMDD',100,500,'y'])

#Data CS
hewl_CS = []
hewl_CS.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'r','hewl-75/33','CS',10,50,'y'])
hewl_CS.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'g','hewl-75/33','CS',30,500,'y'])
hewl_CS.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'b','hewl-75/33','CS',100,500,'y'])
#hewl_CS.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'c','hewl-75/33','CS',10,50,'n'])
#hewl_CS.append([genfromtxt('/home/tlinnet/kte/080716_cpmgDisp_HEWLpH65_normal_CS30_MDD500_NoRandom.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'y','hewl-75/33','CS',30,500,'n'])

acbp_CS = []
acbp_CS.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'r','acbp-17/82','CS',10,50,'y'])
acbp_CS.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS30_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'g','acbp-17/82','CS',30,500,'y'])
acbp_CS.append([genfromtxt('/home/tlinnet/kte/acbp/acbp_cpmg_disp_04MGuHCl_40C_041223_CS100_MDD500.fid/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'b','acbp-17/82','CS',100,500,'y'])

bbl_CS = []
bbl_CS.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'r','bbl-19/112','CS',10,50,'y'])
bbl_CS.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS30_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'g','bbl-19/112','CS',30,500,'y'])
bbl_CS.append([genfromtxt('/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_CS100_MDD500/analysis_FT/int_corr_ft_method_all_awk_full/CS/allplanes_CS.stats'),'b','bbl-19/112','CS',100,500,'y'])


### Fit function
#http://www.swvgs.us/rfisher/mathematics/exponentials_and_logarithms/section_1/exponential_translations_1b.htm
#http://www.scipy.org/Cookbook/FittingData
fitfunc = lambda p, x: p[0]**(x-p[1])+p[2]*x # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

# Get paramaters from data
def getparam(D):
    P = []
    for i in range(len(D)):
        D[i][0][:,1]=1-D[i][0][:,0]/max(D[i][0][:,0])
        datX = D[i][0][:,1]
        datY = D[i][0][:,38]
        pb = [1e8, 1.19, 1.] 
        pa, suc = scipy.optimize.leastsq(errfunc, pb[:], args=(datX, datY))
        fitY = fitfunc(pa,datX)

        fitYi = next(i for i,v in enumerate(fitY) if v > co)
        fitYiy = fitY[fitYi]
        fitYix = datX[fitYi]

        datYi = next(i for i,v in enumerate(datY) if v > co)
        datYiy = datY[datYi]
        datYix = datX[datYi]

        print('%3.3f %s %s %s %s %s'%(datYix,D[i][2],D[i][3],D[i][4],D[i][5],D[i][6]))   
        #P.append([pa,suc,fitYix,fitYiy])
        P.append([pa,suc,datYix,datYiy])
    return P

def plotdata(D,P):
    figA = figure(figsize=(figsize, figsize/1.618))
    figA.suptitle('Tada',fontsize=titfont)
    #figA.gca().invert_xaxis()
    #figA.gca().axes.get_xaxis().set_visible(False)
    #figA.gca().axes.get_yaxis().set_visible(False)
    #Fig 1
    figA1 = figA.add_subplot(111)
    for i in range(len(D)):
        figA1.plot(D[i][0][:,1], D[i][0][:,38],'.-',color=D[i][1],linewidth=0.5,label='%s %s %s %s %s'%(D[i][2],D[i][3],D[i][4],D[i][5],D[i][6]))
        #figA1.plot(D[i][0][:,1], fitfunc(P[i][0],D[i][0][:,1]),'-',color=D[i][1],linewidth=0.5,label='fit')
        figA1.annotate('%3.2f'%(P[i][2]), xy=(P[i][2],P[i][3]),xytext=(P[i][2],0.1*P[i][3]+i*0.1*P[i][3]),arrowprops=dict(facecolor=D[i][1], shrink=0.05, width=2))
    #figA1.set_title("Overlay of all experiments",fontsize=labfont)
    figA1.set_ylim([0,0.01])
    figA1.legend(loc="upper left")
    figA1.set_ylabel('Residual stdev',fontsize=labfont)
    figA1.set_xlabel("Sparseness",fontsize=labfont)

#Param/plot FT
hewl_FTp = getparam(hewl_FT)
acbp_FTp = getparam(acbp_FT)
bbl_FTp = getparam(bbl_FT)

#plotdata(hewl_FT,hewl_FTp)
#plotdata(acbp_FT,acbp_FTp)
#plotdata(bbl_FT,bbl_FTp)

#Param/plot MDD
hewl_MDDp = getparam(hewl_MDD)
acbp_MDDp = getparam(acbp_MDD)
bbl_MDDp = getparam(bbl_MDD)

#plotdata(hewl_MDD,hewl_MDDp)
#plotdata(acbp_MDD,acbp_MDDp)
#plotdata(bbl_MDD,bbl_MDDp)

#Param/plot coMDD
hewl_coMDDp = getparam(hewl_coMDD)
acbp_coMDDp = getparam(acbp_coMDD)
bbl_coMDDp = getparam(bbl_coMDD)

plotdata(hewl_coMDD,hewl_coMDDp)
plotdata(acbp_coMDD,acbp_coMDDp)
plotdata(bbl_coMDD,bbl_coMDDp)

#Param/plot CS
hewl_CSp = getparam(hewl_CS)
acbp_CSp = getparam(acbp_CS)
bbl_CSp = getparam(bbl_CS)

plotdata(hewl_CS,hewl_CSp)
plotdata(acbp_CS,acbp_CSp)
plotdata(bbl_CS,bbl_CSp)

show()
