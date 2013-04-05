from pylab import *
import scipy.optimize
import os
import TB
reload(TB)
dn = os.path.dirname(os.path.realpath(__file__))

#Data BBL
BBL = {}
BBL['NMRpar'] = {'yCAR':117.843,'centerPPM':117.843,'yOBS':76.012,'frq':76.012}
BBL['offset'] = (0,2000,500,5000,1000,20000,0,0,0,1000,2000,1000,2000,1000,2000) # offset/ddof2
BBL['omega1'] = (1903.4,1903.4,1903.4,1903.4,1903.4,1903.4,1379.5,898.1,1113.1,1113.1,1113.1,1379.5,898.1,898.1,1379.5) #omega1/spinlock/slockpwr
BBL['omega1_col'] = {'898.1':'r','1113.1':'g','1379.5':'b','1903.4':'y'} #omega1/spinlock/slockpwr
BBL['time'] = array([0, 0.1, 0.4, 0.04, 0.2])
BBL['guess'] = {'s_R1':1.0,'s_R2':40.0,'s_kEX':10000.0,'s_phi':100000.0,'g_kEX':10000.0}
BBL['desc'] = ('bbl-75/33','FT',10,50,'y')
BBL['NIstop'] = 60
BBL['path'] = os.path.join(dn,'data','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
#BBL['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
BBL['qMDDmet'] = ['CS','coMDD'] # ['FT']
BBL['ser'] = {'pre':'allplanes_','filee':'.ser'}
BBL['stats'] = {'pre':'allplanes_','filee':'.stats'}

TB.getstat(BBL,BBL['qMDDmet'])
TB.getser(BBL,BBL['qMDDmet'])
TB.sortInt(BBL,BBL['qMDDmet'])
TB.getdecay(BBL,BBL['qMDDmet'])
#TB.plotdecays([BBL],BBL['qMDDmet'],fss=range(0,10,5))
p = TB.getrates(BBL,BBL['qMDDmet'])
#TB.plotrates([BBL],BBL['qMDDmet'])
psel = [3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 33]
TB.getglobfit(BBL,BBL['qMDDmet'],psel)
TB.plot_kEX([BBL],BBL['qMDDmet'])
#########################################
#BBL2 = {}
#BBL2['desc'] = ('bbl-75/33','FT',30,500,'y')
#BBL2['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_CS30_MDD500','analysis_FT','int_corr_ft_method_all_awk_full')
#BBL2['qMDDmet'] = ('CS','coMDD')
#BBL2['ser'] = {'pre':'allplanes_','filee':'.ser'}
#BBL2['stats'] = {'pre':'allplanes_','filee':'.stats'}
#TB.getstat(BBL2,BBL2['qMDDmet'])
#TB.getser(BBL2,BBL2['qMDDmet'])
#TB.sortInt(BBL2,BBL2['qMDDmet'])
#TB.plotstats([BBL,BBL2],BBL['qMDDmet'])



