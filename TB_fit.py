from pylab import *
import scipy.optimize
import os, sys
import multiprocessing
import cPickle as pickle
#import pickle as pickle
import yaml
import TB
reload(TB)

dn = os.path.dirname(os.path.realpath(__file__)) #Current directory for script position. See also: os.getcwd() for position of execution.
outdir = os.path.join(os.getcwd(),'pickle')
if not os.path.exists(outdir): os.mkdir(outdir)
data_format = 'pickle' # 'pickle' 'yaml'
load_data = False #True False
save_data = True

multiprocess = True
if multiprocess: jobs = multiprocessing.cpu_count()-1; print "Using Cores=%s for computation"%jobs
else: jobs = 1

#General settings for dataset
ALL = {}
ALL['NMRpar'] = {'yCAR':117.843,'centerPPM':117.843,'yOBS':76.012,'frq':76.012}
ALL['offset'] = (0,2000,500,5000,1000,20000,0,0,0,1000,2000,1000,2000,1000,2000) # offset/ddof2
ALL['omega1'] = (1903.4,1903.4,1903.4,1903.4,1903.4,1903.4,1379.5,898.1,1113.1,1113.1,1113.1,1379.5,898.1,898.1,1379.5) #omega1/spinlock/slockpwr
ALL['omega1_col'] = {'898.1':'r','1113.1':'g','1379.5':'b','1903.4':'y'} #omega1/spinlock/slockpwr
ALL['time'] = array([0, 0.1, 0.4, 0.04, 0.2])
ALL['guess'] = {'s_R1':1.0,'s_R2':40.0,'s_kEX':10000.0,'s_phi':100000.0,'g_kEX':10000.0}
ALL['NIstop'] = 0
ALL['Flags'] = {'multiprocess':multiprocess,'jobs':jobs}
ALL['ser'] = {'pre':'allplanes_','filee':'.ser'}
ALL['stats'] = {'pre':'allplanes_','filee':'.stats'}
ALL['qMDDmet'] = ['coMDD','CS','MDD','FT'] #  ['coMDD','CS','MDD','FT']]

############## Normal #################
if not load_data:
    BBL = {}
    BBL.update(ALL)
    BBL['desc'] = {'name':'bbl-f75-p33','comp_to':'FT','CS_niter':10,'NITER':50,'random_seed':'y'}
    BBL['path'] = os.path.join(dn,'data','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
    #BBL['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
    if not multiprocess: sys.stdout = TB.Logger(os.path.join(outdir,"%s.log"%(BBL['desc']['name']))) #Log to file

    TB.getstat(BBL,BBL['qMDDmet'])
    TB.getser(BBL,BBL['qMDDmet'])
    TB.sortInt(BBL,BBL['qMDDmet'])
    #l = range(94,98,2); BBL['NIarr']['CS'] = l[::-1]; BBL['NIarr']['MDD'] = l[::-1]; BBL['NIarr']['FT'] = l[::-1]
    TB.getdecay(BBL,BBL['qMDDmet'])
    #TB.plotdecays([BBL],BBL['qMDDmet'],peaks=[3],fss=range(0,5,5))
    TB.getrates(BBL,BBL['qMDDmet'])
    #TB.plotrates([BBL],BBL['qMDDmet'], peaks=[1,2,3,4]) 
    ##psel = [3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 33]
    TB.getglobfit(BBL,BBL['qMDDmet'])

    if data_format == 'pickle' and save_data: pickle.dump(BBL,open(os.path.join(outdir,"%s.pickle"%(BBL['desc']['name'])),"wb"))
    elif data_format == 'yaml' and save_data: yaml.dump(BBL,open(os.path.join(outdir,"%s.yaml"%(BBL['desc']['name'])),"w"))

elif load_data:
    if data_format == 'pickle': BBL = pickle.load(open(os.path.join(outdir,"%s.pickle"%(BBL['desc']['name'])),"rb"))
    elif data_format == 'yaml': BBL = yaml.load(open(os.path.join(outdir,"%s.pickle"%(BBL['desc']['name'])),"rb"))

#TB.plot_kEX([BBL],BBL['qMDDmet'])

############## CS30_MDD150 #################
if not load_data:
    BBL2 = {}
    BBL2.update(ALL)
    BBL2['desc'] = {'name':'bbl-f75-p33_CS30_MDD150','comp_to':'FT','CS_niter':30,'NITER':150,'random_seed':'y'}
    BBL2['path'] = os.path.join(dn,'data','bblM_20130104_pH6_5C_0Murea_CS30_MDD150','analysis_FT','int_corr_ft_method_all_awk_full')
    #BBL2['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_CS30_MDD150','analysis_FT','int_corr_ft_method_all_awk_full')
    if not multiprocess: sys.stdout = TB.Logger(os.path.join(outdir,"%s.log"%(BBL2['desc']['name']))) #Log to file

    TB.getstat(BBL2,BBL2['qMDDmet'])
    TB.getser(BBL2,BBL2['qMDDmet'])
    TB.sortInt(BBL2,BBL2['qMDDmet'])
    #l = range(94,98,2); BBL2['NIarr']['CS'] = l[::-1]; BBL2['NIarr']['MDD'] = l[::-1]; BBL2['NIarr']['FT'] = l[::-1]
    TB.getdecay(BBL2,BBL2['qMDDmet'])
    #TB.plotdecays([BBL2],BBL2['qMDDmet'],peaks=[3],fss=range(0,5,5))
    TB.getrates(BBL2,BBL2['qMDDmet'])
    #TB.plotrates([BBL2],BBL2['qMDDmet'], peaks=[1,2,3,4]) 
    ##psel = [3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 33]
    TB.getglobfit(BBL2,BBL2['qMDDmet'])

    if data_format == 'pickle' and save_data: pickle.dump(BBL2,open(os.path.join(outdir,"%s.pickle"%(BBL2['desc']['name'])),"wb"))
    elif data_format == 'yaml' and save_data: yaml.dump(BBL2,open(os.path.join(outdir,"%s.yaml"%(BBL2['desc']['name'])),"w"))

elif load_data:
    if data_format == 'pickle': BBL2 = pickle.load(open(os.path.join(outdir,"%s.pickle"%(BBL2['desc']['name'])),"rb"))
    elif data_format == 'yaml': BBL2 = yaml.load(open(os.path.join(outdir,"%s.pickle"%(BBL2['desc']['name'])),"rb"))

#TB.plot_kEX([BBL],BBL['qMDDmet'])

###########################################


