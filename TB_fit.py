from pylab import *
import scipy.optimize
import os, sys
import multiprocessing
import logging
import cPickle as pickle
import TB

#Logging setup
logger = logging.getLogger("TB")
logger.setLevel(logging.INFO)
logformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

dn = os.path.dirname(os.path.realpath(__file__)) #Current directory for script path position. See also: os.getcwd() for position of execution.
outdir = os.path.join(os.getcwd(),'pickle')
if not os.path.exists(outdir): os.mkdir(outdir)

multiprocess = False  # Only tried multiprocessing in TB.getdecay. It was not faster.
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
ALL['qMDDmet'] = ['coMDD','CS','MDD','FT'] #  ['coMDD','CS','MDD','FT']

############## Normal #################
load_data = True #True False None
save_data = True # True False

BBL = {}
BBL.update(ALL)
BBL['desc'] = {'name':'bbl-f75-p33','comp_to':'FT','CS_niter':10,'NITER':50,'random_seed':'y'}
BBL['path'] = os.path.join(dn,'data','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
#BBL['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full')
logout = os.path.join(outdir,"%s.log"%(BBL['desc']['name']))
tlog=logging.FileHandler(logout); tlog.setFormatter(logformatter); logger.addHandler(tlog); logger.setLevel(logging.INFO)
if load_data==False:
    TB.getstat(BBL,BBL['qMDDmet'])
    TB.getser(BBL,BBL['qMDDmet'])
    #TB.plotstats([BBL],BBL['qMDDmet'])
    TB.sortInt(BBL,BBL['qMDDmet'])
    #l = range(2,10,2); BBL['NIarr']['CS'] = l[::-1]; BBL['NIarr']['MDD'] = l[::-1]; BBL['NIarr']['FT'] = l[::-1]
    TB.getdecay(BBL,BBL['qMDDmet'])
    #TB.plotdecays([BBL],BBL['qMDDmet'],peaks=[3],fss=range(0,5,5))
    TB.getrates(BBL,BBL['qMDDmet'])
    #TB.plotrates([BBL],BBL['qMDDmet'], peaks=[1,2,3,4]) 
    psel = [3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 33]
    TB.getglobfit(BBL,BBL['qMDDmet'],psel)
    #TB.plot_kEX([BBL],BBL['qMDDmet'])

    TB.get_glob_props(BBL,BBL['qMDDmet'])
    #TB.plot_glob_props([BBL],['coMDD','CS'])

    if save_data: 
        pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(BBL['desc']['name'],datetime.date.today())),"wb")
        pickle.dump(BBL,pkl_file)
        pkl_file.close()
elif load_data==True:
    former_date = '2013-04-10'
    pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(BBL['desc']['name'],former_date)),"rb")
    BBL = pickle.load(pkl_file)
    pkl_file.close()

############## CS30_MDD150 #################
load_data = True #True False None
save_data = True # True False

BBL2 = {}
BBL2.update(ALL)
BBL2['desc'] = {'name':'bbl-f75-p33_CS30_MDD150','comp_to':'FT','CS_niter':30,'NITER':150,'random_seed':'y'}
BBL2['path'] = os.path.join(dn,'data','bblM_20130104_pH6_5C_0Murea_CS30_MDD150','analysis_FT','int_corr_ft_method_all_awk_full')
#BBL2['path'] = os.path.join('/','home','tlinnet','kte','t1rho','bblM_20130104_pH6_5C_0Murea_CS30_MDD150','analysis_FT','int_corr_ft_method_all_awk_full')
logout = os.path.join(outdir,"%s.log"%(BBL2['desc']['name']))
tlog=logging.FileHandler(logout); tlog.setFormatter(logformatter); logger.addHandler(tlog); logger.setLevel(logging.INFO)
if load_data==False:
    TB.getstat(BBL2,BBL2['qMDDmet'])
    TB.getser(BBL2,BBL2['qMDDmet'])
    #TB.plotstats([BBL2],BBL2['qMDDmet'])
    TB.sortInt(BBL2,BBL2['qMDDmet'])
    #l = range(94,98,2); BBL2['NIarr']['CS'] = l[::-1]; BBL2['NIarr']['MDD'] = l[::-1]; BBL2['NIarr']['FT'] = l[::-1]
    TB.getdecay(BBL2,BBL2['qMDDmet'])
    #TB.plotdecays([BBL2],BBL2['qMDDmet'],peaks=[3],fss=range(0,5,5))
    TB.getrates(BBL2,BBL2['qMDDmet'])
    #TB.plotrates([BBL2],BBL2['qMDDmet'], peaks=[1,2,3,4]) 
    psel = [3, 5, 6, 7, 8, 9, 11, 14, 15, 16, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 33]
    TB.getglobfit(BBL2,BBL2['qMDDmet'],psel)
    #TB.plot_kEX([BBL2],BBL2['qMDDmet'])

    TB.get_glob_props(BBL2,BBL2['qMDDmet'])
    #TB.plot_glob_props([BBL2],['coMDD','CS'])
                                  
    if save_data: 
        pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(BBL2['desc']['name'],datetime.date.today())),"wb")
        pickle.dump(BBL2,pkl_file)
        pkl_file.close()
elif load_data==True:
    former_date = '2013-04-10'
    pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(BBL2['desc']['name'],former_date)),"rb")
    BBL2 = pickle.load(pkl_file)
    pkl_file.close()
