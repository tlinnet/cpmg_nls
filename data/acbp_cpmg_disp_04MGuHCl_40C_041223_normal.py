from pylab import *
import scipy.optimize
import os, sys
import multiprocessing
import logging
import cPickle as pickle
sys.path.append('/sbinlab2/tlinnet/Desktop/CPMG_code')
#sys.path.append('C:/Users/tlinnet/Desktop/tomat_data/cpmg')
import TB
reload(TB)
import gui_simple

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
ALL['ncyc'] = array([28, 0, 4, 32, 60, 2, 10, 16, 8, 20, 50, 18, 40, 6, 12, 0, 24]) # grep -A 1 "ncyc " procpar
ALL['time_T2'] = 0.06 # grep -A 1 "time_T2 " procpar
ALL['guess'] = {'s_R2':8.0,'s_Domega':3000.0,'s_ka':2.0,'g_ka':2.0}
ALL['NIstop'] = 0
ALL['Flags'] = {'multiprocess':multiprocess,'jobs':jobs}
ALL['ser'] = {'pre':'allplanes_','filee':'.ser'}
ALL['stats'] = {'pre':'allplanes_','filee':'.stats'}
ALL['qMDDmet'] = ['coMDD','CS','MDD','FT'] #  ['coMDD','CS','MDD','FT']

############## Normal #################
load_data = False #True False None
save_data = True # True False

DAT = {}
DAT.update(ALL)
DAT['desc'] = {'name':'acbp-f17-p82','comp_to':'FT','CS_niter':10,'NITER':50,'random_seed':'y'}
DAT['path'] = os.path.join(os.getcwd(),'acbp_cpmg_disp_04MGuHCl_40C_041223_normal.fid','analysis_FT','int_corr_ft_method_all_awk_full')
logout = os.path.join(outdir,"%s.log"%(DAT['desc']['name']))
tlog=logging.FileHandler(logout); tlog.setFormatter(logformatter); logger.addHandler(tlog); logger.setLevel(logging.INFO)
if load_data==False:
    TB.getstat(DAT,DAT['qMDDmet'])
    TB.getser(DAT,DAT['qMDDmet'])
    #TB.plotstats([DAT],DAT['qMDDmet'])
    TB.sortInt(DAT,DAT['qMDDmet'])
    #l = range(128,130,2); DAT['NIarr']['CS'] = l[::-1]; DAT['NIarr']['MDD'] = l[::-1]; DAT['NIarr']['FT'] = l[::-1]
    TB.getrelax(DAT,DAT['qMDDmet']) #x,y,par 
    #TB.plotrelaxation([DAT],DAT['qMDDmet'])
    #TB.plotrelaxation([DAT],['MDD'],NIa=[36],peaks=(['24'])) #TB.plotrelaxation([DAT],DAT['qMDDmet'],peaks=(DAT['peakrange']['MDD'][:10]+['57']))
    psel = [2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 15, 19, 21, 22, 23, 24, 26, 27, 29, 32, 33, 34, 36, 38, 40, 41, 42, 43, 48, 49, 50, 52, 54, 55, 56, 57, 58, 59, 60, 61, 63, 66, 67, 68, 69, 70, 71, 74, 78, 79, 80, 81, 82]
    TB.getglobfit_slow(DAT,DAT['qMDDmet'],psel)
    #TB.plot_globalpar([DAT],DAT['qMDDmet'],globalpar='ka',gkey='relax',ay=[1,5],by=[0,200])

    TB.get_glob_props(DAT,['R2','Domega'],DAT['qMDDmet'],gkey='relax',pkey='R2cpmg_slow')
    #TB.plot_glob_props([DAT],['R2','Domega'],DAT['qMDDmet'],gkey='relax')
    TB.get_glob_pearsons(DAT,['R2','Domega'],DAT['qMDDmet'],gkey='relax',pkey='R2cpmg_slow')
    #TB.plot_single_pearson(DAT,'R2',['coMDD'],NIa=[128,122]) #TB.plot_single_pearson(DAT,'R2',['coMDD','CS'],NIa=[128,122])
    TB.get_glob_pearsons(DAT,['R2','Domega'],DAT['qMDDmet'],Ini=True,gkey='relax',pkey='R2cpmg_slow')
    #TB.plot_single_pearson(DAT,'Domega',mets=['coMDD'],NIa=[98],Ini=True) #TB.plot_single_pearson(DAT,'Domega',['coMDD','CS'],NIa=[128,122],Ini=True)
    #TB.plot_glob_pearsons([DAT],['R2','Domega'],mets=['coMDD'],Ini=True,gkey='relax')

    #warnings.simplefilter('ignore'); gui_simple.loadgui(DAT)
    if save_data:
        pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(DAT['desc']['name'],datetime.date.today())),"wb")
        pickle.dump(DAT,pkl_file)
        pkl_file.close()
elif load_data==True:
    former_date = '2013-04-13'
    pkl_file = open(os.path.join(outdir,"%s_%s.pickle"%(DAT['desc']['name'],former_date)),"rb")
    DAT = pickle.load(pkl_file)
    pkl_file.close()
