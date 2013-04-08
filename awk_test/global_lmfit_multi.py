#-------------------------------------------------------------------------------
# Name:        Test for global fitting with lmfit # http://newville.github.com/lmfit-py/
# Purpose:      To understand how to do global fitting
# Thanks to:   Jonathan, Josef, Charles, Matt Newville and especially Jonathan Helmus
# Reference:   http://mail.scipy.org/pipermail/scipy-user/2013-April/034401.html
# Author:      Troels Emtekaer Linnet
#
# Created:     04-04-2013
# Copyright:   (c) tlinnet 2013
# Licence:     Free
#-------------------------------------------------------------------------------
#
import pylab as pl
import numpy as np
import scipy.optimize
import lmfit
from datetime import datetime
import multiprocessing
from multiprocessing import Pool
import os
import cPickle as pickle
#import pickle as pickle
import yaml
from handythread import foreach

#############################################################################
## Start
#############################################################################
limit = 0.6   # Limit set for chisq test, to select peaks
jobs = multiprocessing.cpu_count()-1
nr_gen_datapeaks = 6000

# Directory for saving
outdir = os.path.join(os.getcwd(),'data')
if not os.path.exists(outdir): os.mkdir(outdir)
data_format = 'pickle' # 'pickle' 'yaml'
load_data = False #True False
save_data = False

multi = 1
if multi==1:
    print "Using Cores=%s for computation"%jobs

############# Fitting functions ################
def sim(pars,x,data=None,eps=None):
    a = pars['a'].value
    b = pars['b'].value
    c = pars['c'].value
    model = a*np.exp(-b*x)+c
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data)/eps

def sim_calc(inp):
    pars,x = inp
    a = pars['a']['val']
    b = pars['b']['val']
    c = pars['c']['val']
    model = a*np.exp(-b*x)+c
    return model

def multi_sim(inp):
    p=inp[0]
    X=inp[1]
    Y=inp[2]
    par = inp[3]
    parlmf = lmfit.Parameters();
    parlmf.add('b', value=par['b']['val'], vary=True, min=par['b']['min']); parlmf.add('a', value=par['a']['val'], vary=True); parlmf.add('c', value=par['c']['val'], vary=True)
    lmf = lmfit.minimize(sim, parlmf, args=(X, Y),method='leastsq')
    dic2 = unpack_sim(parlmf,X,Y)
    return(dic2)

def unpack_sim(par,X,Y):
    dic2 = {}
    Yfit = sim(par,X)
    dic2['Yfit']=Yfit
    dic2['par']={}
    dic2['par']['a_v'] = par['a'].value; dic2['par']['a_e'] = par['a'].stderr
    dic2['par']['b_v'] = par['b'].value; dic2['par']['b_e'] = par['b'].stderr
    dic2['par']['c_v'] = par['c'].value; dic2['par']['c_e'] = par['c'].stderr
    residual = Yfit - Y
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    return(dic2)
#
def err_global_WRONG(pars,x_arr,y_arr,sel_p): #Defining the lmfit.Parameters() each time, takes a LOOONG time. Double the time at least.
    toterr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        par = lmfit.Parameters()
        par.add('b', value=pars['b'].value, vary=True)
        par.add('a', value=pars['a%s'%p].value, vary=True)
        par.add('c', value=pars['c%s'%p].value, vary=True)
        x = x_arr[i]
        y = y_arr[i]
        Yfit = sim(par,x)
        erri = Yfit - y
        toterr = np.concatenate((toterr, erri))
    #print len(toterr), type(toterr)
    return toterr

def err_global(pars,x_arr,y_arr,sel_p):  # Much faster to use a dictionary, and pass to a calc function.
    toterr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        par = {'b':{'val':pars['b'].value,},'a':{'val':pars['a%s'%p].value},'c':{'val':pars['c%s'%p].value}}
        x = x_arr[i]
        y = y_arr[i]
        Yfit = sim_calc([par,x])
        erri = Yfit - y
        toterr = np.concatenate((toterr, erri))
    #print len(toterr), type(toterr)
    #print pars['b'].value
    return toterr

def multi_err_global_WRONG(pars,x_arr,y_arr,sel_p):
    toterr = np.array([])
    inp_arr = []
    for i in range(len(sel_p)):
        p = sel_p[i]
        par = {'b':{'val':pars['b'].value,'min':0.0},'a':{'val':pars['a%s'%p].value},'c':{'val':pars['c%s'%p].value}}
        x = x_arr[i]
        y = y_arr[i]
        inp_arr.append([par,x])
    pool = Pool(processes=jobs)
    res = pool.map(sim_calc, inp_arr)
    for i in range(len(sel_p)):
        Yfit = res[i]
        y = y_arr[i]
        erri = Yfit - y
        toterr = np.concatenate((toterr, erri))
    print pars['b'].value
    return toterr

def plotfunction(dic,pd):
    sel_p = dic['sel_p']
    # Start plotting
    fig = pl.figure()
    sel_p = sel_p[:4]
    for i in range(len(sel_p)):
        p = sel_p[i]
        # Create figure
        ax = fig.add_subplot('%s1%s'%(len(sel_p),i+1))
        X = dic[str(p)]['X']
        Y = dic[str(p)]['Y']
        Ymeas = dic[str(p)]['Yran']
        Yfit = dic[str(p)]['fit']['Yfit']
        Yfit_global = dic[str(p)]['gfit']['Yfit']
        rpar = pd[str(p)]
        fpar = dic[str(p)]['fit']['par']
        gpar = dic[str(p)]['gfit']['par']
        fchisq = dic[str(p)]['fit']['chisqr']
        gchisq = dic[str(p)]['gfit']['chisqr']
        # plot
        ax.plot(X,Y,".-",label='real. Peak: %s'%p)
        ax.plot(X,Ymeas,'o',label='Measured (real+noise)')
        ax.plot(X,Yfit,'.-',label='leastsq fit. chisq:%3.3f'%fchisq)
        ax.plot(X,Yfit_global,'.-',label='global fit. chisq:%3.3f'%gchisq)
        # annotate
        ax.annotate('p%s. real    par: %1.3f %1.3f %1.3f'%(p, rpar['b_v'],rpar['a_v'],rpar['c_v']), xy=(1,1), xycoords='data', xytext=(0.4, 0.8), textcoords='axes fraction')
        ax.annotate('p%s. single  par: %1.3f %1.3f %1.3f'%(p, fpar['b_v'],fpar['a_v'],fpar['c_v']), xy=(1,1), xycoords='data', xytext=(0.4, 0.6), textcoords='axes fraction')
        ax.annotate('p%s. global  par: %1.3f %1.3f %1.3f'%(p, gpar['b_v'],gpar['a_v'],gpar['c_v']), xy=(1,1), xycoords='data', xytext=(0.4, 0.4), textcoords='axes fraction')
        # set title and axis name
        ax.set_title('Fitting for peak %s'%p)
        ax.set_ylabel('Decay')
        # Put legend to the right
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
        ax.grid('on')
    ax.set_xlabel('Time')
    pl.show()
    return()

################ Random peak data generator ###########################
def gendat(nr):
    pd = {}
    for i in range(1,nr+1):
        pd[str(i)] = {}
        b = 0.15
        a = np.random.random_integers(1, 80)/10.
        c = np.random.random_integers(1, 80)/100.
        pd[str(i)]['b_v'] = b; pd[str(i)]['a_v'] = a; pd[str(i)]['c_v'] = c
        #pd[str(i)]
    return(pd)
##########################################################################
def main():
    # set up the data
    data_x = np.linspace(0, 20, 50)
    pd = {} # Parameter dictionary, the "true" values of the data sets
    pd = gendat(nr_gen_datapeaks)  # You can generate a large number of peaks to test

    #Start making a dictionary, which holds all data
    global dic; dic = {}; dic['peaks']=range(1,len(pd)+1)
    for p in dic['peaks']:
        dic['%s'%p] = {}
        dic[str(p)]['X'] = data_x
        par = lmfit.Parameters(); par.add('b', value=pd[str(p)]['b_v'], vary=True); par.add('a', value=pd[str(p)]['a_v'], vary=True); par.add('c', value=pd[str(p)]['c_v'], vary=True)
        dic[str(p)]['Y'] = sim(par,data_x)
        dic[str(p)]['Yran'] = dic[str(p)]['Y'] + np.random.normal(size=len(dic[str(p)]['Y']), scale=0.12)
        dic[str(p)]['fit'] = {}  # Make space for future fit results
        dic[str(p)]['gfit'] = {}  # Male space for future global fit results

    # independent fitting of the trajectories
    print "Fitting single peaks N=%s %s"%(len(pd),(datetime.now()-startTime))
    if multi == 0:
        for p in dic['peaks']:
            X = dic[str(p)]['X']
            Y = dic[str(p)]['Yran']
            par = lmfit.Parameters(); par.add('b', value=2.0, vary=True, min=0.0); par.add('a', value=2.0, vary=True); par.add('c', value=2.0, vary=True)
            lmf = lmfit.minimize(sim, par, args=(X, Y),method='leastsq')
            dic2 = unpack_sim(par,X,Y)
            dic[str(p)]['fit'].update(dic2)
            #b=dic2['par']['b_v'];a=dic2['par']['a_v'];c=dic2['par']['c_v']
            #print "Best fit parameter for peak %s. %3.2f %3.2f %3.2f."%(p,b,a,c),
            #print "Compare to real paramaters. %3.2f %3.2f %3.2f."%(pd[str(p)]['b_v'],pd[str(p)]['a_v'],pd[str(p)]['c_v'])

    elif multi == 1:
        inp_arr = []
        for p in dic['peaks']:
            par = {'b':{'val':2.0,'min':0.0},'a':{'val':2.0},'c':{'val':2.0}}
            inp_arr.append([p,dic[str(p)]['X'],dic[str(p)]['Yran'],par])
        pool = Pool(processes=jobs)
        res = pool.map(multi_sim, inp_arr)
        #res = foreach(multi_sim,inp_arr,threads=jobs,return_=True)
        for i in range(len(dic['peaks'])):
            p = dic['peaks'][i]
            dic2 = res[i]
            dic[str(p)]['fit'].update(dic2)

    print "Done with fitting single peaks %s\n"%(datetime.now()-startTime)
    # Make a selection flag, based on some test. Now a chisq value, but could be a Ftest between a simple and advanced model fit.
    print "Make a test on chisqr %s"%(datetime.now()-startTime)
    sel_p = []
    for p in dic['peaks']:
        chisqr = dic[str(p)]['fit']['chisqr']
        if chisqr < limit:
            dic[str(p)]['Pval'] = 1.0
            #print "Peak %s passed test"%p
            sel_p.append(p)
        else:
            dic[str(p)]['Pval'] = False
    dic['sel_p'] = sel_p
    print 'Done with test on chisqr %s\n'%(datetime.now()-startTime)
    #print sel_p

    # Global fitting
    # Pick up x,y-values and parameters that passed the test
    X_arr = []
    Y_arr = []
    P_arr = lmfit.Parameters(); P_arr.add('b', value=1.0, vary=True, min=0.0)
    dic['gfit'] = {} # Make room for globat fit result
    print "Prepare for global fit %s"%(datetime.now()-startTime)
    for p in sel_p:
        a=dic[str(p)]['fit']['par']['a_v'];c=dic[str(p)]['fit']['par']['c_v']
        X = dic[str(p)]['X']
        Y = dic[str(p)]['Yran']
        X_arr.append(X)
        Y_arr.append(Y)
        P_arr.add('a%s'%p, value=a, vary=True)
        P_arr.add('c%s'%p, value=c, vary=True)
    print "Doing global fit %s"%(datetime.now()-startTime)

    if multi == 0: #or multi == 1
        lmf = lmfit.minimize(err_global, P_arr, args=(X_arr, Y_arr, sel_p),method='leastsq')
        print "Done with global fit %s"%(datetime.now()-startTime)
    elif multi == 1:
        lmf = lmfit.minimize(multi_err_global, P_arr, args=(X_arr, Y_arr, sel_p),method='leastsq')
        #lmf = lmfit.minimize(err_global, P_arr, args=(X_arr, Y_arr, sel_p),method='leastsq')

    #Unpack results
    for i in range(len(sel_p)):
        p = sel_p[i]
        b = P_arr['b']
        a = P_arr['a%s'%p]
        c = P_arr['c%s'%p]
        par_calc = lmfit.Parameters()
        par_calc['b'] = b; par_calc['a'] = a; par_calc['c'] = c
        X = dic[str(p)]['X']
        Y = dic[str(p)]['Yran']
        dic2 = unpack_sim(par_calc,X,Y)
        dic[str(p)]['gfit'].update(dic2)

    print "global fit unpacked %s \n"%(datetime.now()-startTime)

    #Check result
    #for p in sel_p:
        #ip= pd[str(p)]; sp = dic[str(p)]['fit']['par']; gp = dic[str(p)]['gfit']['par']
        #print p, "Single fit. %3.2f %3.2f %3.2f"%(sp['b_v'],sp['a_v'],sp['c_v']),
        #print "Global fit. %3.2f %3.2f %3.2f"%(gp['b_v'],gp['a_v'],gp['c_v'])
        #print p, "Single fit. %3.2f %3.2f %3.2f"%(sp['b_v']-ip['b_v'],sp['a_v']-ip['a_v'],sp['c_v']-ip['c_v']),
        #print "Global fit. %3.2f %3.2f %3.2f"%(gp['b_v']-ip['b_v'],gp['a_v']-ip['a_v'],gp['c_v']-ip['c_v'])##

    #Saving data
    if data_format == 'pickle' and save_data:
        pickle.dump( dic, open( os.path.join(outdir,"test_dic.pickle"), "wb" ) )
        pickle.dump( pd, open( os.path.join(outdir,"test_pd.pickle"), "wb" ) )
    elif data_format == 'yaml' and save_data:
        yaml.dump( dic, open( os.path.join(outdir,"test_dic.yaml"), "w" ) )
        yaml.dump( pd, open( os.path.join(outdir,"test_pd.yaml"), "w" ) )
    print "Making figure %s"%(datetime.now()-startTime)
    #plotfunction(dic,pd)

if __name__ == "__main__":
    startTime = datetime.now()
    if not load_data:
        main()
    elif load_data:
        if data_format == 'pickle':
            dic = pickle.load( open( os.path.join(outdir,"test_dic.pickle"), "rb" ) )
            pd = pickle.load( open( os.path.join(outdir,"test_pd.pickle"), "rb" ) )
        elif data_format == 'yaml':
            dic = pickle.load( open( os.path.join(outdir,"test_dic.pickle"), "rb" ) )
            pd = pickle.load( open( os.path.join(outdir,"test_pd.pickle"), "rb" ) )
        plotfunction(dic,pd)
