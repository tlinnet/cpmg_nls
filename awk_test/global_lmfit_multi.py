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

#############################################################################
## Start
#############################################################################
limit = 0.6   # Limit set for chisq test, to select peaks
jobs = multiprocessing.cpu_count()-1
nr_gen_datapd = 4
multi = 1

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

def multi_sim_fit(inp):
    p=inp[0]
    X=inp[1]
    Y=inp[2]
    par = lmfit.Parameters(); par.add('b', value=2.0, vary=True, min=0.0); par.add('a', value=2.0, vary=True); par.add('c', value=2.0, vary=True)
    lmf = lmfit.minimize(sim, par, args=(X, Y),method='leastsq')
    dic2 = unpack_sim(par,X,Y)
    return(dic2)

def unpack_sim(par,X,Y):
    dic2 = {}
    Yfit = sim(par,X)
    dic2['Yfit']=Yfit
    dic2['a_value'] = par['a'].value; dic2['a_stderr'] = par['a'].stderr
    dic2['b_value'] = par['b'].value; dic2['b_stderr'] = par['b'].stderr
    dic2['c_value'] = par['c'].value; dic2['c_stderr'] = par['c'].stderr
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
def err_global(pars,x_arr,y_arr,sel_p):
    toterr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        par = lmfit.Parameters()
        par.add('b', value=pars['b'].value, vary=True, min=0.0)
        par.add('a', value=pars['a%s'%p].value, vary=True)
        par.add('c', value=pars['c%s'%p].value, vary=True)
        x = x_arr[i]
        y = y_arr[i]
        Yfit = sim(par,x)
        erri = Yfit - y
        toterr = np.concatenate((toterr, erri))
    #print len(toterr), type(toterr)
    return toterr
#
def unpack_global(par, p_list):
    for i in range(len(p_list)):
        p = p_list[i]
        par_calc = lmfit.Parameters()
        b = dic['gfit']['par']['b']
        a = dic['gfit']['par']['a%s'%p]
        c = dic['gfit']['par']['c%s'%p]
        #par_calc['b'] = b; par_calc['a'] = a; par_calc['c'] = c
        #dic[str(p)]['gfit']['par'] = par
        # Calc other parameters for the fit
        #Yfit = sim(par, dic[str(p)]['X'])
        #dic[str(p)]['gfit']['Yfit'] = Yfit
        #residual = Yfit - dic[str(p)]['Yran']
        #dic[str(p)]['gfit']['residual'] = residual
        #chisq = sum(residual**2)
        #dic[str(p)]['gfit']['chisq'] = chisq
        #NDF = len(residual)-len(par)
        #dic[str(p)]['gfit']['NDF'] = NDF
        #dic[str(p)]['gfit']['what_is_this_called'] = np.sqrt(chisq/NDF)
        #dic[str(p)]['gfit']['redchisq'] = chisq/NDF
    return()
################ Random peak data generator ###########################
def gendat(nr):
    pd = {}
    for i in range(1,nr+1):
        b = 0.15
        a = np.random.random_integers(1, 80)/10.
        c = np.random.random_integers(1, 80)/100.
        par = lmfit.Parameters(); par.add('b', value=b, vary=True); par.add('a', value=a, vary=True); par.add('c', value=c, vary=True)
        pd[str(i)] = par
    return(pd)

def main():
    startTime = datetime.now()
    # set up the data
    data_x = np.linspace(0, 20, 50)
    pd = {} # Parameter dictionary, the "true" values of the data sets
    par = lmfit.Parameters(); par.add('b', value=0.15, vary=True); par.add('a', value=2.5, vary=True); par.add('c', value=0.5, vary=True)
    pd['1'] = par # parameters for the first trajectory
    par = lmfit.Parameters(); par.add('b', value=0.15, vary=True); par.add('a', value=4.2, vary=True); par.add('c', value=0.2, vary=True)
    pd['2'] = par       # parameters for the second trajectory, same b
    par = lmfit.Parameters(); par.add('b', value=0.15, vary=True); par.add('a', value=1.2, vary=True); par.add('c', value=0.3, vary=True)
    pd['3'] = par       # parameters for the third trajectory, same b
    pd = gendat(nr_gen_datapd)  # You can generate a large number of peaks to test
    #
    #Start making a dictionary, which holds all data
    global dic; dic = {}; dic['peaks']=range(1,len(pd)+1)
    for p in dic['peaks']:
        dic['%s'%p] = {}
        dic[str(p)]['X'] = data_x
        dic[str(p)]['Y'] = sim(pd[str(p)],data_x)
        dic[str(p)]['Yran'] = dic[str(p)]['Y'] + np.random.normal(size=len(dic[str(p)]['Y']), scale=0.12)
        dic[str(p)]['fit'] = {}  # Make space for future fit results
        dic[str(p)]['gfit'] = {}  # Male space for future global fit results
    #print "keys for start dictionary:", dic.keys()
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
            Yfit = dic2['Yfit']
            #b=dic2['b_value'];a=dic2['a_value'];c=dic2['c_value']
            #print "Best fit parameter for peak %s. %3.2f %3.2f %3.2f."%(p,b,a,c),
            #print "Compare to real paramaters. %3.2f %3.2f %3.2f."%(pd[str(p)]['b'].value,pd[str(p)]['a'].value,pd[str(p)]['c'].value)

    elif multi == 1:
        inp_arr = []
        for p in dic['peaks']:
            inp_arr.append([p,dic[str(p)]['X'],dic[str(p)]['Yran']])
        #res = handythread.foreach(multi_sim_fit,inp_arr,threads=jobs,return_=True)
        pool = Pool(processes=jobs)
        res = pool.map(multi_sim_fit, inp_arr)
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
    print 'Done with test on chisqr %s\n'%(datetime.now()-startTime)
    #print sel_p
    #
    # Global fitting
    # Pick up x,y-values and parameters that passed the test
    X_arr = []
    Y_arr = []
    P_arr = lmfit.Parameters(); P_arr.add('b', value=1.0, vary=True, min=0.0)
    dic['gfit'] = {} # Make room for globat fit result
    print "Prepare for global fit %s"%(datetime.now()-startTime)
    for p in sel_p:
        a=dic[str(p)]['fit']['a_value'];c=dic[str(p)]['fit']['c_value']
        X_arr.append(dic[str(p)]['X'])
        Y_arr.append(dic[str(p)]['Yran'])
        P_arr.add('a%s'%p, value=a, vary=True)
        P_arr.add('c%s'%p, value=c, vary=True)
    print "Doing global fit %s"%(datetime.now()-startTime)
    lmf = lmfit.minimize(err_global, P_arr, args=(X_arr, Y_arr, sel_p),method='leastsq')
    print "Done with global fit %s"%(datetime.now()-startTime)
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
        Yfit = dic2['Yfit']
    print "global fit unpacked %s \n"%(datetime.now()-startTime)
##      #
    # Check result
    #for p in sel_p:
    #s    ip= pd[str(p)]; sp = dic[str(p)]['fit']['par']; gp = dic[str(p)]['gfit']['par']
        #print p, "Single fit. %3.2f %3.2f %3.2f"%(sp['b'].value,sp['a'].value,sp['c'].value),
        #print "Global fit. %3.2f %3.2f %3.2f"%(gp['b'].value,gp['a'].value,gp['c'].value)
        #print p, "Single fit. %3.2f %3.2f %3.2f"%(sp['b'].value-ip['b'].value,sp['a'].value-ip['a'].value,sp['c'].value-ip['c'].value),
        #print "Global fit. %3.2f %3.2f %3.2f"%(gp['b'].value-ip['b'].value,gp['a'].value-ip['a'].value,gp['c'].value-ip['c'].value)##
##    #
##    # Start plotting
    print "Making figure %s"%(datetime.now()-startTime)
    fig = pl.figure()
    sel_p = sel_p[:9]
    for i in range(len(sel_p)):
        p = sel_p[i]
        # Create figure
        ax = fig.add_subplot('%s1%s'%(len(sel_p),i+1))
        X = dic[str(p)]['X']
        Y = dic[str(p)]['Y']
        Ymeas = dic[str(p)]['Yran']
        Yfit = dic[str(p)]['fit']['Yfit']
        Yfit_global = dic[str(p)]['gfit']['Yfit']
##        rpar = pd[str(p)]
##        fpar = dic[str(p)]['fit']['par']
##        gpar = dic[str(p)]['gfit']['par']
##        fchisq = dic[str(p)]['fit']['lmf'].chisqr
##        gchisq = dic[str(p)]['gfit']['chisq']
##        # plot
##        ax.plot(X,Y,".-",label='real. Peak: %s'%p)
##        ax.plot(X,Ymeas,'o',label='Measured (real+noise)')
##        ax.plot(X,Yfit,'.-',label='leastsq fit. chisq:%3.3f'%fchisq)
##        ax.plot(X,Yfit_global,'.-',label='global fit. chisq:%3.3f'%gchisq)
##        # annotate
##        ax.annotate('p%s. real    par: %1.3f %1.3f %1.3f'%(p, rpar['b'].value,rpar['a'].value,rpar['c'].value), xy=(1,1), xycoords='data', xytext=(0.4, 0.8), textcoords='axes fraction')
##        ax.annotate('p%s. single  par: %1.3f %1.3f %1.3f'%(p, fpar['b'].value,fpar['a'].value,fpar['c'].value), xy=(1,1), xycoords='data', xytext=(0.4, 0.6), textcoords='axes fraction')
##        ax.annotate('p%s. global  par: %1.3f %1.3f %1.3f'%(p, gpar['b'].value,gpar['a'].value,gpar['c'].value), xy=(1,1), xycoords='data', xytext=(0.4, 0.4), textcoords='axes fraction')
##        # set title and axis name
##        #ax.set_title('Fitting for peak %s'%p)
##        ax.set_ylabel('Decay')
##        # Put legend to the right
##        box = ax.get_position()
##        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
##        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
##        ax.grid('on')
##    ax.set_xlabel('Time')
##    #
##    print "ready to show figure %s"%(datetime.now()-startTime)
##    #pl.show()
##    return()

if __name__ == "__main__":
    res = main()