#-------------------------------------------------------------------------------
# Name:        Test for global fitting with scipy.optimize.leastsq
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
#
############# Fitting functions ################
def sim(x, p):
    b, a, c  = p # Unpacking of shared variables should come first, then the vary parameters
    return a*np.exp(-b * x) + c
#
def err(p, x, y):
    return sim(x, p) - y
#
def err_global(P_arr, x_arr, y_arr):
    toterr = np.array([])
    s = nr_shared_par # Number of shared parameters. Getting from the set global parameter
    v = nr_vary_par # Number of parameters that vary. Getting from the set global parameter
    for i in range(len(x_arr)):
        par = np.array(P_arr[:s])
        par = np.concatenate((par,P_arr[s+i*v:s+i*v+v]))
        #print p
        x = x_arr[i]
        y = y_arr[i]
        erri = err(par, x, y)
        toterr = np.concatenate((toterr, erri))
    #print len(toterr), type(toterr)
    return toterr
#
def unpack_global(dic, p_list):
    s = nr_shared_par # Number of shared parameters. Getting from the set global parameter
    v = nr_vary_par # Number of parameters that vary. Getting from the set global parameter
    for i in range(len(p_list)):
        p = p_list[i]
        par_shared = dic['gfit']['par'][:s]
        par_vary = dic['gfit']['par'][s+i*v:s+i*v+v]
        par_all = np.concatenate((par_shared, par_vary))
        dic[str(p)]['gfit']['par'] = par_all # Store paramaters
        # Calc other parameters for the fit
        Yfit = sim(dic[str(p)]['X'], par_all)
        dic[str(p)]['gfit']['Yfit'] = Yfit
        residual = Yfit - dic[str(p)]['Yran']
        dic[str(p)]['gfit']['residual'] = residual
        chisq = sum(residual**2)
        dic[str(p)]['gfit']['chisq'] = chisq
        NDF = len(residual)-len(par_all)
        dic[str(p)]['gfit']['NDF'] = NDF
        dic[str(p)]['gfit']['what_is_this_called'] = np.sqrt(chisq/NDF)
        dic[str(p)]['gfit']['redchisq'] = chisq/NDF
    return()
################ Extract parameters from output of global fit ###########################
def getleastsstat(result):
    # http://mail.scipy.org/pipermail/scipy-user/2009-March/020516.html
    dic = {}
    dic['par'], dic['cov_x'], dic['infodict'], dic['mesg'], dic['ier'] = result
    dic['residual'] = dic['infodict']['fvec']
    dic['chisq']=sum(dic['residual']**2) # calculate final chi square
    dic['NDF']=len(dic['residual'])-len(dic['par'])
    dic['what_is_this_called'] = np.sqrt(dic['chisq']/dic['NDF'])
    dic['redchisq'] = dic['chisq']/dic['NDF']
    return(dic)
################ Random peak data generator ###########################
def gendat(nr):
    pd = {}
    for i in range(1,nr+1):
        b = 0.15
        a = np.random.random_integers(1, 80)/10.
        c = np.random.random_integers(1, 80)/100.
        pd[str(i)] = [b,a,c]
    return(pd)
#############################################################################
## Start
#############################################################################
limit = 0.6   # Limit set for chisq test, to select peaks
# Global fitting
global nr_shared_par ; nr_shared_par = 1 # Number of shared parameters
global nr_vary_par ; nr_vary_par = 2 # Number of parameters that vary
#############################################################################
# set up the data
data_x = np.linspace(0, 20, 50)
pd = {} # Parameter dictionary, the "true" values of the data sets
pd['1'] = [0.15, 2.5, 0.5]       # parameters for the first trajectory
pd['2'] = [0.15, 4.2, 0.2]       # parameters for the second trajectory, same b
pd['3'] = [0.15, 1.2, 0.3]       # parameters for the third trajectory, same b
pd = gendat(9)  # You can generate a large number of peaks to test
#
#Start making a dictionary, which holds all data
dic = {}; dic['peaks']=range(1,len(pd)+1)
for p in dic['peaks']:
    dic['%s'%p] = {}
    dic[str(p)]['X'] = data_x
    dic[str(p)]['Y'] = sim(data_x, pd[str(p)])
    dic[str(p)]['Yran'] = dic[str(p)]['Y'] + np.random.normal(size=len(dic[str(p)]['Y']), scale=0.12)
    dic[str(p)]['fit'] = {}  # Make space for future fit results
    dic[str(p)]['gfit'] = {}  # Male space for future global fit results
#print "keys for start dictionary:", dic.keys()
#
# independent fitting of the trajectories
for p in dic['peaks']:
    pguess = [2.0, 2.0, 2.0]
    res = scipy.optimize.leastsq(err, pguess, args=(dic[str(p)]['X'], dic[str(p)]['Yran']), full_output=1)
    res_dic = getleastsstat(res)
    dic[str(p)]['fit'].update(res_dic)
    Yfit = sim(dic[str(p)]['X'], dic[str(p)]['fit']['par'])
    #Yfit2 = dic[str(p)]['Yran']+res_dic['residual']
    #print sum(Yfit-Yfit2), "Test for difference in two ways to get the fitted Y-values "
    dic[str(p)]['fit']['Yfit'] = Yfit
    print "Best fit parameter for peak %s"%p, dic[str(p)]['fit']['par'],
    print "Compare to real paramaters", pd[str(p)]
#
# Make a selection flag, based on some test. Now a chisq value, but could be a Ftest between a simple and advanced model fit.
sel_p = []
for p in dic['peaks']:
    chisq = dic[str(p)]['fit']['chisq']
    if chisq < limit:
        dic[str(p)]['Pval'] = 1.0
        #print "Peak %s passed test"%p
        sel_p.append(p)
    else:
        dic[str(p)]['Pval'] = False
#print sel_p
#
# Global fitting
# Pick up x,y-values and parameters that passed the test
X_arr = []
Y_arr = []
P_arr = [1.0] # Pack guess for shared values in first.
dic['gfit'] = {} # Make room for globat fit result
for p in sel_p:
    par = dic[str(p)]['fit']['par']
    X_arr.append(dic[str(p)]['X'])
    Y_arr.append(dic[str(p)]['Yran'])
    P_arr.append(par[1])
    P_arr.append(par[2])
#print P_arr
res = scipy.optimize.leastsq(err_global, P_arr, args=(X_arr, Y_arr), full_output=1)  # Do the fitting
res_dic = getleastsstat(res) # Extract parameters from result
dic['gfit'].update(res_dic) # Update the data dictionary from the returned parameter
unpack_global(dic, sel_p) # Unpack the paramerts into the selected peaks
#
# Check result
for p in sel_p:
    print p, "Single fit%s"%dic[str(p)]['fit']['par'], "Global fit%s"%dic[str(p)]['gfit']['par'] , "Real par%s"%pd[str(p)]
    #print p, "Single fit%s"%(dic[str(p)]['fit']['par']-pd[str(p)]), "Global fit%s"%(dic[str(p)]['gfit']['par']-pd[str(p)])
#
# Start plotting
fig = pl.figure('Peak')
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
    fchisq = dic[str(p)]['fit']['chisq']
    gchisq = dic[str(p)]['gfit']['chisq']
    # plot
    ax.plot(X,Y,".-",label='real. Peak: %s'%p)
    ax.plot(X,Ymeas,'o',label='Measured (real+noise)')
    ax.plot(X,Yfit,'.-',label='leastsq fit. chisq:%3.3f'%fchisq)
    ax.plot(X,Yfit_global,'.-',label='global fit. chisq:%3.3f'%gchisq)
    # annotate
    ax.annotate('p%s. real    par: %1.3f %1.3f %1.3f'%(p, rpar[0],rpar[1],rpar[2]), xy=(1,1), xycoords='data', xytext=(0.4, 0.8), textcoords='axes fraction')
    ax.annotate('p%s. single  par: %1.3f %1.3f %1.3f'%(p, fpar[0],fpar[1],fpar[2]), xy=(1,1), xycoords='data', xytext=(0.4, 0.6), textcoords='axes fraction')
    ax.annotate('p%s. global  par: %1.3f %1.3f %1.3f'%(p, gpar[0],gpar[1],gpar[2]), xy=(1,1), xycoords='data', xytext=(0.4, 0.4), textcoords='axes fraction')
    # set title and axis name
    #ax.set_title('Fitting for peak %s'%p)
    ax.set_ylabel('Decay')
    # Put legend to the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
    ax.grid('on')
ax.set_xlabel('Time')
#
pl.show()