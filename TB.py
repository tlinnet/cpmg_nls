from pylab import *
import numpy
import scipy.optimize
#import matplotlib as mpl
#mpl.rcParams['text.usetex']=True
import scipy.stats.distributions
import logging
import collections
from datetime import datetime
import os
import lmfit #See http://newville.github.com/lmfit-py/parameters.html
from multiprocessing import Pool
module_logger = logging.getLogger("TB.TB")

#### Standard size of some figures
figsize = 12
titfont = 26
labfont = 12
figfont = 12

##### To handle "RuntimeWarning" from fitting functions to be handled as errors
import warnings
warnings.simplefilter('error')
#warnings.simplefilter('ignore')

####################################### Fit functions ##################################
#Function R2cpmg_fast(w,x) : FitFunc
#	Wave w
#	Variable x
#   //CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
#   //CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
#   //CurveFitDialog/ Equation:
#   //CurveFitDialog/    f(x) = R20+(Psi/kEX)*(1-((4*x)/kEX)*tanh(kEX/(4*x)))
#   //CurveFitDialog/
#   //CurveFitDialog/ End of Equation
#   //CurveFitDialog/ Independent Variables 1
#   //CurveFitDialog/ x
#   //CurveFitDialog/ Coefficients 3
#   //CurveFitDialog/ w[0] = R20
#   //CurveFitDialog/ w[1] = kEX
#   //CurveFitDialog/ w[2] = Psi
#	   return w[0]+(w[2]/w[1])*(1-((4*x)/w[1])*tanh(w[1]/(4*x)))
##############################
def f_R2s(pars,time,data=None):
    R2 = pars['R2'].value
    model = R2+zeros(len(time))
    if data is None:
        return model
    return (model-data)

def f_R2s_calc(par,time):
    R2 = par['R2_v']
    model = R2+zeros(len(time))
    return model

def unpack_f_R2s(par,X,Y,lmf=None):
    dic2 = {}
    Yfit = f_R2s(par,X)
    dic2['Yfit']=Yfit
    dic2['par']={}
    dic2['par']['R2_v'] = par['R2'].value
    dic2['par']['R2_e'] = par['R2'].stderr
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    #print NDF, lmf.nfree
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_fit = array([X,Y,Yfit]).T
    dic2['X_Y_Fit'] = x_y_fit
    return(dic2)
#################################################
#Function R2cpmg(w,x) : FitFunc
#    Wave w
#    Variable x
#
#   //CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
#   //CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
#   //CurveFitDialog/ Equation:
#   //CurveFitDialog/ Variable psi,z,Dm,Dp,nm,np
#   //CurveFitDialog/
#   //CurveFitDialog/ psi = kEX^2-w^2
#   //CurveFitDialog/ z = -2*w*kEX*(2*pA-1)
#   //CurveFitDialog/ Dp = 0.5*(1+(psi+2*w^2)/sqrt(psi^2+z^2))
#   //CurveFitDialog/ Dm = 0.5*(-1+(psi+2*w^2)/sqrt(psi^2+z^2))
#   //CurveFitDialog/ np = 1.0/2/sqrt(2)/x*sqrt(psi+sqrt(psi^2+z^2))
#   //CurveFitDialog/ nm = 1.0/2/sqrt(2)/x*sqrt(-psi+sqrt(psi^2+z^2))
#   //CurveFitDialog/
#   //CurveFitDialog/ f(x) = R20+0.5*(kEX-2*x*acosh(Dp*cosh(np)-Dm*cos(nm)))
#   //CurveFitDialog/ End of Equation
#   //CurveFitDialog/ Independent Variables 1
#   //CurveFitDialog/ x
#   //CurveFitDialog/ Coefficients 4
#   //CurveFitDialog/ w[0] = R20
#   //CurveFitDialog/ w[1] = kEX
#   //CurveFitDialog/ w[2] = w
#   //CurveFitDialog/ w[3] = pA

#    Variable psi,z,Dm,Dp,nm,np
#    psi = w[1]^2-w[2]^2
#    z = -2*w[2]*w[1]*(2*w[3]-1)
#    Dp = 0.5*(1+(psi+2*w[2]^2)/sqrt(psi^2+z^2))
#    Dm = 0.5*(-1+(psi+2*w[2]^2)/sqrt(psi^2+z^2))
#    np = 1.0/2/sqrt(2)/x*sqrt(psi+sqrt(psi^2+z^2))
#    nm = 1.0/2/sqrt(2)/x*sqrt(-psi+sqrt(psi^2+z^2))
#    return w[0]+0.5*(w[1]-2*x*acosh(Dp*cosh(np)-Dm*cos(nm)))
#End

def f_R2cpmg_RC72(pars,nu,data=None):
    #http://dx.doi.org.ep.fjernadgang.kb.dk/10.1016/0022-2364(72)90090-X
    # Richard & Carver 1972
    # A general two-site solution for the chemical exchange produced dependence of T2 upon the carr-Purcell pulse separation
    # Journal of Magnetic Resonance (1972), Volume 6, Issue 1, January 1972, Pages 89-105
    #
    # Variable: R2, kEX, w, pA
    R2 = pars['R2'].value
    kEX = pars['kEX'].value
    w = pars['w'].value
    pA = pars['pA'].value
    #
    psi = kEX**2-w**2
    z = -2*w*kEX*(2*pA-1)
    Dp = 0.5*(1+(psi+2*w**2)/sqrt(psi**2+z**2))
    Dm = 0.5*(-1+(psi+2*w**2)/sqrt(psi**2+z**2))
    np = 1.0/2/sqrt(2)/nu*sqrt(psi+sqrt(psi**2+z**2))
    nm = 1.0/2/sqrt(2)/nu*sqrt(-psi+sqrt(psi**2+z**2))
    #
    model = R2+0.5*(kEX-2*nu*arccosh(Dp*cosh(np)-Dm*cos(nm)))
    if data is None:
        return model
    return (model-data)

def f_R2cpmg_RC72_calc(par,nu):
    R2 = par['R2_v']
    kEX = par['kEX_v']
    w = par['w_v']
    pA = par['pA_v']
    #
    psi = kEX**2-w**2
    z = -2*w*kEX*(2*pA-1)
    Dp = 0.5*(1+(psi+2*w**2)/sqrt(psi**2+z**2))
    Dm = 0.5*(-1+(psi+2*w**2)/sqrt(psi**2+z**2))
    np = 1.0/2/sqrt(2)/nu*sqrt(psi+sqrt(psi**2+z**2))
    nm = 1.0/2/sqrt(2)/nu*sqrt(-psi+sqrt(psi**2+z**2))
    #
    model = R2+0.5*(kEX-2*nu*arccosh(Dp*cosh(np)-Dm*cos(nm)))
    return model

def f_R2cpmg_RC72_fields(pars,fields,nu_a,R2eff_a,R2eff_err_a=None):
    toterr = np.array([])
    #print pars['ka'].value
    for i,field in enumerate(fields):
        nu =nu_a[i];R2eff=R2eff_a[i]
#        if R2eff_err_a is not None:
#            R2eff_err=R2eff_err_a[i]
        # Much faster to use a dictionary, and pass to a calc function.
        wcalc = float(field)/100.0*pars['w'].value
        #par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX_%s'%field].value,'pA_v':pars['pA_%s'%field].value,'w_v':wcalc}
        #par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX_%s'%field].value,'pA_v':pars['pA'],'w_v':wcalc}
        par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX'].value,'pA_v':pars['pA'].value,'w_v':wcalc}
        Yfit = f_R2cpmg_RC72_calc(par,nu)
        if R2eff_err_a is None:
            erri = Yfit - R2eff
        else:
            erri = (Yfit - R2eff)/R2eff_err
        toterr = np.concatenate((toterr, erri))
    return toterr

def unpack_f_R2cpmg_RC72_fields(pars,fields,nu_a,R2eff_a,lmf=None):
    dic2 = {}
    dic2['par']={}
    dic2['par']['kEX_v'] = par['kEX'].value; dic2['par']['kEX_e'] = par['kEX'].stderr
    dic2['par']['pA_v'] = par['pA'].value; dic2['par']['pA_e'] = par['pA'].stderr
    dic2['par']['w_v'] = par['w'].value; dic2['par']['w_e'] = par['w'].stderr

    Yfits = []
    residuals = []
    totresiduals = np.array([])
    for i,field in enumerate(fields):
        nu =nu_a[i];R2eff=R2eff_a[i]
        Yfit = f_R2cpmg_RC72_calc(par,nu)
        Yfits.append(Yfit)
        dic2['par']['R2_v'] = par['R2_%s'%field].value; dic2['par']['R2_e'] = par['R2_%s'%field].stderr
        residual = Yfit - R2eff
        residuals.append(residual)
        totresidual = np.concatenate((totresidual, residual))
    dic2['Yfits']=Yfits
    dic2['residuals'] = residuals
    dic2['totresidual'] = totresidual
    chisqr = sum(totresidual**2)
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    #print NDF, lmf.nfree
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_fit = array([nu,Y,Yfit]).T
    dic2['X_Y_Fit'] = x_y_fit
    return(dic2)

###########################################

def f_R2cpmg_slow(pars,nu,data=None):
    #The inverted chevron plot measured by NMR relaxation reveals a native-like unfolding intermediate in acyl-CoA binding protein.
    #Kaare Teilum, Flemming M Poulsen, Mikael Akke in Proceedings of the National Academy of Sciences of the United States of America (2006)
    # Tollinger_kay http://pubs.acs.org/doi/abs/10.1021/ja011300z Slow Dynamics in Folded and Unfolded States of an SH3 Domain,
    # Martin Tollinger , Nikolai R. Skrynnikov , Frans A. A. Mulder , Julie D. Forman-Kay ,and Lewis E. Kay.  J. Am. Chem. Soc., 2001, 123 (46)
    R2 = pars['R2'].value
    Domega = pars['Domega'].value
    ka = pars['ka'].value
    tau_cpmg = 1.0/(4*nu) #Page 6881
    model = R2+ka*(1.0-sin(Domega*tau_cpmg)/(Domega*tau_cpmg))
    if data is None:
        return model
    return (model-data)

def f_R2cpmg_slow_calc(par,nu):
    R2 = par['R2_v']
    Domega = par['Domega_v']
    ka = par['ka_v']
    tau_cpmg = 1.0/(4*nu) #Page 6881
    model = R2+ka*(1.0-sin(Domega*tau_cpmg)/(Domega*tau_cpmg))
    return model

def unpack_f_R2cpmg_slow(par,nu,Y,lmf=None):
    dic2 = {}
    Yfit = f_R2cpmg_slow(par,nu)
    dic2['Yfit']=Yfit
    dic2['par']={}
    dic2['par']['R2_v'] = par['R2'].value; dic2['par']['R2_e'] = par['R2'].stderr
    dic2['par']['Domega_v'] = par['Domega'].value; dic2['par']['Domega_e'] = par['Domega'].stderr
    dic2['par']['ka_v'] = par['ka'].value; dic2['par']['ka_e'] = par['ka'].stderr
    #print par['ka'].value, par['R2'].value, par['Domega'].value
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    #print NDF, lmf.nfree
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_fit = array([nu,Y,Yfit]).T
    dic2['X_Y_Fit'] = x_y_fit
    return(dic2)

def f_R2cpmg_slow_global(pars,sel_p,nu_a,R2eff_a,R2eff_err_a=None):
    toterr = np.array([])
    #print pars['ka'].value
    for i in range(len(sel_p)):
        p = sel_p[i]
        nu =nu_a[i];R2eff=R2eff_a[i]
        if R2eff_err_a is not None:
            R2eff_err=R2eff_err_a[i]
        # Much faster to use a dictionary, and pass to a calc function.
        par = {'ka_v':pars['ka'].value,'R2_v':pars['R2%s'%p].value,'Domega_v':pars['Domega%s'%p].value}
        Yfit = f_R2cpmg_slow_calc(par,nu)
        if R2eff_err_a is None:
            erri = Yfit - R2eff
        else:
            erri = (Yfit - R2eff)/R2eff_err
        toterr = np.concatenate((toterr, erri))
    return toterr

def unpack_global_slow(P_arr,sel_p,nu_a,R2eff_a,R2eff_err_a=None,lmf=None):
    dic2 = {}
    dic_calc = {}
    residual_arr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        nu =nu_a[i];R2eff=R2eff_a[i]
        if R2eff_err_a is not None:
            R2eff_err=R2eff_err_a[i]
        ka_glob = P_arr['ka']
        R2 = P_arr['R2%s'%p]
        Domega = P_arr['Domega%s'%p]
        dic_calc['par']={}
        dic_calc['par']['ka_v'] = ka_glob.value; dic_calc['par']['ka_e'] = ka_glob.stderr
        dic_calc['par']['R2_v'] = R2.value; dic_calc['par']['R2_e'] = R2.stderr
        dic_calc['par']['Domega_v'] = Domega.value; dic_calc['par']['Domega_e'] = Domega.stderr
        Yfit = f_R2cpmg_slow_calc(dic_calc['par'],nu)
        residual = Yfit - R2eff
        residual_arr = np.concatenate((residual_arr, residual))
    dic2['par'] = {}
    dic2['par']['ka_v'] = ka_glob.value; dic2['par']['ka_e'] = ka_glob.stderr
    #print sum(residual_arr), sum(lmf.residual)
    dic2['residual'] = residual_arr
    chisqr = sum(residual_arr**2)
    #print chisqr, lmf.chisqr, sum(lmf.residual**2)
    dic2['chisqr'] = chisqr
    NDF = len(residual_arr)-len(P_arr)
    #print NDF, lmf.nfree
    dic2['NDF'] = NDF
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    #x_y_sigma_fit = array([X,Y,sigma,Yfit]).T
    #dic2['X_Y_Sigma_Fit'] = x_y_sigma_fit
    return(dic2)

######################################

def f_R2cpmg_fast(pars,time,data=None):
    R2 = pars['R2'].value
    phi = pars['phi'].value
    kEXb = pars['kEXb'].value
    model = R2+(phi/kEXb)*(1-2*tanh(kEXb/4/time)/kEXb*time*2)
    if data is None:
        return model
    return (model-data)

def f_R2cpmg_fast_calc(par,time):
    R2 = par['R2_v']
    phi = par['phi_v']
    kEXb = par['kEXb_v']
    model = R2+(phi/kEXb)*(1-2*tanh(kEXb/4/time)/kEXb*time*2)
    return model

def unpack_f_R2cpmg_fast(par,X,Y,lmf=None):
    dic2 = {}
    Yfit = f_R2cpmg_fast(par,X)
    dic2['Yfit']=Yfit
    dic2['par']={}
    dic2['par']['R2_v'] = par['R2'].value; dic2['par']['R2_e'] = par['R2'].stderr
    dic2['par']['phi_v'] = par['phi'].value; dic2['par']['phi_e'] = par['phi'].stderr
    dic2['par']['kEXb_v'] = par['kEXb'].value; dic2['par']['kEXb_e'] = par['kEXb'].stderr
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    #print NDF, lmf.nfree
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_fit = array([X,Y,Yfit]).T
    dic2['X_Y_Fit'] = x_y_fit
    return(dic2)

#################################################
def f_expdecay(pars,time,data=None): #KTE: extract_sums_to_table.pl. Line 68.
    amp = pars['amp'].value
    decay = pars['decay'].value
    model = amp*exp(-decay*time)
    if data is None:
        return model
    return (model-data)

def f_expdecay_calc(par,time):
    amp = par['amp_v']
    decay = par['decay_v']
    model = amp*exp(-decay*time)
    return model

def multi_f_expdecay(inp):
    par = inp[0]
    X=inp[1]
    Y=inp[2]
    parlmf = lmfit.Parameters();
    parlmf.add('amp', value=par['amp_v'], vary=True)
    parlmf.add('decay', value=par['decay_v'], vary=True, min=par['decay_min'])
    lmf = lmfit.minimize(f_expdecay, parlmf, args=(X, Y),method='leastsq')
    dic2 = unpack_f_expdecay(parlmf,X,Y)
    return(dic2)

def unpack_f_expdecay(par,X,Y,lmf=None):
    dic2 = {}
    Yfit = f_expdecay(par,X)
    dic2['Yfit']=Yfit
    dic2['par']={}
    dic2['par']['amp_v'] = par['amp'].value; dic2['par']['amp_e'] = par['amp'].stderr
    dic2['par']['decay_v'] = par['decay'].value; dic2['par']['decay_e'] = par['decay'].stderr
    dic2['R1r_rates'] = par['decay'].value
    dic2['R1r_err'] = par['decay'].stderr
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    dic2['NDF'] = NDF
    #print NDF, lmf.nfree
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_fit = array([X,Y,Yfit]).T
    dic2['X_Y_Fit'] = x_y_fit
    return(dic2)
#################################################
def f_R1r(pars,tiltAngle,data=None,eps=None): #KTE: R1rhoAnalysis.ipf. Line 144. Use see line 96
    #http://newville.github.com/lmfit-py/fitting.html
    R1 = pars['R1'].value
    R2 = pars['R2'].value
    model = R1*cos(tiltAngle*pi/180.0)**2+R2*sin(tiltAngle*pi/180.0)**2
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data)/eps

def f_R1r_calc(par,tiltAngle):
    R1 = par['R1_v']
    R2 = par['R2_v']
    model = R1*cos(tiltAngle*pi/180.0)**2+R2*sin(tiltAngle*pi/180.0)**2
    return model

def unpack_f_R1r(par,X,Y,sigma,lmf=None):
    dic2 = {}
    dic2['par']={}
    dic2['par']['R1_v'] = par['R1'].value; dic2['par']['R1_e'] = par['R1'].stderr
    dic2['par']['R2_v'] = par['R2'].value; dic2['par']['R2_e'] = par['R2'].stderr
#    Yfit = f_R1r(par,X)
    Yfit = f_R1r_calc(dic2['par'],X) #Should be faster
    dic2['Yfit']=Yfit
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    #print NDF, lmf.nfree
    dic2['NDF'] = NDF
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_sigma_fit = array([X,Y,sigma,Yfit]).T
    dic2['X_Y_Sigma_Fit'] = x_y_sigma_fit
    return(dic2)
#################################################
def f_R1r_exch(pars,inp,data=None,eps=None): #KTE: R1rhoAnalysis.ipf. Line 144. Use see line 96
    #http://newville.github.com/lmfit-py/fitting.html
    tiltAngle,omega1=inp
    R1 = pars['R1'].value
    R2 = pars['R2'].value
    kEX = pars['kEX'].value
    phi = pars['phi'].value
    model = R1*cos(tiltAngle*pi/180)**2+(R2+phi*kEX/((2*pi*omega1/tan(tiltAngle*pi/180))**2+(2*pi*omega1)**2+kEX**2))*sin(tiltAngle*pi/180)**2
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data)/eps

def f_R1r_exch_calc(par,inp):
    tiltAngle,omega1=inp
    R1 = par['R1_v']
    R2 = par['R2_v']
    kEX = par['kEX_v']
    phi = par['phi_v']
    model = R1*cos(tiltAngle*pi/180)**2+(R2+phi*kEX/((2*pi*omega1/tan(tiltAngle*pi/180))**2+(2*pi*omega1)**2+kEX**2))*sin(tiltAngle*pi/180)**2
    return model

def unpack_f_R1r_exch(par,X,Y,sigma,lmf=None):
    dic2 = {}
    dic2['par']={}
    dic2['par']['R1_v'] = par['R1'].value; dic2['par']['R1_e'] = par['R1'].stderr
    dic2['par']['R2_v'] = par['R2'].value; dic2['par']['R2_e'] = par['R2'].stderr
    dic2['par']['kEX_v'] = par['kEX'].value; dic2['par']['kEX_e'] = par['kEX'].stderr
    dic2['par']['phi_v'] = par['phi'].value; dic2['par']['phi_e'] = par['phi'].stderr
#    Yfit = f_R1r_exch(par,X)
    Yfit = f_R1r_exch_calc(dic2['par'],X) #Should be faster
    dic2['Yfit']=Yfit
    residual = Yfit - Y
    #print sum(residual - lmf.residual)
    dic2['residual'] = residual
    chisqr = sum(residual**2)
    #print chisqr - lmf.chisqr
    dic2['chisqr'] = chisqr
    NDF = len(residual)-len(par)
    #print NDF, lmf.nfree
    dic2['NDF'] = NDF
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    x_y_sigma_fit = array([X,Y,sigma,Yfit]).T
    dic2['X_Y_Sigma_Fit'] = x_y_sigma_fit
    return(dic2)
#################################################
def f_R1r_exch_global(pars,sel_p,tilt_a,om1_a,R1rex_a,R1rex_err_a=None):
    toterr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        tilt =tilt_a[i];om1=om1_a[i];R1rex=R1rex_a[i]
        if R1rex_err_a is not None:
            R1rex_err=R1rex_err_a[i]
        # Much faster to use a dictionary, and pass to a calc function.
        par = {'kEX_v':pars['kEX'].value,'R1_v':pars['R1%s'%p].value,'R2_v':pars['R2%s'%p].value,'phi_v':pars['phi%s'%p].value}
        datX = [array(tilt), array(om1)]
        Yfit = f_R1r_exch_calc(par,datX)
        if R1rex_err_a is None:
            erri = Yfit - R1rex
        else:
            erri = (Yfit - R1rex)/R1rex_err
        toterr = np.concatenate((toterr, erri))
    return toterr

def unpack_global(P_arr,sel_p,tilt_a,om1_a,R1rex_a,R1rex_err_a=None,lmf=None):
    dic2 = {}
    dic_calc = {}
    residual_arr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        tilt =tilt_a[i];om1=om1_a[i];R1rex=R1rex_a[i]
        if R1rex_err_a is not None:
            R1rex_err=R1rex_err_a[i]
        #par = lmfit.Parameters()
        kEX_glob = P_arr['kEX']
        R1 = P_arr['R1%s'%p]
        R2 = P_arr['R2%s'%p]
        phi = P_arr['phi%s'%p]
        #par['kEX'] = kEX_glob; par['R1'] = R1; par['R2'] = R2; par['phi'] = phi
        #print "Peak %s .chisqr=%3.2f. kEX= "%(p,chisqr,kEX.value)
        #lmfit.printfuncs.report_errors(par)
        dic_calc['par']={}
        dic_calc['par']['kEX_v'] = kEX_glob.value; dic_calc['par']['kEX_e'] = kEX_glob.stderr
        dic_calc['par']['R1_v'] = R1.value; dic_calc['par']['R1_e'] = R1.stderr
        dic_calc['par']['R2_v'] = R2.value; dic_calc['par']['R2_e'] = R2.stderr
        dic_calc['par']['phi_v'] = phi.value; dic_calc['par']['phi_e'] = phi.stderr
        datX = [array(tilt), array(om1)]
        #Yfit = f_R1r_exch(par,datX)
        Yfit = f_R1r_exch_calc(dic_calc['par'],datX)
        residual = Yfit - R1rex
        residual_arr = np.concatenate((residual_arr, residual))
    dic2['par'] = {}
    dic2['par']['kEX_v'] = kEX_glob.value; dic2['par']['kEX_e'] = kEX_glob.stderr
    #print sum(residual_arr), sum(lmf.residual)
    dic2['residual'] = residual_arr
    chisqr = sum(residual_arr**2)
    #print chisqr, lmf.chisqr, sum(lmf.residual**2)
    dic2['chisqr'] = chisqr
    NDF = len(residual_arr)-len(P_arr)
    #print NDF, lmf.nfree
    dic2['NDF'] = NDF
    dic2['what_is_this_called'] = np.sqrt(chisqr/NDF)
    dic2['redchisqr'] = chisqr/NDF
    #x_y_sigma_fit = array([X,Y,sigma,Yfit]).T
    #dic2['X_Y_Sigma_Fit'] = x_y_sigma_fit
    return(dic2)
#################################################
def Ftest(ss1,df1,ss2,df2):
    ## KTE: R1rhoAnalysis.ipf. Line 364.
    #fRatio = ((ss1-ss2)/(df1-df2))/(ss2/df2) # fRatio = ((ss1-ss2)/(df1-df2))/(ss2/df2)  -  Fval=(WSSR1mWSSR2/P2mP1)/(WSSR2val/NmP2)
    #a = df2/2.0 # a = df2/2
    #b = (df1-df2)/2.0 # b = (df1-df2)/2
    #x = df2/(df2+(df1-df2)*fRatio) # x = df2/(df2+(df1-df2)*fRatio)
    #Pvalue = scipy.stats.stats.betai(a, b, x) # Pvalue = betai(a, b, x) - scipy.stats.stats.betai
    # Possible to set distlimit
    distlimit = 0.95
    Fval=((ss1-ss2)/(df1-df2))/(ss2/df2)
    Fdist= scipy.stats.distributions.f.ppf(distlimit, df1 - df2, df2)
    if Fval > Fdist:
        Pval = (1-scipy.stats.distributions.f.cdf(Fval, df1 - df2, df2))
    else:
        Pval = False
    #print Pvalue, Pval, (Pvalue-Pval)
    return Fval, Fdist, Pval
###############################################
def f5(seq, idfun=None):
   # order preserving. http://www.peterbe.com/plog/uniqifiers-benchmark
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def get_ecolor(omega1_arr,colordic):
    col_arr = []
    for omega1 in omega1_arr:
        col = colordic[str(omega1)]
        col_arr.append(col)
    return(col_arr)

def peakname(dic,peaks=False,mets=False):
    if not mets: mets = dic['qMDDmet']
    for met in mets:
        if not peaks: peaks = dic['peakrange']
        for peak in peaks:
            resn = dic['Int'][met][str(peak)]['resn']
            print peak, resn
################### Data loading
def getstat(dic,dt):
    logger = logging.getLogger("TB.TB.getstat")
    startTime = datetime.now()
    fpath = dic['path']
    dic['stats']['data'] = {}
    dic['NIarr']={}
    dic['NImax']={}
    for met in dt:
        pre = dic['stats']['pre']
        filee = dic['stats']['filee']
        path = os.path.join(fpath,met,pre+met+filee)
        data = genfromtxt(path)
        dic['stats']['data'][met] = data
        NIarr = data[:,0].astype(int)
        dic['NIarr'][met] = NIarr
        dic['NImax'][met] = max(NIarr)
    print "Done getting stat file. It took: %s"%(datetime.now()-startTime)
    logger.info("Done getting stat file. It took: %s"%(datetime.now()-startTime))
    return()

def getser(dic,dt):
    logger = logging.getLogger("TB.TB.getser")
    startTime = datetime.now()
    fpath = dic['path']
    dic['ser']['data'] = {}
    dic['ser']['datas'] = {}
    dic['filenr'] = {}
    dic['peakrange'] = {}
    for met in dt:
        pre = dic['ser']['pre']
        filee = dic['ser']['filee']
        path = os.path.join(fpath,met,pre+met+filee)
        data = genfromtxt(path)
        datas = genfromtxt(path,usecols=(6),converters = {6: lambda s: str(s)})
        dic['ser']['data'][met] = data
        dic['ser']['datas'][met] = datas
        dic['peakrange'][met] = range(int(data[:,0][0]),int(data[:,0][-1]+1))
        dic['filenr'][met] = int(data[:,0].size/len(dic['peakrange'][met]))
    print "Done getting ser file. It took: %s"%(datetime.now()-startTime)
    logger.info("Done getting ser file. It took: %s"%(datetime.now()-startTime))
    return()
##################### Sort data into peaks
def sortInt(dic,dt):
    logger = logging.getLogger("TB.TB.sortInt")
    startTime = datetime.now()
    dic['Int'] = {}
    for met in dt:
        dic['Int'][met] = {}
        NIs = dic['NIarr'][met]
        NIsm = dic['NImax'][met]
        peaki = dic['ser']['data'][met][:,0].astype(int)
        pm = max(peaki)
        datas = dic['ser']['datas']
        data = dic['ser']['data'][met]
        for j in range(pm):
            peakn = datas[met][j]
            peak = str(peaki[j])
            dic['Int'][met][peak] = {}
            dic['Int'][met][peak]['resn'] = peakn
            intens = data[:,5][j::pm]
            dic['Int'][met][peak]['FTInt'] = intens
            CS_H = data[:,3][j::pm]
            dic['Int'][met][peak]['CS_H'] = CS_H
            CS_N = data[:,4][j::pm]
            dic['Int'][met][peak]['CS_N'] = CS_N
            dic['Int'][met][peak]['NIInt'] = {}
            for k in range(len(NIs)):
                NI = str(NIs[k])
                NIint = data[:,8+k][j::pm]
                dic['Int'][met][peak]['NIInt'][NI] = NIint
    print "Done sorting intensities. It took: %s"%(datetime.now()-startTime)
    logger.info("Done sorting intensities. It took: %s"%(datetime.now()-startTime))
    return()
##################### Calculate rates
def getdecay(dic,mets=False,NIstop=False):
    logger = logging.getLogger("TB.TB.getdecay")
    startTime = datetime.now()
    multiprocess = dic['Flags']['multiprocess']
    multi_inp_arr = []
    multi_sort_arr = []
    datX = dic['time']
    centerPPM = dic['NMRpar']['centerPPM']
    frq = dic['NMRpar']['frq']
    slicet = len(datX)
    dic['decay'] = {}
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        Int = dic['Int'][met]
        filenr = dic['filenr'][met]
        dic['decay'][met] = {}
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            print "%s - Getting decay rates for NI=%s"%(met,NI)
            dic['decay'][met][str(NI)] = {}
            for peak in dic['peakrange'][met]:
            #for peak in ['1']:
                dic['decay'][met][str(NI)][str(peak)] = {}
                peakname = Int[str(peak)]['resn']
                dic['decay'][met][str(NI)][str(peak)]['resn'] = peakname
                # Getting Chemical Shifts
                CS_N = Int[str(peak)]['CS_N']
                CS_H = Int[str(peak)]['CS_H']
                i = 0
                for fs in range(0,filenr,slicet):
                ##for fs in range(0,5,slicet):
                    fe = fs+slicet
                    FTInt = Int[str(peak)]['FTInt'][fs:fe]
                    NIInt = Int[str(peak)]['NIInt'][str(NI)][fs:fe]
                    divi = FTInt.argmax()
                    datY = FTInt*NIInt/(FTInt[divi]*NIInt[divi])
                    par = lmfit.Parameters()
                    par.add('amp', value=1.0, vary=True, min=0.0)
                    par.add('decay', value=10.0, vary=True, min=0.0)
                    parmulti = {'amp_v':1.0,'decay_v':10.0,'decay_min':0.0}
                    multi_inp_arr.append([parmulti,datX,datY])
                    multi_sort_arr.append([met,NI,peak,fs])
                    dic['decay'][met][str(NI)][str(peak)][str(fs)] = {}
                    if not multiprocess:
                        try:
                            lmf = lmfit.minimize(f_expdecay, par, args=(datX, datY),method='leastsq')
                            dic2 = unpack_f_expdecay(par,datX,datY,lmf)
                            dic['decay'][met][str(NI)][str(peak)][str(fs)].update(dic2)
                            OK_decay = True
                            dic['decay'][met][str(NI)][str(peak)]['OK_fit'] = OK_decay
                        except (Exception) as e:
                            print "Cannot fit expdecay for %s %s. Reason: %s"%(peak, peakname, e)
                            logger.info("Cannot fit expdecay for %s %s. Reason: %s"%(peak, peakname, e))
                            OK_decay = False
                            dic['decay'][met][str(NI)][str(peak)]['OK_fit'] = OK_decay
                    # Setting keys
                    offset = dic['offset'][i]
                    omega1 = dic['omega1'][i]
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['offset'] = offset
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['omega1'] = omega1
                    # Setting average chemical shifts
                    CS_H_fs = CS_H[fs:fe]
                    CS_H_fs_mean =CS_H_fs.mean()
                    CS_N_fs = CS_N[fs:fe]
                    CS_N_fs_mean = CS_N_fs.mean()
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['CS_H_mean'] = CS_H_fs_mean
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['CS_N_mean'] = CS_N_fs_mean
                    # Calculate NMR properties # KTE: residue_files.pl
                    OMEGA=(centerPPM-CS_N_fs_mean)*frq+offset
                    omegaEFF=sqrt(OMEGA**2+omega1**2)
                    if omega1/OMEGA > 0:
                        theta = 180/pi*abs(arctan(omega1/OMEGA))
                    else:
                        theta = 180- 180/pi*abs(arctan(omega1/OMEGA))
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['OMEGA'] = OMEGA
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['omegaEFF'] = omegaEFF
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['theta'] = theta #theta/tiltAngle
                    i+=1
    if multiprocess:
        print "Doing multiprocess fitting for decay rates"
        logger.info("Doing multiprocess fitting for decay rates")
        pool = Pool(processes=dic['Flags']['multiprocess'])
        multi_res = pool.map(multi_f_expdecay, multi_inp_arr)
        #from handythread import foreach
        #multi_res = foreach(multi_f_expdecay,multi_inp_arr,threads=dic['Flags']['multiprocess'],return_=True)
        for i in range(len(multi_sort_arr)):
            multi_sort = multi_sort_arr[i]
            met = multi_sort[0]; NI = multi_sort[1]
            peak = multi_sort[2]; fs = multi_sort[3]
            multi_sort[1]
            dic2 = multi_res[i]
            dic['decay'][met][str(NI)][str(peak)][str(fs)].update(dic2)
    print "Done with Getting decay rates. It took: %s"%(datetime.now()-startTime)
    logger.info("Done with Getting decay rates. It took: %s"%(datetime.now()-startTime))
    return()

def getrelax(dic,mets=False,NIstop=False):
    startTime = datetime.now()
    logger = logging.getLogger("TB.TB.getrelax")
    ncyc_arr = dic['ncyc']
    time_T2 = dic['time_T2']
    slicet = len(ncyc_arr)
    dic['relax'] = {}
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        Int = dic['Int'][met]
        filenr = dic['filenr'][met]
        dic['relax'][met] = {}
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            print "%s - Getting relaxation rates for NI=%s"%(met,NI)
            dic['relax'][met][str(NI)] = {}
            Pval_peaks = []
            peaks = dic['peakrange'][met]
            #peaks = ['57']
            for peak in peaks:
                dic['relax'][met][str(NI)][str(peak)] = {}
                peakname = Int[str(peak)]['resn']
                dic['relax'][met][str(NI)][str(peak)]['resn'] = peakname
                #print peakname
                MetInt_arr = []
                nu_arr = []
                i = 0
                for fs in range(filenr):
                    dic['relax'][met][str(NI)][str(peak)][str(fs)] = {}
                    FTInt = Int[str(peak)]['FTInt'][fs]
                    NIInt = Int[str(peak)]['NIInt'][str(NI)][fs]
                    MetInt = FTInt*NIInt
                    MetInt_arr.append(MetInt)
                    nu = ncyc_arr[i]/time_T2
                    nu_arr.append(nu)
                    i+=1
                nu_arr,MetInt_arr = zip(*sorted(zip(nu_arr, MetInt_arr)))
                nu_slice = next(x[0] for x in enumerate(nu_arr) if x[1] > 0.001)
                averageZero = average(MetInt_arr[:nu_slice])
                dic['relax'][met][str(NI)][str(peak)]['averageZero'] = averageZero
                #
                CS_N = Int[str(peak)]['CS_N']
                CS_H = Int[str(peak)]['CS_H']
                CS_H_arr = []
                CS_N_arr = []
                R2eff_arr = []
                nu_arr = []
                i = 0
                for fs in range(filenr):
                    FTInt = Int[str(peak)]['FTInt'][fs]
                    NIInt = Int[str(peak)]['NIInt'][str(NI)][fs]
                    MetInt = FTInt*NIInt
                    CS_H_arr.append(CS_H[fs])
                    dic['relax'][met][str(NI)][str(peak)][str(fs)]['CS_H'] = CS_H[fs]
                    CS_N_arr.append(CS_N[fs])
                    dic['relax'][met][str(NI)][str(peak)][str(fs)]['CS_N'] = CS_N[fs]
                    try:
                        R2eff = -1.0/time_T2*log(MetInt/averageZero)
                        R2eff_arr.append(R2eff)
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['R2eff'] = R2eff
                        OK_R2eff = True
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['OK_R2eff'] = OK_R2eff
                        #
                        nu = ncyc_arr[i]/time_T2
                        nu_arr.append(nu)
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['nu'] = nu
                    except (Exception) as e:
                        print "Cannot log calc R2eff for NI=%s Peak=%s %s. nu=%3.2f, MetInt=%3.2f, averageZero=%3.2f, fs=%s Reason: %s"%(NI, peak, peakname, nu, MetInt, averageZero, fs, e)
                        logger.info("Cannot log calc R2eff for NI=%s Peak=%s %s. nu=%3.2f, MetInt=%3.2f, averageZero=%3.2f, fs=%s Reason: %s"%(NI, peak, peakname, nu, MetInt, averageZero, fs, e))
                        OK_R2eff = False
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['OK_R2eff'] = OK_R2eff
                    i+=1
                dic['relax'][met][str(NI)][str(peak)]['CS_H_mean'] = mean(CS_H_arr)
                dic['relax'][met][str(NI)][str(peak)]['CS_N_mean'] = mean(CS_N_arr)
                nu_arr,R2eff_arr = zip(*sorted(zip(nu_arr, R2eff_arr)))
                nu_slice = next(x[0] for x in enumerate(nu_arr) if x[1] > 0.001)
                nu_arr_s = array(nu_arr[nu_slice:])
                R2eff_arr_s = array(R2eff_arr[nu_slice:])
                dic['relax'][met][str(NI)][str(peak)]['nu_arr'] = nu_arr
                dic['relax'][met][str(NI)][str(peak)]['nu_arr_s'] = nu_arr_s
                dic['relax'][met][str(NI)][str(peak)]['R2eff_arr'] = R2eff_arr
                dic['relax'][met][str(NI)][str(peak)]['R2eff_arr_s'] = R2eff_arr_s
                dic['relax'][met][str(NI)][str(peak)]['nu_slice'] = nu_slice
                # Calculate for R2 simple
                Fval, Fdist, Pval = False, False, False
                dic['relax'][met][str(NI)][str(peak)]['R2s'] = {}
                dic['relax'][met][str(NI)][str(peak)]['R2cpmg_slow'] = {}
                try:
                    par_R2s = lmfit.Parameters()
                    par_R2s.add('R2', value=dic['guess']['s_R2'], vary=True, min=0.0)
                    lmf_R2s = lmfit.minimize(f_R2s, par_R2s, args=(nu_arr_s, R2eff_arr_s),method='leastsq')
                    dic_R2s = unpack_f_R2s(par_R2s,nu_arr_s,R2eff_arr_s,lmf_R2s)
                    dic['relax'][met][str(NI)][str(peak)]['R2s'].update(dic_R2s)
                    OK_R2s = True
                    dic['relax'][met][str(NI)][str(peak)]['R2s']['OK_fit'] = OK_R2s
                except (Exception) as e:
                    print "Cannot fit R2s for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
                    logger.info("Cannot fit R2s for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
                    OK_R2s = False
                    dic['relax'][met][str(NI)][str(peak)]['R2s']['OK_fit'] = OK_R2s
                # Calculate for R2 slow
                try:
                    par_R2cpmg_slow = lmfit.Parameters()
                    par_R2cpmg_slow.add('R2', value=dic['guess']['s_R2'], vary=True, min=0.0)
                    par_R2cpmg_slow.add('Domega', value=dic['guess']['s_Domega'], vary=True, min=0.0)
                    par_R2cpmg_slow.add('ka', value=dic['guess']['s_ka'], vary=True, min=0.0)
                    lmf_R2cpmg_slow = lmfit.minimize(f_R2cpmg_slow, par_R2cpmg_slow, args=(nu_arr_s, R2eff_arr_s),method='leastsq')
                    dic_R2cpmg_slow = unpack_f_R2cpmg_slow(par_R2cpmg_slow,nu_arr_s,R2eff_arr_s,lmf_R2cpmg_slow)
                    dic['relax'][met][str(NI)][str(peak)]['R2cpmg_slow'].update(dic_R2cpmg_slow)
                    OK_R2cpmg_slow = True
                    dic['relax'][met][str(NI)][str(peak)]['R2cpmg_slow']['OK_fit'] = OK_R2cpmg_slow
                    #print lmf_R2cpmg_slow.success
                except (Exception) as e:
                    print "Cannot fit R2cpmg_slow for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
                    logger.info("Cannot fit R2cpmg_slow for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
                    OK_R2cpmg_slow = False
                    dic['relax'][met][str(NI)][str(peak)]['R2cpmg_slow']['OK_fit'] = OK_R2cpmg_slow
                if OK_R2s and OK_R2cpmg_slow:
                    Fval, Fdist, Pval = Ftest(dic_R2s['chisqr'],dic_R2s['NDF'],dic_R2cpmg_slow['chisqr'],dic_R2cpmg_slow['NDF'])
                dic['relax'][met][str(NI)][str(peak)]['Fval'] = Fval
                dic['relax'][met][str(NI)][str(peak)]['Fdist'] = Fdist
                dic['relax'][met][str(NI)][str(peak)]['Pval'] = Pval
                if Pval==False:
                    pass
                elif Pval!=False:
                    Pval_peaks.append(peak)

            dic['relax'][met][str(NI)]['Pval_peaks'] = Pval_peaks
            print "Following peak numbers passed the Ftest. %s/%s Met:%s NI=%s."%(len(Pval_peaks),len(peaks),met, NI)
            logger.info("Following peak numbers passed the Ftest. %s/%s Met:%s NI=%s."%(len(Pval_peaks),len(peaks),met, NI))
            print "Peaks:%s"%(Pval_peaks)
            logger.info("Peaks:%s"%(Pval_peaks))
    print "Done relaxation rates. It took: %s"%(datetime.now()-startTime)
    logger.info("Done relaxation rates. It took: %s"%(datetime.now()-startTime))
    return() #nu_arr_s,R2eff_arr_s,dic_R2cpmg_slow['par']

def getrelax_RC72(dics,mets=False,NIstop=False,sguess=False):
    startTime = datetime.now()
    logger = logging.getLogger("TB.TB.getrelax_RC72")
    rdic = collections.OrderedDict()
    for dic in dics:
        ncyc_arr = dic['ncyc']
        time_T2 = dic['time_T2']
        field = dic['field']
        rdic[str(field)] = {}
        slicet = len(ncyc_arr)
        dic['relax'] = {}
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            Int = dic['Int'][met]
            filenr = dic['filenr'][met]
            dic['relax'][met] = {}
            NIarr = dic['NIarr'][met]
            rdic[str(field)][met] = {}
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                print "%s - Getting relaxation rates for NI=%s"%(met,NI)
                dic['relax'][met][str(NI)] = {}
                peaks = dic['peakrange'][met]
                #peaks = ['57']
                rdic[str(field)][met][str(NI)] = {}
                for peak in peaks:
                    dic['relax'][met][str(NI)][str(peak)] = {}
                    peakname = Int[str(peak)]['resn']
                    if peakname == '?-?':
                        #print "Skipping peak: ?-?"
                        continue
                    dic['relax'][met][str(NI)][str(peak)]['resn'] = peakname
                    #print peakname
                    MetInt_arr = []
                    nu_arr = []
                    i = 0
                    for fs in range(filenr):
                        dic['relax'][met][str(NI)][str(peak)][str(fs)] = {}
                        FTInt = Int[str(peak)]['FTInt'][fs]
                        NIInt = Int[str(peak)]['NIInt'][str(NI)][fs]
                        MetInt = FTInt*NIInt
                        MetInt_arr.append(MetInt)
                        nu = ncyc_arr[i]/time_T2
                        nu_arr.append(nu)
                        i+=1
                    nu_arr,MetInt_arr = zip(*sorted(zip(nu_arr, MetInt_arr)))
                    nu_slice = next(x[0] for x in enumerate(nu_arr) if x[1] > 0.001)
                    averageZero = average(MetInt_arr[:nu_slice])
                    dic['relax'][met][str(NI)][str(peak)]['averageZero'] = averageZero
                    #
                    CS_N = Int[str(peak)]['CS_N']
                    CS_H = Int[str(peak)]['CS_H']
                    CS_H_arr = []
                    CS_N_arr = []
                    R2eff_arr = []
                    nu_arr = []
                    i = 0
                    for fs in range(filenr):
                        FTInt = Int[str(peak)]['FTInt'][fs]
                        NIInt = Int[str(peak)]['NIInt'][str(NI)][fs]
                        MetInt = FTInt*NIInt
                        CS_H_arr.append(CS_H[fs])
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['CS_H'] = CS_H[fs]
                        CS_N_arr.append(CS_N[fs])
                        dic['relax'][met][str(NI)][str(peak)][str(fs)]['CS_N'] = CS_N[fs]
                        try:
                            R2eff = -1.0/time_T2*log(MetInt/averageZero)
                            R2eff_arr.append(R2eff)
                            dic['relax'][met][str(NI)][str(peak)][str(fs)]['R2eff'] = R2eff
                            OK_R2eff = True
                            dic['relax'][met][str(NI)][str(peak)][str(fs)]['OK_R2eff'] = OK_R2eff
                            #
                            nu = ncyc_arr[i]/time_T2
                            nu_arr.append(nu)
                            dic['relax'][met][str(NI)][str(peak)][str(fs)]['nu'] = nu
                        except (Exception) as e:
                            print "Cannot log calc R2eff for NI=%s Peak=%s %s. nu=%3.2f, MetInt=%3.2f, averageZero=%3.2f, fs=%s Reason: %s"%(NI, peak, peakname, nu, MetInt, averageZero, fs, e)
                            logger.info("Cannot log calc R2eff for NI=%s Peak=%s %s. nu=%3.2f, MetInt=%3.2f, averageZero=%3.2f, fs=%s Reason: %s"%(NI, peak, peakname, nu, MetInt, averageZero, fs, e))
                            OK_R2eff = False
                            dic['relax'][met][str(NI)][str(peak)][str(fs)]['OK_R2eff'] = OK_R2eff
                        i+=1
                    dic['relax'][met][str(NI)][str(peak)]['CS_H_mean'] = mean(CS_H_arr)
                    dic['relax'][met][str(NI)][str(peak)]['CS_N_mean'] = mean(CS_N_arr)
                    nu_arr,R2eff_arr = zip(*sorted(zip(nu_arr, R2eff_arr)))
                    nu_slice = next(x[0] for x in enumerate(nu_arr) if x[1] > 0.001)
                    nu_arr_s = array(nu_arr[nu_slice:])
                    R2eff_arr_s = array(R2eff_arr[nu_slice:])
                    dic['relax'][met][str(NI)][str(peak)]['nu_arr'] = nu_arr
                    dic['relax'][met][str(NI)][str(peak)]['nu_arr_s'] = nu_arr_s
                    dic['relax'][met][str(NI)][str(peak)]['R2eff_arr'] = R2eff_arr
                    dic['relax'][met][str(NI)][str(peak)]['R2eff_arr_s'] = R2eff_arr_s
                    dic['relax'][met][str(NI)][str(peak)]['nu_slice'] = nu_slice

                    rdic[str(field)][met][str(NI)][peakname]={'peak':peak,'R2eff_arr_s':R2eff_arr_s,'nu_arr_s':nu_arr_s}
    # Collecting nu arrays
    fields = rdic.keys()
    mets = rdic[fields[0]]
    for met in mets:
        NIarr = rdic[fields[0]][met].keys()
        NIarr = map(int, NIarr)
        rdic[met] = {}
        for NI in NIarr:
            Pval_peaks = []
            peaknames = rdic[fields[0]][met][str(NI)].keys()
            rdic[met][str(NI)] = {}
            #peaknames = peaknames[:4]
            for peakname in peaknames:
                rdic[met][str(NI)][peakname] = {}
                R2eff_arr_s_list = []
                nu_arr_s_list = []
                par_R2cpmg_RC72 = lmfit.Parameters()
                for field in fields:
                    peak = rdic[str(field)][met][str(NI)][peakname]['peak']
                    R2eff_arr_s = rdic[str(field)][met][str(NI)][peakname]['R2eff_arr_s']
                    R2eff_arr_s_list.append(R2eff_arr_s)
                    nu_arr_s = rdic[str(field)][met][str(NI)][peakname]['nu_arr_s']
                    nu_arr_s_list.append(nu_arr_s)

                    par_R2cpmg_RC72.add('R2_%s'%field, value=sguess['s_R2'], vary=True, min=5.0, max=50.0)
                    #par_R2cpmg_RC72.add('pA_%s'%field, value=sguess['s_pA'], vary=True, min=0.5, max=1.0)
                    #par_R2cpmg_RC72.add('kEX_%s'%field, value=sguess['s_kEX'], vary=True, min=100.0, max=10000.0)
                par_R2cpmg_RC72.add('pA', value=sguess['s_pA'], vary=True, min=0.5, max=1.0)
                par_R2cpmg_RC72.add('kEX', value=sguess['s_kEX'], vary=True, min=100.0, max=10000.0)
                par_R2cpmg_RC72.add('w', value=sguess['s_w'], vary=True, min=0.0)
                try:
                    lmf_R2cpmg_RC72 = lmfit.minimize(f_R2cpmg_RC72_fields, par_R2cpmg_RC72, args=(fields, nu_arr_s_list, R2eff_arr_s_list),method='leastsq')
                    dic_R2cpmg_RC72 = unpack_f_R2cpmg_RC72_fields(par_R2cpmg_RC72, fields, nu_arr_s_list, R2eff_arr_s_list, lmf_R2cpmg_RC72)
                    #dic['relax'][met][str(NI)][str(peak)]['R2s'].update(dic_R2s)
                    #OK_R2cpmg_RC72 = True
                    #dic['relax'][met][str(NI)][str(peak)]['R2s']['OK_fit'] = OK_R2s


                fig = figure(figsize=(figsize, figsize/1.618))
                ax = fig.add_subplot(111)
                pars = par_R2cpmg_RC72
                print peakname,
                for key in pars:
                    print key, pars[key].value,
                print ""
                for i,field in enumerate(fields):
                    nu =nu_arr_s_list[i];R2eff=R2eff_arr_s_list[i]
                    ax.plot(nu,R2eff,".",label='%s %s %s %s'%(met, NI, peakname, field))

                    wcalc = float(field)/100.0*pars['w'].value
                    #par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX_%s'%field].value,'pA_v':pars['pA_%s'%field].value,'w_v':wcalc}
                    #par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX_%s'%field].value,'pA_v':pars['pA'].value,'w_v':wcalc}
                    par = {'R2_v':pars['R2_%s'%field].value,'kEX_v':pars['kEX'].value,'pA_v':pars['pA'].value,'w_v':wcalc}
                    Yfit = f_R2cpmg_RC72_calc(par,nu)
                    ax.plot(nu,Yfit,"-",color=ax.lines[-1].get_color(), label='R2_%s=%3.2f, kEX=%3.2f, pA=%1.2f, w=%3.2f'%(field,pars['R2_%s'%field].value,pars['kEX'].value,pars['pA'].value,wcalc))
                    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8})
                    ax.legend(loc='best', prop={'size':8})
                    savefig('img/%s_%s_%s.png'%(peakname,NI,met))


    #            # Calculate for R2 simple
    #            #Fval, Fdist, Pval = False, False, False
    #            #dic['relax'][met][str(NI)][str(peak)]['R2s'] = {}
    #            #dic['relax'][met][str(NI)][str(peak)]['R2cpmg_RC72'] = {}
    #            try:
    #                par_R2s = lmfit.Parameters()
    #                par_R2s.add('R2', value=dic['guess']['s_R2'], vary=True, min=0.0)
    #                lmf_R2s = lmfit.minimize(f_R2s, par_R2s, args=(nu_arr_s, R2eff_arr_s),method='leastsq')
    #                dic_R2s = unpack_f_R2s(par_R2s,nu_arr_s,R2eff_arr_s,lmf_R2s)
    #                dic['relax'][met][str(NI)][str(peak)]['R2s'].update(dic_R2s)
    #                OK_R2s = True
    #                dic['relax'][met][str(NI)][str(peak)]['R2s']['OK_fit'] = OK_R2s
    #            except (Exception) as e:
    #               print "Cannot fit R2s for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
    #                logger.info("Cannot fit R2s for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
    #                OK_R2s = False
    #                dic['relax'][met][str(NI)][str(peak)]['R2s']['OK_fit'] = OK_R2s
                # Calculate for R2 slow
                #try:
                #    par_R2cpmg_RC72 = lmfit.Parameters()
                #    par_R2cpmg_RC72.add('R2', value=dic['guess']['s_R2'], vary=True, min=5.0, max=50.0)
                #    par_R2cpmg_RC72.add('kEX', value=dic['guess']['s_kEX'], vary=True, min=500.0, max=10000.0)
                #    par_R2cpmg_RC72.add('w', value=dic['guess']['s_w'], vary=True, min=0.0)
                #    par_R2cpmg_RC72.add('pA', value=dic['guess']['s_pA'], vary=True, min=0.5, max=1.0)
                #    lmf_R2cpmg_RC72 = lmfit.minimize(f_R2cpmg_RC72, par_R2cpmg_RC72, args=(nu_arr_s_list[0], R2eff_arr_s_list[0]),method='leastsq')
                #    #dic_R2cpmg_RC72 = unpack_f_R2cpmg_RC72(par_R2cpmg_RC72,nu_arr_s,R2eff_arr_s,lmf_R2cpmg_RC72)
                #    #dic['relax'][met][str(NI)][str(peak)]['R2cpmg_RC72'].update(dic_R2cpmg_RC72)
                #    OK_R2cpmg_RC72 = True
                #    dic['relax'][met][str(NI)][str(peak)]['R2cpmg_RC72']['OK_fit'] = OK_R2cpmg_RC72
                #    print lmf_R2cpmg_RC72.success
                #except (Exception) as e:
                #    print "Cannot fit R2cpmg_RC72 for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
    #           #     logger.info("Cannot fit R2cpmg_RC72 for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
    #           #     OK_R2cpmg_RC72 = False
    #           #     dic['relax'][met][str(NI)][str(peak)]['R2cpmg_RC72']['OK_fit'] = OK_R2cpmg_RC72
                #fig = figure(figsize=(figsize, figsize/1.618))
                #ax = fig.add_subplot(111)
                #datY = f_R2cpmg_RC72(par_R2cpmg_RC72,nu_arr_s_list[0])
                #ax.plot(nu_arr_s_list[0],R2eff_arr_s_list[0],'.-')
                #ax.plot(nu_arr_s_list[0],datY,'.-')
                #show()
    #            #if OK_R2s and OK_R2cpmg_RC72:
    #            #    Fval, Fdist, Pval = Ftest(dic_R2s['chisqr'],dic_R2s['NDF'],dic_R2cpmg_RC72['chisqr'],dic_R2cpmg_RC72['NDF'])
    #            #dic['relax'][met][str(NI)][str(peak)]['Fval'] = Fval
    #            #dic['relax'][met][str(NI)][str(peak)]['Fdist'] = Fdist
    #            #dic['relax'][met][str(NI)][str(peak)]['Pval'] = Pval
    #            #if Pval==False:
    #            #    pass
    #            #elif Pval!=False:
    #            #    Pval_peaks.append(peak)
    #
    #        #dic['relax'][met][str(NI)]['Pval_peaks'] = Pval_peaks
    #        #print "Following peak numbers passed the Ftest. %s/%s Met:%s NI=%s."%(len(Pval_peaks),len(peaks),met, NI)
    #        #logger.info("Following peak numbers passed the Ftest. %s/%s Met:%s NI=%s."%(len(Pval_peaks),len(peaks),met, NI))
        #print "Peaks:%s"%(Pval_peaks)
        #logger.info("Peaks:%s"%(Pval_peaks))

    #show()
    print "Done relaxation rates. It took: %s"%(datetime.now()-startTime)
    #logger.info("Done relaxation rates. It took: %s"%(datetime.now()-startTime))
    return(rdic,par_R2cpmg_RC72,par_R2cpmg_RC72)

def getrates(dic,mets=False,NIstop=False):
    logger = logging.getLogger("TB.TB.getrates")
    startTime = datetime.now()
    print "Getting exchange rates"
    logger.info("Getting exchange rates")
    dic['rates'] = {}
    slicet = len(dic['time'])
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        dic['rates'][met] = {}
        filenr = dic['filenr'][met]
        Dec = dic['decay'][met]
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            #print "%s - Getting rates for NI=%s"%(met,NI)
            dic['rates'][met][str(NI)] = {}
            Pval_peaks = []
            for peak in dic['peakrange'][met]:
                dic['rates'][met][str(NI)][str(peak)] = {}
                tiltAngle_arr = []
                R1r_rates_arr = []
                R1r_err_arr = []
                peakname = Dec[str(NI)][str(peak)]['resn']
                dic['rates'][met][str(NI)][str(peak)]['resn'] = peakname
                omega1_arr = []
                for fs in range(0,filenr,slicet):
                ##for fs in range(0,5,slicet):
                    theta = Dec[str(NI)][str(peak)][str(fs)]['theta']
                    tiltAngle_arr.append(theta)
                    R1r_rates = Dec[str(NI)][str(peak)][str(fs)]['R1r_rates']
                    R1r_rates_arr.append(R1r_rates)
                    R1r_err = Dec[str(NI)][str(peak)][str(fs)]['R1r_err']
                    R1r_err_arr.append(R1r_err)
                    omega1 = Dec[str(NI)][str(peak)][str(fs)]['omega1']
                    omega1_arr.append(omega1)
                # Prepare data
                datX_f_R1r = array(tiltAngle_arr)
                datX_f_R1r_exch = [array(tiltAngle_arr), array(omega1_arr)]
                datY = array(R1r_rates_arr)
                f_sigma = array(R1r_err_arr) # The std. error on Y
                Fval, Fdist, Pval = False, False, False
                dic['rates'][met][str(NI)][str(peak)]['R1r'] = {}
                dic['rates'][met][str(NI)][str(peak)]['R1r_exch'] = {}
                # Calculate for R1r
                try:
                    par_R1r = lmfit.Parameters()
                    par_R1r.add('R1', value=dic['guess']['s_R1'], vary=True, min=0.0)
                    par_R1r.add('R2', value=dic['guess']['s_R2'], vary=True, min=0.0)
                    lmf_R1r = lmfit.minimize(f_R1r, par_R1r, args=(datX_f_R1r, datY, f_sigma),method='leastsq')
                    dic_R1r = unpack_f_R1r(par_R1r,datX_f_R1r,datY,f_sigma,lmf_R1r)
                    dic['rates'][met][str(NI)][str(peak)]['R1r'].update(dic_R1r)
                    OK_R1r = True
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['OK_fit'] = OK_R1r
                except (Exception) as e:
                    print "Cannot fit R1r for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
                    logger.info("Cannot fit R1r for NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
                    OK_R1r = False
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['OK_fit'] = OK_R1r
                # Calculate for R1r_exch
                try:
                    par_R1r_exch = lmfit.Parameters()
                    par_R1r_exch.add('R1', value=dic['guess']['s_R1'], vary=True, min=0.0)
                    par_R1r_exch.add('R2', value=dic['guess']['s_R2'], vary=True, min=0.0)
                    par_R1r_exch.add('kEX', value=dic['guess']['s_kEX'], vary=True, min=0.0)
                    par_R1r_exch.add('phi', value=dic['guess']['s_phi'], vary=True, min=0.0)
                    lmf_R1r_exch = lmfit.minimize(f_R1r_exch, par_R1r_exch, args=(datX_f_R1r_exch, datY, f_sigma),method='leastsq')
                    dic_R1r_exch = unpack_f_R1r_exch(par_R1r_exch,datX_f_R1r_exch,datY,f_sigma,lmf_R1r_exch)
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch'].update(dic_R1r_exch)
                    OK_R1r_exch = True
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['OK_fit'] = OK_R1r_exch
                except (Exception) as e:
                    print "Cannot fit R1r_exch NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e)
                    logger.info("Cannot fit R1r_exch NI=%s Peak=%s %s. Reason: %s"%(NI, peak, peakname, e))
                    OK_R1r_exch = False
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['OK_fit'] = OK_R1r_exch
                if OK_R1r and OK_R1r_exch:
                    Fval, Fdist, Pval = Ftest(dic_R1r['chisqr'],dic_R1r['NDF'],dic_R1r_exch['chisqr'],dic_R1r_exch['NDF'])
                dic['rates'][met][str(NI)][str(peak)]['Fval'] = Fval
                dic['rates'][met][str(NI)][str(peak)]['Fdist'] = Fdist
                dic['rates'][met][str(NI)][str(peak)]['Pval'] = Pval
                if Pval==False:
                    pass
                elif Pval!=False:
                    Pval_peaks.append(peak)
            dic['rates'][met][str(NI)]['Pval_peaks'] = Pval_peaks
            print "Following peak numbers passed the Ftest. Met:%s NI=%s."%(met, NI)
            logger.info("Following peak numbers passed the Ftest. Met:%s NI=%s."%(met, NI))
            print "Peaks:%s"%(Pval_peaks)
            logger.info("Peaks:%s"%(Pval_peaks))
    print "Done exchange rates. It took: %s"%(datetime.now()-startTime)
    logger.info("Done exchange rates. It took: %s"%(datetime.now()-startTime))
    return()

def getglobfit(dic,mets=False,peaklist=False,NIstop=False):
    dic['gfit'] = {}
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        dic['gfit'][met] = {}
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            if not peaklist:
                peaks = dic['rates'][met][str(NI)]['Pval_peaks']# dic['peakrange'][met]
                testPval = True
            else:
                peaks = peaklist
                testPval = False
            if len(peaks) > 0 :
                dic['gfit'][met][str(NI)] = {}
                tilt = []
                om1 = []
                R1rex = []
                R1rex_err = []
                sel_p = []
                P_arr = lmfit.Parameters(); P_arr.add('kEX', value=dic['guess']['g_kEX'], vary=True, min=0.0)
                for peak in peaks:
                    dc = dic['rates'][met][str(NI)][str(peak)]
                    peakname = dc['resn']
                    if testPval:
                        Pval = dic['rates'][met][str(NI)][str(peak)]['Pval']
                    else:
                        Pval = True
                    OK_fit = dc['R1r_exch']['OK_fit']
                    if Pval!=False and OK_fit:
                    # Get values for R1r_exch
                        dic['gfit'][met][str(NI)][str(peak)] = {}
                        X_Y_Sigma_Fit_exch = dc['R1r_exch']['X_Y_Sigma_Fit']
                        datX_f_R1r_exch = X_Y_Sigma_Fit_exch[0]
                        datY_f_R1r_exch = X_Y_Sigma_Fit_exch[1]
                        datY_f_R1r_exch_err = X_Y_Sigma_Fit_exch[2]
                        calcR1r_exch = X_Y_Sigma_Fit_exch[3]
                        par_R1r_exch = dc['R1r_exch']['par']

                        R1,R2,kEX,phi = par_R1r_exch['R1_v'],par_R1r_exch['R2_v'],par_R1r_exch['kEX_v'],par_R1r_exch['phi_v']
                        #print "p: R1=%3.2f, R2=%3.2f, kEX=%3.2f, phi=%3.2f"%(R1,R2,kEX,phi)
                        tiltAngle_arr_s = datX_f_R1r_exch[0]; omega1_arr_s = datX_f_R1r_exch[1]
                        tilt.append(tiltAngle_arr_s); om1.append(omega1_arr_s)
                        R1rex.append(datY_f_R1r_exch)
                        R1rex_err.append(datY_f_R1r_exch_err)
                        P_arr.add('R1%s'%peak, value=R1, vary=True, min=0.0)
                        P_arr.add('R2%s'%peak, value=R2, vary=True, min=0.0)
                        P_arr.add('phi%s'%peak, value=phi, vary=True, min=0.0)
                        sel_p.append(peak)
                        dic['gfit'][met][str(NI)][str(peak)]['resn'] = peakname
                dic['gfit'][met][str(NI)]['gfit_peaks'] = sel_p
                if len(sel_p) > 0 :
                    # Do global fit
                    startTime = datetime.now()
                    print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                    print "%s Global fit for NI=%s, %s Peaks:%s"%(met,NI,len(sel_p),sel_p)
                    lmf_R1r_exch_glob = lmfit.minimize(f_R1r_exch_global, P_arr, args=(sel_p,tilt,om1,R1rex,R1rex_err),method='leastsq')
                    # Unpack result into each peak
                    for i in range(len(sel_p)):
                        p = sel_p[i]
                        dc = dic['rates'][met][str(NI)][str(p)]
                        kEX_glob = P_arr['kEX']
                        R1 = P_arr['R1%s'%p]
                        R2 = P_arr['R2%s'%p]
                        phi = P_arr['phi%s'%p]
                        par_calc = lmfit.Parameters()
                        par_calc['kEX'] = kEX_glob
                        par_calc['R1'] = R1
                        par_calc['R2'] = R2
                        par_calc['phi'] = phi
                        X_Y_Sigma_Fit_exch = dc['R1r_exch']['X_Y_Sigma_Fit']
                        datX_f_R1r_exch = X_Y_Sigma_Fit_exch[0]
                        datY_f_R1r_exch = X_Y_Sigma_Fit_exch[1]
                        datY_f_R1r_exch_err = X_Y_Sigma_Fit_exch[2]
                        dic_R1r_exch_glob = unpack_f_R1r_exch(par_calc,datX_f_R1r_exch,datY_f_R1r_exch,datY_f_R1r_exch_err,lmf_R1r_exch_glob)
                        dic['gfit'][met][str(NI)][str(p)]['R1r_exch'] = {}
                        dic['gfit'][met][str(NI)][str(p)]['R1r_exch'].update(dic_R1r_exch_glob)
                    dic_glob = unpack_global(P_arr,sel_p,tilt,om1,R1rex,R1rex_err,lmf_R1r_exch_glob)
                    dic['gfit'][met][str(NI)].update(dic_glob)

                    print "Medthod=%s, NI=%s, kEX=%4.4f, chisqr-lmf=%4.4f, chisqr-calc=%4.4f"%(met,NI,kEX_glob.value,lmf_R1r_exch_glob.chisqr, dic_glob['chisqr'])
                    print "It took: %s"%(datetime.now()-startTime)
                else:
                    print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                    print "%s No Global fit for NI=%s, since len(sel_p):%s"%(met,NI,len(peaks))
            else:
                print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                print "%s No Global fit for NI=%s, since len(peaks):%s"%(met,NI,len(peaks))
    return()

def getglobfit_slow(dic,mets=False,peaklist=False,NIstop=False):
    dic['gfit'] = {}
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        dic['gfit'][met] = {}
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            if not peaklist:
                peaks = dic['relax'][met][str(NI)]['Pval_peaks']# dic['peakrange'][met]
                testPval = True
            else:
                peaks = peaklist
                testPval = False
            if len(peaks) > 0 :
                dic['gfit'][met][str(NI)] = {}
                nu_arr = []
                R2eff_arr = []
                sel_p = []
                P_arr = lmfit.Parameters(); P_arr.add('ka', value=dic['guess']['g_ka'], vary=True, min=0.0)
                for peak in peaks:
                    dc = dic['relax'][met][str(NI)][str(peak)]
                    peakname = dc['resn']
                    if testPval:
                        Pval = dic['relax'][met][str(NI)][str(peak)]['Pval']
                    else:
                        Pval = True
                    OK_fit = dc['R2cpmg_slow']['OK_fit']
                    if Pval!=False and OK_fit:
                    # Get values for R2cpmg_slow
                        dic['gfit'][met][str(NI)][str(peak)] = {}
                        X_Y_Fit = dc['R2cpmg_slow']['X_Y_Fit']
                        datX_nu_arr = X_Y_Fit[:,0]
                        datY_f_R2cpmg_slow = X_Y_Fit[:,1]
                        FitR2cpmg_slow = X_Y_Fit[:,2]
                        par_R2cpmg_slow = dc['R2cpmg_slow']['par']

                        R2,Domega= par_R2cpmg_slow['R2_v'],par_R2cpmg_slow['Domega_v']
                        nu_arr.append(datX_nu_arr)
                        R2eff_arr.append(datY_f_R2cpmg_slow)
                        P_arr.add('R2%s'%peak, value=R2, vary=True, min=0.0)
                        P_arr.add('Domega%s'%peak, value=Domega, vary=True, min=0.0)
                        sel_p.append(peak)
                        dic['gfit'][met][str(NI)][str(peak)]['resn'] = peakname
                dic['gfit'][met][str(NI)]['gfit_peaks'] = sel_p
                if len(sel_p) > 0 :
                    # Do global fit
                    startTime = datetime.now()
                    print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                    print "%s Global fit for NI=%s, %s Peaks:%s"%(met,NI,len(sel_p),sel_p)
                    lmf_R2cpmg_slow_glob = lmfit.minimize(f_R2cpmg_slow_global, P_arr, args=(sel_p,nu_arr,R2eff_arr),method='leastsq')
                    # Unpack result into each peak
                    for i in range(len(sel_p)):
                        p = sel_p[i]
                        dc = dic['relax'][met][str(NI)][str(p)]
                        ka_glob = P_arr['ka']
                        R2 = P_arr['R2%s'%p]
                        Domega = P_arr['Domega%s'%p]
                        par_calc = lmfit.Parameters()
                        par_calc['ka'] = ka_glob
                        par_calc['R2'] = R2
                        par_calc['Domega'] = Domega
                        X_Y_Fit = dc['R2cpmg_slow']['X_Y_Fit']
                        datX_nu_arr = X_Y_Fit[:,0]
                        datY_f_R2cpmg_slow = X_Y_Fit[:,1]
                        FitR2cpmg_slow = X_Y_Fit[:,2]
                        #
                        dic_R2cpmg_slow_glob = unpack_f_R2cpmg_slow(par_calc,datX_nu_arr,datY_f_R2cpmg_slow,lmf_R2cpmg_slow_glob)
                        dic['gfit'][met][str(NI)][str(p)]['R2cpmg_slow'] = {}
                        dic['gfit'][met][str(NI)][str(p)]['R2cpmg_slow'].update(dic_R2cpmg_slow_glob)
                    dic_glob = unpack_global_slow(P_arr,sel_p,nu_arr,R2eff_arr,None,lmf=lmf_R2cpmg_slow_glob)
                    dic['gfit'][met][str(NI)].update(dic_glob)

                    print "Medthod=%s, NI=%s, ka=%4.4f, chisqr-lmf=%4.4f, chisqr-calc=%4.4f"%(met,NI,ka_glob.value,lmf_R2cpmg_slow_glob.chisqr, dic_glob['chisqr'])
                    print "It took: %s"%(datetime.now()-startTime)
                else:
                    print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                    print "%s No Global fit for NI=%s, since len(sel_p):%s"%(met,NI,len(peaks))
            else:
                print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
                print "%s No Global fit for NI=%s, since len(peaks):%s"%(met,NI,len(peaks))
    return()

##################################

def del_par_props(dic,pars,mets=False,NIstop=False): #TB.del_chisq_prop(BBL,BBL['qMDDmet'])
    #TB.del_par_props(BBL,['phi','R1','R2','phi_v','R1_v','R2_v'],BBL['qMDDmet'])
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
            if len(peaks) > 0 :
                try:
                    del(dic['gfit'][met][str(NI)]['chisqr_prop'])
                except (KeyError) as e:
                    print "Could not delete key:%s"%e
                for par in pars:
                    try:
                        del(dic['gfit'][met][str(NI)]['%s_prop'%par])
                    except (KeyError) as e:
                        print "Could not delete key:%s"%e
                    try:
                        del(dic['gfit'][met][str(NI)]['par']['%s_prop'%par])
                    except (KeyError) as e:
                        print "Could not delete key:%s"%e
    return()

def get_glob_props(dic,pars,mets=False,NIstop=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    getglob_chisqr_prop(dic,mets,NIstop,gkey,pkey)
    for par in pars:
        getglob_par_prop(dic,par,mets,NIstop,gkey,pkey)
    return()

def getglob_chisqr_prop(dic,mets=False,NIstop=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
            if len(peaks) > 0 :
                NIlist = []
                proplist = []
                for peak in peaks:
                    gf = dic['gfit'][met][str(NI)][str(peak)]['%s'%pkey]
                    sf = dic['%s'%gkey][met][str(NI)][str(peak)]['%s'%pkey]
                    chisqr_glob = gf['chisqr']
                    chisqr_sing = sf['chisqr']
                    prop = chisqr_glob/chisqr_sing
                    NIlist.append(float(NI))
                    proplist.append(prop)
                NI_prop = array([NIlist,proplist]).T
                dic['gfit'][met][str(NI)]['chisqr_prop'] = NI_prop
    return()

def getglob_par_prop(dic,par,mets=False,NIstop=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
            if len(peaks) > 0 :
                NIlist = []
                proplist = []
                for peak in peaks:
                    gf = dic['gfit'][met][str(NI)][str(peak)]['%s'%pkey]
                    sf = dic['%s'%gkey][met][str(NI)][str(peak)]['%s'%pkey]
                    par_glob = gf['par']['%s_v'%par]
                    par_sing = sf['par']['%s_v'%par]
                    prop = par_glob/par_sing
                    NIlist.append(float(NI))
                    proplist.append(prop)
                NI_prop = array([NIlist,proplist]).T
                dic['gfit'][met][str(NI)]['par']['%s_prop'%par] = NI_prop
    return()

def get_glob_pearsons(dic,pars,mets=False,NIstop=False,Ini=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    getglob_chisqr_pearson(dic,mets,NIstop,Ini,gkey,pkey)
    for par in pars:
        a,b,c = getglob_par_pearson(dic,par,mets,NIstop,Ini,gkey,pkey)
    return(a,b,c)

def getglob_chisqr_pearson(dic,mets=False,NIstop=False,Ini=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
            try:
                test = dic['gfit'][met][str(NI)]['par']['Corr']
            except (KeyError) as e:
                dic['gfit'][met][str(NI)]['par']['Corr'] = {}
            if len(peaks) > 0:
                if not Ini: NIsf = NI; t = ''
                else: NIsf = max(NIarr); t = '_ini'
                NIlist = []
                chisqr_sing_list = []
                chisqr_glob_list = []
                for peak in peaks:
                    gf = dic['gfit'][met][str(NI)][str(peak)]['%s'%pkey]
                    if not Ini:
                        sf = dic['%s'%gkey][met][str(NIsf)][str(peak)]['%s'%pkey]
                    else:
                        sf = dic['gfit'][met][str(NIsf)][str(peak)]['%s'%pkey]
                    chisqr_sing = sf['chisqr']
                    chisqr_glob = gf['chisqr']
                    NIlist.append(float(NI))
                    chisqr_sing_list.append(chisqr_sing)
                    chisqr_glob_list.append(chisqr_glob)
                NI_pear = array([NIlist,chisqr_sing_list,chisqr_glob_list]).T
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_data%s'%(t)] = NI_pear
                Pearson_Corr_Coeff,Pearson_Corr_Coeff_tailed_p_value = scipy.stats.pearsonr(NI_pear[:,1], NI_pear[:,2])
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_Pearson_Corr_Coeff%s'%(t)] = Pearson_Corr_Coeff
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_Pearson_Corr_Coeff_tailed_p_value%s'%(t)] = Pearson_Corr_Coeff_tailed_p_value
                lin_slope, lin_inter, lin_r_value, lin_p_value, lin_std_err = scipy.stats.linregress(NI_pear[:,1],NI_pear[:,2])
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_lin_slope%s'%(t)] = lin_slope
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_lin_inter%s'%(t)] = lin_inter
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_lin_r_value%s'%(t)] = lin_r_value
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_lin_p_value%s'%(t)] = lin_p_value
                dic['gfit'][met][str(NI)]['par']['Corr']['chisqr_lin_std_err%s'%(t)] = lin_std_err
    return()

def getglob_par_pearson(dic,par,mets=False,NIstop=False,Ini=False,gkey='rates',pkey='R1r_exch'):
    if not mets: mets = dic['qMDDmet'][0]
    for met in mets:
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
            try:
                test = dic['gfit'][met][str(NI)]['par']['Corr']
            except (KeyError) as e:
                dic['gfit'][met][str(NI)]['par']['Corr'] = {}
            if len(peaks) > 0 :
                if not Ini: NIsf = NI; t = ''
                else: NIsf = max(NIarr); t = '_ini'
                NIlist = []
                par_sing_list = []
                par_sing_list_e = []
                par_glob_list = []
                par_glob_list_e = []
                peaknames = []
                for peak in peaks:
                    resn = dic['Int'][met][str(peak)]['resn']
                    gf = dic['gfit'][met][str(NI)][str(peak)]['%s'%pkey]
                    if not Ini:
                        sf = dic['%s'%gkey][met][str(NIsf)][str(peak)]['%s'%pkey]
                    else:
                        sf = dic['gfit'][met][str(NIsf)][str(peak)]['%s'%pkey]
                    par_sing = sf['par']['%s_v'%par]
                    par_sing_e = sf['par']['%s_e'%par]
                    par_glob = gf['par']['%s_v'%par]
                    par_glob_e = gf['par']['%s_e'%par]
                    NIlist.append(float(NI))
                    par_sing_list.append(par_sing)
                    par_sing_list_e.append(par_sing_e)
                    par_glob_list.append(par_glob)
                    par_glob_list_e.append(par_glob_e)
                    peaknames.append('%s %s'%(peak,resn))
                NI_pear = array([NIlist,par_sing_list,par_sing_list_e,par_glob_list,par_glob_list_e]).T
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_data%s'%(par,t)] = NI_pear
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_resn%s'%(par,t)] = peaknames
                Pearson_Corr_Coeff,Pearson_Corr_Coeff_tailed_p_value = scipy.stats.pearsonr(NI_pear[:,1], NI_pear[:,3])
                print NI, Pearson_Corr_Coeff
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_Pearson_Corr_Coeff%s'%(par,t)] = Pearson_Corr_Coeff
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_Pearson_Corr_Coeff_tailed_p_value%s'%(par,t)] = Pearson_Corr_Coeff_tailed_p_value
                lin_slope, lin_inter, lin_r_value, lin_p_value, lin_std_err = scipy.stats.linregress(NI_pear[:,1],NI_pear[:,3])
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_slope%s'%(par,t)] = lin_slope
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_inter%s'%(par,t)] = lin_inter
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_r_value%s'%(par,t)] = lin_r_value
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_p_value%s'%(par,t)] = lin_p_value
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_std_err%s'%(par,t)] = lin_std_err
                single_x = NI_pear[:,1][:,np.newaxis]
                glob_x = NI_pear[:,3][:,np.newaxis]
                slope, _, _, _ = np.linalg.lstsq(single_x, glob_x)
                slope = float(slope)
                dic['gfit'][met][str(NI)]['par']['Corr']['%s_slope%s'%(par,t)] = slope
    return(NI_pear[:,1],NI_pear[:,3],NI)

################ Plot functions ###########################
def plotstats(dics,mets):
    for dic in dics:
        desc_name = dic['desc']['name']
        fig = figure('Stats_%s'%(desc_name),figsize=(figsize, figsize/1.618))
        fig.suptitle('Tada',fontsize=titfont)
        ax = fig.add_subplot(111)
        for met in mets:
            data = dic['stats']['data'][met]
            NIs = data[:,0]
            datX = 1 - NIs/dic['NImax'][met]
            datY = data[:,38]
            dc = dic['desc']
            ax.plot(datX,datY,'.-', label='%s : %s'%(met, dc))
        ax.set_ylim([0,0.1])
        ax.set_ylabel('Residual stdev',fontsize=labfont)
        ax.set_xlabel("Sparseness",fontsize=labfont)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
        ax.grid('on')
        print "Write: show() , to see graphs"
        #show()
    return()

def plotdecays(dics,mets=False,peaks=False,NIa=False,fss=[0]):
    #TB.plotdecays([BBL],['coMDD'])
    #TB.plotdecays([BBL],['coMDD','CS'],fss=range(0,35,5))
    for dic in dics:
        desc_name = dic['desc']['name']
        fig = figure('Decay_%s'%(desc_name),figsize=(figsize, figsize/1.618))
        ax = fig.add_subplot(111)
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            if not NIa: NIarr = list(dic['NIarr'][met][:2]) # [dic['NImax'][met]] #
            else: NIarr = NIa
            for NI in NIarr:
                if not peaks: peaks = ['1']
                #print peaks, type(peaks)
                for peak in peaks:
                    peakname = dic['decay'][met][str(NI)][str(peak)]['resn']
                    for fs in fss:
                        X_Y_Fit = dic['decay'][met][str(NI)][str(peak)][str(fs)]['X_Y_Fit']
                        datX = X_Y_Fit[:,0]
                        datY = X_Y_Fit[:,1]
                        fitY = X_Y_Fit[:,2]
                        par = dic['decay'][met][str(NI)][str(peak)][str(fs)]['par']
                        a = par['amp_v']
                        a_e = par['amp_e']
                        b = par['decay_v']
                        b_e = par['decay_e']
                        datXs = sort(datX)
                        fitXlin = linspace(datXs[0],datXs[-1],100)
                        fitYlin = f_expdecay_calc(par,fitXlin)
                        ax.plot(datX,datY,".",label='%s %s %s %s. fit [a*exp(-b*x]: a=%3.2f +-%3.2f, a=%3.2f +-%3.2f '%(met, peakname, NI, fs, a,a_e,b,b_e))
                        ax.plot(fitXlin,fitYlin,"-",color=ax.lines[-1].get_color())
                        ax.plot(datX,fitY,"o",mfc='none',mec=ax.lines[-1].get_color())
        title('Decay')
        ax.set_xlabel('Time')
        ax.set_ylabel('Normalized intensity')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
        ax.grid('on')
        print "Write: show() , to see graphs"
        #show()
    return()

def plotrelaxation(dics,mets=False,peaks=False,NIa=False):
    for dic in dics:
        desc_name = dic['desc']['name']
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            if not peaks: peaks = dic['peakrange'][met]
            if not NIa: NIarr = list(dic['NIarr'][met][:1]) #
            else: NIarr = NIa
            for NI in NIarr:
                figR2s = figure('R2_%s_%s'%(NI,desc_name),figsize=(figsize, figsize/1.618))
                ax = figR2s.add_subplot(111)
                for peak in peaks:
                    dc = dic['relax'][met][str(NI)][str(peak)]
                    peakname = dc['resn']
                    Pval = dic['relax'][met][str(NI)][str(peak)]['Pval']
                    if Pval==False:
                        # Get values for R2
                        X_Y_Fit = dc['R2s']['X_Y_Fit']
                        #print peakname
                        datX_f_R2s = X_Y_Fit[:,0]
                        datY_f_R2s = X_Y_Fit[:,1]
                        calcR2s = X_Y_Fit[:,2]
                        par_R2s = dc['R2s']['par']
                        datXs_f_R2s_lin = linspace(datX_f_R2s[0],datX_f_R2s[-1],500)
                        datYs_f_R2s_lin = f_R2s_calc(par_R2s,datXs_f_R2s_lin)
                        ax.errorbar(datX_f_R2s, datY_f_R2s,fmt='.',label='%s %s-%s %s. R2=%3.2f+-%3.2f '%(met, peak, peakname, NI, par_R2s['R2_v'], par_R2s['R2_e']))
                        ax.plot(datX_f_R2s,datY_f_R2s,"-")
                        ax.plot(datXs_f_R2s_lin,datYs_f_R2s_lin,"-",color=ax.lines[-1].get_color())
                        ax.plot(datX_f_R2s,calcR2s,"o",mfc='none',mec=ax.lines[-1].get_color())#
                    if Pval!=False:
                        figR2cpmg_slow = figure('R2cpmg_slow_%s_%s_%s'%(NI,peak,desc_name),figsize=(figsize, figsize/1.618))
                        bx = figR2cpmg_slow.add_subplot(111)
                        # Get values for R2cpmg_slow
                        X_Y_Fit_exch = dc['R2cpmg_slow']['X_Y_Fit']
                        datX_f_R2cpmg_slow = X_Y_Fit_exch[:,0]
                        datY_f_R2cpmg_slow = X_Y_Fit_exch[:,1]
                        calcR2cpmg_slow = X_Y_Fit_exch[:,2]
                        par_R2cpmg_slow = dc['R2cpmg_slow']['par']
                        datXs_f_R2cpmg_slow_lin = linspace(datX_f_R2cpmg_slow[0],datX_f_R2cpmg_slow[-1],500)
                        datYs_f_R2cpmg_slow_lin = f_R2cpmg_slow_calc(par_R2cpmg_slow,datXs_f_R2cpmg_slow_lin )
                        bx.errorbar(datX_f_R2cpmg_slow, datY_f_R2cpmg_slow,fmt='.',
                        label='%s %s-%s %s. R2=%3.2f ka=%3.2f Domega=%3.2f'%(met, peak, peakname, NI,par_R2cpmg_slow['R2_v'],par_R2cpmg_slow['ka_v'],par_R2cpmg_slow['Domega_v']))
                        #bx.plot(datX_f_R2cpmg_slow,datY_f_R2cpmg_slow,"-")
                        bx.plot(datXs_f_R2cpmg_slow_lin,datYs_f_R2cpmg_slow_lin,"-",color=bx.lines[-1].get_color())
                        bx.plot(datX_f_R2cpmg_slow,calcR2cpmg_slow,"o",mfc='none',mec=bx.lines[-1].get_color())#
                        #R2cpmg_slow graph
                        bx.set_title('R2cpmg slow fitting')
                        #bx.set_xlabel('Time')
                        bx.set_xlabel(r'$\nu$ CPMG')
                        bx.set_ylabel('R2cpmg slow coef')
                        box = bx.get_position()
                        bx.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
                        bx.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                        bx.grid('on')
                #R1r graph
                ax.set_title('R2 fitting')
                #ax.set_xlabel('Time')
                ax.set_xlabel(r'$\nu$ CPMG')
                ax.set_ylabel('R2 coef')
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                ax.grid('on')
    print "Write: show() , to see graphs"
#    #show()
    return()

def plotrates(dics,mets=False,peaks=False,NIa=False,glob=False):
    for dic in dics:
        desc_name = dic['desc']['name']
        unique_omega1 = sort(f5(dic['omega1']))
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            if not peaks: peaks = dic['peakrange'][met]
            if not NIa: NIarr = list(dic['NIarr'][met][:2]) #
            else: NIarr = NIa
            print NIarr
            for NI in NIarr:
                print NI
                figR1r = figure('R1r_%s_%s'%(NI,desc_name),figsize=(figsize, figsize/1.618))
                ax = figR1r.add_subplot(111)
                for peak in peaks:
                    dc = dic['rates'][met][str(NI)][str(peak)]
                    peakname = dc['resn']
                    Pval = dic['rates'][met][str(NI)][str(peak)]['Pval']
                    if Pval==False:
                        # Get values for R1r
                        X_Y_Sigma_Fit = dc['R1r']['X_Y_Sigma_Fit']
                        datX_f_R1r = X_Y_Sigma_Fit[:,0]
                        datY_f_R1r = X_Y_Sigma_Fit[:,1]
                        datY_f_R1r_err = X_Y_Sigma_Fit[:,2]
                        calcR1r = X_Y_Sigma_Fit[:,3]
                        par_R1r = dc['R1r']['par']
                        datXs_f_R1r = sort(datX_f_R1r)
                        datYs_f_R1r = f_R1r_calc(par_R1r,datXs_f_R1r)
                        datXs_f_R1r_lin = linspace(datXs_f_R1r[0],datXs_f_R1r[-1],100)
                        datYs_f_R1r_lin = f_R1r_calc(par_R1r,datXs_f_R1r_lin)
                        ax.errorbar(datX_f_R1r, datY_f_R1r,yerr=datY_f_R1r_err,fmt='.',label='%s %s-%s %s. R1=%3.2f+-%3.2f '%(met, peak, peakname, NI, par_R1r['R1_v'], par_R1r['R1_e']))
                        #ax.plot(datXs_f_R1r,datYs_f_R1r,"-",color=ax.lines[-1].get_color())
                        ax.plot(datXs_f_R1r_lin,datYs_f_R1r_lin,"-",color=ax.lines[-1].get_color())
                        ax.plot(datX_f_R1r,calcR1r,"o",mfc='none',mec=ax.lines[-1].get_color())#
                    if Pval!=False:
                        figR1r_exch = figure('R1r_exch_%s_%s_%s'%(NI,peak,desc_name),figsize=(figsize, figsize/1.618))
                        bx = figR1r_exch.add_subplot(111)
                        # Get values for R1r_exch
                        if glob:
                            X_Y_Sigma_Fit_exch = dic['gfit'][met][str(NI)][str(peak)]['R1r_exch']['X_Y_Sigma_Fit']
                            print "plotting global fit"
                        else:
                            X_Y_Sigma_Fit_exch = dc['R1r_exch']['X_Y_Sigma_Fit']
                        datX_f_R1r_exch = X_Y_Sigma_Fit_exch[0]
                        datY_f_R1r_exch = X_Y_Sigma_Fit_exch[1]
                        datY_f_R1r_exch_err = X_Y_Sigma_Fit_exch[2]
                        calcR1r_exch = X_Y_Sigma_Fit_exch[3]
                        if glob:
                            par_R1r_exch = dic['gfit'][met][str(NI)][str(peak)]['R1r_exch']['par']
                        else:
                            par_R1r_exch = dc['R1r_exch']['par']
                        print par_R1r_exch
                        tiltAngle_arr_s, omega1_arr_s = zip(*sorted(zip(datX_f_R1r_exch[0], datX_f_R1r_exch[1])))
                        datXs_f_R1r_exch = [array(tiltAngle_arr_s), array(omega1_arr_s)]
                        datYs_f_R1r_exch = f_R1r_exch_calc(par_R1r_exch,datXs_f_R1r_exch)

                        omega_plots = []
                        for om1 in unique_omega1:
                            tiltAngle_arr_s_lin = linspace(tiltAngle_arr_s[0],tiltAngle_arr_s[-1],100)
                            om1_arr = ones(len(tiltAngle_arr_s_lin))*om1
                            datXs_f_R1r_exch_lin = [array(tiltAngle_arr_s_lin), array(om1_arr)]
                            datYs_f_R1r_exch_lin = f_R1r_exch_calc(par_R1r_exch,datXs_f_R1r_exch_lin)
                            omega_plots.append([datXs_f_R1r_exch_lin,datYs_f_R1r_exch_lin])
                        ecolor_arr=get_ecolor(datX_f_R1r_exch[1],dic['omega1_col'])
                        bx.errorbar(datX_f_R1r_exch[0], datY_f_R1r_exch, yerr=datY_f_R1r_exch_err,fmt=None,ecolor='k',zorder=1,label='%s %s %s'%(met, peakname, NI)) #  fit f_R1r_exch: %s, dc['R1r_exch']['p']
                        bx.scatter(datX_f_R1r_exch[0],datY_f_R1r_exch,marker='s',c=ecolor_arr,edgecolors='None',zorder=4)
                        for omp in omega_plots:
                            colval = dic['omega1_col'][str(omp[0][1][0])]
                            bx.plot(omp[0][0],omp[1],"-", color=colval, zorder=2) # ,color=bx.lines[-1].get_color()
                        bx.scatter(datX_f_R1r_exch[0],calcR1r_exch,marker='x',c=ecolor_arr,zorder=3)#,color=bx.lines[-1].get_color()
                        #R1r_exch graph
                        bx.set_title('R1_exch fitting')
                        bx.set_xlabel('Tilt angle')
                        bx.set_ylabel('R1r_coef')
                        box = bx.get_position()
                        bx.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
                        bx.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                        bx.grid('on')
                #R1r graph
                ax.set_title('R1 fitting')
                ax.set_xlabel('Tilt angle')
                ax.set_ylabel('R1r_coef')
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                ax.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()

def plot_globalpar(dics,mets=False,NIstop=False,globalpar='kEX',gkey='rates',ay=False,by=False):
    for dic in dics:
        desc_name = dic['desc']['name']
        fig = figure(desc_name,figsize=(figsize, figsize/1.618))
        ax = fig.add_subplot(211)
        bx = fig.add_subplot(212)
        if not mets: mets = dic['qMDDmet'][0]
        imets = 0
        for met in mets:
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            globalpar_Y = []
            globalpar_Y_e= []
            chisqr_Y = []
            chisqr_prop_list = []
            chisqr_prop_peak_list_NI = array([])
            chisqr_prop_peak_list = array([])
            peak_plot_list = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]

                globalpar_v = gf['par']['%s_v'%globalpar]
                globalpar_e = gf['par']['%s_e'%globalpar]
                chisqr = gf['chisqr']
                gfit_peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
                Pval_peaks = dic['%s'%gkey][met][str(NI)]['Pval_peaks']
                gfit_peaks_in_Pval = [x for x in gfit_peaks if x in Pval_peaks]
                NI_X.append(NI)
                globalpar_Y.append(globalpar_v)
                globalpar_Y_e.append(globalpar_e)
                chisqr_Y.append(chisqr)
            NI_X = 1.0-array(NI_X)/float(NImax)
            globalpar_Y = array(globalpar_Y)
            globalpar_Y_e= array(globalpar_Y_e)
            chisqr_Y = array(chisqr_Y)
            ax.errorbar(NI_X, globalpar_Y,yerr=globalpar_Y_e,fmt='.-',label='%s'%(met))
            bx.plot(NI_X,chisqr_Y,'-',label='%s'%(met))
    ax.plot(NI_X,ones(len(NI_X))*globalpar_Y[0],c='k')
    ax.annotate('%s=%3.1f'%(globalpar,globalpar_Y[0]), xy=(NI_X[0],globalpar_Y[0]), xycoords='data', xytext=(0.5*NI_X[1],globalpar_Y[0]+0.05*globalpar_Y[0]), textcoords='data') #'axes fraction'
    bx.plot(NI_X,ones(len(NI_X))*chisqr_Y[0],c='k')
    ax.set_title('Global %s fitting'%globalpar)
    #
    ax.set_ylabel('%s'%globalpar)
    if ay: ax.set_ylim(ay[0],ay[1])
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
    ax.legend(loc='lower left')
    ax.grid('on')
    #
    bx.set_ylabel('chisqr')
    if by: bx.set_ylim(by[0],by[1])
    bx.set_xlabel('Sparseness')
    #box = bx.get_position()
    #bx.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    #bx.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
    bx.legend(loc='best')
    bx.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()

def plot_glob_props(dic,pars,mets=False,NIstop=False,gkey='rates'):
    if not mets: mets = dic['qMDDmet'][0]
    plot_chisqr_glob_prop(dic,mets,NIstop,gkey)
    for par in pars:
        plot_par_glob_prop(dic,par,mets,NIstop,gkey)
    return()

def plot_chisqr_glob_prop(dics,mets=False,NIstop=False,gkey='rates'):
    for dic in dics:
        desc_name = dic['desc']['name']
        fig = figure('chisqr_%s'%desc_name,figsize=(figsize, figsize/1.618))
        if not mets: mets = dic['qMDDmet'][0]
        imets = 0
        for met in mets:
            ax = fig.add_subplot('%s1%s'%(len(mets),imets+1))
            ax2 = ax.twinx()
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            chisqr_prop_list = []
            chisqr_prop_peak_list_NI = array([])
            chisqr_prop_peak_list = array([])
            peak_plot_list = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]
                gfit_peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
                Pval_peaks = dic['%s'%gkey][met][str(NI)]['Pval_peaks']
                gfit_peaks_in_Pval = [x for x in gfit_peaks if x in Pval_peaks]
                NI_X.append(NI)
                #
                chisqr_prop_dic = gf['chisqr_prop']
                chisqr_prop_NI = 1.0-chisqr_prop_dic[:,0]/float(NImax)
                chisqr_prop = chisqr_prop_dic[:,1]
                chisqr_mean = np.mean(chisqr_prop)
                chisqr_min = chisqr_prop.min()
                chisqr_max = chisqr_prop.max()
                chisqr_rmsd = np.sqrt(np.mean((chisqr_prop-chisqr_mean)**2))
                chisqr_prop_list.append([chisqr_mean,chisqr_min,chisqr_max,chisqr_rmsd])
                chisqr_prop_peak_list_NI = np.concatenate((chisqr_prop_peak_list_NI, chisqr_prop_NI))
                chisqr_prop_peak_list = np.concatenate((chisqr_prop_peak_list, chisqr_prop))
                peak_plot_list.append([len(gfit_peaks),len(Pval_peaks),len(gfit_peaks_in_Pval)])

            NI_X = 1.0-array(NI_X)/float(NImax)
            #
            #ax.plot(chisqr_prop_peak_list_NI,chisqr_prop_peak_list,'.',c='b')
            peak_plot_arr = array(peak_plot_list).T
            ax2.plot(NI_X, peak_plot_arr[0],'o',c='b',label='Peaks used for global fit')
            ax2.plot(NI_X, peak_plot_arr[1],'o',c='r',label='Ftest peaks')
            ax2.plot(NI_X, peak_plot_arr[2],'s',markersize=8,c='g',alpha=0.5,label='Fitted peaks in Ftest')
            #
            chisqr_prop_arr = array(chisqr_prop_list).T
            ax.errorbar(NI_X,chisqr_prop_arr[0],yerr=chisqr_prop_arr[3],marker='s',c='b',label='%s \nChisqr proportion \nglobal fit/single fit'%(met)) #ax.lines[-1].get_color()
            ax.plot(NI_X,chisqr_prop_arr[1],c='b',alpha=0.1)
            ax.plot(NI_X,chisqr_prop_arr[2],c='b',alpha=0.1)
            imets+=1
            #
            ax.set_ylabel('chisqr prop')
            ax.set_ylim(0,10)
            ax.set_xlabel('Sparseness')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.95),prop={'size':8}) # Put a legend to the right of the current axis
            #ax.legend(loc='best')
            ax.grid('on')
            #
            peaks = dic['peakrange'][met]
            ax2.set_ylabel('Number of peaks')
            ax2.set_ylim(0,len(peaks))
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.65),prop={'size':8}) # Put a legend to the right of the current axis
            #ax2.legend(loc='best')
            ax2.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()

def plot_par_glob_prop(dics,par,mets=False,NIstop=False,gkey='rates'):
    for dic in dics:
        desc_name = dic['desc']['name']
        fig = figure('%s_%s'%(par,desc_name),figsize=(figsize, figsize/1.618))
        if not mets: mets = dic['qMDDmet'][0]
        imets = 0
        for met in mets:
            ax = fig.add_subplot('%s1%s'%(len(mets),imets+1))
            ax2 = ax.twinx()
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            par_prop_list = []
            par_prop_peak_list_NI = array([])
            par_prop_peak_list = array([])
            peak_plot_list = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]
                gfit_peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
                Pval_peaks = dic['%s'%gkey][met][str(NI)]['Pval_peaks']
                gfit_peaks_in_Pval = [x for x in gfit_peaks if x in Pval_peaks]
                NI_X.append(NI)
                #
                par_prop_dic = gf['par']['%s_prop'%par]
                par_prop_NI = 1.0-par_prop_dic[:,0]/float(NImax)
                par_prop = par_prop_dic[:,1]
                par_mean = np.mean(par_prop)
                par_min = par_prop.min()
                par_max = par_prop.max()
                par_rmsd = np.sqrt(np.mean((par_prop-par_mean)**2))
                par_prop_list.append([par_mean,par_min,par_max,par_rmsd])
                par_prop_peak_list_NI = np.concatenate((par_prop_peak_list_NI, par_prop_NI))
                par_prop_peak_list = np.concatenate((par_prop_peak_list, par_prop))
                peak_plot_list.append([len(gfit_peaks),len(Pval_peaks),len(gfit_peaks_in_Pval)])

            NI_X = 1.0-array(NI_X)/float(NImax)
            #
            #ax.plot(par_prop_peak_list_NI,par_prop_peak_list,'.',c='b') # All peaks proportion
            peak_plot_arr = array(peak_plot_list).T
            ax2.plot(NI_X, peak_plot_arr[0],'o',c='b',label='Peaks used for global fit')
            ax2.plot(NI_X, peak_plot_arr[1],'o',c='r',label='Ftest peaks')
            ax2.plot(NI_X, peak_plot_arr[2],'s',markersize=8,c='g',alpha=0.5,label='Fitted peaks in Ftest')
            #
            par_prop_arr = array(par_prop_list).T
            ax.errorbar(NI_X,par_prop_arr[0],yerr=par_prop_arr[3],marker='s',c='b',label='%s \n%s proportion \nglobal fit/single fit'%(met,par)) #ax.lines[-1].get_color() #Rmsd for peaks proportion
            ax.plot(NI_X,par_prop_arr[1],c='b',alpha=0.1) #Upper graph
            ax.plot(NI_X,par_prop_arr[2],c='b',alpha=0.1) #Lower graph
            imets+=1
            #
            ax.set_ylabel('%s prop'%par)
            ax.set_ylim(0,10)
            ax.set_xlabel('Sparseness')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.95),prop={'size':8}) # Put a legend to the right of the current axis
            #ax.legend(loc='best')
            ax.grid('on')
            #
            peaks = dic['peakrange'][met]
            ax2.set_ylabel('Number of peaks')
            ax2.set_ylim(0,len(peaks))
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.65),prop={'size':8}) # Put a legend to the right of the current axis
            #ax2.legend(loc='best')
            ax2.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()

def plot_single_pearson(dic,par,mets=False,NIa=False,Ini=False):
    if not mets: mets = dic['qMDDmet'][0]
    if not NIa: NIarr = dic['NIarr'][met][:2] #
    else: NIarr = NIa
    print NIarr
    for NI in NIarr:
        if not Ini: t = ''
        else: t = '_ini'
        #fig = figure('Pearson_NI=%s, par=%s %s'%(NI,par,t)) #figsize=(figsize, figsize)#/1.618
        imet = 1
        for met in mets:
            print 'Pearson_NI=%s, par=%s met=%s %s'%(NI,par,met,t)
            fig = figure('Pearson_NI=%s, par=%s met=%s %s'%(NI,par,met,t),figsize=(figsize, figsize)) #figsize=(figsize, figsize)#/1.618
            #fig = figure()
            ax = fig.add_subplot('111')
            #ax = fig.add_subplot('%s1%s'%(len(mets),imet))
            data = dic['gfit'][met][str(NI)]['par']['Corr']['%s_data%s'%(par,t)]
            peaknames = dic['gfit'][met][str(NI)]['par']['Corr']['%s_resn%s'%(par,t)]
            Pearson_Corr_Coeff = dic['gfit'][met][str(NI)]['par']['Corr']['%s_Pearson_Corr_Coeff%s'%(par,t)]
            Pearson_Corr_Coeff_tailed_p_value = dic['gfit'][met][str(NI)]['par']['Corr']['%s_Pearson_Corr_Coeff_tailed_p_value%s'%(par,t)]
            lin_slope = dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_slope%s'%(par,t)]
            lin_inter = dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_inter%s'%(par,t)]
            lin_r_value = dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_r_value%s'%(par,t)]
            lin_p_value = dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_p_value%s'%(par,t)]
            lin_std_err = dic['gfit'][met][str(NI)]['par']['Corr']['%s_lin_std_err%s'%(par,t)]
            slope = dic['gfit'][met][str(NI)]['par']['Corr']['%s_slope'%par]
            ax.errorbar(data[:,1],data[:,3],xerr=data[:,2],yerr=data[:,4],fmt='.',label='%s %s %s'%(met,par,t)) #ax.lines[-1].get_color()
            ax.plot(np.concatenate(([0], data[:,1])),np.concatenate(([0], data[:,1])),'-',color=ax.lines[-1].get_color(),label='Slope=1')
            ax.plot(np.concatenate(([0], data[:,1])),np.concatenate(([0], data[:,1]))*slope,'-',label='a=%3.2f, pearson=%3.2f'%(slope,Pearson_Corr_Coeff))
            for i in range(len(data[:,1])):
                ax.annotate('%s'%(peaknames[i]), xy=(data[:,1][i],data[:,3][i]), xycoords='data', xytext=(data[:,1][i],data[:,3][i]), textcoords='data', rotation=-45, size=8) #'axes fraction'
            ax.set_ylabel('Global fit')
            if not Ini: ax.set_xlabel('Single fit')
            else: ax.set_xlabel('Initial global fit')
            #ax.set_ylim(0,max(data[:,1]*1.05))
            #ax.set_xlim(0,max(data[:,1]*1.05))
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.95),prop={'size':8}) # Put a legend to the right of the current axis
            ax.legend(loc='best')
            ax.grid('on')
            imet+=1
        print "Write: show() , to see graphs"
    return()

def plot_glob_pearsons(dic,pars,mets=False,NIstop=False,Ini=False,gkey='rates'):
    if not mets: mets = dic['qMDDmet'][0]
    plot_chisqr_glob_pearson(dic,mets,NIstop,Ini,gkey)
    for par in pars:
        plot_glob_pearson(dic,par,mets,NIstop,Ini,gkey)
    return()

def plot_chisqr_glob_pearson(dics,mets=False,NIstop=False,Ini=False,gkey='rates'):
    for dic in dics:
        desc_name = dic['desc']['name']
        if not Ini: t = ''
        else: t = '_ini'
        fig = figure('chisqr_%s%s'%(desc_name,t),figsize=(figsize, figsize/1.618))
        if not mets: mets = dic['qMDDmet'][0]
        imets = 0
        for met in mets:
            ax = fig.add_subplot('%s1%s'%(len(mets),imets+1))
            ax2 = ax.twinx()
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            chisqr_prop_list = []
            peak_plot_list = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]['par']['Corr']
                gfit_peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
                Pval_peaks = dic['%s'%gkey][met][str(NI)]['Pval_peaks']
                gfit_peaks_in_Pval = [x for x in gfit_peaks if x in Pval_peaks]
                NI_X.append(NI)
                #
                Pearson_Corr_Coeff = gf['chisqr_Pearson_Corr_Coeff%s'%(t)]
                #print Pearson_Corr_Coeff, NI, met
                chisqr_prop_list.append(Pearson_Corr_Coeff)
                peak_plot_list.append([len(gfit_peaks),len(Pval_peaks),len(gfit_peaks_in_Pval)])

            NI_X = 1.0-array(NI_X)/float(NImax)
            # The thing we want to see
            ax.plot(NI_X, chisqr_prop_list,'.-',c='b',label='chisqr %s %s'%(met,t))
            # Peak numbers
            peak_plot_arr = array(peak_plot_list).T
            ax2.plot(NI_X, peak_plot_arr[0],'o',c='b',label='Peaks used for global fit')
            ax2.plot(NI_X, peak_plot_arr[1],'o',c='r',label='Ftest peaks')
            ax2.plot(NI_X, peak_plot_arr[2],'s',markersize=8,c='g',alpha=0.5,label='Fitted peaks in Ftest')
            imets+=1
            #
            ax.set_ylabel('chisqr pearson')
            ax.set_ylim(0.0,1.2)
            ax.set_xlabel('Sparseness')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.95),prop={'size':8}) # Put a legend to the right of the current axis
            #ax.legend(loc='best')
            ax.grid('on')
            #
            peaks = dic['peakrange'][met]
            ax2.set_ylabel('Number of peaks')
            ax2.set_ylim(0,len(peaks))
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.65),prop={'size':8}) # Put a legend to the right of the current axis
            #ax2.legend(loc='best')
            ax2.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()

def plot_glob_pearson(dics,par,mets=False,NIstop=False,Ini=False,gkey='rates'):
    for dic in dics:
        desc_name = dic['desc']['name']
        if not Ini: t = ''
        else: t = '_ini'
        fig = figure('%s_%s%s'%(par,desc_name,t),figsize=(figsize, figsize/1.618))
        if not mets: mets = dic['qMDDmet'][0]
        imets = 0
        for met in mets:
            ax = fig.add_subplot('%s1%s'%(len(mets),imets+1))
            ax2 = ax.twinx()
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            par_prop_list = []
#            par_prop_peak_list_NI = array([])
#            par_prop_peak_list = array([])
            peak_plot_list = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]['par']['Corr']
                gfit_peaks = dic['gfit'][met][str(NI)]['gfit_peaks']
                Pval_peaks = dic['%s'%gkey][met][str(NI)]['Pval_peaks']
                gfit_peaks_in_Pval = [x for x in gfit_peaks if x in Pval_peaks]
                NI_X.append(NI)
                #
                Pearson_Corr_Coeff = gf['%s_Pearson_Corr_Coeff%s'%(par,t)]
                #print Pearson_Corr_Coeff, NI, met
                par_prop_list.append(Pearson_Corr_Coeff)
                peak_plot_list.append([len(gfit_peaks),len(Pval_peaks),len(gfit_peaks_in_Pval)])

            NI_X = 1.0-array(NI_X)/float(NImax)
            # The thing we want to see
            ax.plot(NI_X, par_prop_list,'.-',c='b',label='%s %s %s'%(par,met,t))
            # Peak numbers
            peak_plot_arr = array(peak_plot_list).T
            ax2.plot(NI_X, peak_plot_arr[0],'o',c='b',label='Peaks used for global fit')
            ax2.plot(NI_X, peak_plot_arr[1],'o',c='r',label='Ftest peaks')
            ax2.plot(NI_X, peak_plot_arr[2],'s',markersize=8,c='g',alpha=0.5,label='Fitted peaks in Ftest')
            imets+=1
            #
            ax.set_ylabel('%s pearson'%par)
            ax.set_ylim(0.0,1.2)
            ax.set_xlabel('Sparseness')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.95),prop={'size':8}) # Put a legend to the right of the current axis
            #ax.legend(loc='best')
            ax.grid('on')
#            #
            peaks = dic['peakrange'][met]
            ax2.set_ylabel('Number of peaks')
            ax2.set_ylim(0,len(peaks))
            box = ax2.get_position()
            ax2.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # Shink current axis by 20%
            ax2.legend(loc='center left', bbox_to_anchor=(1.05, 0.65),prop={'size':8}) # Put a legend to the right of the current axis
            #ax2.legend(loc='best')
            ax2.grid('on')
    print "Write: show() , to see graphs"
    #show()
    return()


