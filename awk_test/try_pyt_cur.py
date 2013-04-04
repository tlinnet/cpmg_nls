from pylab import *
import scipy.optimize
import lmfit

fitfunc = lambda x,a,b,c:a*np.exp(-b*x)+c # Target fitfunction
errfitfunc = lambda p, x, y: fitfunc(x,*p) - y # Distance to the target fitfunction
def lmfitfunc(pars, x, data=None,eps=None):
    amp = pars['amp'].value
    decay = pars['decay'].value
    const = pars['const'].value
    model = amp*np.exp(-decay*x)+const
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data)/eps

datX = np.linspace(0,4,50)
pguess = [2.5, 1.3, 0.5]
datY = fitfunc(datX,*pguess)
datYran = datY + 0.2*np.random.normal(size=len(datX))

####### Try least squares
lea = {}
lea['par'], lea['cov_x'], lea['infodict'], lea['mesg'], lea['ier'] = scipy.optimize.leastsq(errfitfunc, pguess, args=(datX, datYran), full_output=1)
print lea['par'], lea['ier']
datY_lea = fitfunc(datX,*lea['par'])

####### Try curve_fit
cur = {}
cur['par'], cur['pcov'], cur['infodict'], cur['mesg'], cur['ier'] = scipy.optimize.curve_fit(fitfunc, datX, datYran, p0=pguess, full_output=1)
datY_cur=fitfunc(datX,*cur['par'])
cur['par_variance'] = diagonal(cur['pcov']); cur['par_stderr'] = sqrt(cur['par_variance'])
# Read this: http://mail.scipy.org/pipermail/scipy-user/2009-March/020516.html
cur['chisq']=sum(cur['infodict']['fvec']*cur['infodict']['fvec']) # calculate final chi square
cur['NDF']=len(datY)-len(cur['par'])
cur['RMS_residuals'] = sqrt(cur['chisq']/cur['NDF'])
print cur['par'], cur['ier'], cur['chisq'], cur['par_stderr']

####### Try lmfit
par = lmfit.Parameters()
par.add('amp', value=2.5, vary=True)
par.add('decay', value=1.3, vary=True)
par.add('const', value=0.5, vary=True)
lmf = lmfit.minimize(lmfitfunc, par, args=(datX, datYran),method='leastsq')
#datY_lmfit =datYran+lmf.residual
datY_lmfit = lmfitfunc(par,datX)
# See http://cars9.uchicago.edu/software/python/lmfit/fitting.html#fit-results-label
print par['amp'].value, par['amp'].stderr, par['amp'].correl
print lmf.nfev, lmf.success, lmf.errorbars, lmf.nvarys, lmf.ndata, lmf.nfree, lmf.chisqr, lmf.redchi
lmfit.printfuncs.report_errors(par) #lmf.params

#####################  This part is just to explore lmfit
#ci, trace = lmfit.conf_interval(lmf,sigmas=[0.68,0.95],trace=True, verbose=0)
#lmfit.printfuncs.report_ci(ci)
#x, y, grid=lmfit.conf_interval2d(lmf,'amp','decay',30,30)
#x1,y1,prob1=trace['amp']['amp'], trace['amp']['decay'],trace['amp']['prob']
#x2,y2,prob2=trace['decay']['decay'], trace['decay']['amp'],trace['decay']['prob']

#figure(1)
#contourf(x,y,grid)
#scatter(x1,y1,c=prob1,s=30)
#scatter(x2,y2,c=prob2,s=30)
#xlabel('amp');
#colorbar();
#ylabel('decay');
######################

figure(2)
subplot(3,1,1)
plot(datX,datY,".-",label='real')
plot(datX,datYran,'o',label='random')
plot(datX,datY_lea,'.-',label='leastsq fit')
legend(loc="best")
subplot(3,1,2)
plot(datX,datY,".-",label='real')
plot(datX,datYran,'o',label='random')
plot(datX,datY_cur,'.-',label='curve fit')
legend(loc="best")
subplot(3,1,3)
plot(datX,datY,".-",label='real')
plot(datX,datYran,'o',label='random')
plot(datX,datY_lmfit,'.-',label='lmfit')
legend(loc="best")

show()
