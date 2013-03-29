from pylab import *
from scipy import optimize

fitfunc = lambda x,a,b,c:a*np.exp(-b*x)+c # Target fitfunction
errfitfunc = lambda p, x, y: fitfunc(x,*p) - y # Distance to the target fitfunction

datX = np.linspace(0,4,50)
pguess = [2.5, 1.3, 0.5]
datY = fitfunc(datX,*pguess)
datYran = datY + 0.2*np.random.normal(size=len(datX))



lea = {}
lea['par'], lea['cov_x'], lea['infodict'], lea['mesg'], lea['ier'] = optimize.leastsq(errfitfunc, pguess, args=(datX, datYran), full_output=1)
print lea['par'], lea['ier']
datY_lea = fitfunc(datX,*lea['par'])

cur = {}
cur['par'], cur['pcov'], cur['infodict'], cur['mesg'], cur['ier'] = optimize.curve_fit(fitfunc, datX, datYran, p0=pguess, full_output=1)
#print datX, type(datX)
datY_cur=fitfunc(datX,*cur['par'])
cur['par_variance'] = diagonal(cur['pcov']); cur['par_stderr'] = sqrt(cur['par_variance'])
# Read this: http://mail.scipy.org/pipermail/scipy-user/2009-March/020516.html
cur['chisq']=sum(cur['infodict']['fvec']*cur['infodict']['fvec']) # calculate final chi square
cur['NDF']=len(datY)-len(cur['par'])
cur['RMS_residuals'] = sqrt(cur['chisq']/cur['NDF'])
print cur['par'], cur['ier'], cur['chisq']

subplot(2,1,1)
plot(datX,datY,".-",label='real')
plot(datX,datYran,'o',label='random')
plot(datX,datY_lea,'.-',label='leastsq fit')
legend(loc="best")
subplot(2,1,2)
plot(datX,datY,".-",label='real')
plot(datX,datYran,'o',label='random')
plot(datX,datY_cur,'.-',label='curve fit')
legend(loc="best")

show()
