from pylab import *
from scipy import optimize

def fitfunc(x,a,b):
    return a*x + b # Target function
def errfitfunc(p, x, y):
    return fitfunc(x,*p) - y # Distance to the target function

dat = genfromtxt('testdataP225.txt')
datX = dat[:,0]
datY = dat[:,1]
pguess = (2.0, 2.0)

lea = {}
lea['par'], lea['cov_x'], lea['infodict'], lea['mesg'], lea['ier'] = optimize.leastsq(errfitfunc, pguess, args=(datX, datY), full_output=1)
print lea['par'], lea['ier']
datY_lea = fitfunc(datX,*lea['par'])

cur = {}
cur['par'], cur['pcov'], cur['infodict'], cur['mesg'], cur['ier'] = optimize.curve_fit(fitfunc, datX, datY, p0=pguess, full_output=1)
print datX, type(datX)
datY_cur=fitfunc(datX,*cur['par'])
cur['par_variance'] = diagonal(cur['pcov']); cur['par_stderr'] = sqrt(cur['par_variance'])
# Read this: http://mail.scipy.org/pipermail/scipy-user/2009-March/020516.html
cur['chisq']=sum(cur['infodict']['fvec']*cur['infodict']['fvec']) # calculate final chi square
cur['NDF']=len(datY)-len(cur['par'])
cur['RMS_residuals'] = sqrt(cur['chisq']/cur['NDF'])
print cur['par'], cur['ier'], cur['chisq']

subplot(2,1,1)
plot(datX,datY,"o",label='data')
plot(datX,datY_lea, ".-",label='least')
legend(loc="upper left")
subplot(2,1,2)
plot(datX,datY,"o",label='data')
plot(datX,datY_cur, ".-",label='cur')
show()
