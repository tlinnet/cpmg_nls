from pylab import *
from scipy import optimize

def fitfunc(x,a,b):
    return a*x + b # Target function
def errfunc(p, x, y):
    return fitfunc(x,*p) - y # Distance to the target function

dat = genfromtxt('testdataP225.txt')
datX = dat[:,0]
datY = dat[:,1]
print datX

p = (2.0, 2.0) 
pa_lea,suc = optimize.leastsq(errfunc, p, args=(datX, datY))
print suc
fitYlea = fitfunc(datX,*pa_lea)

pa_cur, pcov = optimize.curve_fit(fitfunc, datX, datY,p0=p)
print pcov
fitYcur = fitfunc(datX,*pa_cur)
variance = diagonal(pcov) #Refer [3]
StdErr = sqrt(variance)
print StdErr

plot(datX, datY, "o",label='data') # Plot of the data
plot(datX, fitYlea, ".-",label='least') # Plot of the data
plot(datX, fitYcur, ".-",label='cur') # Plot of the data
legend(loc="upper left")

show()
