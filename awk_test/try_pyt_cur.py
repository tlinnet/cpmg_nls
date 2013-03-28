from pylab import *
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

func = lambda x,a,b,c:a*np.exp(-b*x)+c # Target function
errfunc = lambda p, x, y: func(x,*p) - y # Distance to the target function

x = np.linspace(0,4,50)
p = [2.5, 1.3, 0.5]
yreal = func(x,*p)
yran = yreal + 0.2*np.random.normal(size=len(x))

p_lea,suc = leastsq(errfunc, p, args=(x, yreal))
print suc
yfit_lea = func(x,*p_lea)

p_cur, pc_cur = curve_fit(func, x, yran)
print pc_cur
yfit_cur=func(x,*p_cur)

plot(x,yreal,label='real')
plot(x,yran,'o',label='random')
plot(x,yfit_lea,'.-',label='leastsq fit')
plot(x,yfit_cur,'.-',label='curve fit')
legend(loc="best")

show()
