from pylab import *
import numpy
import scipy.optimize
import os

figsize = 16
titfont = 26
labfont = 12
figfont = 12

def R1r(x,R1,R2):
    return R1*cos(x*pi/180)**2+R2*sin(x*pi/180)**2

def R1r_exch(x,y,R1,R2,phi,kEX):
    return R1*cos(x*pi/180)**2+(R2+phi*kEX/((2*pi*y/tan(x*pi/180))**2+(2*pi*y)**2+kEX**2))*sin(x*pi/180)**2

def insdat(dat,dt):
    p = dat[1]
    for i in range(len(dat[2])):
        mf=dat[2][i]
        for j in range(len(dt)):
            ft = dt[j]
            path = os.path.join(p,mf[0],ft[0]+mf[0]+ft[1])
            if ft[-1] == '.ser':
                data = genfromtxt(path)
                datas = genfromtxt(path,usecols=(6),converters = {6: lambda s: str(s)})
            else: 
                data = genfromtxt(path)
                datas=[]
            dat[2][i][1].append([data,datas])
    return(dat)

def plotstats(dat,ft):
    dc = dat[0]
    figA = figure(figsize=(figsize, figsize/1.618))
    figA.suptitle('Tada',fontsize=titfont)
    figA1 = figA.add_subplot(111)
    for i in range(len(dat[2])):
        if dat[2][i][0] in ft:
            NIs= dat[2][i][1][1][0][:,0]
            datX = 1 - NIs/max(NIs)
            datY = dat[2][i][1][1][0][:,38]
            figA1.plot(datX,datY,'.-',label='%s : %s %s %s %s'%(dat[2][i][0], dc[0],dc[1],dc[2],dc[3]))
    figA1.set_ylim([0,0.01])
    figA1.legend(loc="upper left")
    figA1.set_ylabel('Residual stdev',fontsize=labfont)
    figA1.set_xlabel("Sparseness",fontsize=labfont)
    show()

def arrdata(dat,ft):
    collect = []
    for i in range(len(dat[2])):
        if dat[2][i][0] in ft:
            met = dat[2][i][0]
            NIs = dat[2][i][1][1][0][:,0]
            NIs = NIs.astype(int)
            NIsm = max(NIs)
            peaki = dat[2][i][1][0][0][:,0]
            peaki = peaki.astype(int)
            pm = max(peaki)
            peaks = []
            intens = []
            for j in range(pm):
                peak = peaki[j]
                peakn = dat[2][i][1][0][1][j]
                intens = dat[2][i][1][0][0][:,5][j::pm]
                NIintens = []
                for k in range(len(NIs)):
                    NIint = [NIs[k],dat[2][i][1][0][0][:,8+k][j::pm]]
                    NIintens.append(NIint)
                peaks.append([peak,peakn,['Full FT int',intens],NIintens])
            collect.append([met,NIsm,peaks])
    return(collect)
