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

def getstat(dic,dt):
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
    return()

def getser(dic,dt):
    fpath = dic['path']
    dic['ser']['data'] = {}
    dic['ser']['datas'] = {}
    dic['filenr'] = {}
    for met in dt:
        pre = dic['ser']['pre']
        filee = dic['ser']['filee']
        path = os.path.join(fpath,met,pre+met+filee)
        data = genfromtxt(path)
        datas = genfromtxt(path,usecols=(6),converters = {6: lambda s: str(s)})
        dic['ser']['data'][met] = data
        dic['ser']['datas'][met] = datas
        dic['filenr'][met] = int(data[:,0].size/data[:,0][-1])
    return()

def sortInt(dic,dt):
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
            dic['Int'][met][peak]['NIInt'] = {}
            for k in range(len(NIs)):
                NI = str(NIs[k])
                NIint = data[:,8+k][j::pm]
                dic['Int'][met][peak]['NIInt'][NI] = NIint
    return()

def plotstats(dics,dt):
    figA = figure(figsize=(figsize, figsize/1.618))
    figA.suptitle('Tada',fontsize=titfont)
    figA1 = figA.add_subplot(111)
    for dic in dics:
        for met in dt:
            data = dic['stats']['data'][met]
            NIs = data[:,0]        
            datX = 1 - NIs/dic['NImax'][met]
            datY = data[:,38]
            dc = dic['desc']
            figA1.plot(datX,datY,'.-', label='%s : %s'%(met, dc))
    figA1.set_ylim([0,0.01])
    figA1.legend(loc="upper left")
    figA1.set_ylabel('Residual stdev',fontsize=labfont)
    figA1.set_xlabel("Sparseness",fontsize=labfont)
    show()
    return()

def expdecay(x,a,b):
    return a*exp(-b*x) # Target function

def getdecay(dic,dt,NIrem=10):
    datX = dic['time']
    slicet = len(datX)
    dic['decay'] = {}
    for met in dt:
        Int = dic['Int'][met]
        NIarr = dic['NIarr'][met]
        dic['decay'][met] = {}
        for peak in Int.iterkeys(): #for peak in ['1']:
            dic['decay'][met][peak] = {}
            for NI in NIarr[:-NIrem]:
                dic['decay'][met][peak][str(NI)] = {}
                for fs in range(0,dic['filenr'][met],slicet): #for fs in range(45,65,slicet):
                    dic['decay'][met][peak][str(NI)][str(fs)] = {}
                    fe = fs+slicet  

                    peakname = Int[peak]['resn']
                    FTInt = Int[peak]['FTInt'][fs:fe]
                
                    NIInt = Int[peak]['NIInt'][str(NI)][fs:fe]
                    divi = FTInt.argmax()
                    datY = FTInt*NIInt/(FTInt[divi]*NIInt[divi])

                    pguess = (1.0,10.0)
                    par, parcov = scipy.optimize.curve_fit(expdecay, datX, datY,p0=pguess)
                    var = diagonal(parcov)
                    SE = sqrt(var)

                    fitY = expdecay(datX,*par)
                    residuals = datY-fitY

                    x_y_fit_resi = array([datX,datY,fitY,residuals]).T
                
                    dic['decay'][met][peak][str(NI)][str(fs)]['data'] = x_y_fit_resi
                    dic['decay'][met][peak][str(NI)][str(fs)]['a'] = [par[0],SE[0]]
                    dic['decay'][met][peak][str(NI)][str(fs)]['b'] = [par[1],SE[1]]
    return()

def plotdecays(dics,mets=['CS'],peaks=[1],NIs=[96],fss=[0]):
    #TB.plotdecays([BBL],['coMDD'])
    #TB.plotdecays([BBL],['coMDD','CS'],fss=range(0,35,5))
    fig = figure()
    ax = fig.add_subplot(111)
    for dic in dics:
        for met in mets:
            for peak in peaks:
                peakname = dic['Int'][met][str(peak)]['resn']
                for NI in NIs:
                    for fs in fss:
                        datX = dic['decay'][met][str(peak)][str(NI)][str(fs)]['data'][:,0]
                        datY = dic['decay'][met][str(peak)][str(NI)][str(fs)]['data'][:,1]
                        fitY = dic['decay'][met][str(peak)][str(NI)][str(fs)]['data'][:,2]
                        resi = dic['decay'][met][str(peak)][str(NI)][str(fs)]['data'][:,3]
                        a = dic['decay'][met][str(peak)][str(NI)][str(fs)]['a'][0]
                        a_e = dic['decay'][met][str(peak)][str(NI)][str(fs)]['a'][1]
                        b = dic['decay'][met][str(peak)][str(NI)][str(fs)]['b'][0]
                        b_e = dic['decay'][met][str(peak)][str(NI)][str(fs)]['b'][1]

                        datXs = sort(datX)
                        fitXlin = linspace(datXs[0],datXs[-1],100)
                        fitYlin = expdecay(fitXlin,a,b)
                    
                        ax.plot(datX, datY, "o",label='%s %s %s %s'%(met, peakname, NI, fs)) 
                        ax.plot(fitXlin, fitYlin, "r-",label='fit [a*exp(-b*x]: a=%3.2f +-%3.2f, a=%3.2f +-%3.2f'%(a,a_e,b,b_e))
                        #ax.plot(datX,fitY,"rp")
    title('Decay')
    ax.set_xlabel('Time')
    ax.set_ylabel('Normalized intensity')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':6}) # Put a legend to the right of the current axis
    ax.grid('on')
    show()
    return()




#			$OMEGA=($centerPPM-$chemShift{$peakName})*$frq+$offset{$i};
#			$omegaEFF=sqrt($OMEGA**2+$omega1{$i}**2);
#			if (($omega1{$i}/$OMEGA) > 0){
#				$theta=180/$PI*abs(atan($omega1{$i}/$OMEGA));
#			}else{
#				$theta=180-180/$PI*abs(atan($omega1{$i}/$OMEGA));
#			$chemShift{$process[6]}=$process[4]; chemshift is N
