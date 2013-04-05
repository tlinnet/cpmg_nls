from pylab import *
import numpy
import scipy.optimize
import scipy.stats.distributions
import os
import lmfit #See http://newville.github.com/lmfit-py/parameters.html

import sys
class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "w")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
sys.stdout = Logger("Log_file.txt")

figsize = 16
titfont = 26
labfont = 12
figfont = 12

############ Fit functions
def f_expdecay_lmfit(pars,time,data=None): #KTE: extract_sums_to_table.pl. Line 68.
    amp = pars['amp'].value
    decay = pars['decay'].value
    model = amp*exp(-decay*time)
    if data is None:
        return model
    return (model-data)

def f_R1r_lmfit(pars,tiltAngle,data=None,eps=None): #KTE: R1rhoAnalysis.ipf. Line 144. Use see line 96
    #http://newville.github.com/lmfit-py/fitting.html
    R1 = pars['R1'].value
    R2 = pars['R2'].value
    model = R1*cos(tiltAngle*pi/180.0)**2+R2*sin(tiltAngle*pi/180.0)**2
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data)/eps

def f_R1r_exch_lmfit(pars,inp,data=None,eps=None): #KTE: R1rhoAnalysis.ipf. Line 144. Use see line 96
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

def f_R1r_exch_global(pars,sel_p,tilt_a,om1_a,R1rex_a,R1rex_err_a=None):
    toterr = np.array([])
    for i in range(len(sel_p)):
        p = sel_p[i]
        tilt =tilt_a[i];om1=om1_a[i];R1rex=R1rex_a[i];R1rex_err=R1rex_err_a[i]
        par = lmfit.Parameters()
        par.add('kEX', value=pars['kEX'].value, vary=True, min=0.0)
        par.add('R1', value=pars['R1%s'%p].value, vary=True)
        par.add('R2', value=pars['R2%s'%p].value, vary=True)
        par.add('phi', value=pars['phi%s'%p].value, vary=True, min=0.0)
        datX = [array(tilt), array(om1)]
        Yfit = f_R1r_exch_lmfit(par,datX)
        if R1rex_err_a is None:
            erri = Yfit - R1rex
        else:
            erri = (Yfit - R1rex)/R1rex_err
        toterr = np.concatenate((toterr, erri))
    #print par['kEX']
    return (toterr)

def unpack_global(dic, p_list, met, NI):
    for i in range(len(p_list)):
        p = p_list[i]
        par = lmfit.Parameters()
        dic['gfit'][met][str(NI)][str(p)] = {}
        kEX = dic['gfit'][met][str(NI)]['par']['kEX']
        R1 = dic['gfit'][met][str(NI)]['par']['R1%s'%p]
        R2 = dic['gfit'][met][str(NI)]['par']['R2%s'%p]
        phi = dic['gfit'][met][str(NI)]['par']['phi%s'%p]
        par['kEX'] = kEX; par['R1'] = R1; par['R2'] = R2; par['phi'] = phi
        dic['gfit'][met][str(NI)][str(p)]['par'] = par
        # Calc other parameters for the fit
        datX_f_R1r_exch = dic['rates'][met][str(NI)][str(p)]['R1r_exch']['data'][0]
        datY_f_R1r_exch = dic['rates'][met][str(NI)][str(p)]['R1r_exch']['data'][1]
        Yfit = f_R1r_exch_lmfit(par,datX_f_R1r_exch)
        dic['gfit'][met][str(NI)][str(p)]['Yfit'] = Yfit
        residual = Yfit - datY_f_R1r_exch
        dic['gfit'][met][str(NI)][str(p)]['residual'] = residual
        chisqr = sum(residual**2)
        dic['gfit'][met][str(NI)][str(p)]['chisqr'] = chisqr
        NDF = len(residual)-len(par)
        dic['gfit'][met][str(NI)][str(p)]['NDF'] = NDF
        dic['gfit'][met][str(NI)][str(p)]['what_is_this_called'] = sqrt(chisqr/NDF)
        dic['gfit'][met][str(NI)][str(p)]['redchisqr'] = chisqr/NDF
        #print "#####################  PEAK %s ##################################"%p
        #print "Peak %s .chisqr=%3.2f. kEX= "%(p,chisqr,kEX.value)
        #lmfit.printfuncs.report_errors(par)
    return()

############ Stat functions
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

############ Misc. functions
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

############# Data import
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
    return()
#############################################
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
            CS_H = data[:,3][j::pm]
            dic['Int'][met][peak]['CS_H'] = CS_H
            CS_N = data[:,4][j::pm]
            dic['Int'][met][peak]['CS_N'] = CS_N
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

def getdecay(dic,mets,NIstop=False):
    datX = dic['time']
    centerPPM = dic['NMRpar']['centerPPM']
    frq = dic['NMRpar']['frq']
    slicet = len(datX)
    dic['decay'] = {}
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
            print "%s - Getting decay for NI=%s"%(met,NI)
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
                    dic['decay'][met][str(NI)][str(peak)][str(fs)] = {}
                    fe = fs+slicet
                    FTInt = Int[str(peak)]['FTInt'][fs:fe]
                    NIInt = Int[str(peak)]['NIInt'][str(NI)][fs:fe]
                    divi = FTInt.argmax()
                    datY = FTInt*NIInt/(FTInt[divi]*NIInt[divi])
                    par = lmfit.Parameters()
                    par.add('amp', value=1.0, vary=True)
                    par.add('decay', value=10.0, vary=True)
                    lmf = lmfit.minimize(f_expdecay_lmfit, par, args=(datX, datY),method='leastsq')
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['par'] = par
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['lmf'] = lmf
                    #fitY = datY+lmf.residual
                    fitY = f_expdecay_lmfit(par,datX)
                    R1r_rates = par['decay'].value
                    R1r_err = par['decay'].stderr
                    x_y_fit_resi = array([datX,datY,fitY]).T
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'] = x_y_fit_resi
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['R1r_rates'] = R1r_rates
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['R1r_err'] = R1r_err
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
    return()

def plotdecays(dics,mets=['CS'],peaks=False,NIarr=False,fss=[0]):
    #TB.plotdecays([BBL],['coMDD'])
    #TB.plotdecays([BBL],['coMDD','CS'],fss=range(0,35,5))
    fig = figure()
    ax = fig.add_subplot(111)
    for dic in dics:
        for met in mets:
            if not NIarr: NIarr = [dic['NImax'][met]]
            for NI in NIarr:
                if not peaks: peaks = ['1']
                print peaks, type(peaks)
                for peak in peaks:
                    peakname = dic['decay'][met][str(NI)][str(peak)]['resn']
                    for fs in fss:
                        datX = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,0]
                        datY = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,1]
                        fitY = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,2]
                        par = dic['decay'][met][str(NI)][str(peak)][str(fs)]['par']
                        a = par['amp'].value
                        a_e = par['amp'].stderr
                        b = par['decay'].value
                        b_e = par['decay'].stderr
                        datXs = sort(datX)
                        fitXlin = linspace(datXs[0],datXs[-1],100)
                        fitYlin = f_expdecay_lmfit(par,fitXlin)
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
    show()
    return()

def getrates(dic,mets=['CS'],NIstop=False):
    dic['rates'] = {}
    slicet = len(dic['time'])
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
            print "%s - Getting rates for NI=%s"%(met,NI)
            dic['rates'][met][str(NI)] = {}
            Pval_peaks = []
            for peak in dic['peakrange'][met]:
                dic['rates'][met][str(NI)][str(peak)] = {}
                dic['rates'][met][str(NI)][str(peak)]['R1r'] = {}
                dic['rates'][met][str(NI)][str(peak)]['R1r_exch'] = {}
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
                # Calculate for R1r
                try:
                    par_R1r = lmfit.Parameters()
                    par_R1r.add('R1', value=2.0, vary=True)
                    par_R1r.add('R2', value=23.0, vary=True)
                    lmf_R1r = lmfit.minimize(f_R1r_lmfit, par_R1r, args=(datX_f_R1r, datY, f_sigma),method='leastsq')
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['par'] = par_R1r
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['lmf'] = lmf_R1r
                    #calcR1r = datY+lmf_R1r.residual    # NOT WORKING?
                    calcR1r = f_R1r_lmfit(par_R1r,datX_f_R1r)
                    x_y_fit_resi = array([datX_f_R1r,datY,f_sigma,calcR1r]).T
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['data'] = x_y_fit_resi
                except RuntimeError as e:
                    print "Cannot fit R1r for %s %s. Reason: %s"%(peak, peakname, e)
                # Calculate for R1r_exch
                try:
                    par_R1r_exch = lmfit.Parameters()
                    par_R1r_exch.add('R1', value=dic['guess']['s_R1'], vary=True)
                    par_R1r_exch.add('R2', value=dic['guess']['s_R2'], vary=True)
                    par_R1r_exch.add('kEX', value=dic['guess']['s_kEX'], vary=True, min=0.0)
                    par_R1r_exch.add('phi', value=dic['guess']['s_phi'], vary=True, min=0.0)
                    #print par_R1r_exch['kEX'].value, par_R1r_exch['kEX'].min
                    lmf_R1r_exch = lmfit.minimize(f_R1r_exch_lmfit, par_R1r_exch, args=(datX_f_R1r_exch, datY, f_sigma),method='leastsq')
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['par'] = par_R1r_exch
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['lmf'] = lmf_R1r_exch
                    #calcR1r_exch = datY+lmf_R1r_exch.residual    # NOT WORKING?
                    calcR1r_exch = f_R1r_exch_lmfit(par_R1r_exch,datX_f_R1r_exch)
                    x_y_fit_resi = array([datX_f_R1r_exch,datY,f_sigma,calcR1r_exch]).T
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['data'] = x_y_fit_resi
                    Fval, Fdist, Pval = Ftest(lmf_R1r.chisqr,lmf_R1r.nfree,lmf_R1r_exch.chisqr,lmf_R1r_exch.nfree)
                except RuntimeError as e:
                    print "Cannot fit R1r_exch for %s %s. Reason: %s"%(peak, peakname, e)
                dic['rates'][met][str(NI)][str(peak)]['Fval'] = Pval
                dic['rates'][met][str(NI)][str(peak)]['Fdist'] = Fdist
                dic['rates'][met][str(NI)][str(peak)]['Pval'] = Pval
                if Pval==False:
                    pass
                elif Pval!=False:
                    Pval_peaks.append(peak)
            dic['rates'][met][str(NI)]['Pval_peaks'] = Pval_peaks
    print "Following peak numbers passed the Ftest"
    print Pval_peaks
    return(Pval_peaks)

def plotrates(dics,mets=False,peaks=False,NIarr=False):
    for dic in dics:
        unique_omega1 = sort(f5(dic['omega1']))
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            if not peaks: peaks = dic['peakrange'][met]
            if not NIarr: NIarr = [dic['NImax'][met]]
            for NI in NIarr:
                figR1r = figure('R1r %s'%NI)
                ax = figR1r.add_subplot(111)
                for peak in peaks:
                    dc = dic['rates'][met][str(NI)][str(peak)]
                    peakname = dc['resn']
                    Pval = dic['rates'][met][str(NI)][str(peak)]['Pval']
                    if Pval==False:
                    # Get values for R1r
                        datX_f_R1r = dc['R1r']['data'][:,0]
                        datY_f_R1r = dc['R1r']['data'][:,1]
                        datY_f_R1r_err = dc['R1r']['data'][:,2]
                        calcR1r = dc['R1r']['data'][:,3]
                        par_R1r = dc['R1r']['par']
                        datXs_f_R1r = sort(datX_f_R1r)
                        datYs_f_R1r = f_R1r_lmfit(par_R1r,datXs_f_R1r)
                        datXs_f_R1r_lin = linspace(datXs_f_R1r[0],datXs_f_R1r[-1],100)
                        datYs_f_R1r_lin = f_R1r_lmfit(par_R1r,datXs_f_R1r_lin)
                        #plot
                        ax.errorbar(datX_f_R1r, datY_f_R1r,yerr=datY_f_R1r_err,fmt='.',label='%s %s %s. R1=%3.2f+-%3.2f '%(met, peakname, NI, par_R1r['R1'].value, par_R1r['R1'].stderr))
                        ax.plot(datXs_f_R1r,datYs_f_R1r,"-",color=ax.lines[-1].get_color())
                        ax.plot(datX_f_R1r,calcR1r,"o",mfc='none',mec=ax.lines[-1].get_color())#
                        #ax.plot(datXs_f_R1r_lin,datYs_f_R1r_lin,"-",color=ax.lines[-1].get_color())
                    if Pval!=False:
                        figR1r_exch = figure('R1r_exch %s %s'%(NI,peak))
                        bx = figR1r_exch.add_subplot(111)
                    # Get values for R1r_exch
                        datX_f_R1r_exch = dc['R1r_exch']['data'][0]
                        datY_f_R1r_exch = dc['R1r_exch']['data'][1]
                        datY_f_R1r_exch_err = dc['R1r_exch']['data'][2]
                        calcR1r_exch = dc['R1r_exch']['data'][3]
                        par_R1r_exch = dc['R1r_exch']['par']
                        #plot
                        tiltAngle_arr_s, omega1_arr_s = zip(*sorted(zip(datX_f_R1r_exch[0], datX_f_R1r_exch[1])))
                        datXs_f_R1r_exch = [array(tiltAngle_arr_s), array(omega1_arr_s)]
                        datYs_f_R1r_exch = f_R1r_exch_lmfit(par_R1r_exch,datXs_f_R1r_exch)
                        omega_plots = []
                        for om1 in unique_omega1:
                            tiltAngle_arr_s_lin = linspace(tiltAngle_arr_s[0],tiltAngle_arr_s[-1],100)
                            om1_arr = ones(len(tiltAngle_arr_s_lin))*om1
                            datXs_f_R1r_exch_lin = [array(tiltAngle_arr_s_lin), array(om1_arr)]
                            datYs_f_R1r_exch_lin = f_R1r_exch_lmfit(par_R1r_exch,datXs_f_R1r_exch_lin)
                            omega_plots.append([datXs_f_R1r_exch_lin,datYs_f_R1r_exch_lin])
                        bx.errorbar(datX_f_R1r_exch[0], datY_f_R1r_exch, yerr=datY_f_R1r_exch_err,fmt=".",label='%s %s %s'%(met, peakname, NI)) #  fit f_R1r_exch: %s, dc['R1r_exch']['p']
                        for omp in omega_plots:
                            colval = dic['omega1_col'][str(omp[0][1][0])]
                            #print colval
                            bx.plot(omp[0][0],omp[1],"-", color=colval) # ,color=bx.lines[-1].get_color()
                        #bx.plot(datXs_f_R1r_exch[0],datYs_f_R1r_exch,"-",color=bx.lines[-1].get_color())
                        bx.plot(datX_f_R1r_exch[0],calcR1r_exch,"o",mfc='none',mec=bx.lines[-1].get_color())#,color=bx.lines[-1].get_color()
                        #R1r_exch graph
                        bx.set_title('R1_exch fitting')
                        bx.set_xlabel('Tilt angle')
                        bx.set_ylabel('R1r_coef')
                        box = bx.get_position()
                        bx.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
                        bx.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                        bx.grid('on')
                #R1r graph
                ax.set_title('R1 fitting')
                ax.set_xlabel('Tilt angle')
                ax.set_ylabel('R1r_coef')
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
                ax.grid('on')
    show()
    return()

def getglobfit(dic,mets=['coMDD'],peaks=False,NIstop=False):
    dic['gfit'] = {}
    for met in mets:
        dic['gfit'][met] = {}
        if not peaks:
            peaks = dic['rates'][met][str(NI)]['Pval_peaks']# dic['peakrange'][met]
            testPval = True
        else: testPval = False
        NIarr = dic['NIarr'][met]
        for NI in NIarr:
            if NIstop:
                if NI <= NIstop: break
                else: pass
            else:
                if NI <= dic['NIstop']: break
                else: pass
            print "################### METHOD %s########Global Fit NI=%s##############"%(met,NI)
            print "Getting Global fit for NI=%s"%NI
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
                if Pval!=False:
                # Get values for R1r_exch
                    dic['gfit'][met][str(NI)][str(peak)] = {}
                    datX_f_R1r_exch = dc['R1r_exch']['data'][0]
                    datY_f_R1r_exch = dc['R1r_exch']['data'][1]
                    datY_f_R1r_exch_err = dc['R1r_exch']['data'][2]
                    calcR1r_exch = dc['R1r_exch']['data'][3]
                    par_R1r_exch = dc['R1r_exch']['par']
                    R1,R2,kEX,phi = par_R1r_exch['R1'].value,par_R1r_exch['R2'].value,par_R1r_exch['kEX'].value,par_R1r_exch['phi'].value
                    #print "p: R1=%3.2f, R2=%3.2f, kEX=%3.2f, phi=%3.2f"%(R1,R2,kEX,phi)
                    tiltAngle_arr_s = datX_f_R1r_exch[0]; omega1_arr_s = datX_f_R1r_exch[1]
                    tilt.append(tiltAngle_arr_s); om1.append(omega1_arr_s)
                    R1rex.append(datY_f_R1r_exch)
                    R1rex_err.append(datY_f_R1r_exch_err)
                    P_arr.add('R1%s'%peak, value=par_R1r_exch['R1'].value, vary=True)
                    P_arr.add('R2%s'%peak, value=par_R1r_exch['R2'].value, vary=True)
                    P_arr.add('phi%s'%peak, value=par_R1r_exch['phi'].value, vary=True)
                    sel_p.append(peak)
                    dic['gfit'][met][str(NI)][str(peak)]['resn'] = peakname
            dic['gfit'][met][str(NI)]['gfit_peaks'] = sel_p
            lmf = lmfit.minimize(f_R1r_exch_global, P_arr, args=(sel_p,tilt,om1,R1rex,R1rex_err),method='leastsq')
            dic['gfit'][met][str(NI)]['par']= P_arr
            dic['gfit'][met][str(NI)]['lmf']= lmf
            unpack_global(dic, sel_p, met, NI) # Unpack the paramerts into the selected peaks
            print "Medthod=%s, NI=%s, kEX=%4.4f, chisqr=%4.4f"%(met,NI,P_arr['kEX'].value,lmf.chisqr)
    return()

def plot_kEX(dics,mets=False,peaks=False,NIstop=False):
    fig = figure()
    ax = fig.add_subplot(211)
    bx = fig.add_subplot(212)
    for dic in dics:
        if not mets: mets = dic['qMDDmet'][0]
        for met in mets:
            NIarr = dic['NIarr'][met]
            NImax = dic['NImax'][met]
            NI_X = []
            kEX_Y = []
            kEX_Y_e= []
            chisqr_Y = []
            for NI in NIarr:
                if NIstop:
                    if NI <= NIstop: break
                    else: pass
                else:
                    if NI <= dic['NIstop']: break
                    else: pass
                gf = dic['gfit'][met][str(NI)]
                kEX = gf['par']['kEX'].value
                kEX_e = gf['par']['kEX'].stderr
                chisqr = gf['lmf'].chisqr
                NI_X.append(NI)
                kEX_Y.append(kEX)
                kEX_Y_e.append(kEX_e)
                chisqr_Y.append(chisqr)
            NI_X = 1.0-array(NI_X)/float(NImax)
            kEX_Y = array(kEX_Y)
            kEX_Y_e= array(kEX_Y_e)
            chisqr_Y = array(chisqr_Y)
            ax.errorbar(NI_X, kEX_Y,yerr=kEX_Y_e,fmt='.-',label='%s'%(met))
            ax.plot(NI_X, ones(len(NI_X))*kEX_Y[0])
            bx.plot(NI_X,chisqr_Y,'.-',label='%s'%(met))
    ax.set_title('Global kEX fitting')
    ax.set_ylabel('kEX')
    ax.set_ylim(10000,20000)
    ax.legend(loc='upper left')
    bx.set_ylabel('chisqr')
    bx.set_xlabel('Sparseness')
    bx.legend(loc='upper left')
    show()
    return()