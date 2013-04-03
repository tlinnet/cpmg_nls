from pylab import *
import numpy
import scipy.optimize
import scipy.stats.distributions
import os
import lmfit

figsize = 16
titfont = 26
labfont = 12
figfont = 12

############ Fit functions
def f_expdecay(time,a,b): #KTE: extract_sums_to_table.pl. Line 68
    return a*exp(-b*time)

def f_R1r(tiltAngle,R1,R2): #KTE: R1rhoAnalysis.ipf. Line 144. Use see line 96
    #R1=R1r_coef, R2=R1r_rates
    return R1*cos(tiltAngle*pi/180.0)**2+R2*sin(tiltAngle*pi/180.0)**2

def f_R1r_exch(inp,R1,R2,kEX,phi): #KTE:  R1rhoAnalysis.ipf. Line 115. Use see line 102
    #R1=R1r_exch_coef, R2=R1r_rates
    tiltAngle,omega1=inp
    return R1*cos(tiltAngle*pi/180)**2+(R2+phi*kEX/((2*pi*omega1/tan(tiltAngle*pi/180))**2+(2*pi*omega1)**2+kEX**2))*sin(tiltAngle*pi/180)**2

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

############# Data modifications
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

def getdecay(dic,dt,NIrem=10):
    datX = dic['time']
    centerPPM = dic['NMRpar']['centerPPM']
    frq = dic['NMRpar']['frq']
    slicet = len(datX)
    dic['decay'] = {}
    for met in dt:
        Int = dic['Int'][met]
        NIarr = dic['NIarr'][met]
        filenr = dic['filenr'][met]
        dic['decay'][met] = {}
        ### !!!!
        ## for NI in NIarr[:-NIrem]:
        for NI in NIarr[:1]:
            dic['decay'][met][str(NI)] = {}
            for peak in dic['peakrange'][met]:
            ##for peak in ['1']:
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

                    pguess = (1.0,10.0)
                    res = scipy.optimize.curve_fit(f_expdecay, datX, datY,p0=pguess, full_output=1)
                    dec_dic = getstatpar(res, datY)

                    fitY = f_expdecay(datX,*dec_dic['p'])

                    R1r_rates = dec_dic['p'][1]
                    R1r_err = dec_dic['psterr'][1]
                    x_y_fit_resi = array([datX,datY,fitY,dec_dic['resi']]).T

                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'] = x_y_fit_resi
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['a'] = [dec_dic['p'][0],dec_dic['psterr'][0]]
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['b'] = [dec_dic['p'][1],dec_dic['psterr'][1]]
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['R1r_rates'] = R1r_rates
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['R1r_err'] = R1r_err
                    dic['decay'][met][str(NI)][str(peak)][str(fs)].update(dec_dic)
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
                    # Calculate NMR properties
                    ########## KTE: residue_files.pl
                    #OMEGA=($centerPPM-$chemShift{$peakName})*$frq+$offset{$i};
                    #$omegaEFF=sqrt($OMEGA**2+$omega1{$i}**2);
                    #if (($omega1{$i}/$OMEGA) > 0){
                    #$theta=180/$PI*abs(atan($omega1{$i}/$OMEGA));
                    #}else{
                    #$theta=180-180/$PI*abs(atan($omega1{$i}/$OMEGA));
                    ################################
                    OMEGA=(centerPPM-CS_N_fs_mean)*frq+offset
                    omegaEFF=sqrt(OMEGA**2+omega1**2)
                    if omega1/OMEGA > 0:
                        theta = 180/pi*abs(arctan(omega1/OMEGA))
                    else:
                        theta = 180- 180/pi*abs(arctan(omega1/OMEGA))
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['OMEGA'] = OMEGA
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['omegaEFF'] = omegaEFF
                    dic['decay'][met][str(NI)][str(peak)][str(fs)]['theta'] = theta #theta/tiltAngle
                    # tail residue_files.pl
                    #print OMEGA, omega1, omegaEFF, theta, R1r, R1r_err
                    # R1rhoAnalysis.ipf. Line 311.
                    # dof, slockpwr, w_eff, tiltAngle, R1r_rates, R1r_err
                    i+=1
    return()

def plotdecays(dics,mets=['CS'],peaks=[1],NIs=[96],fss=[0]):
    #TB.plotdecays([BBL],['coMDD'])
    #TB.plotdecays([BBL],['coMDD','CS'],fss=range(0,35,5))
    fig = figure()
    ax = fig.add_subplot(111)
    for dic in dics:
        for met in mets:
            for NI in NIs:
                for peak in peaks:
                    peakname = dic['decay'][met][str(NI)][str(peak)]['resn']
                    for fs in fss:
                        datX = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,0]
                        datY = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,1]
                        fitY = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,2]
                        resi = dic['decay'][met][str(NI)][str(peak)][str(fs)]['data'][:,3]
                        a = dic['decay'][met][str(NI)][str(peak)][str(fs)]['a'][0]
                        a_e = dic['decay'][met][str(NI)][str(peak)][str(fs)]['a'][1]
                        b = dic['decay'][met][str(NI)][str(peak)][str(fs)]['b'][0]
                        b_e = dic['decay'][met][str(NI)][str(peak)][str(fs)]['b'][1]

                        datXs = sort(datX)
                        fitXlin = linspace(datXs[0],datXs[-1],100)
                        fitYlin = f_expdecay(fitXlin,a,b)

                        ax.plot(datX,datY,".",label='%s %s %s %s. fit [a*exp(-b*x]: a=%3.2f +-%3.2f, a=%3.2f +-%3.2f '%(met, peakname, NI, fs, a,a_e,b,b_e))
                        ax.plot(fitXlin,fitYlin,"-",color=ax.lines[-1].get_color())
                        #ax.plot(,label='fit [a*exp(-b*x]: a=%3.2f +-%3.2f, a=%3.2f +-%3.2f'%(a,a_e,b,b_e))
                        #ax.plot(datX,fitY,"rp")
    title('Decay')
    ax.set_xlabel('Time')
    ax.set_ylabel('Normalized intensity')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
    ax.grid('on')
    show()
    return()

def getstatpar(result,datY):
    # http://mail.scipy.org/pipermail/scipy-user/2009-March/020516.html
    dic = {}
    dic['p'], dic['pcov'], dic['infodict'], dic['mesg'], dic['ier'] = result
    dic['pvar'] = diagonal(dic['pcov'])
    dic['psterr'] = sqrt(dic['pvar'])
    dic['resi'] = dic['infodict']['fvec']
    dic['chisq']=sum(dic['infodict']['fvec']*dic['infodict']['fvec']) # calculate final chi square
    dic['NDF']=len(datY)-len(dic['p'])
    dic['RMSresi'] = sqrt(dic['chisq']/dic['NDF'])
    dic['redchisq'] = dic['chisq']/dic['NDF']
    return(dic)

def f5(seq, idfun=None):
   # order preserving. http://www.peterbe.com/plog/uniqifiers-benchmark
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def getrates(dic,mets=['CS'],NIs=[96]):
    R1r_coef=[2,23]
    R1r_exch_coef=[1,40,10000,100000]
    dic['rates'] = {}
    slicet = len(dic['time'])
    for met in mets:
        dic['rates'][met] = {}
        filenr = dic['filenr'][met]
        Dec = dic['decay'][met]
        for NI in NIs:
            dic['rates'][met][str(NI)] = {}
            Pval_peaks = []
            for peak in dic['peakrange'][met]:
            #for peak in [21]: #21,22,27
            #for peak in [1,5]: #21,22,27
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
                    result_R1r = scipy.optimize.curve_fit(f_R1r, datX_f_R1r, datY, p0=R1r_coef, sigma=f_sigma, full_output=1)
                    R1r_dic = getstatpar(result_R1r, datY)
                    dic['rates'][met][str(NI)][str(peak)]['R1r'].update(R1r_dic)
                    calcR1r = f_R1r(datX_f_R1r,*R1r_dic['p'])
                    x_y_fit_resi = array([datX_f_R1r,datY,f_sigma,calcR1r,R1r_dic['resi']]).T
                    dic['rates'][met][str(NI)][str(peak)]['R1r']['data'] = x_y_fit_resi
                except RuntimeError as e:
                    print "Cannot fit R1r_exch for %s. Reason: %s"%(peakname, e)

                # Calculate for R1r_exch
                try:
                    result_R1r_exch = scipy.optimize.curve_fit(f_R1r_exch, datX_f_R1r_exch, datY, p0=R1r_exch_coef, sigma=f_sigma, full_output=1)
                    R1r_exch_dic = getstatpar(result_R1r_exch, datY)
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch'].update(R1r_exch_dic)
                    calcR1r_exch = f_R1r_exch(datX_f_R1r_exch,*R1r_exch_dic['p'])
                    x_y_fit_resi = array([datX_f_R1r_exch,datY,f_sigma,calcR1r_exch,R1r_exch_dic['resi']]).T
                    dic['rates'][met][str(NI)][str(peak)]['R1r_exch']['data'] = x_y_fit_resi

                    Fval, Fdist, Pval = Ftest(R1r_dic['chisq'],R1r_dic['NDF'],R1r_exch_dic['chisq'],R1r_exch_dic['NDF'])
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
    return()

def plotrates(dics,mets=False,peaks=False,NIs=[96]):
    for dic in dics:
        unique_omega1 = sort(f5(dic['omega1']))
        if not mets:
            mets = dic['qMDDmet'][0]
        for met in mets:
            if not peaks:
                peaks = dic['peakrange'][met]
            for NI in NIs:
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
                        datXs_f_R1r = sort(datX_f_R1r)
                        datYs_f_R1r = f_R1r(datXs_f_R1r,*dc['R1r']['p'])
                        datXs_f_R1r_lin = linspace(datXs_f_R1r[0],datXs_f_R1r[-1],100)
                        datYs_f_R1r_lin = f_R1r(datXs_f_R1r_lin,*dc['R1r']['p'])
                        #plot
                        ax.errorbar(datX_f_R1r, datY_f_R1r,yerr=datY_f_R1r_err,fmt='.',label='%s %s %s. fit f_R1r: %s '%(met, peakname, NI, dc['R1r']['p']))
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
                        #plot
                        tiltAngle_arr_s, omega1_arr_s = zip(*sorted(zip(datX_f_R1r_exch[0], datX_f_R1r_exch[1])))
                        datXs_f_R1r_exch = [array(tiltAngle_arr_s), array(omega1_arr_s)]
                        datYs_f_R1r_exch = f_R1r_exch(datXs_f_R1r_exch,*dc['R1r_exch']['p'])

                        omega_plots = []
                        for om1 in unique_omega1:
                            tiltAngle_arr_s_lin = linspace(tiltAngle_arr_s[0],tiltAngle_arr_s[-1],100)
                            om1_arr = ones(len(tiltAngle_arr_s_lin))*om1
                            datXs_f_R1r_exch_lin = [array(tiltAngle_arr_s_lin), array(om1_arr)]
                            datYs_f_R1r_exch_lin = f_R1r_exch(datXs_f_R1r_exch_lin,*dc['R1r_exch']['p'])
                            omega_plots.append([datXs_f_R1r_exch_lin,datYs_f_R1r_exch_lin])
                        bx.errorbar(datX_f_R1r_exch[0], datY_f_R1r_exch, yerr=datY_f_R1r_exch_err,fmt=".",label='%s %s %s'%(met, peakname, NI)) #  fit f_R1r_exch: %s, dc['R1r_exch']['p']
                        for omp in omega_plots:
                            colval = dic['omega1_col'][str(omp[0][1][0])]
                            #print colval
                            bx.plot(omp[0][0],omp[1],"-", color=colval) # ,color=bx.lines[-1].get_color()
                        #bx.plot(datXs_f_R1r_exch[0],datYs_f_R1r_exch,"-",color=bx.lines[-1].get_color())
                        bx.plot(datX_f_R1r_exch[0],calcR1r_exch,"o",mfc='none',mec=bx.lines[-1].get_color())#,color=bx.lines[-1].get_color()
                        #R1r_exch graph
                        bx.set_title('R1_exch fitting: %s'%met)
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

def getglobfit(dic,mets=['coMDD'],peaks=False,NIs=[96]):
    for met in mets:
        if not peaks:
            peaks = dic['peakrange'][met]
        for NI in NIs:
            for peak in peaks:
                dc = dic['rates'][met][str(NI)][str(peak)]
                peakname = dc['resn']
                Pval = dic['rates'][met][str(NI)][str(peak)]['Pval']
                if Pval!=False:
                # Get values for R1r_exch
                    datX_f_R1r_exch = dc['R1r_exch']['data'][0]
                    datY_f_R1r_exch = dc['R1r_exch']['data'][1]
                    datY_f_R1r_exch_err = dc['R1r_exch']['data'][2]
                    calcR1r_exch = dc['R1r_exch']['data'][3]
                    R1,R2,kEX,phi = dc['R1r_exch']['p']
                    print "p: R1=%3.2f, R2=%3.2f, kEX=%3.2f, phi=%3.2f"%(R1,R2,kEX,phi)
    return()
