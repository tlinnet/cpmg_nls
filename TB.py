from pylab import *
import numpy
import scipy.optimize
import os

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
        for peak in dic['peakrange'][met]:
        ##for peak in ['1']:
            dic['decay'][met][str(peak)] = {}
            peakname = Int[str(peak)]['resn']
            dic['decay'][met][str(peak)]['resn'] = peakname
            # Getting Chemical Shifts
            CS_N = Int[str(peak)]['CS_N']
            CS_H = Int[str(peak)]['CS_H']
            ### !!!!
            ## for NI in NIarr[:-NIrem]:
            for NI in NIarr[:1]:
                dic['decay'][met][str(peak)][str(NI)] = {}
                i = 0
                for fs in range(0,filenr,slicet):
                ##for fs in range(0,5,slicet):
                    dic['decay'][met][str(peak)][str(NI)][str(fs)] = {}
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

                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['data'] = x_y_fit_resi
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['a'] = [dec_dic['p'][0],dec_dic['psterr'][0]]
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['b'] = [dec_dic['p'][1],dec_dic['psterr'][1]]
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['R1r_rates'] = R1r_rates
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['R1r_err'] = R1r_err
                    dic['decay'][met][str(peak)][str(NI)][str(fs)].update(dec_dic)
                    # Setting keys
                    offset = dic['offset'][i]
                    omega1 = dic['omega1'][i]
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['offset'] = offset
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['omega1'] = omega1
                    # Setting average chemical shifts
                    CS_H_fs = CS_H[fs:fe]
                    CS_H_fs_mean =CS_H_fs.mean()
                    CS_N_fs = CS_N[fs:fe]
                    CS_N_fs_mean = CS_N_fs.mean()
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['CS_H_mean'] = CS_H_fs_mean
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['CS_N_mean'] = CS_N_fs_mean
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
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['OMEGA'] = OMEGA
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['omegaEFF'] = omegaEFF
                    dic['decay'][met][str(peak)][str(NI)][str(fs)]['theta'] = theta #theta/tiltAngle
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
            for peak in peaks:
                peakname = dic['decay'][met][str(peak)]['resn']
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

def plotR1(dics,mets=['CS'],NIs=[96]):
    fig = figure()
    ax = fig.add_subplot(211)
    bx = fig.add_subplot(212)
    R1r_coef=[2,23]
    R1r_exch_coef=[1,40,10000,100000]
    for dic in dics:
        slicet = len(dic['time'])
        for met in mets:
            filenr = dic['filenr'][met]
            Dec = dic['decay'][met]
            for NI in NIs:
                #for peak in dic['peakrange'][met]:
                for peak in range(1,2):
                    tiltAngle_arr = []
                    R1r_rates_arr = []
                    R1r_err_arr = []
                    peakname = Dec[str(peak)]['resn']
                    for fs in range(0,filenr,slicet):
                    ##for fs in range(0,5,slicet):
                        theta = Dec[str(peak)][str(NI)][str(fs)]['theta']
                        tiltAngle_arr.append(theta)
                        R1r_rates = Dec[str(peak)][str(NI)][str(fs)]['R1r_rates']
                        R1r_rates_arr.append(R1r_rates)
                        R1r_err = Dec[str(peak)][str(NI)][str(fs)]['R1r_err']
                        R1r_err_arr.append(R1r_err)
                    datX = array(tiltAngle_arr)
                    datY = array(R1r_rates_arr)
                    f_sigma = array(R1r_err_arr)
                    print "before", R1r_coef,
                    resa = scipy.optimize.curve_fit(f_R1r, datX, datY, p0=R1r_coef, sigma=f_sigma, full_output=1)
                    R1r_dic = getstatpar(resa, datY)
                    R1r_coef = R1r_dic['p']
                    print "after", R1r_coef
                    calcR1r = f_R1r(array(tiltAngle_arr),*R1r_coef)

                    datXs = sort(datX)
                    fitXlin = linspace(datXs[0],datXs[-1],100)
                    fitYlina = f_R1r(fitXlin,*R1r_coef)
                    ax.plot(datX, datY, ".",label='%s %s %s. fit f_R1r: %s '%(met, peakname, NI, R1r_coef))
                    ax.plot(fitXlin,fitYlina,"-",color=ax.lines[-1].get_color())
    title('R1 fitting')
    ax.set_xlabel('Tilt angle')
    ax.set_ylabel('R1r_coef')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height]) # Shink current axis by 20%
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':8}) # Put a legend to the right of the current axis
    ax.grid('on')
    show()
    return()
