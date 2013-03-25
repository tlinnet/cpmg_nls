from pylab import *
import scipy.optimize
import os
import TB

dt = [('allplanes_','.ser'),('allplanes_','.stats')]
#Data BBL
#BBLn = [('bbl-75/33','FT',10,50,'y'),'/home/tlinnet/kte/t1rho/bblM_20130104_pH6_5C_0Murea_normal/analysis_FT/int_corr_ft_method_all_awk_full/',[['CS',[]],['coMDD',[]]]]
BBLn = [('bbl-75/33','FT',10,50,'y'),os.path.join('.','data','bblM_20130104_pH6_5C_0Murea_normal','analysis_FT','int_corr_ft_method_all_awk_full'),[['CS',[]],['coMDD',[]]]]
BBLn=TB.insdat(BBLn,dt)

#TB.plotstats(BBLn,('CS','coMDD'))

test = TB.arrdata(BBLn,('CS'))
#tr = test[0][2][0]

