#!/sbinlab2/software/python-enthought-dis/epd-7.3-2-rh5-x86_64/bin/python2.7

import sys, os
import numpy as np

if len(sys.argv) > 1:
    filename = sys.argv[1]
    filebase, fileext = os.path.splitext(filename)
else:
    print "Please execute with arguments: table_plots/table_simple_complex.fit 0.95 1 50"
    sys.exit()

outdata=open(filebase+".sherekhan","w")
outdata.write("60.12\n0.040000\n")
outdata.write("#nu_cpmg(Hz) \t R2(1/s) \t Esd(R2)\n")

tf = np.recfromtxt(filename,names=True)
datanames= tf.dtype.names
resis = datanames[1:]
curfreq = -1
for i,resi in enumerate(resis):
    val = []
    prevfreq = tf[0][0]
    outdata.write("# %s\n"%resi)
    for j,line in enumerate(tf):
        freq = line[0]
        pint = line[i+1]
        if freq != prevfreq and prevfreq != 0.00:
            val2 = np.array(val)
            #print resi, freq, prevfreq, val2, val2.std(), val2.mean()
            outdata.write("%3.2f \t %3.6f \t %3.6f\n"%(prevfreq,val2.mean(),val2.std()))
        if freq != curfreq:
            curfreq = freq
            val = [pint]
        else:
            val.append(pint)
        prevfreq=freq
    val2 = np.array(val)
    #print resi, freq, prevfreq, val2, val2.std(), val2.mean()
    outdata.write("%3.2f \t %3.6f \t %3.6f\n"%(freq, val2.mean(),val2.std()))

outdata.close()            
