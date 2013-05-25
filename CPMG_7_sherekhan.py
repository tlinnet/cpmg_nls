#!/sbinlab2/software/python-enthought-dis/epd-7.3-2-rh5-x86_64/bin/python2.7

import sys, os
import numpy as np
import re
import math
import collections

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

if len(sys.argv) > 3:
    filename = sys.argv[1]
    filebase, fileext = os.path.splitext(filename)
    TIMET2 = sys.argv[2]
    YOBS = sys.argv[3]
else:
    print "Please execute with arguments: ./CPMG_7_sherekhan.py table.txt $TIMET2 $YOBS"
    print "./CPMG_7_sherekhan.py table.txt 0.05 76.012 "
    sys.exit()

outdata=open(filebase+".sherekhan","w")
outdata.write("%s\n%s\n"%(YOBS,TIMET2))
outdata.write("#nu_cpmg(Hz) \t R2(1/s) \t Esd(R2)\n")

tf = np.recfromtxt(filename,names=True)
datanames= tf.dtype.names
resis = datanames[1:]
all_resis = []
non_assigned = []

for i,resi in enumerate(resis):
    resinr = re.findall(r"\d+", resi)[0]
    resit = resi[0]
    if resit.isdigit():
        non_assigned.append(int(resinr))
    else:
        all_resis.append(int(resinr))

max_resis = max(all_resis)
nearest = int(math.ceil(max_resis / 100.0)) * 100
print "not assigned peaks with line nr in SPARKY file: %s" % non_assigned

curfreq = -1
dic = collections.OrderedDict()
freqs = []
for i,resi in enumerate(resis):
    val = []
    prevfreq = tf[0][0]
    resinr = re.findall(r"\d+", resi)[0]
    resit = resi[0]
    if resit.isdigit():
        resinr = nearest + int(resinr)
        resn = 'X'
    else:
        resn = resi[0]
    dic[str(resinr)] = collections.OrderedDict({'resn':resn,'resi':resinr})
    #outdata.write("# %s\n"%resi)

    for j,line in enumerate(tf):
        freq = line[0]
        pint = line[i+1]
        if freq != prevfreq and prevfreq != 0.00:
            val = np.array(val)
            #print resi, freq, prevfreq, val, val.std(), val.mean()
            valstd = val.std()
            dic[str(resinr)][str(prevfreq)]={'freq':prevfreq,'mean':val.mean(),'std':valstd}
            freqs.append(prevfreq)
        if freq != curfreq:
            curfreq = freq
            val = [pint]
        else:
            val.append(pint)
        prevfreq=freq
    val = np.array(val)
    valstd = val.std()
    dic[str(resinr)][str(prevfreq)]={'freq':prevfreq,'mean':val.mean(),'std':valstd}
    freqs.append(prevfreq)
    #print resi, freq, prevfreq, val, val.std(), val.mean()

dic = collections.OrderedDict(sorted(dic.items(), key=lambda t: int(t[0])))

freqs = f7(freqs)
for key, value in dic.iteritems():
    #print key, value
    resi = value['resi']
    resn = value['resn']
    outdata.write("# %s%s\n"%(resn,resi))
    maxvalstd = 0
    for freq in freqs:
        valstd = value[str(freq)]['std']
        if valstd > maxvalstd:
            maxvalstd = valstd
    for freq in freqs:
        valstd = value[str(freq)]['std']
        if valstd == 0.0:
            valstd = maxvalstd+maxvalstd*0.01
        valmean = value[str(freq)]['mean']
        outdata.write("%3.2f \t %3.7f \t %3.7f \n"%(freq,valmean,valstd))

print "Now use out ", filebase+".sherekhan", "to upload to sherekhan.bionmr.org/app/calculation"
outdata.close()
