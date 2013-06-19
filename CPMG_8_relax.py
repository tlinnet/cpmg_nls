#!/sbinlab2/software/python-enthought-dis/epd-7.3-2-rh5-x86_64/bin/python2.7

import sys, os
import numpy as np
import re
import math
import collections
from collections import defaultdict

skiplinesser = 13
mol = 'protein'

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>1)

if len(sys.argv) > 4:
    timerecords = sys.argv[1]
    if not os.path.exists(timerecords):
        print "%s does not exists"%timerecords
        sys.exit()
    else:
        timebase, timeext = os.path.splitext(timerecords)
        times = open(timerecords, 'r').readlines()
    sernames = sys.argv[2]
    if not os.path.exists(sernames):
        print "%s does not exists"%sernames
        sys.exit()
    else:
        serbase, serext = os.path.splitext(sernames)
        serfiles = open(sernames, 'r').readlines()
    if len(serfiles) != len(times):
        print "Number of files listed in %s(%s), does not match number in %s(%s)"%(serfiles,len(serfiles),times,len(times))
    TIMET2 = float(sys.argv[3])
    SFRQ = float(sys.argv[4])

else:
    print "Please execute with arguments: ./CPMG_8_relax.py time.dat time.dat $SEL $TIMET2 $SFRQ"
    print "Please execute with arguments: ./CPMG_8_relax.py time.dat table_ser_files.txt 0.05 750.0614444"
    sys.exit()

outmodelname = serbase+"_model.txt"
outmodel=open(outmodelname,"w")
outrelaxname = serbase+"_relax.py"
outrelax=open(outrelaxname,"w")
iniser = open(serfiles[0].rstrip(), 'r')
iniserlines = iniser.readlines()

all_resis = []
non_assigned = []
for i,line in enumerate(iniserlines[skiplinesser:]):
    INDEX, X_AXIS, Y_AXIS, X_PPM, Y_PPM, VOL, ASS, Z_A0 = line.split()
    resinr = re.findall(r"\d+", ASS)[0]
    REHN = re.findall(r"HN", ASS)
    resitest = ASS[0]
    if resitest.isdigit():
        non_assigned.append(int(resinr))
    elif not resitest.isdigit() and len(REHN) == 0:
        non_assigned.append(int(resinr))
    else:
        all_resis.append(int(resinr))    
iniser.close()

max_resis = max(all_resis)
nearest = int(math.ceil(max_resis / 100.0)) * 100
print "not assigned peaks with line nr in SPARKY file: %s" % non_assigned

dic = collections.OrderedDict()
serfile = open(serfiles[0].rstrip(), 'r')
serfilelines = serfile.readlines()
for j,line in enumerate(serfilelines[skiplinesser:]):
    INDEX, X_AXIS, Y_AXIS, X_PPM, Y_PPM, VOL, ASS, Z_A0 = line.split()
    resinr = re.findall(r"\d+", ASS)[0]
    REHN = re.findall(r"HN", ASS)
    resitest = ASS[0]
    if resitest.isdigit():
        resinr = nearest + int(resinr)
        resn = 'X'
    elif not resitest.isdigit() and len(REHN) == 0:
        resinr = nearest + int(resinr)
        resn = ASS
    else:
        resn = ASS[0]
    dic[str(resinr)] = collections.OrderedDict([('resn', resn), ('resi', resinr),('resn', resn), ('ASS', ASS), ('X_PPM', X_PPM), ('Y_PPM', Y_PPM), ('VOL', [])])
serfile.close()

timelist = []
for i,serfilei in enumerate(serfiles):
    serfile = open(serfilei.rstrip(), 'r')
    time = int(times[i].split()[1])
    timelist.append(time)
    serfilelines = serfile.readlines()
    for j,line in enumerate(serfilelines[skiplinesser:]):
        INDEX, X_AXIS, Y_AXIS, X_PPM, Y_PPM, VOL, ASS, Z_A0 = line.split()
        VOL = float(VOL)
        resinr = re.findall(r"\d+", ASS)[0]
        REHN = re.findall(r"HN", ASS)
        resitest = ASS[0]
        if resitest.isdigit():
            resinr = nearest + int(resinr)
        elif not resitest.isdigit() and len(REHN) == 0:
            resinr = nearest + int(resinr)

        vollist = dic[str(resinr)]['VOL']
        vollist.append(VOL)
        dic[str(resinr)]['VOL'] = vollist
    serfile.close()

for key, value in dic.iteritems():
    resinr = key
    resi = value['resi']
    resn = value['resn']
    ASS = value['ASS']
    X_PPM = value['X_PPM']
    Y_PPM = value['Y_PPM']
    VOL = value['VOL']
    string = "%s %s %s %s N "%(mol,resi,resn,resi)
    for intens in VOL:
        string = string + "%1.6e "%intens
    outmodel.write(string + "\n")
outmodel.close()

intensity_cols = range(6,6+len(serfiles))
print "The spectrum ID string:"
print timelist
print "The frequence list:"
print np.array(timelist)/TIMET2
print "The intensity columns:"
print intensity_cols, "\n"


outrelax.write("# Create the 'rx' data pipe." + "\n")
outrelax.write("pipe.create(pipe_name='origin rx', pipe_type='relax_disp', bundle='rx')" + "\n")
outrelax.write("" + "\n")
outrelax.write("# The type of experiment" + "\n")
outrelax.write("relax_disp.exp_type(exp_type='cpmg fixed')" + "\n")
outrelax.write("" + "\n")
outrelax.write("# Read the sequence from file" + "\n")
#outrelax.write("sequence.read(file='%s',"%(os.path.join(os.getcwd(),outmodelname)) + "\n")
outrelax.write("sequence.read(file='%s',"%(os.path.join(outmodelname)) + "\n")
outrelax.write("dir=None, spin_id_col=None, mol_name_col=1, res_num_col=2, res_name_col=3, spin_num_col=4, spin_name_col=5, sep=None, spin_id=None)" + "\n")
outrelax.write("" + "\n")

outrelax.write("# Read the intensities from columns" + "\n")

# Settings
int_method='point sum' # height
integration_points = 1

for i,time in enumerate(timelist):
    spec_id = i
    freq = time/TIMET2
    int_col = intensity_cols[i]
    if freq == 0:
        CPMGfreq = None
    else:
        CPMGfreq = freq
    #outrelax.write("spectrum.read_intensities(file='%s',"%(os.path.join(os.getcwd(),outmodelname)) + "\n")
    outrelax.write("spectrum.read_intensities(file='%s',"%(os.path.join(outmodelname)) + "\n")
    outrelax.write("dir=None, spectrum_id='%s_%s', heteronuc='N', proton='HN', int_method='%s', int_col=(%s), spin_id_col=None, mol_name_col=1,"%(spec_id,time,int_method,int_col) + "\n")
    outrelax.write("res_num_col=2, res_name_col=3, spin_num_col=4, spin_name_col=5, sep=None, spin_id=None, ncproc=None)" + "\n")
    outrelax.write("relax_disp.cpmg_frq(spectrum_id='%s_%s', cpmg_frq=%s)"%(spec_id,time,CPMGfreq) + "\n")
    outrelax.write("relax_disp.relax_time(spectrum_id='%s_%s', time=%s)"%(spec_id,time,TIMET2) + "\n")
    outrelax.write("spectrometer.frequency(id='%s_%s', frq=%s, units='MHz')"%(spec_id,time,SFRQ) + "\n")
    #outrelax.write("spectrum.integration_points(N=%s, spectrum_id='0_2', spin_id=None)"%integration_points + "\n")
    outrelax.write("\n")

outrelax.write("# Replicated spectrums" + "\n")
print "Replicated spectrums"
for dup in sorted(list_duplicates(timelist)):
    string = '['
    for dupi in dup[1]:
        string += "'%s_%s',"%(dupi,dup[0])
    string += ']'
    #print dup[0],dup[1],
    print string
    outrelax.write("spectrum.replicated(spectrum_ids=%s)"%string + "\n")

outrelax.write("# Spin isotope" + "\n")
outrelax.write("spin.isotope(isotope='15N', spin_id='@N*', force=True)" + "\n")

outrelax.close()
print "Now run:"
print "relax_disp -g -l LOGFILE.txt %s"%outrelaxname
