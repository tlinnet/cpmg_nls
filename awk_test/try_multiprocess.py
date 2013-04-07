import numpy as np
import multiprocessing
from multiprocessing import Pool
from datetime import datetime
from joblib import Parallel, delayed
#http://www.scipy.org/Cookbook/Multithreading?action=AttachFile&do=view&target=test_handythread.py
from handythread import foreach

def getsqrt(n):
    res = np.sqrt(n**2)
    return(res)

def main():
    jobs = multiprocessing.cpu_count()-1
    a = range(10000)
    for method in ['normal','multi Pool','joblib delayed','handythread']:
        startTime = datetime.now()
        sprint=True
        if method=='normal':
            res = []
            for i in a:
                b = getsqrt(i)
                res.append(b)
        elif method=='multi Pool':
            pool = Pool(processes=jobs)
            res = pool.map(getsqrt, a)
        elif method=='joblib delayed':
            res = Parallel(n_jobs=jobs)(delayed(getsqrt)(i) for i in a)
        elif method=='handythread':
            res = foreach(getsqrt,a,threads=jobs,return_=True)
        else:
            sprint=False
        if sprint:
            print "Method was %s"%method
            print "Done :%s"%(datetime.now()-startTime)
            print res[-10:], type(res[-1])
    return(res)

if __name__ == "__main__":
    res = main()