# Fit Cepheids
# Finds best fits of periodic function to data from Cepheid Variables experiment

from matplotlib import pyplot as plt
from math import cos, pi, log, inf
import json

def ff(x, period):
    #function must be linear and have magnitude 1 at x=period,
    #or periodic over period and have magnitude 1 at x=period.
    return NotImplemented

def ff1(x, period):
    return x/period
    
def ff2(x, period):
    const = 2
    return (x % period)**const / period**const

def ff4(x, period, const):
    cutoff = period / (const + 1)
    m1 = 1 / (2 * cutoff)
    m2 = 1 / (2 * cutoff * const)
    if x % period <= cutoff:
        return x % period * m1
    if x % period > cutoff:
        return 0.5 + (x % period - cutoff) * m2

def ff42(x, period):
    return ff2(ff4(x, period) * period, period)

def yf(f, x, amplitude, period, phase, yoff, const):
    y = amplitude * cos(2*pi* f(x + phase, period, const)) + yoff
    return y

def chi2f(Y, fY, dY):
    return sum(((a-b)/berr)**2 for a, b, berr in zip(Y, fY, dY) if Y != None)

#import data from prev.txt

with open("prev.txt", "r+") as f:
    try:
        load = json.load(f)
    except:
        load = dict()

#calculate magnitudes, errors etc.

#Title (S#R# where S stands for Star and R stands for Region. Add E on the end if points have been excluded e.g. S1R2E)
Title = "S1R2"

#imported data
c_data = "62.7	78.365	39.189	54.968	28.087	40.186	41.511	65.193	40.55	36.602	49.442	75.656"
cerr_data = "5.03	4.878	4.36	4.578	4.162	4.367	4.392	4.767	4.456	4.518	4.713	4.8"

#python data
k = 1.067e12

D = dict()
force_save = False
re_eval = False

if Title in load and not re_eval:
    existing = True
    print("\nExisting data for {}\n".format(Title))
    D = load[Title]
else:
    existing = False
    print("\nNo existing data for {}\n".format(Title))
    D["chi2"] = inf
    
    omit = () #days to omit data
    
    D["c"] = [float(s) if not i+1 in omit else None for i, s in enumerate(c_data.split(sep = "\t"))]
    D["cerr"] = [float(s) if not i+1 in omit else None for i, s in enumerate(cerr_data.split(sep = "\t"))]
    D["c_omit"] = [float(s) if i+1 in omit else None for i, s in enumerate(c_data.split(sep = "\t"))]
    
    D["m"] = [-2.5*log(C/k, 10) for C in D["c"] if C != None]
    D["m_omit"] = [-2.5*log(C/k, 10) for C in D["c_omit"] if C != None]
    D["merr"] = [Cerr/C for Cerr, C in zip(D["cerr"],D["c"]) if C != None]

#line fitting variables
use_saved = True #Recall saved data, including the fit
auto_find = False #Attempt to find an improvement on the fit, automatically

period = 5.700000000000005
amplitude = -0.4
phase = -1.0999999999999999
yoff = 25.599999999999998
const = 0.30000000000000004

if use_saved and existing:
    amplitude, period, phase, yoff, const = D["A"], D["T"], D["dx"], D["dy"], D["C"]
    print("Using best known values:\nperiod = {1}\namplitude = {0}\nphase = {2}\nyoff = {3}\nconst = {4}\n".format(amplitude, period, phase, yoff, const))    

#construct datapoints
X=[i for i in range(1,len(D["c"])+1) if D["c"][i-1] != None]
X_omit=[i for i in range(1,len(D["c"])+1) if D["c"][i-1] == None]
fit=[yf(ff4, x, amplitude, period, phase, yoff, const) for x in X]

#auto fitfinder
changes = list()
for i in range(3**5):
    changes.append((i//(3**4)-1, (i//(3**3))%3-1, (i//(3**2))%3-1, (i//3)%3-1, i%3-1))
    #print(changes[-1])

testchi2 = D["chi2"]
improved = True
if auto_find and (not use_saved) and (existing):
    amplitude, period, phase, yoff, const = D["A"], D["T"], D["dx"], D["dy"], D["C"]
    print("AUTO")
    while improved:
        improved = False
        for change in changes:
            testfit = [yf(ff4, x, amplitude+change[0]*0.1, period+change[1]*0.1, phase+change[2]*0.1, yoff+change[3]*0.1, const+change[4]*0.1) for x in X]
            
            if chi2f(D["m"], testfit, D["merr"]) < testchi2:
                testchi2 = chi2f(D["m"], testfit, D["merr"])
                amplitude += change[0]*0.1
                period += change[1]*0.1
                phase += change[2]*0.1
                yoff += change[3]*0.1
                const += change[4]*0.1
                improved = True
                fit = testfit
                break

#visual datapoints
Xvis=[i*0.01 for i in range(0,1401)]
Yvis=[yf(ff4, x, amplitude, period, phase, yoff, const) for x in Xvis]
Ym=[i for i in D["m"]]
Ym_omit=[i for i in D["m_omit"]]

#fit calculations and data saving
newchi2 = chi2f(D["m"], fit, D["merr"])
print("old chi2: {}\nnew chi2: {}".format(D["chi2"], newchi2))

if newchi2 < D["chi2"] or force_save:
    linecol = "g"
    D["chi2"], D["A"], D["T"], D["dx"], D["dy"], D["C"] = newchi2, amplitude, period, phase, yoff, const
    load[Title] = D
    with open("prev.txt", "w") as f:
        json.dump(load, f)
        
elif newchi2 == D["chi2"]:
    linecol = "b"
    
else:
    linecol = "r"

#plot
plt.figure(figsize=(20,10))
    
plt.title(Title)
plt.xlabel("Day")
plt.ylabel("Apparent Magnitude")

plt.plot(Xvis, Yvis, linecol+"-")
#plt.plot(X, Ym, "r-")
plt.errorbar(X, Ym, D["merr"], fmt = "ko")
#plt.plot(X2, Yfit, linecol+"+")
#plt.plot(X_omit, Ym_omit, "kx")

ymin, ymax = plt.ylim()
plt.ylim((ymax, ymin))