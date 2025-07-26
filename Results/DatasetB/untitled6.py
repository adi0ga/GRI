import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
fil=open("intersect.txt","r")
num=[]
x=[]
y=[]
z=[]
val=[]
#index=0
try:
    while True:
        string=fil.readline()
        lis=string.split("-",maxsplit=1)
        num.append(eval(lis[0]))
        val.append(1)
        lis1=lis[1].split(",")
        lis2=lis1[2].split("--")
        znew=eval(lis2[0])
        xnew=eval(lis1[0])
        ynew=eval(lis1[1])
        if (xnew>=0.0)and(xnew<=50.0)and(ynew>=0.0)and(ynew<=50.0)and(znew>=0.0)and(znew<=400.0):
            x.append(xnew)
            y.append(ynew)
            z.append(znew)    
except:
    fil.close()


def func3(xq, yq):
    zq=[]
    for xi in range(len(xq)):
        zq.append(len(yq)*[0])
    for xi in range(len(xq)):
        for yi in range(len(yq)):
            for i in range(len(x)):
                if (((x[i]+0.5)>=xq[xi])and((x[i]-0.5)<=xq[xi]))and(((y[i]+0.5)>=yq[yi])and((y[i]-0.5)<=yq[yi])):
                    zq[xi][49-yi]=zq[xi][49-yi]+1
    return zq

def fit_func(x,a,mu,sigma,c):
    """gaussian function used for the fit"""
    return a * norm.pdf(x,loc=mu,scale=sigma) + c

#make up some normally distributed data and do a histogram
hist,left = np.histogram(x,50,(0.0,50.0))
centers = left[:-1] + (left[1] - left[0])/2

#fit the histogram
p0 = [10,15,15,5] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(x,50,(0.0,50.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
ax.plot(ket,y_fit,'r-')
plt.xlabel("X-Coordinate")
plt.ylabel("Net Lambda(b)")
plt.title(f"X-Coordinate for intersection data with mu={p1[1]:.6}and sigma={p1[2]:.6}")
plt.show()
hist,left = np.histogram(y,50,(0.0,50.0))
centers = left[:-1] + (left[1] - left[0])/2

#fit the histogram
p0 = [100,25,15,0] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(y,50,(0.0,50.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
ax.plot(ket,y_fit,'r-')
plt.xlabel("Y-Coordinate")
plt.ylabel("Net Lambda(b)")

plt.title(f"Y-Coordinate for intersection data with mu={p1[1]:.6} and sigma={p1[2]:.6}")

#yh=plt.hist(y,200,weights=val)
plt.show()
hist,left = np.histogram(z,400,(0,400.0))
centers = left[:-1] + (left[1] - left[0])/2

#fit the histogram
p0 = [10,200,50,0] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(z,400,(0.0,400.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
ax.plot(ket,y_fit,'r-')
plt.xlabel("Z-Coordinate")
plt.ylabel("Net Lambda(b)")

plt.title(f"Z-Coordinate for intersection data with mu={p1[1]:.6}and sigma={p1[2]:.6}")

#zh=plt.hist(z,200,weights=val)
plt.show()