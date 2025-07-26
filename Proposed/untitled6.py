import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import norm
import matplotlib.mlab as mlab
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
        if (xnew>=-25.0)and(xnew<=25.0)and(ynew>=-134.0)and(ynew<=134.0)and(znew>=-134.0)and(znew<=134.0):
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


# make these smaller to increase the resolution
dx, dy = 1,  1

x1 = np.arange(0.5, 50.5, dx)
y1= np.arange(0.5, 50.5, dy)
X, Y = np.meshgrid(x1, y1)
x2=list(x1)
y2=list(y1)
dicx={}
for i in range(len(x2)):
    dicx[str(x2[i])]=0
    for j in  range(len(x)):
        if x2[i]==x[j]:
            dicx[str(x2[i])]=dicx[str(x2[i])]+val[j]
dicy={}
for i in range(len(y2)):
    dicy[str(y2[i])]=0
    for j in  range(len(y)):
        if y2[i]==y[j]:
            dicy[str(y2[i])]=dicy[str(y2[i])]+val[j]

# when layering multiple images, the images need to have the same
# extent.  This does not mean they need to have the same shape, but
# they both need to render to the same coordinate system determined by
# xmin, xmax, ymin, ymax.  Note if you use different interpolations
# for the images their apparent extent could be different due to
# interpolation edge effects

extent = np.min(x1), np.max(x1), np.min(y1), np.max(y1)
fig = plt.figure(frameon=False)
"""
Z1 = np.add.outer(range(8), range(8)) % 2  # chessboard
im1 = plt.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
                 extent=extent)
"""
Z2 = func3(x2, y2)

im2 = plt.imshow(Z2, cmap=plt.cm.viridis, alpha=1, interpolation='bilinear',
                 extent=extent)

plt.show()
sns.heatmap(Z2)
plt.show()
#xh=plt.hist(x,200,weights=val)
def fit_func(x,a,mu,sigma,c):
    """gaussian function used for the fit"""
    return a * norm.pdf(x,loc=mu,scale=sigma) + c

#make up some normally distributed data and do a histogram
hist,left = np.histogram(x,200,(-25.0,25.0))
centers = left[:-1] + (left[1] - left[0])

#fit the histogram
p0 = [10,0,15,0] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(x,200,(-25.0,25.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
#ax.plot(ket,y_fit,'r-')

plt.title(f"X-Coordinate for intersection data with mu={p1[1]:.6}and sigma={p1[2]:.6}")
plt.show()
hist,left = np.histogram(y,200,(-134.0,134.0))
centers = left[:-1] + (left[1] - left[0])

#fit the histogram
p0 = [100,0,15,0] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(y,268,(-134.0,134.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
#ax.plot(ket,y_fit,'r-')

plt.title(f"Y-Coordinate for intersection data with mu={p1[1]:.6} and sigma={p1[2]:.6}")

#yh=plt.hist(y,200,weights=val)
plt.show()
hist,left = np.histogram(z,268,(-134.0,134.0))
centers = left[:-1] + (left[1] - left[0])

#fit the histogram
p0 = [100,0,50,0] #starting values for the fit
p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)

#plot the histogram and fit together
fig,ax = plt.subplots()
ax.hist(z,268,(-134.0,134.0))
ket = np.linspace(left[0],left[-1],1000)
y_fit = fit_func(ket, *p1)
#ax.plot(ket,y_fit,'r-')

plt.title(f"Z-Coordinate for intersection data with mu={p1[1]:.6}and sigma={p1[2]:.6}")

#zh=plt.hist(z,200,weights=val)
plt.show()