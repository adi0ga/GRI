import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
for iterno in range(5):
    fil=open(f"iterate{iterno}.txt","r")
    num=[]
    x=[]
    y=[]
    z=[]
    val=[]
    #index=0
    try:
        while True:
            string=fil.readline()
            lis=string.split("&")
            num.append(eval(lis[0]))
            val.append(eval(lis[2]))
            lis1=lis[1].split(",")
            x.append(eval(lis1[0]))
            y.append(eval(lis1[1]))
            z.append(eval(lis1[2]))    
    except:
        fil.close()
    def fit_func(x,a,mu,sigma,c):
        """gaussian function used for the fit"""
        return a * norm.pdf(x,loc=mu,scale=sigma) + c
    data_x=np.multiply(x,val)
    mu=np.mean(data_x)
    sigma=np.std(data_x)
    #make up some normally distributed data and do a histogram
    hist,left = np.histogram(x,50,(-25.0,25.0),weights=val)
    centers = left[:-1] + (left[1] - left[0])/2
    
    #fit the histogram
    p0 = [10,0,15,5] #starting values for the fit
    p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)
    
    #plot the histogram and fit together
    fig,ax = plt.subplots()
    ax.hist(x,50,(-25.0,25.0),weights=val,)
    ket = np.linspace(left[0],left[-1],1000)
    y_fit = fit_func(ket, *p1)
    ax.plot(ket,y_fit,'r-')
    plt.xlabel("X-Coordinate")
    plt.ylabel("Net Lambda(b)")
    plt.title(f"X-Coordinate for iteration {iterno+1} with mu={p1[1]:.6}and sigma={p1[2]:.6}")
    plt.show()
    hist,left = np.histogram(y,268,(-134.0,134.0),weights=val)
    centers = left[:-1] + (left[1] - left[0])/2
    
    #fit the histogram
    p0 = [12000,0,10,1] #starting values for the fit
    p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)
    data_y=np.multiply(y,val)
    mu=np.mean(data_y)
    sigma=np.std(data_y)
    #plot the histogram and fit together
    fig,ax = plt.subplots()
    ax.hist(y,268,(-134.0,134.0),weights=val,)
    ket = np.linspace(left[0],left[-1],1000)
    y_fit = fit_func(ket, *p1)
    #ax.plot(ket,y_fit,'r-')
    plt.xlabel("Y-Coordinate")
    plt.ylabel("Net Lambda(b)")
    
    plt.title(f"Y-Coordinate for iteration {iterno+1} ")
    
    #yh=plt.hist(y,200,weights=val)
    plt.show()
    hist,left = np.histogram(z,268,(-134.0,134.0),weights=val)
    centers = left[:-1] + (left[1] - left[0])/2
    
    #fit the histogram
    p0 = [12000,0,10,1] #starting values for the fit
    p1,_ = curve_fit(fit_func,centers,hist,p0,maxfev=10000)
    data_z=np.multiply(z,val)
    mu=np.mean(data_z)
    sigma=np.std(data_z)   
    #plot the histogram and fit together
    fig,ax = plt.subplots()
    ax.hist(z,268,(-134.0,134.0),weights=val,)
    ket = np.linspace(left[0],left[-1],1000)
    y_fit = fit_func(ket, *p1)
    #ax.plot(ket,y_fit,'r-')
    plt.xlabel("Z-Coordinate")
    plt.ylabel("Net Lambda(b)")
    
    plt.title(f"Z-Coordinate for iteration {iterno+1} ")
    
    #zh=plt.hist(z,200,weights=val)
    plt.show()


