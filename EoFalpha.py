# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 13:44:55 2015

@author: Fermi
"""

#This is the alpha test version of the integrated python EoF exercise script.
#Expect things that do not work


#here begins the main programs
#--------------------------------------------------------------

# Set parameters for the run
# load the data
#---------------

yearly=1; # if yearly==0 load monthly data from 1976 to 2005
          # if yearly==1 load yearly data from 1870 to present
          

# filtering parameter
#--------------------
          
trend=1;  # if trend==0 remove the linear trend, if trend==1 keep the trend.
seas=1;   # if seas==0  remove the seasonal cycle, if seas==1 keep the seasonal cycle

nym=0;    # width in years of the time filtering hanning window.

# reduce size 
#------------

pow2=0; # power of 2 of spatial reduction
        # if your computer is slow and/or you have low memory, you can perform the analysis 
        # on a reduce grid  of 2^pow2 degrees


# EOF analysis parameter
#------------------------

neof=8;     # number of EOF to be computed

            # region where EOF will be computed
# latmin=-5;latmax=5;lonmin=120;lonmax=280; # pacific equatorial
# latmin=0;latmax=60;lonmin=290;lonmax=360; # north atlantic
latmin=-90;latmax=90;lonmin=0;lonmax=360; # world

# added plot
#-----------

mean_pl=0;  # if mean_pl==1 plot the mean temperature 
trend_pl=0; # if trend_pl==1 plot hte trend
seas_pl=0;  # if seas_pl==1 plot the seasonal cycle

hoev_pl=0;  # if hoev_pl==1 plot time-longitude diagram at the equator
            # for both sst and sst_anomalies

##################################################################   
#       PARAMETER SETTING ENDS HERE
##################################################################

# import packages and functions that are used
import numpy as np
import scipy as sp
from numpy import transpose, reshape, shape, linalg
from scipy import signal, interpolate, io, sparse
import matplotlib.pyplot as plt

# Here begin the definition of functions
#-----------------------------------------------------------------

def sst_tmfilt(sst,timesst,yearly,trend,seas,nym):
    
    #annual average sst
    if yearly:
        timesst2=timesst[:,np.arange(12,np.size(timesst),12,dtype='int16')];
        facm=1;
    if not(yearly):
        timesst2=timesst;
        facm=12;
    # remove trend
    # the detrend function is used as scipy.signal.detrend
    tmp=sst*0;
    for i in range(np.size(sst,0)):
        for j in range(np.size(sst,1)):
            if not(trend):
                tmp[i,j,:]=signal.detrend(sst[i,j,:]); 
            if trend:
                tmp[i,j,:]=signal.detrend(sst[i,j,:],type='constant');
                
    
    if trend:
        print ('sst detrended in constant type')
            
        
    if not(trend):
        print('sst detrended')
        
    
    trnd=sst-tmp;
    sst=tmp;
    sstmean=np.mean(trnd,2); 
    ssttrend=(trnd[:,:,-1]-trnd[:,:,0])/(timesst2[0,-1]-timesst2[0,0]);
    
    #Remove seasonal cycle
    sstseas=sst[:,:,np.arange(12)]*0; 
    
    if not(yearly) and not(seas):
        for i in range(np.size(sst,0)):
            for j in range(np.size(sst,1)):
                for n in range(12):
                    sstseas[i,j,n]=np.mean(sst[i,j,np.arange(n,np.size(sst,2),12)]);
                    sst[i,j,np.arange(n,np.size(sst,2),12)]=sst[i,j,np.arange(n,np.size(sst,2),12)]-sstseas[i,j,n];
        
        print('seasonal cycle removed')
        
    
    #apply hanning window filter
    nmm=(nym)*facm/2;
    
    if nmm>0:
        coef=np.hanning(2*nmm-1)/nmm; #numpy has a hanning function
        #coef=ones(2*nmm-1,1)/(2*nmm-1);
        tmp=sst;
        for i in range(0,np.size(sst,2)):
            for n in range(nmm,np.size(sst,0)-nmm+1):
                tmp[n,:,i]=np.squeeze(sst[np.arange(int(n-nmm+1),int(n+nmm,1)),:,i])*coef;
                
            
        sst=tmp[np.arange(int(nmm),int(np.size(sst,0)-nmm+2)),:,:];
        timesst2=timesst2[np.arange(int(nmm),int(end-nmm+1))];

    print('time hanning filter completed')
    
    return sst,timesst2,facm,ssttrend,sstmean,sstseas
    

#define eof function for eigenvalue analysis
def sst_eof(sst,xsst,ysst,latmin,latmax,lonmin,lonmax,neof):
    
    # indices of the spatial area (cound be refined)
    condj=np.zeros(np.size(ysst))
    indj=ysst*0
    for ii in np.arange(np.size(ysst)):
        condj[ii]=(ysst[ii]>=latmin and ysst[ii] <=latmax)
        if condj[ii]==True:
            indj[ii]=ii
    
    condi=np.zeros(np.size(xsst))
    indi=xsst*0
    if lonmax>lonmin:
        for ii in np.arange(np.size(xsst)):
            condi[ii]=(xsst[ii]<=lonmax and xsst[ii]>=lonmin)
            if condi[ii]==True:
                indi[ii]=ii
    if not (lonmax>lonmin):
        for ii in np.arange(np.size(xsst)):
            condi[ii]=(xsst[ii]<=lonmax or xsst[ii]>=lonmin)
            if condi[ii]==True:
                indi[ii]=ii
    
    indi=np.compress(condi,indi).astype(int) # to make them integers as indices
    indj=np.compress(condj,indj).astype(int)
    
    #here the brackets are changed from () to []. label 1-0 is not required here
    
    #use extract to get falues for eof from True/Falce index condition function
    
    xeof = np.extract(condi,xsst)
    yeof = np.extract(condj,ysst)
    
    # indices of the time span
    indsst=np.arange(np.size(sst,2)) 
    nsst=np.size(indsst,0)
    
    # initialization of the matrice 2D space-time
    sstsqueezing=sst[indi,:,:]
    sstsqueezing=sstsqueezing[:,indj,:]
    sstsqueezing=sstsqueezing[:,:,indsst]
    data=sstsqueezing
    A=np.size(data,0);
    B=np.size(data,1);
    data=transpose(data,[2,0,1]);
    data=reshape(data,[nsst,-1])
    
    # eof
    pc, eig, eof = sparse.linalg.svds(data[:,:],neof) 
    #the scipy.sparse.linalg.svds function is not the same as scipy.linalg.svd
    #also notice that the eigenvalues are given as a 1Xn matrix instead of a nXn diagonal matrix, so the following line is needed:
    
    eig=np.diag(eig)
    # notice that eigenvalues are in increasing order for this function, which could be a problem if neof>1 and needs rearrangement.
    
    eof=reshape(eof,(neof, xeof.shape[0], yeof.shape[0])); 
    
    eofreshape=transpose(eof,[0,2,1]) #get the orientation of the data correct
    eofreshape=eofreshape[::-1,::-1,:] 
    
    
    eof=eofreshape #ordreing in this data matrix would be [# of eigenvector,ydata,xdata]. Refer to this for future steps.
    
    # normalization of variance
    tracedata2=np.trace(np.dot(data,transpose(data))); #here matrix multiplication is used as dot(*,*)
    
    var=np.zeros(neof) #define var before assign its values
    for n in range(0,neof):
        var[n]=eig[n,n]**2/tracedata2;
    
    var=var[::-1] #correct the order of the variations
    
    # normalization of principal component
    for n in range(0,neof):
        pc[:,n]=pc[:,n]/np.std(pc[:,n], ddof=1);    
    
    return eof,pc,var,xeof,yeof


# this defines the regression function
# here some careful tests needs to be done to see if pc and eof are transposed in the right way.

def sst_regression(sst,pc,eof,xsst,ysst):
    
    #regression
    neof=np.size(pc,1);
    reg=np.zeros([neof,np.size(sst,1),np.size(sst,0)]);
    
    #sst data is in the index form of [x,y,t]
    
    for n in range(0,neof):
        for i in range(0,np.size(sst,0)):
            for j in range(0,np.size(sst,1)):
                #linalg.lstsq is used to give the same value as matlab backslash operator do for two column matrices
                reg[n,j,i]=linalg.lstsq(pc[:,n].reshape(-1,1),np.squeeze(sst[i,j,:]))[0]; #L not fixed for indices
            
            
    reg=transpose(reg,[0,2,1])  #This is to change reg in the format from [t,y,x] into [t,x,y] for plotting purpose
    reg=reg[::-1,:,:]  #This is to order the eigenvalues from largest to smallest

    #check the sigh for eigenfunctions (does not work properlly yet)
    condj=(ysst>=-4)*(ysst<=4)
    condj=np.squeeze(condj)
    indj=np.arange(np.size(condj))[condj]
    
    condi=(xsst>=251)*(xsst<=279)
    condi=np.squeeze(condi)
    indi=np.arange(np.size(condi))[condi]
    
    for n in range(0,neof):
        sumchecker=sum((sum(sum(reg))))
        print (sumchecker,n)
        if sumchecker<=0: 
            pc[n,:]=-pc[n,:];
            reg[n,:,:]=-reg[n,:,:];
            eof[n,:,:]=-eof[n,:,:];
            
            
    return reg,pc,eof

# define plotting function for pc
def plotpc(pc,timesst):
    
    neof=np.size(pc,1);
    plt.figure()
    plt.plot(transpose(timesst),pc)
    plt.grid()
    plt.legend((np.arange(neof)+1).astype(str))
    plt.title('pc time series')
    plt.axis('image')
    plt.show()

# define plotting for eigenfunctions
def ploteof(reg,masksst,xsst,ysst,pc,timesst,var,latmin,latmax,lonmin,lonmax):
    
    
    neof=np.size(reg,0);
    print (shape(reg),'regshape')
    
    pcplot=pc[:,::-1]
    
    # modify the dimension of xsst and ysst
    xplot=np.squeeze(xsst)
    yplot=np.squeeze(ysst)
    mm=np.arange(neof)*0.0
    for n in range(0,neof):
        plt.figure()
        plt.subplot(2,1,1)
        mm[n]=max(np.amax(reg[n,:,:]),np.amax(-reg[n,:,:])) #This is to extract the region that will be desplayed
        tmp=reg[n,:,:]; #L

        
        #this is to maskout land areas in plotting
        for i in range (np.size(tmp,0)):
            for j in range (np.size(tmp,1)):
                if masksst[i,j]==1:
                    tmp[i,j]=+100
        
        
        plt.pcolormesh(xplot,yplot,tmp.T,vmin=-1.1*mm[n],vmax=1.1*mm[n]) #this is the function used for pcolor
        
        plt.xlim(np.amin(xplot),np.amax(xplot))
        plt.ylim(np.amin(yplot),np.amax(yplot))
        plt.axis('image')
        
        plt.colorbar()
        plt.title('eof var=%2.0f'%(var[n]*100))
        
        plt.subplot(2,1,2)
        plt.plot(transpose(timesst),pcplot[:,n]) 
        plt.grid()
        plt.axis('image')
        plt.title('pc%1.0f'%(n+1))
        

#plot the trend part (linear of constant) of SST
def plottrend(ssttrend,sstmean,xsst,ysst,masksst):
    
    xplot=np.squeeze(xsst)
    yplot=np.squeeze(ysst)
    
    plt.figure()
    tmp=ssttrend
    mm=max(np.amax(tmp),np.amax(-tmp))
    plt.pcolormesh(xplot,yplot,tmp.transpose(),vmin=-1.1*mm,vmax=1.1*mm)
    plt.xlim(np.amin(xplot),np.amax(xplot))
    plt.ylim(np.amin(yplot),np.amax(yplot))
    plt.axis('image')
    plt.title('trend in deg/year')
    plt.colorbar()
    

#plot time-longitude diagram at the equator for both sst and sst_anomalies
def plothoev(sst,timesst,xsst,ysst,masksst,anom,yearly):
    
    Xgrid=np.arange(2,362,2)
    data=np.zeros([np.size(sst,2),np.size(Xgrid)])
    Ygrid=[0]
    
    if yearly:
        if anom==0:
            timesstplot=timesst[:,np.arange(12,np.size(timesst),12)];
        else:
            timesstplot=timesst
    else:
        timesstplot=timesst;
    
    for n in range(np.size(sst,2)):
        datafunction=interpolate.interp2d(xsst,ysst,sst[:,:,n].transpose())   #2D interpolation from scipy. Mean value may need transpose.
        data[n,:]=datafunction(Xgrid,Ygrid)
        if anom:
            cax=[-2,2]
        else:
            cax=[23,30]
    
    plt.figure()
    plt.pcolormesh(Xgrid,np.squeeze(timesstplot),data,vmin=cax[0],vmax=cax[1])
    plt.xlim(np.amin(Xgrid),np.amax(Xgrid))
    plt.ylim(np.amin(timesstplot),np.amax(timesstplot))
    plt.title(r'lat=%3.2f$^o$'%(Ygrid))
    plt.axis('image')
    plt.colorbar()


#plot the average value for sst
def plotmean(sstmean,xsst,ysst,masksst):
    
    xplot=np.squeeze(xsst)
    yplot=np.squeeze(ysst)
    
    plt.figure()
    tmp=sstmean
    mma=np.amax(tmp)
    mmi=np.amin(tmp)
    
    plt.pcolormesh(xplot,yplot,tmp.transpose(),vmin=mmi,vmax=mma)
    plt.xlim(np.amin(xplot),np.amax(xplot))
    plt.ylim(np.amin(yplot),np.amax(yplot))
    plt.title('mean temperature')
    plt.axis('image')
    plt.colorbar()


#plot the seasonal shape of SST field
def plotseas(sstseas,xsst,ysst,masksst):
    
    xplot=np.squeeze(xsst)
    yplot=np.squeeze(ysst)
    
    plt.figure()
    for n in range(6):
        plt.subplot(3,2,n+1)
        tmp=sstseas[:,:,2*n]
        mm=7.5
        plt.pcolormesh(xplot,yplot,tmp.transpose(),vmin=-1.1*mm,vmax=1.1*mm)
        plt.xlim(np.amin(xplot),np.amax(xplot))
        plt.ylim(np.amin(yplot),np.amax(yplot))
        plt.title('month %d - seasonal cycle'%(2*n+1))
        plt.axis('image')
        plt.colorbar()


########################################################
# Here begins the main part of the program
########################################################

# load the data in python way with selections
if yearly:
    matdata = io.loadmat('SST_1870-present.mat')
if not (yearly):
    matdata = io.loadmat('SST_1976-2005.mat')

print ('data loaded')

masksst=matdata['masksst']
sst=matdata['sst']
timesst=matdata['timesst']
xsst=matdata['xsst']
ysst=matdata['ysst']
# reduce size 
#[sst,xsst,ysst,masksst]=sst_redsize(sst,xsst,ysst,masksst,pow2);

print shape(timesst)

if hoev_pl:
    anom=0;
    plothoev(sst,timesst,xsst,ysst,masksst,anom,yearly)
    anom=1;


# Time-filtering
#---------------
sst,timesst,facms,ssttrend,sstmean,sstseas=sst_tmfilt(sst,timesst,yearly,trend,seas,nym);

# added plot
#-----------------------------

if hoev_pl:
    plothoev(sst,timesst,xsst,ysst,masksst,anom,yearly)

if mean_pl:
    plotmean(sstmean,xsst,ysst,masksst);

if trend_pl:
    plottrend(ssttrend,xsst,ysst,masksst);

if seas_pl:
    plotseas(sstseas,xsst,ysst,masksst);


# compute EOFs
#-------------

[eof,pc,var,xeof,yeof]=sst_eof(sst,xsst,ysst,latmin,latmax,lonmin,lonmax,neof);
[reg,pc,eof]=sst_regression(sst,pc,eof,xsst,ysst);


# plot eofs and pc
#%------------------

ploteof(reg,masksst,xsst,ysst,pc,timesst,var,latmin,latmax,lonmin,lonmax)
plotpc(pc,timesst)

print('end of the run')

plt.show()