#!/usr/bin/env python

"""
Usage notes: for the sake of locating everything in one file I have put all of the Python code needed to generate figures
3, 4, 5, 7, 8, and 9 in this file. For the other figures see `figures.jl`.
It's a little messy but should run (and generate all figures) if you run the whole file, 
and many parts will work in isolation if you copy and paste from between the comments for the relevant section.
This file also includes code to generate the table published in the manuscript with the best fit parameters and their uncertainties.
"""

################ GENERAL SETUP ################
from PTemceeFit import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys,math
import IPython

font = {'family' : 'DejaVu Serif',
    'weight' : 'normal',
    'size'   : 16}
plt.rc('font', **font) #set all plot attribute defaults

def pltFormatter(fig,axList,**kwargs):
    for ax in axList:
        ax.minorticks_on()
        ax.grid(b=True,which="major",alpha=0.5)
        ax.grid(b=True,which="minor",alpha=0.3)
        ax.set_xticks([2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22])
        legend=kwargs.get("legend")
        if legend != None:
            legend.get_frame().set_edgecolor('black') 

def trackPercent(place,totalLength,strLen): #percent output tracker
    percent = place/totalLength*100
    string="{:.2f} % complete".format(percent)
    sys.stdout.write("\r") #this "moves the cursor" to the beginning of the I0 line
    sys.stdout.write(" "*strLen) #this "clears" whatever was on the line last time by writing whitespace
    sys.stdout.write("\r") #move the cursor back to the start again
    sys.stdout.write(string) #display the current percent we are at
    sys.stdout.flush() #flush finishes call to print() (this is like what's under the hood of print function)
    strLen=len(string) #return the new string length for next function call
    return strLen

getProfiles = DiskWind.getProfiles

################ READ IN FIT RESULTS ################

summit0 = readPickle('jPyPTEmceeVar0.p'); summit1 = readPickle('jPyPTEmceeVar1.p'); summit2 = readPickle('jPyPTEmceeVar2.p'); summit3 = readPickle('jPyPTEmceeVar3.p')
summit4 = readPickle('jPyPTEmceeVar4.p'); summit5 = readPickle('jPyPTEmceeVar5.p'); summit6 = readPickle('jPyPTEmceeVar6.p'); summit7 = readPickle('jPyPTEmceeVar7.p')
summit8 = readPickle('jPyPTEmceeVar8.p'); summit9 = readPickle('jPyPTEmceeVar9.p'); summit10 = readPickle('jPyPTEmceeVar10.p'); summit11 = readPickle('jPyPTEmceeVar11.p')
summit12 = readPickle('jPyPTEmceeVar12.p'); summit13 = readPickle('jPyPTEmceeVar13.p')
allResults = [np.concatenate((summit0[i],summit1[i],summit2[i],summit3[i],summit4[i],
                              summit5[i],summit6[i],summit7[i],summit8[i],summit9[i],
                             summit6[i],summit7[i],summit8[i],summit9[i],summit10[i],
                             summit11[i],summit12[i],summit13[i]),axis=1) for i in range(len(summit0))]
flat_samples,pos,prob,lhood,acor = allResults

################ GET FIT RESULTS TABLE ################

from IPython.display import display, Math
labels=["i","rBar","MFac","rFac","f1","f2","f3","f4","PA","scale","cenShift"] #units are degrees, r_s, 3e8*Msun, rMax/rMin, unitless (0<f1<1), unitless (0<f2<1), unitless (0<f3<1), unitless (0<f4<1), degrees (but +90 offset from standard PA convention), maxFlux(model)/maxFlux(data), microns
lowiMask=flat_samples[0,:,0] < 45 #low inclination only
avgParams = []
for i in range(len(labels)):
    mcmc = np.percentile(flat_samples[0,:, i], [16, 68, 84]) #PT way
    #mcmc = np.percentile(flat_samples[0,lowiMask, i], [16, 68, 84]) #PT way -- un comment to see low inclination only results 
    avgParams.append(mcmc[1])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.5f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    display(Math(txt))

################ PLOT ALL PHASES (Figure 7) ################

def plotPhases(data,θList,mα=0.1):
    indx=[0,1,2,6,7,8,12,13,14,18,19,20]; oindx=[3,4,5,9,10,11,15,16,17,21,22,23]
    fig,axs = plt.subplots(nrows=4,ncols=6,figsize=(20,12),sharex=True,sharey=True,facecolor="white")
    place = 0; strLen = 0; N = len(θList)
    for θ in θList:
        i,rMin,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift = θ
        λCen=2.172; ν = (data[0]-λCen)/λCen*3e5; ν = data[0] #not really ν anymore but testing wavelength plotting
        ν,line,phaseList = getProfiles(np.array(θ,dtype=float),data)
        ind = 0
        for ax in axs.reshape(-1):
            onoff = "on" if ind in indx else "off"
            ax.errorbar(ν,data[4][ind],yerr=data[5][ind],marker="o",ms=3,markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
            ax.fill_between(ν,data[4][ind]-data[5][ind],data[4][ind]+data[5][ind],color="dodgerblue",alpha=0.5)
            ax.plot(ν,phaseList[ind],lw=2,c="crimson",alpha=mα)
            if ind in indx:
                ax.set_title(r"index = {0} ($\textbf{{off}}$)".format(ind)) #off axis, changing naming scheme
            else:
                ax.set_title(r"index = {0} (on)".format(ind))
#             ax.set_xlabel("velocity [km/s]"); ax.set_ylabel("phase [deg]")
            ind += 1
            ax.grid(visible=True,which="major",alpha=0.4)
            ax.minorticks_on()
        strLen = trackPercent(place,N,strLen)
    fig.supxlabel("$\lambda$ [$\mu$m]"); fig.supylabel(r"Differential phase [$^\circ$]")
    plt.tight_layout()
    return fig,axs


fac = 2
SMALL_SIZE = 8*fac
MEDIUM_SIZE = 10*fac
BIGGER_SIZE = 12*fac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
mpl.rcParams['text.usetex']=True
data = readPickle("3c273_juljanmarmay_append_gilles_specirf_wide_v6.p")
fig,axs=plotPhases(data,[avgParams],1); #avgParams from a few cells down
fig.savefig("avgParamsPhases.png")

################ MAKE CORNER PLOT (Figure 9) ################

import corner
import matplotlib as mpl
import pickle
fac = 3
SMALL_SIZE = 6*fac
MEDIUM_SIZE = 8*fac
BIGGER_SIZE = 10*fac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
mpl.rcParams['xtick.major.pad']='2'
mpl.rcParams['ytick.major.pad']='2'
mpl.rcParams["axes.titlepad"] = 10
mpl.rcParams['axes.labelpad']=10
mpl.rcParams['font.family']="DejaVu Serif"
mpl.rcParams['text.usetex']=True

flat_samples_T0 = np.copy(flat_samples[0,:,:])
flat_samples_T0[:,8] -= 90. #make units match paper
flat_samples_T0[:,2] *= 30. #make units match paper

n = len(flat_samples[0,:,:])
samples_4_corner = flat_samples_T0[:,0:-2] #get rid of lambda and n because they are the least "physical" of the parameters
labels = [r"$i$", r"$\bar{r}$", r"$M_{BH}$", r"$r_{\textrm{fac}}$",r"$f_1$",r"$f_2$",r"$f_3$",r"$f_4$",r"$\theta_{\textrm{PA}}$",r"$n$",r"$\Delta \lambda_c$"]
labels_4_corner = labels[:-2]
figure = corner.corner(samples_4_corner, labels=labels_4_corner,
                       quantiles=[0.16, 0.5, 0.84],size=(40,40),fill_contours=True,top_ticks=False,
                       show_titles=True,labelpad=0.1)

lhood_rescale = lhood[0,:] 
maxInd = np.argmax(lhood_rescale)
flat_samples_T0 = flat_samples[0,:,:]
θBest = flat_samples_T0[maxInd]

ndim = len(labels_4_corner)
axes = np.array(figure.axes).reshape((ndim, ndim))
# Loop over the diagonal
for i in range(ndim):
    ax = axes[i, i]
    ax.axvline(θBest[i], color="dodgerblue")
    
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
        ax.axvline(θBest[xi], color="dodgerblue")
        ax.axhline(θBest[yi], color="dodgerblue")
        ax.plot(θBest[xi], θBest[yi], "darkblue")
        ax.xaxis.label.set_size(20)
        ax.yaxis.label.set_size(20)

figure.savefig("corner.svg",bbox_inches="tight")

################ SHOW MODEL AND DATA CENTROIDS (Figure 4) ################

from scipy.optimize import curve_fit
import matplotlib as mpl
plt.close("all")
plt.rcParams['text.usetex']=False
fontFac=2
SMALL_SIZE = 8*fontFac
MEDIUM_SIZE = 10*fontFac
BIGGER_SIZE = 12*fontFac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

from scipy.optimize import minimize, least_squares
λData = data[0]
UData = data[1]; VData = data[2]
ϕData = data[4]
IData = data[3]
x0 = np.array([0,0],dtype=int)
centroids=[]; strLen = 0

def baselineFit(x,y,UData,VData,IData,ϕData,λi,onlyOn=True):
    z = 0
    indx=[0,1,2,6,7,8,12,13,14,18,19,20]
    for i in range(len(UData)):
        if onlyOn:
            if i in indx:
                Bi = i
                z += f2d([UData[Bi],VData[Bi],IData[λi]],x,y,ϕData[Bi,λi]) #wavelength,baseline indices
        else:
            Bi = i
            z += f2d([UData[Bi],VData[Bi],IData[λi]],x,y,ϕData[Bi,λi]) #wavelength,baseline indices
    return z

def f2min(x0,args=np.zeros(6)):
    UData,VData,IData,ϕData,λi,onlyOn = args
    z = baselineFit(x0[0],x0[1],UData,VData,IData,ϕData,λi,onlyOn)
    return z


for i in range(len(λData)):
    ϕ = ϕData[:,i]
    I = IData[i]
    tmp = [0.,0.]
    res = minimize(f2min,x0,args=[UData,VData,IData,ϕData,i,False],method="Nelder-Mead",options={'xatol': 1e-12})
    if res.success == True:
        centroids.append(res.x)
    else:
        print("problem at i = {0}".format(i))
        break
    strLen = trackPercent(i+1,len(λData),strLen)

def IWrapper(params,data,nr=1024,nϕ=2048): #this is ~3x as fast as python version!
    i,r̄,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift = params; scale = 1.; windWeight = 1.; γ = 1.; A0 = 1.; τ = 10. #some parameters fixed for now
    α,β,r,νloc,ϕ,sini,cosi,dA,rMin,rMax = DiskWind.setup(i,nr,nϕ,r̄,rFac,γ)
    I,γ,A0,τ = DiskWind.getIntensity(r,ϕ,sini,cosi,rMin,rMax,γ,A0,τ,f1=f1,f2=f2,f3=f3,f4=f4)
    return I,α,β,νloc,dA
I,α,β,ν,dA = IWrapper(avgParams,data)

BLRAng = θBest[2]*3e8*2e33*6.67e-8/9e20/548/3.09e24
fig,axs=plt.subplots(nrows=1,ncols=2,figsize=(20,8),facecolor="white")
channelCenters = λData[mask]
channels = []
for i in range(len(channelCenters)):
    if i == 0:
        Δ = np.abs(channelCenters[1]-channelCenters[0])
        channels.append(channelCenters[0] - Δ/2) #left edge
    elif i == len(channelCenters)-1:
        Δ = np.abs(channelCenters[-1]-channelCenters[-2])
        channels.append(channelCenters[-1] - Δ/2) #left edge
        channels.append(channelCenters[-1] + Δ/2) #right edge
    else:
        Δ = np.abs(channelCenters[i]-channelCenters[i-1])
        channels.append(channelCenters[i]-Δ/2) #left edge
        
        
λmod = (2.172+avgParams[-1])/ν

rot = -(avgParams[-3])/180*np.pi

xCD = np.array([centroid[0] for centroid in centroids]); yCD = np.array([centroid[1] for centroid in centroids])
mask = ((λData > 2.166) & (λData < 2.18))
conv = 206264806719 #rad to μas

for i in range(2):
    ax = axs.flatten()[i]; Iloop = I if i==0 else IData
    cmap = plt.cm.RdBu_r
    bounds = np.linspace(2.166,2.180,8)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    if i == 0:
        xC,yC = np.zeros(len(channels)-1),np.zeros(len(channels)-1)
        for j in range(len(channels)-1):
            if j == 0 or j == len(channels)-2:
                maskMod = np.where((λmod>=channels[j]) & (λmod <= channels[j+1]))
            else:
                maskMod = np.where((λmod>channels[j]) & (λmod <= channels[j+1]))

            X = np.sum(Iloop[maskMod]*α[maskMod]*dA[maskMod])/np.sum(Iloop[maskMod]*dA[maskMod]); Y = np.sum(Iloop[maskMod]*β[maskMod]*dA[maskMod])/np.sum(Iloop[maskMod]*dA[maskMod])
            xC[j],yC[j] = X,Y
    
        xCRot = (np.cos(rot)*xC-np.sin(rot)*yC); yCRot = np.sin(rot)*xC+np.cos(rot)*yC
        s=ax.scatter(-xCRot*BLRAng/4.848e-12,yCRot*BLRAng/4.848e-12,c=channelCenters,cmap=cmap,norm=norm,
                edgecolors="k",s=150)#,c=channels,cmap="RdBu_r")#,c=channels*1e3,cmap="RdBu")
        def f(x,m,b):
            return m*x + b
        fit,cov = curve_fit(f,-xCRot*BLRAng/4.848e-12,yCRot*BLRAng/4.848e-12)
        xline=-xCRot*2*BLRAng/4.848e-12
        
        ax.plot(xline,f(xline,*fit),lw=2,ls="--",c="k")
        
    else:
        ax.plot(xCD[mask]*conv,yCD[mask]*conv,c="k",ls="-",lw=0.5)
        ax.annotate("PA Jet = 222°",(15,20))
        s=ax.scatter(xCD[mask]*conv,yCD[mask]*conv,c=λData[mask],cmap=cmap,norm=norm,
                edgecolors="k",s=200,zorder=10)#,c=channels,cmap="RdBu_r")#,c=channels*1e3,cmap="RdBu")
        fit,cov = curve_fit(f,-xCRot*BLRAng/4.848e-12,yCRot*BLRAng/4.848e-12)
        ax.plot(-xCRot*BLRAng/4.848e-12,f(-xCRot*BLRAng/4.848e-12,*fit),lw=2,ls="--",c="k")
        
    
    ax.set_xlim(-31,31); ax.set_ylim(-31,31)
    ax.set_aspect("equal")
    jetx,jety = np.zeros(100),np.linspace(-45,45,100)
    jetθ = -222/180*np.pi #rotating from N into + RA
    jetX = -np.sin(jetθ)*jetx + np.cos(jetθ)*jety; jetY = np.sin(jetθ) + np.cos(jetθ)*jety
    if i>0:
        ax.plot(jetX,jetY,color='k',lw=2)
    ax.invert_xaxis()
    
    fig.colorbar(s, ax=ax,fraction=0.05,pad=0.,label="wavelength (μm)",spacing='proportional',ticks=bounds,boundaries=bounds)
    title = "Data".format(i,fList[i-1]) if i>0 else "Model result"
    ax.set_title(title)
    ax.set_facecolor("white")
    ax.minorticks_on()
    ax.tick_params(which="both",direction="in",top=True,right=True,labeltop=False,labelright=False)
    ax.set_xlabel("ΔRA (μas)")
    ax.set_ylabel("ΔDec (μas)")

fig.savefig("centroid_comparison.pdf")

################ SHOW BASELINES (FIGURE 8) ################

fontFac=3
SMALL_SIZE = 8*fontFac
MEDIUM_SIZE = 10*fontFac
BIGGER_SIZE = 12*fontFac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
fig,ax=plt.subplots(figsize=(14,14),facecolor='white')
indx=[0,1,2,6,7,8,12,13,14,18,19,20]; oindx=[3,4,5,9,10,11,15,16,17,21,22,23]
for i in range(len(UData)):
    if i in indx:
        ax.scatter(UData[i],VData[i],color='k',s=60*fontFac)
        ax.scatter(-UData[i],-VData[i],color='k',s=60*fontFac)
    else:
        ax.scatter(UData[i],VData[i],color='grey',s=60*fontFac)
        ax.scatter(-UData[i],-VData[i],color='grey',s=60*fontFac)

jetx,jety = np.zeros(100),np.linspace(-100,100,100)
jetθ = -222/180*np.pi #rotating from N into + RA
jetX = -np.sin(jetθ)*jetx + np.cos(jetθ)*jety; jetY = np.sin(jetθ) + np.cos(jetθ)*jety
ax.plot(jetX,jetY,color='k',lw=2)
ax.annotate("PA Jet = 222°",(40,50))
ax.set_xlim(60,-60)
ax.set_ylim(-60,60)
ax.tick_params(which="both",direction="in",labeltop=False,labelright=False)

ax2 = ax.twiny()
ax3 = ax.twinx()
mticks = np.array([100,50,0,-50,-100])
new_tick_locations = mticks/2.172

ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels([str(tick) for tick in mticks])
ax3.set_ylim(ax.get_ylim())
ax3.set_yticks(new_tick_locations)
ax3.set_yticklabels([str(tick) for tick in mticks])
ax2.tick_params(which="both",direction="in",labeltop=True,labelright=False)
ax3.tick_params(which="both",direction="in",labeltop=False,labelright=True)
ax3.set_ylabel("projected baseline length [m]",labelpad=25)
ax2.set_xlabel("projected baseline length [m]",labelpad=25)
ax.set_xlabel(r"u [M$\lambda$]")
ax.set_ylabel(r"v [M$\lambda$]")
ax.plot([70,-70],[0,0],color='grey',alpha=0.5)
ax.plot([0,0],[70,-70],color='grey',alpha=0.5)
ax.minorticks_on()
ax2.minorticks_on()
ax3.minorticks_on()
fig.tight_layout()
fig.savefig("baselines.svg")

################ SHOW MODEL LP AND PHASES FROM FITS (FIGURE 3) ################

def plotParamsStacked(data,θList,mα=0.1,cList=[],labelList=[],lwList=[],lsList=[],last=True):
    λCen=2.172; ν = (data[0]-λCen)/λCen*3e5; ν = data[0] #not really ν anymore but testing wavelength space
    indx=[0,1,2,6,7,8,12,13,14,18,19,20]; oindx=[3,4,5,9,10,11,15,16,17,21,22,23]
    fig,ax1 = plt.subplots(figsize=(14,12),facecolor="white")
    ax2 = ax1.twinx()
    
    meanOn = np.mean(np.array(data[4])[indx],axis=0)
    onErr = np.sqrt(1/np.sum((np.array(data[5])[indx])**(-2),axis=0)) #see weighted average here: http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
    
    dodgerBlue=(.12,0.56,1.00)
    ax1.errorbar(ν,data[3]+1.,yerr=data[6],marker="o",ms=4,label="3C 273",markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
    ax2.errorbar(ν,meanOn,yerr=onErr,label="3C 273",marker="o",ms=4,markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
    ax1.fill_between(ν,data[3]-data[6]+1.,data[3]+data[6]+1.,color=dodgerBlue,alpha=0.5)
    ax2.fill_between(ν,meanOn-onErr,meanOn+onErr,color=dodgerBlue,alpha=0.5)
    strLen = 0; place = 1; N = len(θList)
    for θ in θList:
        try:
            i,rMin,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift = θ
            λCen=2.172+cenShift; #ν = (data[0]-λCen)/λCen*3e5
        except:
            i,rMin,Mfac,rFac,f1,f2,f3,pa,scale,cenShift = θ
            f4 = np.copy(f3); f3 = np.copy(f2); f2 = np.copy(f1)
            θ = np.array([i,rMin,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift])
        try:
            ν,line,phaseList = getProfiles(np.array(θ,dtype=float),data)
        except:
            print("""problem with: i = {0:.2f}, rMin = {1:.2f}, MFac = {2:.2f}, rFac = {3:.2f}, f1 (sin^2 1) = {4:.2f}, f2 (sin^2 2) = {5:.2f}, f3 (sin*cos)= {6:.2f}, f4 (cos^2)= {7:.2f} \
pa = {8:.2f}, scale = {9:.2f}, cenShift = {10:.4f}""".format(*θ))
        phase = np.mean(np.array(phaseList)[indx],axis=0); #phaseo = np.mean(np.array(phaseList)[oindx],axis=0)
        label = "Model fits".format(N) if place == N else ""
        label = labelList[place-1] if len(labelList)>0 else label
        c = cList[place-1] if len(cList)>0 else "crimson"
        mαLoc = mα[place-1] if type(mα) is list else mα
        lw = lwList[place-1] if len(lwList)>0 else 2; 
        ls = lsList[place-1] if len(lsList)>0 else "-"

        ax1.plot(ν,line+1.,label=label,lw=lw,c=c,alpha=mαLoc,ls=ls)
        ax2.plot(ν,phase,label=label,lw=lw,c=c,alpha=mαLoc,ls=ls)
        if last == True and place >= N-1:
            label = "Line center" 
            ax1.vlines(λCen,0.2,1.8,label=label,colors=c,ls=ls,lw=lw,alpha=mαLoc)
        strLen = trackPercent(place,N,strLen); place+=1
    ax1.set_title("Line and phase profile comparison")
    ax1.set_ylim(0.2,1.8)
    ax2.set_ylim(-0.4,1.2)
    #ax2.set_xlabel("Velocity [km/s]")
    ax1.set_xlabel("λ [μm]")
    ax1.set_ylabel("Normalized flux")
    ax2.set_ylabel("Δϕ [deg]")
    l = ax2.legend(loc='upper left')
    pltFormatter(fig,[ax1],legend=l)
    ax2.set_xticks([2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22])
    ax2.set_xlim(2.126,2.21)
    ax1.set_xlim(2.126,2.21)
    fig.tight_layout()
    return fig,ax1,ax2

import matplotlib as mpl
mpl.rcParams['xtick.major.pad']='2'
mpl.rcParams['ytick.major.pad']='2'
mpl.rcParams['axes.labelpad']=10
mpl.rcParams['font.family']="DejaVu Serif"
mpl.rcParams['text.usetex']=False
fontFac=3
SMALL_SIZE = 8*fontFac
MEDIUM_SIZE = 10*fontFac
BIGGER_SIZE = 12*fontFac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

data = readPickle("3c273_juljanmarmay_append_gilles_specirf_wide_v6.p")
labelList=["" for i in range(100)];labelList.append("Best fit"); #labelList.append("Best low inclination fit")
cList = ["crimson" for i in range(100)];cList.append("crimson"); #cList.append("blueviolet")
mαList = [0.05 for i in range(100)]; mαList.append(1.); #mαList.append(0.6)
lwList = np.ones(100); lwList = np.append(lwList,np.array([4,4]))
lsList = ["-" for i in range(100)]; lsList.append("--",)#;lsList.append("--")

allResults = [np.concatenate((summit0[i],summit1[i],summit2[i],summit3[i],summit4[i],
                              summit5[i],summit6[i],summit7[i],summit8[i],summit9[i],
                             summit6[i],summit7[i],summit8[i],summit9[i],summit10[i],
                             summit11[i],summit12[i],summit13[i]),axis=1) for i in range(len(summit0))]
flat_samples,pos,prob,lhood,acor = allResults
maxInd = np.argmax(lhood[0,:])
θBest = flat_samples[0,maxInd,:] #prob matches with pos, NOT with flat_samples (idk why we have that one seems redundant to save both?)
best_low_i = test[np.argmax(lhood[0,lowiMask]),:]

θList = flat_samples[0,np.random.randint(len(flat_samples[0]),size=100),:] #looks good!
θList=np.vstack((θList,θBest))
fig,ax1,ax2=plotParamsStacked(data,θList,mαList,cList,labelList,lwList,lsList)
fig.savefig('best_params_summit.pdf')

################ SHOW FINE-TUNING PROBLEM AT LOW i (FIGURE 5) ################

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'DejaVu Serif'
mpl.rcParams['mathtext.it'] = 'DejaVu Serif:italic'
mpl.rcParams['mathtext.bf'] = 'DejaVu Serif:bold'

mpl.rcParams['xtick.major.pad']='2'
mpl.rcParams['ytick.major.pad']='2'
mpl.rcParams['axes.labelpad']=10
mpl.rcParams['font.family']="DejaVu Serif"
mpl.rcParams['text.usetex']=False
fontFac=3
SMALL_SIZE = 8*fontFac
MEDIUM_SIZE = 10*fontFac
BIGGER_SIZE = 12*fontFac

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#avgParams[8] += 90. #only run once if correcting for offset in calculating avgParams!!!
θ_finetune = np.copy(best_low_i) 
θ_finetune[7] = 0.2
θ_finetune = [best_low_i,θ_finetune] 
labelList = [r"""$f_1$ = {0:.2f}, $f_2$ = {1:.2f},
$f_3$ = {2:.2f}, $f_4$ = {3:.2f}""".format(*best_low_i[4:8]),r"""$f_1$ = {0:.2f}, $f_2$ = {1:.2f},
$f_3$ = {2:.2f}, $f_4$ = {3:.2f}""".format(*θ_finetune[1][4:8])]
cList = ["crimson","purple"]
lwList = [4,4]
lsList = ["--","-"]
mαList = [0.8,0.8]
θ_finetune.reverse();lsList.reverse();cList.reverse();labelList.reverse();
fig,ax1,ax2=plotParamsStacked(data,θ_finetune,mαList,cList,labelList,lwList,lsList)
fig.savefig('finetune.pdf')