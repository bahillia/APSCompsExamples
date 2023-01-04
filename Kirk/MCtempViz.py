#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys,math,os
#os.chdir("/home/kirk/Documents/research/Dexter/Sp21")
from emceeFit import *
global getProfiles #= DiskWind.getProfiles
getProfiles = DiskWind.getProfiles



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
    if math.floor(percent)==69:
        string="{:.2f} % complete -- nice".format(percent)
    else:
        string="{:.2f} % complete".format(percent)
    sys.stdout.write("\r") #this "moves the cursor" to the beginning of the I0 line
    sys.stdout.write(" "*strLen) #this "clears" whatever was on the line last time by writing whitespace
    sys.stdout.write("\r") #move the cursor back to the start again
    sys.stdout.write(string) #display the current percent we are at
    sys.stdout.flush() #flush finishes call to print() (this is like what's under the hood of print function)
    strLen=len(string) #return the new string length for next function call
    return strLen

def plotParams(data,θList,mα=0.1,cList=[],labelList=[]):
    λCen=2.172; ν = (data[0]-λCen)/λCen*3e5; ν = data[0] #not really ν anymore but testing wavelength space
    indx=[0,1,2,6,7,8,12,13,14,18,19,20]; oindx=[3,4,5,9,10,11,15,16,17,21,22,23]
    fig,axd = plt.subplot_mosaic([['a','b','c']],figsize=(24,6),facecolor="white")
    ax1 = axd["a"]; ax2 = axd["b"]; ax3 = axd["c"]
    ax1.get_shared_x_axes().join(ax1,ax2); ax1.get_shared_x_axes().join(ax1,ax3)
    ax2.get_shared_y_axes().join(ax2,ax3); ax3.set_yticklabels([])
    dodgerBlue=(.12,0.56,1.00)
    ax1.errorbar(ν,data[3],yerr=data[6],marker="o",ms=4,label="3C 273",markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
    meanOn = np.mean(np.array(data[4])[indx],axis=0); meanOff = np.mean(np.array(data[4])[oindx],axis=0)
    onErr = np.sqrt(1/np.sum((np.array(data[5])[indx])**(-2),axis=0)); offErr = np.sqrt(1/np.sum((np.array(data[5])[oindx])**(-2),axis=0)) #see weighted average here: http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
    ax2.errorbar(ν,meanOn,yerr=onErr,label="3C 273",marker="o",ms=4,markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
    ax3.errorbar(ν,meanOff,yerr=offErr,label="3C 273",marker="o",ms=4,markerfacecolor="darkblue",markeredgecolor="darkblue",linewidth=0,elinewidth=0.5,capsize=1.5,capthick=0.5,color='darkblue')
    ax1.fill_between(ν,data[3]-data[6],data[3]+data[6],color=dodgerBlue,alpha=0.5)
    ax2.fill_between(ν,meanOn-onErr,meanOn+onErr,color=dodgerBlue,alpha=0.5)
    ax3.fill_between(ν,meanOff-offErr,meanOff+offErr,color=dodgerBlue,alpha=0.5)
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
        phase = np.mean(np.array(phaseList)[indx],axis=0); phaseo = np.mean(np.array(phaseList)[oindx],axis=0)
        label = "Disk wind model ({} samples)".format(N) if place == 1 else ""
        label = labelList[place-1] if len(labelList)>0 else label
        c = cList[place-1] if len(cList)>0 else "crimson"
        mαLoc = mα[place-1] if type(mα) is list else mα
        ax1.plot(ν,line,label=label,lw=3,c=c,alpha=mαLoc)
        ax2.plot(ν,phase,label=label,lw=3,c=c,alpha=mαLoc)
        ax3.plot(ν,phaseo,label=label,lw=3,c=c,alpha=mαLoc)
        label = "Line center" if place == 1 else ""
        ax1.vlines(λCen,0.,0.6,label=label,colors=c,ls="--",lw=2,alpha=mαLoc)
        ax2.vlines(λCen,-0.4,0.4,label="",colors=c,ls="--",lw=2,alpha=mαLoc)
        ax3.vlines(λCen,-0.4,0.4,label="",colors=c,ls="--",lw=2,alpha=mαLoc)
        strLen = trackPercent(place,N,strLen); place+=1
    ax1.set_title("Line profile comparison")
    ax2.set_title("Phase profile (mean off) comparison")
    ax3.set_title("Phase profile (mean on) comparison")
    #ax2.set_xlabel("Velocity [km/s]")
    ax2.set_xlabel("λ [μm]")
    ax1.set_ylabel("Normalized flux")
    ax2.set_ylabel("Phase [deg]")
    l = ax2.legend(loc='upper left')
    ax1.set_ylim(-0.02,0.7)
    ax2.set_ylim(-0.45,0.45)
    ax3.set_ylim(-0.45,0.45)
    pltFormatter(fig,[ax1,ax2,ax3],legend=l)
    fig.tight_layout()
    return fig,ax1,ax2,ax3

def main():
    font = {'family' : 'DejaVu Serif',
    'weight' : 'normal',
    'size'   : 16}
    plt.rc('font', **font) #set all plot attribute defaults

    SummitResults = readPickle('jPyPTEmceeVar.p')
    flat_samples,pos,prob,lhood,acor = SummitResults
    getProfiles = DiskWind.getProfiles
    summit0 = readPickle('jPyPTEmceeVar0.p'); summit1 = readPickle('jPyPTEmceeVar1.p'); summit2 = readPickle('jPyPTEmceeVar2.p'); summit3 = readPickle('jPyPTEmceeVar3.p')
    summit4 = readPickle('jPyPTEmceeVar4.p'); summit5 = readPickle('jPyPTEmceeVar5.p'); summit6 = readPickle('jPyPTEmceeVar6.p'); summit7 = readPickle('jPyPTEmceeVar7.p')
    summit8 = readPickle('jPyPTEmceeVar8.p'); summit9 = readPickle('jPyPTEmceeVar9.p')

    allResults = [np.concatenate((summit0[i],summit1[i],summit2[i],summit3[i],summit4[i],summit5[i],summit6[i],summit7[i],summit8[i],summit9[i]),axis=1) for i in range(len(summit0))]

    flat_samples,pos,prob,lhood,acor = allResults
    fontFac=1.5
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

    T1Walkers=flat_samples[0].reshape((25000,24,11))
    T3Walkers=flat_samples[2].reshape((25000,24,11))
    T6Walkers=flat_samples[3].reshape((25000,24,11))
    cList=["cyan","grey","crimson"]; labelList=["Hottest","Medium","Coolest"]; alphaList=[1.,0.8,0.6]
    for i in range(2500): #make 100 frame test
        if i%25==0:
            print("at i = {}/2500".format(i))
        fig,ax1,ax2,ax3=plotParams(data,[T6Walkers[i,0,:],T3Walkers[i,0,:],T1Walkers[i,0,:]],alphaList,cList,labelList)
        fig.savefig('tmpPlots3/test_{:04d}.svg'.format(i))
        plt.close("all")

main()
