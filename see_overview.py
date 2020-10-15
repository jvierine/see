#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import scipy.signal as ss
import stuffr
import re
import time
import os

#dnames=["/data0/2020.03.26/test1e6_4.04e6","/data1/2020.03.26/test1e6_4.04e6"]
dnames=["/data0/2020.10.15/test1e6_4.04e6","/data1/2020.10.15/test1e6_4.04e6"]
#dnames=["/data0/2019.11.25/test1e6_2.7e6","/data1/2019.11.25/test1e6_2.7e6"]
f=re.search(".*_(.*)",dnames[0]).group(1)
print(f)
cf=float(f)

ch=["cha","chb","chc","chd","che","chf","chg","chh"]

def plot_cc(ch0="cha",
            ch1="chc",
            n_avg=1,
            fft_len=1024*2,
            dt=10,
            rfi_rem=False,
            plot_phase=True,
            realtime=False,
            flim=[0,0],
            now=False):
    
#    fft_len=512*2
    n_fft=500*2
    sr=1e6
#    cf=4.2e6
    
    step_size=int(n.floor((dt*sr-fft_len)/n_fft))
    
    #print(step_size)
    
    d=drf.DigitalRFReader(dnames)
    print(d.get_properties(ch0))
    print(d.get_channels())
    b=d.get_bounds(ch0)
    print(b)
    print((b[1]-b[0])/1e6)
    
    S=n.zeros([fft_len,n_fft],dtype=n.complex64)
    
    b=d.get_bounds(ch0)
    
    cycle_num=int(n.floor(b[1]/sr/dt))-1
    print("%s %s %d cycles since 1970"%(ch0,ch1,cycle_num))

    if now:
        i0=(long(b[1]/(dt*sr))-1)*long(dt)*long(sr)
    else:
        i0=int(cycle_num*dt*sr)
    
    
#    wf=ss.hann(fft_len)
    wf=ss.chebwin(fft_len,100.0)    
    for si in range(n_fft):
        print("%d/%d"%(si,n_fft))
        for ai in range(n_avg):
            z0=d.read_vector_c81d(i0+step_size*si+fft_len*ai,fft_len,ch0)
            z1=d.read_vector_c81d(i0+step_size*si+fft_len*ai,fft_len,ch1)

            if rfi_rem:
                if ai==0:
                    S[:,si]=n.fft.fftshift(n.fft.fft(wf*z0)*n.conj(n.fft.fft(wf*z1)))
                else:
                    p0=n.abs(S[:,si])
                    p1=n.fft.fftshift(n.fft.fft(wf*z0)*n.conj(n.fft.fft(wf*z1)))
                    gidx=n.where(n.abs(p1) < p0)[0]
                    S[gidx,si]=p1[gidx]
            else:
                S[:,si]+=n.fft.fftshift(n.fft.fft(wf*z0)*n.conj(n.fft.fft(wf*z1)))
                
        
    tv=n.arange(n_fft)*step_size/sr
    fv=n.fft.fftshift(n.fft.fftfreq(fft_len,d=1/sr))+cf

    


#    S2=n.copy(n.abs(S))
 #   for si in range(S.shape[1]):
  #      S2[:,si]=(S2[:,si]-n.median(S2[:,si]))/n.median(n.abs(S2[:,si]-n.median(S2[:,si])))
    dB=10.0*n.log10(S)

    for si in range(S.shape[1]):
        dB[:,si]=dB[:,si]-n.nanmedian(dB[:,si])

    plt.close()
    plt.figure(figsize=(15,10))
    plt.clf()
    plt.pcolormesh(fv/1e6,tv,n.transpose(dB),vmin=-3,cmap="inferno")
    plt.title("%s - %s\n%s-%s"%(stuffr.unix2datestr(i0/sr),stuffr.unix2datestr(i0/sr+dt),ch0,ch1))
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Time (s)")    
    plt.colorbar()
    if flim[0] != 0:
        plt.xlim(flim)
    else:
        plt.xlim([ (cf-sr/2)/1e6, (cf+sr/2)/1e6])
    plt.tight_layout()
    if now:
        ofname="img/see-%1.2f-%s-%s.png"%(time.time(),ch0,ch1)
        plt.savefig(ofname)
        os.system("cp %s img/0latest.png"%(ofname))
        plt.pause(10)
        
    else:
        plt.savefig("img/xcp-%d-%s-%s.png"%(i0,ch0,ch1))
        plt.clf()
        plt.close()

    if plot_phase:
        plt.figure(figsize=(15,10))    
        plt.pcolormesh(fv/1e6,tv,n.transpose(n.angle(S)))
        plt.title("%s - %s\n%s-%s"%(stuffr.unix2datestr(i0/sr),stuffr.unix2datestr(i0/sr+dt),ch0,ch1))
        plt.colorbar()
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Time (s)")    
    
        plt.xlim([ (cf-sr/2)/1e6, (cf+sr/2)/1e6])
        plt.tight_layout()
        if now:
#            plt.savefig("img/see-%1.2f.png"%(time.time()))
            plt.pause(10.0)
#            plt.savefig("img/see-%1.2f.png"%(time.time()))
        else:
            plt.savefig("img/xca-%d-%s-%s.png"%(i0,ch0,ch1))
            plt.clf()
            plt.close()
    
def xc_pngs():
    for c in ch:
        plot_cc(ch0="cha",ch1=c)
        
def xc_rt(dt=10,n_avg=1,plot_phase=False,realtime=False,flim=[0,0],fft_len=1024,ch0="chc",ch1="chc"):
    plot_cc(ch0=ch0,ch1=ch1,dt=dt,n_avg=n_avg,now=True,rfi_rem=False,fft_len=fft_len,plot_phase=plot_phase,realtime=realtime,flim=flim)    

def realtime_updates():
    # sweep
    plt.figure(figsize=(15,10))
#    while True:
#       xc_rt(dt=10,n_avg=4,plot_phase=False,fft_len=1024,realtime=True,flim=[3.9,4.4])
    xc_rt(dt=5*60,n_avg=8,fft_len=2*2*2*1024,plot_phase=False,realtime=True,flim=[3.7,4.4])        
    #time.sleep(10)
        
realtime_updates()

# one-off sweep
#xc_rt(dt=5*60,n_avg=8)
# pulse
#xc_rt(dt=10,n_avg=1)
