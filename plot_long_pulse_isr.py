#!/usr/bin/env python


import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import stuffr
import scipy.signal as ss
import h5py
import time

import lpf_filter as lpf

# Juha is the best, best at radar and radio! Everybody knows it

def decode_all(dir_name=["/data0/2019.10.24/test200e3_7.953e6","/data1/2019.10.24/test200e3_7.953e6"],
               out_dir_name="/data0/magrad_24.10",
               idx0=None,
               ch="cha",
               n_ipp=1000,
               ipp=20000,
               n_fft=100,
               offset=160,     # tx on offset
               plot=True,
               sr=200e3):

    d=drf.DigitalRFReader(dir_name)
    b=d.get_bounds(ch)
    n_figs=int(n.floor((b[1]-b[0])/(ipp*n_ipp)))
    for i in range(n_figs):
        try:
            decode_channel(idx0=b[0]+i*n_ipp*ipp,n_ipp=n_ipp)
        except:
            print("missing data")

def decode_channel(dir_name=["/data0/2019.10.24/test200e3_7.953e6","/data1/2019.10.24/test200e3_7.953e6"],
                   out_dir_name="/data0/magrad_24.10",
                   idx0=None,
                   ch="cha",
                   n_ipp=500,
                   ipp=20000,
                   n_fft=100,
                   offset=160,     # tx on offset
                   plot=True,
                   sr=200e3):

    sr_factor=1e6/sr
    
    d=drf.DigitalRFReader(dir_name)
    b=d.get_bounds(ch)    
    md=d.read_metadata(b[0],b[0]+1,ch)
    for mk in md:
        sr=md[mk]["sample_rate_numerator"]

    if idx0 == None:
        i0=b[0]
    else:
        i0=idx0
    
    print(i0)
    n_rg=ipp/n_fft
    S=n.zeros([n_rg,n_fft])

    WS=n.zeros([n_rg,n_fft])    

    z=d.read_vector_c81d(i0,n_ipp*ipp,ch)
    wf=ss.hann(n_fft)
    
    fvec=n.fft.fftshift(n.fft.fftfreq(n_fft,d=1/float(sr)))/1e3
    rvec=n.arange(n_rg,dtype=n.float64)*n_fft*3e8/sr/1e3/2.0
    ipp_idx=0
    for i in range(n_ipp/10):
        S0=n.zeros([n_rg,n_fft])
        W0=n.zeros([n_rg,n_fft])
        
        for j in range(10):
            
            for ri in range(n_rg):
                zin=z[ipp_idx*ipp + ri*n_fft + n.arange(n_fft)]
                S00=n.abs(n.fft.fftshift(n.fft.fft(wf*zin)))**2.0
                W0[ri,:]+=S00**2.0
                S0[ri,:]+=S00
            ipp_idx+=1                

        W0=W0-S0**2.0
        S+=S0/W0
        WS+=1.0/W0
        

    S=S/WS

    for si in range(S.shape[1]):
        S[:,si]=S[:,si]/n.median(S[:,si])
    plt.figure(figsize=(15,10))
    dB=10.0*n.log10(S)
    nfloor=n.nanmedian(dB)
    plt.pcolormesh(fvec,rvec,dB,vmin=nfloor,vmax=nfloor+20)
    plt.colorbar()
    plt.title("HF ISR spectrum\n%s"%(stuffr.unix2datestr(i0/float(sr))))
    plt.xlabel("Doppler (kHz)")
    plt.ylabel("Range (km)")
    plt.tight_layout()
    plt.ylim([n.min(rvec),n.max(rvec)])
    plt.xlim([n.min(fvec),n.max(fvec)])
    
    plt.savefig("%s/isr-%d.png"%(out_dir_name,i0))
    plt.clf()
    plt.close()
    ho=h5py.File("%s/isr-%d.h5"%(out_dir_name,i0),"w")
    ho["i0"]=i0
    ho["S"]=S
    ho["rvec"]=rvec
    ho["fvec"]=fvec
    ho.close()
    

# real-time plot, just plot one integration period at last data in directory.
import sys
ch="cha"
if len(sys.argv) > 1:
    print("using channel %s"%(sys.argv[1]))
    ch=sys.argv[1]

print(ch)


#decode_channel(plot=True,ch=ch)
decode_all(plot=True,ch=ch)



