#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import scipy.signal as s
import digital_rf as drf

import sys
import see_config as sc
import stuffr
import os

def debug_start(conf,d):
    for i in range(30):
        z=d.read_vector_c81d(conf.t0*conf.sample_rate+conf.offset+i*conf.step_len,conf.step_len*conf.sample_rate,conf.ch[0])
        plt.plot(stuffr.decimate(n.abs(z)**2.0,dec=100))
        plt.title(i)
        plt.show()
        plt.plot(n.fft.fftshift(n.fft.fftfreq(1000000,d=1/conf.sample_rate))/1e6,10.0*n.log10(n.fft.fftshift(n.abs(n.fft.fft(z[0:1000000]))**2.0)))
        plt.show()

def calculate_sweep(conf,d,i0):
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1.0/conf.sample_rate))
    fidx=n.where( (fvec>conf.fmin)&(fvec<conf.fmax))[0]

    carrier_fidx=n.where( (fvec>-10.0)&(fvec<10.0))[0]
    
    n_freq=len(fidx)
    
    nfft=conf.nfft
    step_len=conf.step_len*conf.sample_rate
    t=n.arange(nfft,dtype=n.float32)/conf.sample_rate

    wfun=s.chebwin(nfft,150)

    S=n.zeros([conf.nsteps,n_freq])


    overlap=nfft/2
    nmax_avg=int(n.floor((conf.step_len*conf.sample_rate-nfft-conf.offset)/overlap))
    if conf.fast:
        n_avg=n.min([conf.n_avg,nmax_avg])
    else:
        n_avg=nmax_avg
    
    tvec=n.zeros(conf.nsteps)
    carrier = n.zeros(conf.nsteps,dtype=n.complex64)
    for step_idx in range(conf.nsteps):
        fnow = conf.f0 + step_idx*conf.fstep
        fshift = conf.center_freq-fnow
        dshift=n.exp(1j*2.0*n.pi*fshift*t)

        for ch_i in range(len(conf.ch)):
            for avg_i in range(n_avg):
                inow=i0+step_idx*step_len + avg_i*overlap
                tvec[step_idx]=step_idx*step_len/conf.sample_rate
                print("%s n_avg %d/%d f0 %1.2f"%(stuffr.unix2datestr(inow/conf.sample_rate),n_avg,nmax_avg,fnow/1e6))
                
                z=dshift*d.read_vector_c81d(inow,nfft,conf.ch[ch_i])
                X=n.fft.fftshift(n.fft.fft(wfun*z))
                if conf.debug:
                    plt.plot(fvec,10.0*n.log10(X))
                    plt.show()
                S[step_idx,:]+=n.abs(X[fidx])**2.0

                carrier[step_idx]+=n.sum(n.abs(X[carrier_fidx])**2.0)
        #plt.plot(z.real)
        #plt.plot(z.imag)
        #plt.show()
        #plt.plot(10.0*n.log10(n.abs(X)**2.0))
        #plt.show()
    dB=10.0*n.log10(S)
    dB=dB-n.nanmedian(dB)
    
    if conf.fscale == "kHz":
        fvec=fvec/1e3
        
    plt.pcolormesh(tvec,fvec[fidx],n.transpose(dB),vmin=-3,vmax=60.0)
    plt.xlabel("Time (s)")
    if conf.fscale == "kHz":
        plt.ylabel("Frequency (kHz)")
        plt.ylim([conf.fmin/1e3,conf.fmax/1e3])        
    else:
        plt.ylabel("Frequency (Hz)")
        plt.ylim([conf.fmin,conf.fmax])        
        
    cb=plt.colorbar()
    cb.set_label("dB")

    plt.title("Cycle start %s"%(stuffr.unix2datestr(i0/conf.sample_rate)))
    plt.savefig("img/wb_sweep_%1.2f.png"%(i0/conf.sample_rate))
    if conf.show_plot:
        plt.show()
    else:
        plt.clf()
        plt.close()
    
    plt.plot(tvec,10.0*n.log10(carrier))
    plt.xlabel("Time (s)")
    plt.ylabel("Carrier power (dB)")
    plt.savefig("img/wb_pwr_%1.2f.png"%(i0/conf.sample_rate))
    if conf.show_plot:
        plt.show()
    else:
        plt.clf()
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        conf=sc.see_config(sys.argv[1])
    else:
        conf=sc.see_config()

    print(conf)

    d=drf.DigitalRFReader(conf.data_dirs)
    print("Channels")
    chs=d.get_channels()
    print(chs)
    print("Data extent:")
    b=d.get_bounds(chs[0])
    print("%s-%s"%(stuffr.unix2datestr(b[0]/conf.sample_rate),stuffr.unix2datestr(b[1]/conf.sample_rate)))

    if conf.debug:
        debug_start(conf,d)


    i0=conf.t0*conf.sample_rate + conf.offset
    
    for i in range(conf.n_cycles):
        calculate_sweep(conf,d,i0+i*conf.nsteps*conf.step_len*conf.sample_rate)
