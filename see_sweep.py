#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import scipy.signal as s
import digital_rf as drf

import sys
import see_config as sc
import stuffr
import os


def calculate_sweep(conf,d):
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1.0/conf.sample_rate))
    fidx=n.where( (fvec>conf.fmin)&(fvec<conf.fmax))[0]

    carrier_fidx=n.where( (fvec>-10.0)&(fvec<10.0))[0]
    
    n_freq=len(fidx)
    i0=conf.t0*conf.sample_rate + conf.offset
    nfft=conf.nfft
    step_len=conf.step_len*conf.sample_rate
    t=n.arange(nfft,dtype=n.float32)/conf.sample_rate

    wfun=s.chebwin(nfft,200)

    S=n.zeros([conf.nsteps,n_freq])

    overlap=nfft/2
    n_avg=n.min([10,int(n.floor((conf.step_len*conf.sample_rate-nfft-conf.offset)/overlap))])
    
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
                print(stuffr.unix2datestr(inow/conf.sample_rate))

                z=dshift*d.read_vector_c81d(inow,nfft,conf.ch[ch_i])
 #               df=n.angle(n.mean(z[0:(len(z)-1000)]*n.conj(z[1000:len(z)])))
#                print(df)
#                z=z*n.exp(1j*n.arange(len(z))*df/1000.0)
      #          zd=stuffr.decimate(z,dec=1000)
     #           plt.plot(zd.real)
       #         plt.plot(zd.imag)
        #        plt.show()
#                z=z-n.mean(z)
                X=n.fft.fftshift(n.fft.fft(wfun*z))
                print(len(X))
                print(len(fidx))
                S[step_idx,:]+=n.abs(X[fidx])**2.0

                carrier[step_idx]+=n.sum(n.abs(X[carrier_fidx])**2.0)
        #plt.plot(z.real)
        #plt.plot(z.imag)
        #plt.show()
        #plt.plot(10.0*n.log10(n.abs(X)**2.0))
        #plt.show()
    dB=10.0*n.log10(S)
    dB=dB-n.nanmedian(dB)
    plt.pcolormesh(tvec,fvec[fidx],n.transpose(dB),vmin=-3,vmax=70.0)
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    cb=plt.colorbar()
    cb.set_label("dB")
    plt.show()
    plt.plot(tvec,10.0*n.log10(carrier))
    plt.xlabel("Time (s)")
    plt.ylabel("Carrier power (dB)")
    plt.show()

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

    calculate_sweep(conf,d)
