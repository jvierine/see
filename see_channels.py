#!/usr/bin/env python

import sys
import see_config as sc
import stuffr
import numpy as n
import matplotlib.pyplot as plt
import scipy.signal as s
import digital_rf as drf


def show_specs(conf,d):
    specs=[]
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1/conf.sample_rate))
    chs=d.get_channels()
    b=d.get_bounds(chs[0])
    N=100
    w=s.hann(conf.nfft)
    
    for ci,ch in enumerate(chs):
        spec=n.zeros(conf.nfft)
        for i in range(N):
            spec+=n.fft.fftshift(n.abs(n.fft.fft(w*d.read_vector_c81d(b[1]-N*conf.nfft+i*conf.nfft,conf.nfft,ch)))**2.0)
        specs.append(spec)
        plt.plot(fvec,ci*20+10.0*n.log10(spec),label=ch)
    plt.legend()
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
    show_specs(conf,d)
