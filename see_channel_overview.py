#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import scipy.signal as ss
import stuffr
import re
import time
import os

dnames=["/data0/2020.10.12/test1e6_4.04e6","/data1/2020.10.12/test1e6_4.04e6"]
ch=["cha","chb","chc","chd","che","chf","chg","chh"]
d=drf.DigitalRFReader(dnames)

p=d.get_properties("cha")
print(p)
print("Sample-rate %f "%(float(p["sample_rate_numerator"])/float(p["sample_rate_denominator"])))
print(d.get_channels())
b=d.get_bounds("cha")
print(b)
print((b[1]-b[0])/1e6)




nw=10
nfft=1024
ns=nw*nfft
sr=1e6

f0=4.04e6
fvec=n.fft.fftshift(n.fft.fftfreq(nfft,d=1.0/1e6))+f0
wf=ss.hann(nfft)
for c in ch:
    print(c)
    
    S=n.zeros(nfft)
    for i in range(nw):
        z=d.read_vector_c81d(b[1]-ns+nfft*i,nfft,c)
        S+=n.abs(n.fft.fftshift(n.fft.fft(wf*z)))**2.0
    plt.plot(fvec/1e6,10.0*n.log10(S),label=c)
plt.legend()
plt.show()

ntime=100
S=n.zeros([nfft,ntime])



for c in ch:
    tvec=[]
    b=d.get_bounds(c)
    idx0=b[0]
    step=int(n.floor((b[1]-b[0])/ntime))
    for i in range(ntime):
        z=d.read_vector_c81d(idx0+step*i,nfft,c)
        S[:,i]=n.abs(n.fft.fftshift(n.fft.fft(wf*z)))**2.0
        tvec.append(idx0+step*i)
    tvec=n.array(tvec)
    plt.pcolormesh(tvec,fvec/1e6,10.0*n.log10(S),cmap="plasma")
    plt.colorbar()
    plt.title(c)
    plt.show()
    
