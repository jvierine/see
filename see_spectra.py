#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import digital_rf as drf
import scipy.signal as ss
import stuffr
import re
import time
import os

def create_spectra(ch0="cha",
                   ch1="chb",
                   i0=0,
                   time_span=25*60,                   
                   f0=4.3,
                   n_fft=1000,
                   fft_len=4096,
                   flim=[4.15-0.5,4.15+0.5],
                   sample_rate=1e6,
                   dnames=["/data0/2020.09.08/test1e6_4.15e6"]):

    S=n.zeros([n_fft,fft_len],dtype=n.float64)
    
    
    
                   


d=drf.DigitalRFReader(dnames)
print(d.get_properties(ch0))
print(d.get_channels())
b=d.get_bounds(ch0)
print(b)
print((b[1]-b[0])/1e6)

S=n.zeros([fft_len,n_fft],dtype=n.complex64)

b=d.get_bounds(ch0)

    

realtime_updates(ch0="cha",ch1="chb",dt=25*60,f=4.15e6,flim=[4.15-0.5,4.15+0.5],dnames=["/data0/2020.09.08/test1e6_4.15e6"])
