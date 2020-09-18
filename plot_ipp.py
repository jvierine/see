#!/usr/bin/env python


import digital_rf as drf
import matplotlib.pyplot as plt
import numpy as n
import stuffr

# real-time plot, just plot one integration period at last data in directory.
import sys
ch="ch1"
if len(sys.argv) > 1:
    print("using channel %s"%(sys.argv[1]))
    ch=sys.argv[1]

d=drf.DigitalRFReader("/data0/test200e3_7.953e6")

print("Recorded channels")
print(d.get_channels())
prop=d.get_properties(ch)
b=d.get_bounds(ch)

sr=prop["samples_per_second"]

i0=long(n.floor(b[1]/sr))*sr


ipp=100000/5
z=d.read_vector_c81d(i0,ipp,ch)
plt.plot(z.real)
plt.plot(z.imag)
plt.title("t=%s (UTC)"%(stuffr.unix2datestr(i0/sr)))
plt.show()


