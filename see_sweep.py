#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import scipy.signal as s
import digital_rf as drf

import sys
import see_config as sc
import stuffr
import os
import itertools as ito
import h5py
import clean_image 

def debug_start(conf,d):
    for i in range(30):
        z=d.read_vector_c81d(conf.t0*conf.sample_rate+conf.offset+i*conf.step_len,conf.step_len*conf.sample_rate,conf.ch[0])
        plt.plot(stuffr.decimate(n.abs(z)**2.0,dec=100))
        plt.title(i)
        plt.show()
        plt.plot(n.fft.fftshift(n.fft.fftfreq(1000000,d=1/conf.sample_rate))/1e6,10.0*n.log10(n.fft.fftshift(n.abs(n.fft.fft(z[0:1000000]))**2.0)))
        plt.show()

def phase_channels(conf,d,i0,carrier_width=10.0,cphases=None,camps=None,use_cphases=False):
    """
    simple calculation of phase differences between ch0 and other channels
    as well as channel amplitude.
    can be used to beamform all channels together.
    """
    if use_cphases:
        print("Using calibrated phases")
    else:
        cphases=n.zeros(len(conf.ch))
        camps=n.ones(len(conf.ch))        
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1.0/conf.sample_rate))
    fidx=n.where( (fvec>conf.fmin)&(fvec<conf.fmax))[0]
    carrier_fidx=n.where( (fvec>-10.0)&(fvec<10.0))[0]
    n_freq=len(fidx)
    nfft=conf.nfft
    step_len=conf.step_len*conf.sample_rate
    t=n.arange(nfft,dtype=n.float32)/conf.sample_rate
    wfun=s.chebwin(nfft,150)
    n_chan=len(conf.ch)
    ch_pairs=[]
    n_pair=0
    for i in range(1,n_chan):
        ch_pairs.append((0,i))
        n_pair+=1
    
    S=n.zeros([n_pair,conf.nsteps,n_freq],dtype=n.complex64)
    overlap=int(nfft/conf.overlap_fraction)
    nmax_avg=int(n.floor((conf.step_len*conf.sample_rate-nfft-conf.offset)/overlap))
    n_avg=1
    tvec=n.zeros(conf.nsteps)
    pwrs=n.zeros(len(conf.ch))
    npwrs=n.zeros(len(conf.ch))    
    for step_idx in range(conf.nsteps):
        fnow = conf.f0 + step_idx*conf.fstep
        fshift = conf.center_freq-fnow
        dshift=n.exp(1j*2.0*n.pi*fshift*t)
        
        for pi,ci in enumerate(ch_pairs):
            for avg_i in range(n_avg):
                inow=i0+step_idx*step_len + avg_i*overlap
                print("%s n_avg %d/%d f0 %1.2f"%(stuffr.unix2datestr(inow/conf.sample_rate),n_avg,nmax_avg,fnow/1e6))
                tvec[step_idx]=step_idx*step_len/conf.sample_rate
                
                z0=dshift*d.read_vector_c81d(inow,nfft,conf.ch[ci[0]])*camps[ci[0]]*n.exp(1j*cphases[ci[0]])
                pwrs[ci[0]]+=n.mean(n.abs(z0)**2.0)
                npwrs[ci[0]]+=1.0
                z1=dshift*d.read_vector_c81d(inow,nfft,conf.ch[ci[1]])*camps[ci[1]]*n.exp(1j*cphases[ci[1]])
                pwrs[ci[1]]+=n.mean(n.abs(z1)**2.0)
                npwrs[ci[1]]+=1.0                
                X0=n.fft.fftshift(n.fft.fft(wfun*z0))
                X1=n.fft.fftshift(n.fft.fft(wfun*z1))                
                S[pi,step_idx,:]+=X0[fidx]*n.conj(X1[fidx])
                
    gfidx=n.where( n.abs((fvec[fidx])>carrier_width))[0]    
    phase_diffs=[]
    
    for i in range(n_pair):
        xc=n.copy(S[i,:,gfidx])
        xc=xc.flatten()
        xc_idx=n.argsort(xc)
        phase_diff=n.angle(n.sum(xc[xc_idx[int(len(xc_idx)*0.3):int(len(xc_idx)*0.9)]]))
        phase_diffs.append(phase_diff)
        if use_cphases:
            plt.pcolormesh(n.angle(S[i,:,:]))
            plt.colorbar()
            plt.show()
            plt.plot(n.angle(n.sum(S[i,:,:],axis=0)))
            plt.axhline(phase_diff)
            plt.show()
        
    ch_pwrs=pwrs/npwrs

    phases=n.zeros(n_chan)
    # simple phaseup with reference to channel 0
    # ph0-ph1
    # ph0-ph2
    # ph0-ph3
    # ph0-ph4
    for i in range(n_pair):
        phases[i+1]=phase_diffs[i]
#    plt.plot(phases)
 #   plt.show()
    # pha = 0
    # pha - phb = mab
    # pha - phc = 
    return(1.0/n.sqrt(ch_pwrs),phases)
                
    
        
def calculate_sweep(conf,d,i0,use_cphases=False,cphases=None,camps=None):

    if use_cphases:
        print("Using calcibrated phases")
    else:
        cphases=n.zeros(len(conf.ch))
        camps=n.ones(len(conf.ch))
        
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1.0/conf.sample_rate))
    fidx=n.where( (fvec>conf.fmin)&(fvec<conf.fmax))[0]

    carrier_fidx=n.where( (fvec>-10.0)&(fvec<10.0))[0]
    
    n_freq=len(fidx)
    
    nfft=conf.nfft
    step_len=conf.step_len*conf.sample_rate
    t=n.arange(nfft,dtype=n.float32)/conf.sample_rate

    wfun=s.chebwin(nfft,150)

    S=n.zeros([conf.nsteps*conf.nsubsteps,n_freq])
    N_avgd=n.zeros([conf.nsteps*conf.nsubsteps,n_freq])    

    overlap=int(nfft/conf.overlap_fraction)
    nmax_avg=int(n.floor((conf.step_len*conf.sample_rate-conf.trim_end)/overlap/conf.nsubsteps))+1
    sub_len=int(n.floor((step_len-conf.trim_end)/conf.nsubsteps))
    
    if conf.fast:
        n_avg=n.min([conf.n_avg,nmax_avg])
    else:
        n_avg=nmax_avg
    
    tvec=n.zeros(conf.nsteps*conf.nsubsteps)
    carrier = n.zeros(conf.nsteps*conf.nsubsteps,dtype=n.complex64)
    f0s=n.zeros(conf.nsteps*conf.nsubsteps)
    
    for step_idx in range(conf.nsteps):
        fnow = conf.f0 + step_idx*conf.fstep
        fshift = conf.center_freq-fnow
        dshift=n.exp(1j*2.0*n.pi*fshift*t)
        step_i0=i0+step_idx*step_len
        
        for sub_idx in range(conf.nsubsteps):
            
            f0s[step_idx*conf.nsubsteps+sub_idx]=fnow
            tvec[step_idx*conf.nsubsteps+sub_idx]=(step_idx*step_len+sub_idx*sub_len)/conf.sample_rate

            
            
            for avg_i in range(n_avg):
                inow=i0 + step_idx*step_len + sub_idx*sub_len + avg_i*overlap
                
                # make sure this fits the step
                if (inow+nfft-step_i0+conf.trim_end) < step_len:
                    print("%s n_avg %d/%d f0 %1.2f"%(stuffr.unix2datestr(inow/conf.sample_rate),avg_i,nmax_avg,fnow/1e6))
                    z=n.zeros(nfft,dtype=n.complex128)
                    # beamform signals
                    for ch_i in range(len(conf.ch)):
                        z+=dshift*d.read_vector_c81d(inow,nfft,conf.ch[ch_i])*camps[ch_i]*n.exp(1j*cphases[ch_i])
                    
                    X=n.fft.fftshift(n.fft.fft(wfun*z))
                    if conf.debug:
                        plt.plot(fvec,10.0*n.log10(X))
                        plt.show()
                
                    S[step_idx*conf.nsubsteps+sub_idx,:]+=n.abs(X[fidx])**2.0
                    N_avgd[step_idx*conf.nsubsteps+sub_idx,:]+=1.0
                    carrier[step_idx*conf.nsubsteps+sub_idx]+=n.sum(n.abs(X[carrier_fidx])**2.0)
                
    S=S/N_avgd
    ofname="img/%s_sweep_%1.2f.h5"%(conf.prefix,i0/conf.sample_rate)
    print("Saving %s"%(ofname))
    ho=h5py.File(ofname,"w")
    ho["S"]=S
    ho["phases"]=cphases
    ho["amps"]=camps
    ho["time_vec"]=tvec
    ho["freq_vec"]=fvec[fidx]
    ho["f0s"]=f0s
    ho["carrier_pwr"]=carrier
    ho["center_freq"]=conf.center_freq
    ho["fscale"]=str(conf.fscale)
    ho["t0"]=i0/conf.sample_rate
    ho["date"]=stuffr.unix2datestr(i0/conf.sample_rate)
    ho.close()

    clean_image.plot_file(ofname,show_plot=conf.show_plot)
        
def calculate_sweep_xc(conf,d,i0,use_cphases=False,cphases=None,camps=None):

    if use_cphases:
        print("Using calcibrated phases")
    else:
        cphases=n.zeros(len(conf.ch))
        camps=n.ones(len(conf.ch))
        
    fvec=n.fft.fftshift(n.fft.fftfreq(conf.nfft,d=1.0/conf.sample_rate))
    fidx=n.where( (fvec>conf.fmin)&(fvec<conf.fmax))[0]

    carrier_fidx=n.where( (fvec>-10.0)&(fvec<10.0))[0]
    
    n_freq=len(fidx)
    
    nfft=conf.nfft
    step_len=conf.step_len*conf.sample_rate
    t=n.arange(nfft,dtype=n.float32)/conf.sample_rate

    wfun=s.chebwin(nfft,150)

    n_chan=len(conf.ch)
    cpix=ito.combinations(n.arange(n_chan,dtype=n.int),2)
    ch_pairs=[]
    for i in range(n_chan):
        ch_pairs.append((i,i))
    for cpi in cpix:
        ch_pairs.append((cpi[0],cpi[1]))
        
    n_pairs=len(ch_pairs)
        
    S=n.zeros([n_pairs,conf.nsteps,n_freq],dtype=n.complex64)

    overlap=int(nfft/conf.overlap_fraction)
    
    nmax_avg=int(n.floor((conf.step_len*conf.sample_rate-nfft-conf.offset)/overlap))+1
    
    if conf.fast:
        n_avg=n.min([conf.n_avg,nmax_avg])
    else:
        n_avg=nmax_avg
    
    tvec=n.zeros(conf.nsteps)
    carrier = n.zeros(conf.nsteps,dtype=n.complex64)
    f0s=n.zeros(conf.nsteps)
    
    for step_idx in range(conf.nsteps):
        fnow = conf.f0 + step_idx*conf.fstep
        f0s[step_idx]=fnow
        fshift = conf.center_freq-fnow
        dshift=n.exp(1j*2.0*n.pi*fshift*t)

        tvec[step_idx]=step_idx*step_len/conf.sample_rate
        
        for avg_i in range(n_avg):
            inow=i0+step_idx*step_len + avg_i*overlap
            print("%s n_avg %d/%d f0 %1.2f"%(stuffr.unix2datestr(inow/conf.sample_rate),n_avg,nmax_avg,fnow/1e6))
            
            for pi,chp_i in enumerate(ch_pairs):                                        
                
                z0=dshift*d.read_vector_c81d(inow,nfft,conf.ch[chp_i[0]])*camps[chp_i[0]]*n.exp(1j*cphases[chp_i[0]])
                z1=dshift*d.read_vector_c81d(inow,nfft,conf.ch[chp_i[1]])*camps[chp_i[1]]*n.exp(1j*cphases[chp_i[1]])
                X0=n.fft.fftshift(n.fft.fft(wfun*z0))
                X1=n.fft.fftshift(n.fft.fft(wfun*z1))                
                S[pi,step_idx,:]+=X0[fidx]*n.conj(X1[fidx])
                
    
    if conf.fscale == "kHz":
        fvec=fvec/1e3

    ho=h5py.File("img/%s_sweep_xc_%1.2f.h5"%(conf.prefix,i0/conf.sample_rate),"w")
    ho["S"]=S
    ho["ch_pairs"]=ch_pairs
    ho["phases"]=cphases
    ho["amps"]=camps
    ho["time_vec"]=tvec
    ho["freq_vec"]=fvec[fidx]
    ho["f0s"]=f0s
    ho["center_freq"]=conf.center_freq
    ho["t0"]=i0/conf.sample_rate
    ho["date"]=stuffr.unix2datestr(i0/conf.sample_rate)
    ho.close()



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

    if conf.debug_timing:
        debug_start(conf,d)


    i0=conf.t0*conf.sample_rate + conf.offset

    if conf.xc:
        for i in range(conf.n_cycles):
            calculate_sweep_xc(conf,d,i0+i*conf.cycle_len*conf.sample_rate)
    
    
    if len(conf.ch)>1:
        amps,phases=phase_channels(conf,d,i0)
        # debug phasing
        #phase_channels(conf,d,i0,cphases=phases,camps=amps,use_cphases=True)    
    
    for i in range(conf.cycle_num,conf.n_cycles):
        if len(conf.ch)>1:
            calculate_sweep(conf,d,i0+i*conf.cycle_len*conf.sample_rate,use_cphases=True,cphases=phases,camps=amps)
        else:
            calculate_sweep(conf,d,i0+i*conf.cycle_len*conf.sample_rate,use_cphases=False)
