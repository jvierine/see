#!/usr/bin/env python
import numpy as n
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as s
import scipy.io as sio
import h5py
import sys
import stuffr
import os

def t2(x,sigma,alpha=2.0):
    # outlier remove
    A=n.zeros([len(x)+len(x),len(x)])
    m=[]
    n_meas=0
    for i in range(len(x)):
        if not n.isnan(x[i]):
            A[n_meas,i]=1.0/sigma[i]
            m.append(x[i])
            n_meas+=1
    for i in range(len(x)):
        m.append(0.0)
        if i>1 and i<(len(x)-1):
            A[i+n_meas,i]=-2.0*alpha
        else:
            A[i+n_meas,i]=-1.0*alpha
        if i > 1:
            A[i+n_meas,i-1]=1.0*alpha
        if i < len(x)-1:
            A[i+n_meas,i+1]=1.0*alpha
    m=n.array(m)
    xhat=n.linalg.lstsq(A[0:(n_meas+len(x))],m)[0]
    return(xhat)

def plot_file(fname,show_plot=False,clean_plot=False,vmin=-6,vmax=50.0):
    h=h5py.File(fname,"r+")
    I=h["S"][()]
    fscale=h["fscale"][()]

    dB=10.0*n.log10(I)
    dB=dB-n.nanmedian(10.0*n.log10(I[I>0]))
    fig=plt.figure(figsize=(8*1.5,6*1.5))
    ax1 = fig.add_subplot(111)
    #plt.title("%s"%(stuffr.unix2datestr(h["t0"][()])),y=1.06)


    freq_p=n.copy(h["freq_vec"][()])
    if fscale == "kHz":
        freq_p=freq_p/1e3

    tvec=h["time_vec"][()]
    dBT=n.transpose(dB)
    plt.title("%s $f_0=%1.2f$ (MHz)"%(stuffr.unix2datestr(h["t0"][()]),n.min(h["f0s"][()])/1e6),y=1.06)
        
    cm=ax1.pcolormesh(tvec,freq_p,dBT,vmin=vmin,vmax=vmax,cmap="plasma")
    ax1.set_xlabel("Time (s)")

    ax1.set_ylabel("Frequency offset (%s)"%(fscale))
    cb=fig.colorbar(cm,ax=ax1)
    cb.set_label("dB")
    ax1.set_ylim([n.min(freq_p),n.max(freq_p)])
    ax1.set_xlim([n.min(tvec),n.max(tvec)])    

    if n.max(h["f0s"][()])-n.min(h["f0s"][()]) > 10.0:
        ax2 = ax1.twiny()    
        ax2.plot(h["f0s"][()]/1e6, n.ones(len(h["f0s"][()])),alpha=0)
        ax2.set_xlabel("Heating frequency (MHz)")
        ax2.set_ylim([n.min(freq_p),n.max(freq_p)])
        ax2.set_xlim([n.min(h["f0s"][()]/1e6),n.max(h["f0s"][()]/1e6)])

    plt.tight_layout()
    imgfname="%s.png"%(fname)
    plt.savefig(imgfname)
    os.system("cp %s /tmp/0latest.png"%(imgfname))
    os.system("rsync -av --progress /tmp/0latest.png j@soppel.org:/var/www/soppel.org/see/")
    os.system("rsync -av --progress img/*.png j@soppel.org:/var/www/soppel.org/see/")

    if show_plot:
        plt.show()
    plt.clf()
    plt.close()

    if clean_plot:
        dBc=n.copy(dB)
        s0=dBc.shape[0]
        std_est=n.median(n.abs(dBc[0:(s0-1),:]-dBc[1:(s0),:]))

        diff=n.abs(dBc[0:(s0-1),:]-dBc[1:(s0),:])
        #  plt.pcolormesh(diff)
        # plt.show()
        idx=n.where(diff>20.0*std_est)
#        print(idx)
        dBc[idx[0]+1,idx[1]]=n.nan
        dBc[idx[0],idx[1]]=n.nan
 #       plt.pcolormesh(n.transpose(dBc),vmin=-3,vmax=40,cmap="plasma")
  #      plt.show()
    
        for i in range(dB.shape[1]):
            print("%d/%d"%(i,dB.shape[1]))
            m=dBc[:,i]
            sigma=n.ones(len(m))
            xhat=t2(m,sigma=sigma,alpha=0.1)
            shit_idx=n.where(n.isnan(dBc[:,i]))[0]
            dBc[shit_idx,i]=xhat[shit_idx]

        dBc=dBc-n.nanmedian(dBc)
        if "dB_clean" in h.keys():
            del h["dB_clean"]
            h["dB_clean"]=dBc


        fig=plt.figure(figsize=(8*1.5,6*1.5))
        ax1 = fig.add_subplot(111)
        plt.title("%s $f_0=%1.2f$ (MHz)"%(stuffr.unix2datestr(h["t0"][()]),n.min(h["f0s"][()])/1e6),y=1.06)
        cm=ax1.pcolormesh(h["time_vec"][()],freq_p,n.transpose(dBc),vmin=vmin,vmax=vmax,cmap="plasma")
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Frequency offset (%s)"%(fscale))
        ax1.set_xlim([n.min(tvec),n.max(tvec)])            
        cb=fig.colorbar(cm,ax=ax1)
        cb.set_label("dB")
        ax1.set_ylim([n.min(freq_p),n.max(freq_p)])
        print(h["f0s"][()])
        if n.max(h["f0s"][()])-n.min(h["f0s"][()]) > 10.0:
            ax2 = ax1.twiny()    
            print(h["f0s"][()])
            ax2.plot(h["f0s"][()]/1e6, n.ones(len(h["f0s"][()])),alpha=0)
            ax2.set_xlabel("Heating frequency (MHz)")
            ax2.set_ylim([n.min(freq_p),n.max(freq_p)])
            ax2.set_xlim([n.min(h["f0s"][()]/1e6),n.max(h["f0s"][()]/1e6)])

        plt.tight_layout()
        
        plt.savefig("%s.c.png"%(sys.argv[1]))
        if show_plot:
            plt.show()
        plt.clf()
        plt.close()
        h.close()
    
            
# example
if __name__ == "__main__":
    plot_file(sys.argv[1],show_plot=True,clean_plot=True,vmin=-6,vmax=20)
