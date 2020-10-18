#!/usr/bin/env python
import numpy as n
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as s
import scipy.io as sio
import h5py
import sys
import stuffr

# theory matrix with individual weights for each pixel
def create_theory_matrix(I,sigma,alpha_x=0.1,alpha_y=0.1):
    n_x=I.shape[0]
    n_y=I.shape[1]
    n_pa=n_x*n_y
    row_idx=0

    rows = []
    cols = []
    data = []
    d = []
    
    # measurements
    for i in range(n_x):
        for j in range(n_x):
            rows.append(row_idx)
            cols.append(i*n_y + j)
            data.append(1.0/sigma[i,j])
            d.append(I[i,j]/sigma[i,j])
            row_idx+=1
            
    # laplacian regularization
    for i in range(1,n_x-1):
        for j in range(1,n_y-1):
            rows.append(row_idx)
            cols.append(i*n_y + j-1)
            data.append(alpha_x)
            
            rows.append(row_idx)
            cols.append(i*n_y + j)
            data.append(-2.0*alpha_x)
            
            rows.append(row_idx)
            cols.append(i*n_y + j+1)
            data.append(alpha_x)
            d.append(0.0)
            row_idx += 1

            rows.append(row_idx)
            cols.append((i-1)*n_y + j)
            data.append(alpha_y)
            
            rows.append(row_idx)
            cols.append(i*n_y + j)
            data.append(-2.0*alpha_y)
            
            rows.append(row_idx)
            cols.append((i+1)*n_y + j)
            data.append(alpha_y)
            d.append(0.0)
            row_idx += 1

    n_row=row_idx
    n_col=n_x*n_y
    # create a sparse theory matrix
    G=coo_matrix((data, (rows, cols)), shape=(n_row, n_col), dtype=n.float32)
    # return matrix and measurement
    return(G,d)


def get_clean_image(I,alpha_x=1.0,alpha_y=1.0):

    n_x=I.shape[0]
    n_y=I.shape[1]
    sigma=n.zeros([n_x,n_y])


    # iterative renormalized least-squares
    sigma[:,:]=1.0
    for i in range(5):
        print("i %d"%(i))
        G,d=create_theory_matrix(I,sigma,alpha_x=0.1,alpha_y=1.0)
        i0=s.lsqr(G,d)[0]
        i0.shape=(n_x,n_y)
        sigma=(n.abs(I-i0)+1.0)*0.1
    

    
    # star-field is the residuals
    rfi=(I-i0)
    rfi[rfi<0]=0
    
    plt.subplot(121)
    plt.pcolormesh(rfi)
    plt.colorbar()
    plt.subplot(122)
    plt.pcolormesh(i0)
    plt.colorbar()    
    plt.show()


def t2(x,sigma,alpha=2.0):
    A=n.zeros([len(x)+len(x),len(x)])
#    sigma=n.zeros(len(x))
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
#    plt.pcolormesh(A)
 #   plt.show()
    xhat=n.linalg.lstsq(A[0:(n_meas+len(x))],m)[0]

    return(xhat)
#    plt.plot(xhat)
 #   plt.show()
            
            
# example
if __name__ == "__main__":
    h=h5py.File(sys.argv[1],"r+")
    print(h.keys())
    I=h["S"].value

    dB=10.0*n.log10(I)
    dB=dB-n.nanmedian(dB)
    fig=plt.figure(figsize=(8*1.5,6*1.5))
    ax1 = fig.add_subplot(111)
    plt.title("%s"%(stuffr.unix2datestr(h["t0"].value)),y=1.05)
    cm=ax1.pcolormesh(h["time_vec"].value,h["freq_vec"].value,n.transpose(dB),vmin=-6,vmax=40,cmap="plasma")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Frequency offset (kHz)")
    cb=fig.colorbar(cm,ax=ax1)
    cb.set_label("dB")
    ax1.set_ylim([n.min(h["freq_vec"].value),n.max(h["freq_vec"].value)])

    ax2 = ax1.twiny()    
    print(h["f0s"].value)
    ax2.plot(h["f0s"].value/1e6, n.ones(len(h["f0s"].value)),alpha=0)
    ax2.set_xlabel("Heating frequency (MHz)")
 #   ax2.set_xlim([-1,1])
#    ax2.cla()
    plt.tight_layout()
    plt.savefig("%s.png"%(sys.argv[1]))
    plt.show()
    plt.clf()
    plt.close()

    dBc=n.copy(dB)
    s0=dBc.shape[0]
    std_est=n.median(n.abs(dBc[0:(s0-1),:]-dBc[1:(s0),:]))

    diff=n.abs(dBc[0:(s0-1),:]-dBc[1:(s0),:])
  #  plt.pcolormesh(diff)
   # plt.show()
    idx=n.where(diff>20.0*std_est)
    print(idx)
    dBc[idx[0]+1,idx[1]]=n.nan
    dBc[idx[0],idx[1]]=n.nan
    plt.pcolormesh(n.transpose(dBc),vmin=-3,vmax=40,cmap="plasma")
    plt.show()
    
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
    plt.title("%s"%(stuffr.unix2datestr(h["t0"].value)),y=1.05)
    cm=ax1.pcolormesh(h["time_vec"].value,h["freq_vec"].value,n.transpose(dBc),vmin=-6,vmax=40,cmap="plasma")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Frequency offset (kHz)")
    cb=fig.colorbar(cm,ax=ax1)
    cb.set_label("dB")
    ax1.set_ylim([n.min(h["freq_vec"].value),n.max(h["freq_vec"].value)])
    ax2 = ax1.twiny()    
    print(h["f0s"].value)
    ax2.plot(h["f0s"].value/1e6, n.ones(len(h["f0s"].value)),alpha=0)
    ax2.set_xlabel("Heating frequency (MHz)")
    plt.tight_layout()
    plt.savefig("%s.c.png"%(sys.argv[1]))
    plt.show()
    plt.clf()
    plt.close()
    h.close()
