#!/usr/bin/env python3

import numpy as np
from scipy import signal

def medianfilter2d(matrix,kernel_size=5):

    m,n = matrix.shape

    ksize = np.max(kernel_size)

    if ksize % 2 == 0:
        ksize = ksize + 1

    half = ksize//2

    nullval = 1000

    # pad nullval
    padded = np.zeros(shape=(m+half*2,n+half*2)) + nullval

    padded[half:m+half, half:n+half] = matrix

    # set np.nan value to be nulval
    padded[np.isnan(padded)] = nullval
   
    # median filter
    filt = signal.medfilt(padded,kernel_size=ksize)

    # set nullval to np.nan
    filt[filt==nullval] = np.nan

    result = filt[half:m+half, half:n+half]

    return result


def medianfilter2d_temp(matrix,kernel_size=5):

    m,n = matrix.shape

    #return np.copy(matrix)

    filtered = np.zeros(shape=matrix.shape)

    if isinstance(kernel_size,int):
        h_size0 = kernel_size // 2
        h_size1 = kernel_size // 2
    else:
        h_size0 = kernel_size[0] // 2
        h_size1 = kernel_size[1] // 2
        

    for i in range(m):
        for j in range(n):
            if matrix[i,j] is not np.nan and ( i>=h_size0 or i<=(m-h_size0) or j>=h_size1 or j<=(n-h_size1) ):
                chip = matrix[i-h_size0:i+h_size0, j-h_size1:j+h_size1]
                filtered[i,j] = np.nanmedian(chip)

            else:
                filtered[i,j] = np.nan

    return matrix
