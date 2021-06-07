#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import numpy as np
import pickle
import os
import pathlib

pkl_name = sys.argv[1]
with open(pkl_name, "rb") as f:
    gc_filtered_offset = pickle.load(f)

azOff_re            = gc_filtered_offset['refer_azOffset']       
rngOff_re           = gc_filtered_offset['refer_rngOffset']      

gc_azOff_cmp_std       = gc_filtered_offset['gc_filtered_azOffsetStd'] 
gc_rngOff_cmp_std      = gc_filtered_offset['gc_filtered_rngOffsetStd']

gc_azOff               = gc_filtered_offset['gc_filtered_azOffset']    
gc_rngOff              = gc_filtered_offset['gc_filtered_rngOffset']

gc_hgt                  = gc_filtered_offset['gc_hgt_offset'] 

figdir              = gc_filtered_offset['figdir']
title               = gc_filtered_offset['title']

#(azOff, rngOff, azOff_re, rngOff_re, azOff_cmp, rngOff_cmp, figdir, title= None):

pad = 0.2
tick_sca = 3
frac = 0.08
padbar = 0.05
shrink=0.7

fig, axs = plt.subplots(1,5, sharey=True, figsize=(23,15))
ax = axs[0]
ax.set_title("hgt")
im = ax.imshow(gc_hgt, cmap=cm.jet)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',label='m')
ax.set_ylabel("offset index")
ax.set_xlabel("offset index")


### Plot azimuth offset. ###
if "az_vmin" in gc_filtered_offset:
    vmin10 = gc_filtered_offset["az_vmin"]
    vmax10 = gc_filtered_offset["az_vmax"]
else:
    vmin = np.nanmin(azOff_re)
    vmax = np.nanmax(azOff_re)
    
    vmax = max(abs(vmin), abs(vmax))
    vmin = - vmax
    
    vmin10 = np.floor(vmin*10)/10 - pad
    vmax10 = np.ceil(vmax*10)/10 + pad

tickstep=(vmax10-vmin10)/4

ax = axs[1]
ax.set_title('azimuth offset') 
im = ax.imshow(gc_azOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10, vmax10+1e-5, tickstep),label='m')
ax.set_ylabel("offset index")
ax.set_xlabel("offset index")

ax=axs[3]
ax.set_title('azimuth offset std') 
im = ax.imshow(gc_azOff_cmp_std, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax, fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10,vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

### Plot range offset. ###

if "rng_vmin" in gc_filtered_offset:
    vmin10 = gc_filtered_offset["rng_vmin"]
    vmax10 = gc_filtered_offset["rng_vmax"]
else:
    vmin = np.nanmin(rngOff_re)
    vmax = np.nanmax(rngOff_re)
    
    vmax = max(abs(vmin), abs(vmax))
    vmin = - vmax
    
    vmin10 = np.floor(vmin*10)/10 - pad
    vmax10 = np.ceil(vmax*10)/10 + pad

tickstep=(vmax10-vmin10)/4

ax=axs[2]
ax.set_title('range offset')
im = ax.imshow(gc_rngOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',label='m')
ax.set_xlabel("offset index")

ax=axs[4]
ax.set_title('range offset std')
im = ax.imshow(gc_rngOff_cmp_std, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10, vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

# Save it
fig.savefig(os.path.join(figdir,'gc_offset_' + title + ".pdf"), format='pdf',bbox_inches='tight')
fig.savefig(os.path.join(figdir,'gc_offset_' + title + ".png"), format='png',bbox_inches='tight')

plt.close(fig)

