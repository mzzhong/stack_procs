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
    filtered_offset = pickle.load(f)

azOff   = filtered_offset['azOffset_filtered']    
rngOff  = filtered_offset['rngOffset_filtered']   
azOff_re      = filtered_offset['refer_azOffset']       
rngOff_re     = filtered_offset['refer_rngOffset']      
azOff_cmp            = filtered_offset['azOffset']             
rngOff_cmp           = filtered_offset['rngOffset']            
figdir              = filtered_offset['figdir']
title               = filtered_offset['title']

#(azOff, rngOff, azOff_re, rngOff_re, azOff_cmp, rngOff_cmp, figdir, title= None):

pad = 0.2

# Ranges of values.

### Plot azimuth offset. ###
vmin = np.nanmin(azOff_re)
vmax = np.nanmax(azOff_re)

vmin = - max(abs(vmin), abs(vmax))
vmax =   max(abs(vmin), abs(vmax))

vmin10 = np.floor(vmin*10)/10- pad
vmax10 = np.ceil(vmax*10)/10 + pad

fig, axs = plt.subplots(1,6, sharey=True, figsize=(24,15))

frac = 0.08
padbar = 0.05
shrink=0.7
tickstep=(vmax10-vmin10)/4

#ax = fig.add_subplot(161)
ax = axs[0]
ax.set_title('Predicted azimuth offset') 
im = ax.imshow(azOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10, vmax10+1e-5, tickstep),label='m')
ax.set_ylabel("offset index")
ax.set_xlabel("offset index")

#ax = fig.add_subplot(162)
ax=axs[1]
ax.set_title('Raw azimuth offset') 
im = ax.imshow(azOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax, fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10,vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

#ax = fig.add_subplot(163)
ax=axs[2]
ax.set_title('Filtered azimuth offset')
im = ax.imshow(azOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10, vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

### Plot range offset. ###
vmin = np.nanmin(rngOff_re)
vmax = np.nanmax(rngOff_re)

vmin = - max(abs(vmin), abs(vmax))
vmax =   max(abs(vmin), abs(vmax))

vmin10 = np.floor(vmin*10)/10-pad
vmax10 = np.ceil(vmax*10)/10+pad
tickstep=(vmax10-vmin10)/4

#ax = fig.add_subplot(164)
ax=axs[3]
ax.set_title('Predicted range offset')
im = ax.imshow(rngOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10,vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

#ax = fig.add_subplot(165)
ax=axs[4]
ax.set_title('Raw range offset') 
im = ax.imshow(rngOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10,vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

#ax = fig.add_subplot(166)
ax=axs[5]
ax.set_title('Filtered range offset')          
im = ax.imshow(rngOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
fig.colorbar(im,ax=ax,fraction=frac, pad=padbar, shrink=shrink, orientation='horizontal',ticks=np.arange(vmin10,vmax10+1e-5,tickstep),label='m')
ax.set_xlabel("offset index")

pathlib.Path(figdir).mkdir(parents=True, exist_ok=True)

#plt.tight_layout()
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

fig.savefig(os.path.join(figdir,'offset_' + title + ".pdf"), format='pdf',bbox_inches='tight')
fig.savefig(os.path.join(figdir,'offset_' + title + ".png"), format='png',bbox_inches='tight')

# Cut a line through to find the boundary of Swath, for Ridgecrest S1 only.
#fig  = plt.figure(2, figsize=(10,10))
#ax = fig.add_subplot(111)
#line = azOff[2000,:]
#ax.plot(line)
#fig.savefig(figdir + '/'+'line.png')
#print(np.arange(len(line))[line<-8])

plt.close(fig)

