##########################  Depreacated methods #################################
        def _run_offset_filter(self,data,data_snr,mask=None, label=None, refer=None):

            doc = self.doc
           
            med1_data = np.copy(data)

            # Step 0: mask out using snr
            do_step_0 = True
            
            if do_step_0:
                med1_data[data_snr<4] = np.nan

            # Step 1:
            do_step_1 = True
            if do_step_1:
                # Mask out where the reference is nan.
                med1_data[np.isnan(refer)] = np.nan
            
                # Using the mask
                if mask is not None: 
                    med1_data[mask] = np.nan

            # Step 2: normal median filter
            do_step_2 = True
            if do_step_2:
                med2_data = np.copy(data)
                med2_kernel_size = (7,7)

                # Median filters
                #med2_data = ndimage.median_filter(input=med2_data,size=med2_kernel_size,mode='nearest') 
                #med2_data = signal.medfilt(med1_data,kernel_size=med2_kernel_size)
                med2_data = mdf(med1_data, kernel_size=med2_kernel_size)
                #med2_data = cv2.medianBlur(med1_data,med2_kernel_size[0])
            else:
                med2_data = np.copy(med1_data)

            # Step 3: iteratively fill in small holes, fill in + median filter
            do_step_3 = True

            if do_step_3:
                iteration = 3
                hole_size = 7
                med3_data = np.copy(med2_data)
                for i in range(iteration):
                    med3_data = self.fill_in_holes(data=med3_data, hole_size = hole_size)
                    #med3_data = self.fill(med3_data)

            else:
                med3_data = np.copy(med2_data)


            ## Step 4: set neighboring values of nans as nan.
            do_step_4 = True
            if do_step_4:
    
                med4_kernel_size = 7
                med4_data = np.copy(med3_data)
                ind = np.where(np.isnan(med4_data))
    
                for i in range(med4_kernel_size):
                    for j in range(med4_kernel_size):
                        ind_i = ind[0] + i - med4_kernel_size//2
                        ind_i[ind_i < 0] = 0
                        ind_i[ind_i >= doc.numWinDown]=doc.numWinDown-1
    
                        ind_j = ind[1] + j - med4_kernel_size//2
                        ind_j[ind_j < 0] = 0
                        ind_j[ind_j >= doc.numWinAcross]=doc.numWinAcross-1
    
                        med4_data[(ind_i,ind_j)] = np.nan
            else:
                med4_data = np.copy(med3_data)

            # Step 5: remove the values on the margin, if interpolation is performed.
            if do_step_3:
                bb = (iteration + 1) * (hole_size // 2)
                med4_data[:iteration,:] = np.nan
                med4_data[-iteration:,:] = np.nan
                med4_data[:,:iteration] = np.nan
                med4_data[:,-iteration:] = np.nan

            return med4_data


        def _offset_filter(self, offsetfile, snrfile, offsetLosFile, title):

            # Filter the offset fields.
            # The intermediate unit is meter/day.
            # The final unit is still in pixel.

            doc = self.doc

            ds = gdal.Open(offsetfile)
            azOffset = ds.GetRasterBand(1).ReadAsArray()
            rngOffset = ds.GetRasterBand(2).ReadAsArray()

            ds = gdal.Open(snrfile)
            data_snr = ds.GetRasterBand(1).ReadAsArray()

            ds = gdal.Open(offsetLosFile)
            inc = ds.GetRasterBand(1).ReadAsArray()

            # Generatre a mask for invalid values at margins (useful for S1, GPU ampcor gives (-4, -4) to invalid cross-correlation))
            # Alternative way is to use offsetLosFile, but could be wrong if offsetLosFile doesn't match offset field
            # True: invalid, False: valid
            mask_of_invalid = np.logical_or(np.logical_and(azOffset==-4, rngOffset==-4), inc==0)

            # Reference
            # For glacier applications, load the predicted offset fields.
            if self.stack in ["tops","stripmap"]:
                refer_azOffset = doc.grossDown
                refer_rngOffset = doc.grossAcross

            # For Ridgecrest Earthquake, the reference is simply zero.
            else:
                refer_azOffset = np.zeros(shape=azOffset.shape)
                refer_rngOffset = np.zeros(shape=rngOffset.shape)

            # Convert reference from pixels to meter per day.
            refer_azOffset = refer_azOffset * doc.azPixelSize
            refer_rngOffset = refer_rngOffset * doc.rngPixelSize
 
            # Save the reference (already saved as grossDown/Across)
            #doc.refer_azOffset = refer_azOffset
            #doc.refer_rngOffset = refer_rngOffset

            # Convert observation from pixel to meter per day.
            azOffset = azOffset * doc.azPixelSize / doc.interim
            rngOffset = rngOffset * doc.rngPixelSize / doc.interim

            # Get vmin and vmax
            vmin_az = np.nanmin(refer_azOffset)
            vmax_az = np.nanmax(refer_azOffset)

            vmin_az = - max(abs(vmin_az), abs(vmax_az))
            vmax_az =   max(abs(vmin_az), abs(vmax_az))

            vmin_rng = np.nanmin(refer_rngOffset)
            vmax_rng = np.nanmax(refer_rngOffset)

            vmin_rng = - max(abs(vmin_rng), abs(vmax_rng))
            vmax_rng =   max(abs(vmin_rng), abs(vmax_rng))


            ###########     Filtering
            # Generate mask by reference. Values deviate too much from reference set as True. Use for deriving mis-coreg.
            mask_by_reference = self._get_refer_mask(azOffset, refer_azOffset, rngOffset, refer_rngOffset)
            # Combine the mask of both invalid values and erroreous estimations
            mask = np.logical_or(mask_of_invalid, mask_by_reference)

            if postproc_verbose:
                # Show the mask based on deviation from reference
                plt.figure(2, figsize=(10,10))
                plt.imshow(mask_by_reference)
                plt.title("The first mask: mask based on deviation from reference offset fields")
                plt.savefig('2.png')
     
                plt.figure(3, figsize=(10,10))
                plt.imshow(mask)
                plt.title("Combined mask of invalid values (-4, -4) and values of large deviations")
                plt.savefig('3.png')
    
                # Show the original unfiltered and uncoregistered azimuth offset
                fig = plt.figure(4, figsize=(10,10))
                ax1 = fig.add_subplot(121)
                im1 = ax1.imshow(azOffset,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                fig.colorbar(im1)
    
                ax2 = fig.add_subplot(122)
                im2 = ax2.imshow(rngOffset,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                fig.colorbar(im2)
    
                plt.title("Original Ampcor azimuth offset and range offset")
                plt.savefig('4.png')

            ############ Run filtering ###############################
            # Version 1: Currently deprecated
            #azOffset_filtered = self._run_offset_filter(azOffset, data_snr, mask=mask, label='az',refer=refer_azOffset)
            #rngOffset_filtered = self._run_offset_filter(rngOffset, data_snr, mask=mask, label='rng',refer=refer_rngOffset)

            # Version 2: 1) set the masked value to np.nan and 2) run 2-D median filter across the fields
            azOffset_filtered = self._run_offset_filter_v2(azOffset, data_snr, mask=mask, label='az',refer=refer_azOffset)
            rngOffset_filtered = self._run_offset_filter_v2(rngOffset, data_snr, mask=mask, label='rng',refer=refer_rngOffset)

            # Manually remove and first and second swath (ad hoc to Ridgecreast Earthquake)
            if self.stack == "tops_RC":
                if self.trackname == "track_64":
                    #azOffset_filtered[:,:644] = np.nan
                    #azOffset_filtered[:,1369:] = np.nan
                    #rngOffset_filtered[:,:644] = np.nan
                    #rngOffset_filtered[:,1369:] = np.nan
                    pass

                elif self.trackname == "track_7101":
                    pass


            if postproc_verbose:
                # Show the filtered range and azimuth offset
                fig = plt.figure(5, figsize=(10,10))
                ax1 = fig.add_subplot(121)
                im1 = ax1.imshow(azOffset_filtered,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                fig.colorbar(im1)
    
                ax2 = fig.add_subplot(122)
                im2 = ax2.imshow(rngOffset_filtered,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                fig.colorbar(im2)
                
                plt.title('Filtered but uncoregistered')
                plt.savefig('5.png')

            ############    Miscoregistration correction ###############
            # 1) update mask
            # Option 1
            # Generate mask by reference. Values deviate too much from reference set as True. Use for deriving mis-coreg.
            # May get nothing, because all bad values have been filtered out
            # mask_by_reference = self._get_refer_mask(azOffset_filtered, refer_azOffset, rngOffset_filtered, refer_rngOffset)

            # Option 2
            # Now, invalid values are NaN in filtered offset fields. Include them into mask_by_reference
            mask_by_reference[np.isnan(azOffset_filtered)] = True
            
            mask = mask_by_reference

            if postproc_verbose:
                # Show the mask corresponding the Nan value before miscoregistration correction
                plt.figure(figsize=(10,10))
                im = plt.imshow(mask)
                plt.title("NaN value before miscoregistration correction")
                plt.savefig('9.png')

            # 2) Get miscoregistration
            # Only true values in mask should be used for calculating miscoregistration
            az_mis, rng_mis = self._get_misCoreg(azOffset_filtered, refer_azOffset, rngOffset_filtered, refer_rngOffset, mask)
            print('Estimated miscoregistration (azimuth & range): ', az_mis, rng_mis)

            # Mis-coregistration correction.
            # 2019.03.04
            # Close the offset correction, because it is unnecessary for S1ab
            # Need to open up this for CSK
            miscoreg_correction = True

            if miscoreg_correction:
                
                azOffset_filtered = azOffset_filtered - az_mis
                rngOffset_filtered = rngOffset_filtered - rng_mis

                # Additional correction for Ridgecrest Earthquake
                if self.stack == "tops_RC":
                    if self.trackname == "track_64":
                        # SW3
                        #print(azOffset_filtered[:,1369:1379].tolist())
                        #print(np.nanmean(azOffset_filtered[:,1369:1379]))
                        #azOffset_filtered[:,1369:] -= np.nanmean(azOffset_filtered[:,1369:1379])
                        azOffset_filtered[:,1369:] -= 0.02
                        #print(np.nanmean(azOffset_filtered[:,1369:1379]))
                        #rngOffset_filtered[:,:644] = np.nan
                        pass
    
                    elif self.trackname == "track_7101":
                        pass

            if postproc_verbose:
                # Show the filtered and coregistered offset fields
                fig = plt.figure(10, figsize=(10,10))
                ax1 = fig.add_subplot(121)
                im1 = ax1.imshow(azOffset_filtered,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                fig.colorbar(im1)
    
                ax2 = fig.add_subplot(122)
                im2 = ax2.imshow(rngOffset_filtered,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                fig.colorbar(im2)
                
                plt.title('Filtered and coregistered')
                plt.savefig('10.png')

            # Display (Save the figures)
            # Fail to work, X server issue
            # option 1
            #self._display_az_rng(azOffset_filtered, rngOffset_filtered, refer_azOffset, refer_rngOffset, azOffset, rngOffset, title)
            # option 2
            #figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'radar'))
            #display_az_rng_global(azOffset_filtered, rngOffset_filtered, refer_azOffset, refer_rngOffset, azOffset, rngOffset, figdir, title)
            
            # Option 3
            # Save the file as objects to disk and plot them later
            filtered_offset_folder = disp_temp_folder
            #filtered_offset_folder = disp_temp_folder + '/' + title
            #if not os.path.exists(filtered_offset_folder):
            #    os.mkdir(filtered_offset_folder)
            
            filtered_offset = {}
            filtered_offset['azOffset_filtered'] = azOffset_filtered
            filtered_offset['rngOffset_filtered'] = rngOffset_filtered
            filtered_offset['refer_azOffset'] = refer_azOffset
            filtered_offset['refer_rngOffset'] = refer_rngOffset
            filtered_offset['azOffset'] = azOffset
            filtered_offset['rngOffset'] = rngOffset
            
            figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'radar'))
            filtered_offset['figdir'] = figdir
            
            filtered_offset['title'] = title + '_' + str(doc.runid)

            pkl_name = os.path.join(filtered_offset_folder, title +'.pkl')
            with open(pkl_name, "wb") as f:
                pickle.dump(filtered_offset, f)

            os.system("./display_az_rng.py " + pkl_name)

            # Convert the unit from meter per day to pixel.
            # Save the offset fields.
            self._save_az_rng(azOffset_filtered, rngOffset_filtered, mask_of_invalid)
            
            return 0



        def display_az_rng_global(azOff, rngOff, azOff_re, rngOff_re, azOff_cmp, rngOff_cmp, figdir, title= None):
        
            pad = 0.2
        
            # Ranges of values.
            # Plot azimuth offset.
            vmin = np.nanmin(azOff_re)
            vmax = np.nanmax(azOff_re)
        
            vmin = - max(abs(vmin), abs(vmax))
            vmax =   max(abs(vmin), abs(vmax))
        
            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad
        
            fig = plt.figure(figsize=(18,6))
        
            frac = 0.07
            padbar = 0.1
            tickstep=0.4
        
            ax = fig.add_subplot(161)
            ax.set_title('Predicted azimuth offset') 
            im = ax.imshow(azOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
            #fig.colorbar(im,fraction=0.07, orientation='horizontal',label='meter')
        
            ax = fig.add_subplot(162)
            ax.set_title('Raw azimuth offset') 
            im = ax.imshow(azOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
        
            ax = fig.add_subplot(163)
            ax.set_title('Filtered azimuth offset')
            im = ax.imshow(azOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
        
            # Plot range offset.
            vmin = np.nanmin(rngOff_re)
            vmax = np.nanmax(rngOff_re)
        
            vmin = - max(abs(vmin), abs(vmax))
            vmax =   max(abs(vmin), abs(vmax))
        
            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad
        
            ax = fig.add_subplot(164)
            ax.set_title('Predicted range offset')
            im = ax.imshow(rngOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
        
            ax = fig.add_subplot(165)
            ax.set_title('Raw range offset') 
            im = ax.imshow(rngOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
        
            ax = fig.add_subplot(166)
            ax.set_title('Filtered range offset')          
            im = ax.imshow(rngOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
        
            try:
                if not os.path.exists(figdir):
                    os.makedirs(figdir)
            except:
                pass
        
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
        
            return 0

