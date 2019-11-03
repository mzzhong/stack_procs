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

