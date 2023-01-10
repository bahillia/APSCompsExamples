# -*- coding: utf-8 -*-
"""
Created on Fri May 14 13:31:17 2021

@author: patbe
"""
import numpy as np
from scipy import interpolate
from scipy import integrate
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel, Box1DKernel
from astropy.io import fits
from IPython.display import clear_output

def avg_error(error,fit_mask):
    fit_error = error[fit_mask]
    low_error = np.where(fit_error < 1e-16)[0]
    means = np.mean(
        [fit_error[low_error-2],fit_error[low_error-1],fit_error[low_error+1],fit_error[low_error+2]],axis=0)
    fit_error[low_error]=means
    error[fit_mask] = fit_error
    return error

# Chop takes a wavelength array and a flux array and returns an array with only the values between [min_val,max_val]
def chop(data,drange):
    min_ind = np.nonzero(data[0] >= drange[0])[0][0]
    max_ind = np.nonzero(data[0] >= drange[1])[0][0]
    
    data = np.asarray([data[0][min_ind:max_ind],
            data[1][min_ind:max_ind],
            data[2][min_ind:max_ind]])
    
    return data

def combine(data1,data2):
    #The combined spectrum should be interpolated to a new array
    # Length = combined length of individual arrays
    wave1, flux1, err1,flag1 = data1
    wave2, flux2, err2,flag2 = data2
    
    # If there is no overlap region, return the concatenated data and fill in with 
    #      proxy spectrum later
    if(wave1[-1] < wave2[0]):
        comb_data = np.concatenate((data1,data2),axis=1)
        return comb_data
    
    comb_wave = np.zeros(len(wave1)+len(wave2))
    comb_flux = np.zeros(len(flux1)+len(flux2))
    comb_err = np.zeros(len(err1)+len(err2))
    comb_flag = np.zeros(len(flag1)+len(flag2))
    
    min_ind, max_ind = get_overlap(wave1,wave2)
    
    #Keep first array up to the overlap point
    comb_wave[0:min_ind] = wave1[0:min_ind]
    comb_flux[0:min_ind] = flux1[0:min_ind]
    comb_err[0:min_ind] = err1[0:min_ind]
    comb_flag[0:min_ind] = flag1[0:min_ind]
    
    combined_index = min_ind
    wave1_index = min_ind
    wave2_index = 0
    # Iterate until we reach the end of wave1 (the end of the overlap)
    # in each iteration, add the wave & flux value from whichever array comes next in wavelength
    while(wave1_index <= len(wave1)):
        clear_output(wait=True)
        print('Index:\t%d/%d'%(wave2_index,max_ind))
        #If the current value of wave2 <= current value of wave1, insert wave2 into combined array
        #move on to next index of combined array
        #move on to next index of wave2
        #print('Wave 1\tindex: %d\tvalue: %0.1f'%(wave1_index,wave1[wave1_index]))
        #print('Wave 2\tindex: %d\tvalue: %0.1f'%(wave2_index,wave2[wave2_index]))
        if(wave2[wave2_index] < wave1[wave1_index]):
            # Debugging println
            #print('inserting '+str(wave2[wave2_index])+' at index '+str(combined_index))
            comb_wave[combined_index] = wave2[wave2_index]
            comb_flux[combined_index] = flux2[wave2_index]
            comb_err[combined_index] = err2[wave2_index]
            comb_flag[combined_index] = flag2[wave2_index]
            combined_index += 1
            wave2_index +=1
            #print('Selected wave 2')
        #If the current value of wave2 > current value of wave1, insert wave1 at current combined array location
        #increase current combined array location
        #increase wave1_index
        else:
            #Debugging println
            #print('inserting '+str(wave1[wave1_index])+' at index '+str(combined_index))
            comb_wave[combined_index] = wave1[wave1_index]
            comb_flux[combined_index] = flux1[wave1_index]
            comb_err[combined_index] = err1[wave1_index]
            comb_flag[combined_index] = flag1[wave1_index]
            combined_index += 1
            wave1_index += 1
            #print('selected wave 1')
            if(wave1_index == len(wave1)):
                break
    
    #Keep second array from end of overlap region onward
    comb_wave[combined_index::] = wave2[max_ind::]
    comb_flux[combined_index::] = flux2[max_ind::]
    comb_err[combined_index::] = err2[max_ind::]
    comb_flag[combined_index::] = flag2[max_ind::]
    
    return np.array([comb_wave,comb_flux,comb_err,comb_flag])

def combine_all(data):
    combined_data = data[0]
    for n in range(1,len(data)):
        combined_data = combine(combined_data,data[n])
    return combined_data

def convolve_gaussian(wave,flux,target_width,current_res=None):
    # Interpolate data to linear grid for convolution (eliminates strange behavior at edges)
    x_orig = wave
    y_orig = flux
    f_linear = interpolate.interp1d(x_orig,y_orig)

    x_new = np.linspace(min(wave),max(wave),num=len(wave))
    y_new = f_linear(x_new)

    # Convolve with Gaussian
    # If current_res given, use it. Preferred ( look up average in STIS )
    # otherwise, current_res = (max-min)/len
    if(current_res==None):
        current_res = (max(x_new)-min(x_new))/len(x_new)
    
    FWHM = target_width/current_res # target FWHM of Gaussian in current grating pixels

    std = FWHM/2.355 # convert to std for astropy Gaussian1DKernel
    
    flux_conv = convolve_fft(y_new,Gaussian1DKernel(std))
   
    return [x_new,flux_conv]

# Low should be the array [wave,flux] for the most redward data, mid the middle data, high the bluest data
# settl = model spectrum
# data MUST overlap
def create_combined_spectra(low,mid,high,settl):
    comb_low_wave,comb_low_flux = combine(low[0],mid[0],low[1],mid[1])
    comb_med_wave,comb_med_flux = combine(comb_low_wave,high[0],comb_low_flux,high[1])
    combined_wave,combined_flux = combine(comb_med_wave,settl[0],comb_med_flux,settl[1])
    
    return [combined_wave,combined_flux]

def get_overlap(array1,array2):
    #nonzero(logical argument)[0][0] returns the location of the first nonzero value for the given bool argument
    #ie: first index where array > value is true
    min_ind = np.nonzero(array1 >= min(array2))[0][0]
    max_ind = np.nonzero(array2 >= max(array1))[0][0]
    return min_ind,max_ind


def get_spectra(path,muscles=False,wave_range=None,get_error=True):
    HDU = fits.getdata(path,1)
    if(muscles==True):
        wave = HDU['WAVELENGTH']
        flux = HDU['FLUX']
        if(get_error):
            error = HDU['ERROR']
    else:
        wave = HDU['WAVELENGTH'][0]
        flux = HDU['FLUX'][0]
        if(get_error):
            error = HDU['ERROR'][0]
    if(wave_range != None):
        ind = (wave > wave_range[0]) & (wave < wave_range[1])
        if(get_error):
            return wave[ind],flux[ind],error[ind]
        return wave[ind],flux[ind]
    else:
        if(get_error):
            return wave,flux,error
        return wave,flux
    
    
# Replace values in array (orig) with given values (new)
def replace(orig,new):
    # Get indices where the line region begins & ends
    low_ind = np.where(orig[0] > min(new[0]))[0][0]
    hi_ind = np.where(orig[0] > max(new[0]))[0][0]
    
    # Split the spectrum around the line region
    low_arr = orig[:,:low_ind]
    hi_arr = orig[:,hi_ind::]
    
    # Concatenate with the line region in the middle
    new_spec = np.concatenate((low_arr,new,hi_arr),axis=1)
    return new_spec



