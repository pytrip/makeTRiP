#!/usr/bin/env python
"""
=============================
Creating ddd files for TRiP98
=============================
This script takes scoring files with the name convention: XXX
from SHIELD-HIT10A and converts them into TRiP98 readable .ddd-format.
It strips the input file to reduce the data size and fits the radial files
to determine the FWHMs of the model Gaussians and the factor between the amplitudes.

The first Gaussian is intended to model the primary particle, the second for the
secondaries.


The fitting is done by:
leastsqbound.py
from:
http://pydoc.net/nmrglue/0.2/nmrglue.analysis.leastsqbound
which does a constrained multivariant Levenberg-Marquardt fit, since the SciPy package does not offer one.

WARINING: the input scoring file, must have linear bin spacing!

# Group: Aarhus Particle Therapy Group
# Date: August 2015
"""
import os
import sys
import glob
import math
import time
import datetime

from numpy import *
from pylab import * # for plotting
from scipy.optimize import leastsq # For fitting with Levenberg-Marquardt
from scipy import optimize
from copy import deepcopy

from leastsqbound import *


__authors__ = "Ricky Teiwes, Niels Bassler, Toke Printz Ringbaek"
__copyright__ = "Copyright 2012, Aarhus Particle Therapy Group"
__credits__ = ["APTG"]
__license__ = ""
__version__ = "0.2"
__maintainer__ = ""
__email__ = ""
__status__ = "Production"


def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well

        Parameter
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)


########################### Strip the input files ###########################
# Strips the input file to reduce the amout of data for TRiP98

# It adds the 1. datapoint
# It adaptively chooses the datapoint in between
# it adds the last datapoint ->TODO Why choose the last datapoint? It's often uninteresting?! /TPR

# STRIPPING ALGORITHM #
'''
This algorithm chooses the most important datapoints
The algortihm works by summing the total differences, then choosing points where the difference
is over a certain threshhold
'''
def _strip_data(data, z_offset, z_binsize):
    """
    Strips the z data

    Parameters:
    * data
    * offset (z_offset)
    * spacing (z_binsize)
    
    pointsChoosen is the index of the points choosen
    
    Returns: (pointsChoosen,{cov_x,infodict,mesg},ier)
    """
        ## Start of algorithm ##
    stripped_x = [z_offset] # add the first datapoint to the stripped x data TODO: make it variable!
    stripped_y = [data[0]] # add the first datapoint to the stripped y data
    pointsChoosen = [0] # array that collect the pointnumber of the data points that were choosen for the lateral fits
    strippedData = []
    y_old = data[0] # set the old_y
    sum = 0 # set the difference to zero
    # sum over the differences
    for i in range(len(data)):
        sum = sum + math.fabs(data[i] - y_old)
        y_old = data[i]
    if sum == 0: # if there are no differences between the datapoint, add the last point and end
        #stripped_y.append("%1.5e" % float(data[len(data)-1])) #Better not use scientific notation since this is not used in TRiP /TPR
        stripped_y.append(float(data[len(data)-1]))
        stripped_x.append(z_offset + z_binsize * float(len(data)-1))
        pointsChoosen.append(len(data)-1)
    else:
        y_old = data[0] # set the first datapoint as y_old
        current_sum = 0
        for i in range(len(data)-1):# minus 1 so the last datapoint can not be double inside list
            current_sum = current_sum + math.fabs( float(data[i]) - float(y_old) )
            if ( current_sum >=  sum/79 ): # is the current sum of differences is bigger, add a datapoint and reset the current sum
                #stripped_y.append("%1.5e" % float(data[i])) #Better not use scientific notation since this is not used in TRiP /TPR
                stripped_y.append(float(data[i])) 
                stripped_x.append(float(z_offset + z_binsize * i))
                pointsChoosen.append(i)
                current_sum = 0
                y_old = float(data[i])
        #stripped_y.append("%1.5e" % float(data[len(data)-1])) #Better not use scientific notation since this is not used in TRiP /TPR
        stripped_y.append(float(data[len(data)-1]))
        stripped_x.append(z_offset + z_binsize * float(len(data)-1))
        pointsChoosen.append(len(data)-1)
    stripped_y = [round(x, 3) for x in stripped_y] #truncate dose data instead of scientific notation /TPR 
    #TODO: truncation method above a bit slow but for some reason a "%.3f" % stripped_y approach didn't work! /TPR  
    # find location of Bragg peak
    BP_index = next(x[0] for x in enumerate(stripped_y) if x[1] == max(stripped_y))
    ##
    ## Write the stripped data to file: stripped
    stripped = open('stripped', "w")
    for i in range(len(stripped_y)):
        stripped.write(' ' + str(stripped_x[i] ) + ' ')#
        stripped.write( str(stripped_y[i]) + '\n' )#
        strippedData.append(str(stripped_x[i]) + " " + str(stripped_y[i]))# TODO: Not very nice, but it works -> what does this TODO relates to so I can fix it?!? /TPR
    stripped.close()
    ##
    return pointsChoosen, strippedData, BP_index #added BP_index in output /TPR

def _write_ddd_file(energy, strippedData, lateralFits, BP_index): 
    
    ### Get ion species from config file: /NBassler 2013
    fc = open("../config.sh",'r')
    fl = fc.readlines()
    for line in fl:
        if line[0:3] == "ion":
            ion = line.split("=")[1].strip()
    fc.close() 

    ## Find the cutoff index for the FWHM values ## /TPR 2015
    FWHM_cutoff = len(strippedData) # set cutoff value to last index, which will be used if the cutoff requirement below is never reached /TPR
    for k in range(len(strippedData)):
        if (lateralFits[k][2] >= 4*lateralFits[k][0]) and (k > BP_index): #up until this point, we assume the FWHM2 values make physical sense /TPR
            break #after that point the FWHM2 values tend to become nonsensically large and should no longer be added to the list /TPR
        FWHM_cutoff=k #define index for cutoff value to be used below /TPR     
    
    # make header with current time
    t=  datetime.datetime.now()
    now = datetime.datetime.fromtimestamp(time.mktime(t.timetuple()))
    ''' Make the header '''
    header1 = '!filetype    ddd\n'
    header2 = '!fileversion    '+t.strftime("%Y%m%d")+"\n"
    header3 = '!filedate    ' + now.ctime() + '\n'
    header4 = '!projectile    ' + ion + '\n'
    header5 = '!material      H2O\n'
    header6 = '!composition   H2O\n'
    header7 = '!density 1\n'
    header8 = '!energy ' + str(energy) + '\n'
    header9 = '#   z[g/cm**2] dE/dz[MeV/(g/cm**2)] FWHM1[mm] factor FWHM2[mm]\n' #changed units for FWHM values to mm /TPR
    header10 = '!ddd\n'
    #
    ## Write the DDD file##
    output_path = "./DDD/"
    _mkdir(output_path)
    ddd_filename = output_path  + ion + '.H2O.MeV' + str(energy) + '00.ddd'
    output_file = open(ddd_filename, "w")
    output_file.write(header1 + header2 + header3 + header4 + header5 + header6 + header7 + header8 + header9 + header10)
    for k in range(len(strippedData)):
	if (lateralFits[k][1] == lateralFits[k][1] and lateralFits[k][1] > 0): # Test if the factor is not equal to NaN and non-zero! Skip the current line otherwise.
            if k < FWHM_cutoff: #up til the cutoff value, we add FWHM values normally from the respective lists /TPR
                output_file.write(str(strippedData[k]) + " " + str(lateralFits[k][0]) + " " + str(lateralFits[k][1])+ " " + str(lateralFits[k][2]) + "\n")
            elif k >= FWHM_cutoff: #from this point, we overwrite the FWHM1, factor and FWHM2 values with the last ones before the cutoff instead of writing out the nonsense values given by the fits /TPR 
                output_file.write(str(strippedData[k]) + " " + str(lateralFits[FWHM_cutoff][0]) + " " + str(lateralFits[FWHM_cutoff][1])+ " " + str(lateralFits[FWHM_cutoff][2]) + "\n")
    output_file.close() 
    #
    ## Write file with Bragg peak indices to file: BPs.txt ## /TPR 2015
    BPs_file = open('BPs.txt', "a")
    previous = open('BPs.txt', "r").read()
    if previous == '': # Add a header only if the file is empty
        BPs_file.write(header2 + header3 + '!Bragg peak location and stopping power value' + "\n" + '# energy[MeV] z[g/cm**2] dE/dz[MeV/(g/cm**2)]' + "\n")
    BPs_file.write(str(energy) + " " + str(strippedData[BP_index]) + "\n" )  
    BPs_file.close()  
    #designer note: the BPs.txt data can be used to plot z-E relation to edit SIS tables to match the generated DDD files!

# Method fitting the lateral dose distribution
def _lateral_fit(left_radii, radial_dose_scaled, depth, energy):
    """
    Fitting for lateral dose distribution
    Fits the lateral dose scaled with the radius using a double Gaussian by the use of Levenberg-Marquardt optimization

    Parameters:
    * radii        radius values corresponding to dose
    * dose         dose D(r)
    * dose_scaled  dose times the radii D(r)*r
    * depth        only used for graphing
    * energy       only used for graphing
    
    Returns: (x,{cov_x,infodict,mesg},ier)

    """
    def residuals(p, y, x): # objective function for minimization (vector function)
        err = y - peval(x,mu,p)
        return err
    
    def peval(x, mu, p): # The model function: Two gaussians (sum), the mean mu is locked
        #In guideline with a similar routine for TRiP by U. Weber, the function is scaled with the radius to fit D(r)*r /TPR
        return ( p[0]*exp(- ( (x - mu)**2 )/( 2*p[1]**2) ) + p[2]*exp(- ( (x - mu)**2 )/( 2*p[3]**2) ) ) * x

    x = left_radii 
    y = radial_dose_scaled
    ## Calculate initial guesses ##
    a1 = y.max() # guess for the amplitude of the first gaussian
    mu = 0 # set the mean to zero - a guess would be: sum(x*y)/sum(y) 
    if (sum(y) == 0): # check if denominator is 0
        sigma1 = 1e5
    else:
        sigma1 = sqrt(abs(sum((x-mu)**2*y)/sum(y))) # guess for the deviation
    a2 = deepcopy( a1 ) * 0.05 # guess for the  amplitude of the second gaussian
    sigma2 = 2 # guess for the deviation of the second gaussian
    sigma2 = deepcopy( sigma1 ) * 5.0
    # Collect the guesses in an array
    pname = (['a1','sigma1','a2','sigma2'])
    p0 = array([a1,sigma1,a2,sigma2])
    bounds = [(0, 1e6),(0, 1e6),(0, 1e6),(0, 1e6)] # set bound for the parameters, they should be non-negative
    plsq  = leastsqbound(residuals, p0,bounds, args=(y,x), maxfev=5000) # Calculate the fit with Levenberg-Marquardt
    pfinal = plsq[0] # final parameters
    covar = plsq[1] # covariance matrix
 
    ''' Plot the fit '''
    _mkdir('./pics/')
    figure
    xnum = len(x) 
    xplot = linspace(x[0],x[xnum-1],num=xnum*10)
    plot(x,y,'r.',xplot,peval(xplot,mu,plsq[0]),'b-') # plot the functional fit on a finer grid than data points
    title('Energy: ' + str(energy) + ' , Depth: ' + str(depth))
    xlabel('radius [g/cm**2]')
    ylabel('dose [MeV/(g/cm**2)]')
    print depth
    savefig('./pics/ddd_e_' + str(energy) + '_depth_' + str(depth) + '.png')
    if energy == 40 and depth == 251:
        show()
        close()
    else:
        close()
    close()

    # Calculate the output parameters
    # FWHM = sqrt(8*log(2)) * Gaussian_Sigma
    FWHM1 = 2.354820045 * pfinal[1] # Full Width at Half Maximum for the first Gaussian
    FWHM2 = 2.354820045 * pfinal[3] # Full Width at Half Maximum for the second Gaussian
    # Scale from cm to mm according to Gheorghe from Marburg
    FWHM1 *= 10.0
    FWHM2 *= 10.0
    #Get the amplitudes of the normalized form of the double Gaussian /TPR
    A1 = pfinal[0] * 2*math.pi*FWHM1**2 
    A2 = pfinal[2] * 2*math.pi*FWHM2**2 
    normfactor = A1+A2 #changed to use of the normalized amplitudes /TPR
    if (normfactor == 0): # check if denominator is 0
#        factor = 1.0
	print "normfactor = ", normfactor, energy, depth
        factor = -1e6   # Calculate the factor between the two amplitudes TODO: CHECK THAT THIS IS THE RIGHT ->#what's going on here?!?! /TPR
    else:
        factor = A2/normfactor # Calculate the factor between the two normalized amplitudes: factor = A2 / (A1+A2), This definition was given by Gheorghe and Uli from Marburg!

    return FWHM1, factor, FWHM2


def _compute_ddd(inpFile,density):
    # 1D filename density scoring_radius (x_length) z_offset z_binsize number_of_z_datapoints number_of_radial_datapoints
    input_file_z = inpFile
    energy_string = input_file_z[8:11]
    input_file_radial = 'DDD_Cr_' + energy_string + '.dat'
    energy = int(energy_string)
    ##### Read the datafiles ##### 
    # 1D data z-direction
    file = open(input_file_z, "r")
    zRaw = file.readlines()
    file.close()
    # 2D data radial
    file = open(input_file_radial, "r")
    radialRaw = file.readlines()
    file.close()    
    # read the parameters (done from radial file)
    line2 = radialRaw[1] # line 2 - contains the BIN NR
    x_bin_nr = int( line2.split()[3] ) # x bin nr
    z_bin_nr = int( line2.split()[9] )
    line5 = radialRaw[4] # line 5 - contains the END BIN  
    x_length = float( line5.split()[4].replace('D','E') ) # length of x [cm] replace D with E in data #changed index from 3 to 4 /TPR
    z_length = float( line5.split()[12].replace('D','E') ) # length of z [cm] replace D with E in data #changed index from 9 to 12 /TPR 
    # Calculate the parameters
    x_binsize = x_length/x_bin_nr # size of x bins [cm]
    x_offset = 0.5 * x_binsize # middle of x_bin
    z_binsize = z_length/z_bin_nr # size of z bins [cm]
    z_offset = 0.5 * z_binsize # middle of z_bin   
    number_of_z_datapoints = z_bin_nr
    number_of_radial_datapoints = x_bin_nr
    ########################### Convert z datafiles ##########################
    area = x_length * x_length * math.pi # TODO: Important to be right! #x_length is the scoring radius, #area is a factor used to convert the input files
    orig_dd = zRaw[6:len(zRaw)] # not len -1 here since it should start in line 7
    new_orig_dd = [] #empty list to append below /TPR
    for i in range(number_of_z_datapoints): #making list of tuple with all 4 data columns into list of integers with only relevant column /TPR
        value = orig_dd[i].split()[3] 
        new_orig_dd.append(value)
    depthDoseData = [(float(x) * area) for x in new_orig_dd] # Store the depth dose data as floats and convert to the right unit #changed to new_orig_dd /TPR
    pointsChoosen, strippedData, BP_index = _strip_data(depthDoseData, z_offset, z_binsize); # choose the points, #added BP_index /TPR
    ########################### Convert radial datafiles ##########################
    radialData = radialRaw[6:len(radialRaw)] # not len -1 here # CHECKED!!! start in line 7 of radial input file
    new_radialData = [] #empty list to append below /TPR
    for i in range(number_of_z_datapoints*number_of_radial_datapoints): #making list of tuple with all 4 data columns into list of integers with only relevant column like above /TPR
        value = radialData[i].split()[3] 
        new_radialData.append(value)
    ##### Put data in nice format #####
    data = zeros([number_of_z_datapoints,  number_of_radial_datapoints], 'f') # float array that holds the data
    ''' Create a 2D-array with i-index for depth, j-index for radius '''
    for i in range(number_of_z_datapoints):   # iterate through depths
        for j in range( number_of_radial_datapoints ): # iterate through width
            data[i][j] = float( new_radialData[i*number_of_radial_datapoints+j] ) #changed to new_radialData /TPR
    lateralFits =[]
    for i in pointsChoosen: # where i is the index of the point choosen
        data1 = data[i][:] # right side of the data
        radial_dose = density*data1
        #radial_dose = density * concatenate((data1[::-1],data1)) # no concatenate needed when fitting D(r)*r functions /TPR
        left_radii = density * arange(x_offset,x_length,x_binsize) # construct a vector containing the positive radii 
        #radii = concatenate((-1 * left_radii[::-1],left_radii)) # arange(1,11,3) # no concatenate needed when fitting D(r)*r functions /TPR
        radial_dose_scaled = radial_dose * left_radii #get the D(r)*r function to be fitted /TPR
        fitParameters = _lateral_fit(left_radii, radial_dose_scaled, i, energy) # Fit the lateral data
        lateralFits.append(fitParameters)
    # FakeItTillYouMakeIt
    fitparams = open('fitparams', "w")
    for i in range(len(lateralFits)):
        fitparams.write( str(lateralFits[i][0]) + " " + str(lateralFits[i][2]) + " " + str(lateralFits[i][1]) + '\n' )#
    fitparams.close()
    
    _write_ddd_file(energy_string, strippedData, lateralFits, BP_index); # Write the ddd-file, #added BP_index /TPR

################### SCRIPT START ####################
print '========================================================'
print '========     ddd.py - converts to TRiP98        ========'
print '========================================================'
print 'Copyright (c) 2012'
print 'Dept. of Experimental Clinical Oncology,'
print 'Aarhus University Hospital, Aarhus, Denmark.'
print ''
print 'Build date: March 2015   XXX'
print ''
# Parameters
density = float(1) # density of target material
input_path = '' # 
for input_file in glob.glob( os.path.join(input_path, 'DDD_Cr1_*.dat') ):
    print "current file is: " + input_file
    _compute_ddd(input_file,density)
print '========================================================'
