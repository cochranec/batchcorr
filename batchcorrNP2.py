#!/usr/bin/python
# batchcorrNP
# Reproduces the functionality of batchcorr.exe for dark-correcting files produced
# using the a-Si detector at 1-ID at APS, using Python.
# Current version will run the correction on ALL *.ge2 files in the current working directory.
# Written by Chris Cochrane (cochranec@gmail.com), Dec. 2012
# Current version - Mar. 10, 2013

# Functionality presently requires the NumPy library - doing the file loading and
# corrections would take up to 100 times longer using the native Python functions.
# NumPy can be downloaded from: http://www.numpy.org/
# Written using NumPy version 1.5.1, on Python 2.7.1+
# Not tested on Python 3.x; send any feedback to cochranec@gmail.com

import numpy
import os
import re
import glob
import sys
import argparse

outDir = './'

# USER MUST DEFINE PATH TO BAD PIXEL INFORMATION
# Typically 'C:\DetectorData\1339.6\Full\1339.6Full_BadPixel_d.txt.img'
# Hard-coding this location in is inadvisable, as it makes the code less portable.
# Either a relative or an absolute path can be used.    It's probably safer to use an absolute path.
# (NB: the r prior to the string indicates a raw string, and must be included)
badPixFile = r'/home/chris/Python/batchCorr/GE3Bad.img'

# Command-line parser arguments - make everything more user friendly
parser = argparse.ArgumentParser(
    description='Dark correction and summing of GE2 files.',
    epilog='Written by Chris Cochrane, Dec. 2012. E-mail: cochranec@gmail.com')
parser.add_argument('--lo', type=int, nargs=1, default=0,help='Lower bound of run numbers. NOT YET IMPLEMENTED')
parser.add_argument('--hi', type=int, nargs=1, default=1,help='Upper bound of run numbers. NOT YET IMPLEMENTED')
parser.add_argument('--all','-a', action='store_true', default=True,help='Flag to perform dark correction on all GE2 files in present directory.')
parser.add_argument('--ndel', action='store_true', default=False, help='Print out .cor files.    Will produce a dark corrected output file for each frame in each GE2 file, as well as the sum files.    (WARNING: May use a LOT of space.)')
parser.add_argument('--drk', type=str, nargs=1, default='dark', help='Dark stub.    Some string that is unique to dark files.    Need not be the ENTIRE stub.    Default = "dark"')
clargs = parser.parse_args()

num_X = 2048
num_Y = 2048

lo = clargs.lo
hi = clargs.hi
ndel = clargs.ndel

allfiles = glob.glob('*[0-9].ge[1-4]')

#Recast the file reading function, which reduces runtime
fread = numpy.fromfile

#Preallocate memory space for arrays
sumvalues = numpy.zeros(num_X*num_Y,numpy.float32)
binvalues = numpy.zeros(num_X*num_Y,numpy.float32)
corrected = numpy.array(num_X*num_Y,numpy.float32)
darkvalues= numpy.array(num_X*num_Y,numpy.float32)
badPixels = numpy.array(num_X*num_Y,numpy.float32)

#Read in bad pixel data
try:
    with open(badPixFile, mode='rb') as badPxobj:
        badPxobj.seek(8192)
        badPixels = fread(badPxobj, numpy.uint16, num_X * num_Y)
except IOError as e:
    print '\nUnable to access bad pixel information at ' + badPixFile
    print 'Ensure that the file exists, or change the "badPixFile" variable on line 26 to direct to the file location.\n'
    sys.exit()

# Find dark files.
# The 'dark' stub can be anywhere in the filename (not necessarily at the start)
darks = [x for x in allfiles if clargs.drk.lower() in x.lower()]

print len(darks), "candidate(s) for dark file found."
if len(darks) == 0:
    print "Double-check that dark file is properly located."
    sys.exit()
for i in range(len(darks)):
    print "     (" + str(i+1) + ") " + darks[i]
print "Which dark file should we use?"
dkI = 0
while dkI < 1 or dkI > len(darks):
    try:
        dkI = int(raw_input().strip())
    except:
        dkI = 0
        print('Choose one of the available options. [1-' + str(len(darks)) + ']')
darkfile = darks[dkI-1]
print "Using", darkfile


# Read in dark file
# Average over all the exposures in the file
# This reduces the number of 'over reduced' pixels.
darkfile = darks[0]
with open(darkfile, mode='rb') as darkobj:
    statinfo = os.stat(darkfile)
    nFrames = (statinfo.st_size - 8192) / (2 * num_X * num_Y)
    darkobj.seek(8192)
    for i in range(nFrames):
        binvalues = fread(darkobj, numpy.uint16, num_X * num_Y).astype('float32')
        sumvalues = sumvalues + binvalues

darkvalues = sumvalues / nFrames
sumvalues[:] = 0

# Pixel data is stored as 0, 1, 2, 3
# Any pixel with a non-zero value is deemed 'bad'
badInd = numpy.array(numpy.where(badPixels == 2))
badInd1= numpy.array(numpy.where(badPixels%2== 1))

print "Dark file and bad pixel data read successfully."
#print "High val is",hi,"; low val is",lo

# Produce list of files to be dark corrected
if not clargs.all:
    print [int(re.findall('_([0-9]*).ge2',x)[0]) for x in allfiles]
    files = [x for x in allfiles if lo <= int(re.findall('_([0-9]*).ge([1-4])',x)[0]) <= hi and clargs.drk.lower() not in x.lower()]
else:
    files = [x for x in allfiles if clargs.drk.lower() not in x.lower()]

print len(files), "of", len(allfiles), "GE files in directory are in range. ", len([x for x in allfiles if clargs.drk.lower() in x.lower()]), "dark files ignored."

# Maybe use a try: construct here
c = raw_input('Perform dark correction on all available files? ([y]/n)').strip()
if c.lower() == 'n':
    print "No dark correction will be performed.    Terminating script."
    sys.exit()
else:
    print "Proceeding with dark correction."

#Perform a loop over all files
for f in files:
    statinfo = os.stat(f)
    nFrames = (statinfo.st_size - 8192) / (2 * num_X * num_Y)
    print "\nReading:",f, "\nFile contains", nFrames,"frames.    Summing and dark correcting."

    # Sum all values in this file
    with open(f, mode='rb') as fileobj:
        fileobj.seek(8192)
        for i in range(nFrames):
            binvalues = fread(fileobj, numpy.uint16,num_X * num_Y).astype('float32')
            sumvalues = sumvalues + binvalues
            if ndel:
                corName = outDir + f[:-3] + str(i) + '.cor'
                corSlice = binvalues - darkvalues
                # Correct for bad pixels by taking an average of nearest neighbours
                corSlice[badInd] = (corSlice[badInd + 1] + corSlice[badInd - 1] + corSlice[badInd + num_X] + corSlice[badInd - num_X]) / 4

                # Set border region and negative pixels to 0
                corSlice[badInd1] = 0
                corSlice[numpy.where(corrected<0)] = 0
                with open(corName, mode='wb') as outSlice:
                    corSlice.tofile(outSlice)
                print i, ',',
    # Remove the equivalent dark frame value
    corrected = sumvalues - darkvalues * nFrames

    # Correct for bad pixels by taking an average of nearest neighbours
    corrected[badInd] = (corrected[badInd + 1] + corrected[badInd - 1] + corrected[badInd + num_X] + corrected[badInd - num_X]) / 4

    # Set border region and negative pixels to 0
    corrected[badInd1] = 0
    corrected[numpy.where(corrected<0)] = 0

    sumName = outDir + f[:-3] + 'sum'
    print "Output sum to " + sumName
    with open(sumName, mode='wb') as outFile:
        corrected.tofile(outFile)

    sumvalues[:] = 0

print("Done")

