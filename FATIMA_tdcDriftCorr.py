#!/usr/bin/python2
"""
Creates a file containing offset parameters for FATIMA from ref time
histogram (TH2D).
Needs a root file as input which has a TH2D histogram with run time in
minutes on the X axis and reference time spectrum in ps on the Y axis.

Arguments:
  (1) Root file name
  (2) Name of the TH2D
  (3) Value to which the cuts are aligned to using the offset
  (4) width of the timestamp gate on X axis (in minutes)
  (5) First time stamp (in minutes, optional! default: first timestamp)
  (6) Last time stamp (in minutes, optional! default: last timestamp)

The Cut range in the projected histograms can be adjusted via the
dTRANGE_L and dTRANGE_R parameters below.

"""

import ROOT
from time import sleep
import os.path as path
import os
import sys
import re   #import regex lib
import datetime
import matplotlib.pyplot as plt
import math
import numpy as np

dTRANGE_L = -2000  #low limit on projected dt hist in ps
dTRANGE_R = 2000   #upper limit on projected dt hist in ps

def getoffset(inth2d, startbin, rangebins, alignT = 0):
  """
  Extracts 1D histogram from 2D histogram by projection on Y axis

  Args: inth2d - TH2D root histogram
        startbin  - lower limit of the cut on X axis (in ps)
        rangebins - width of the cut window for the projection (in ps)

  Returns: the resulting histogram as X and Y arrays, and the offset in ps
  """
  maxbin = inth2d.GetXaxis().GetBinLowEdge(inth2d.GetNbinsX()) + inth2d.GetXaxis().GetBinWidth(0)
  if (startbin < 0 or startbin > (maxbin - 1)):
    print "Trying cut, but is out of range. Check input parameter of reference run/ts!"
    exit(0)
  if (startbin + rangebins > maxbin):
    endbin = maxbin
  else:
    endbin = startbin + rangebins
  inth2d.GetXaxis().SetRangeUser(startbin, endbin)
  hdt = inth2d.ProjectionY()
  hdt.GetXaxis().SetRangeUser(dTRANGE_L, dTRANGE_R)
  offset = hdt.GetMean() - alignT
  X = []
  Y = []
  for i in range(0, hdt.GetNbinsX()):
    X.append(hdt.GetBinLowEdge(i) - offset)
    Y.append(hdt.GetBinContent(i))
  return X, Y, offset

def get_int (xarray, yarray):
  tot_int = 0
  for i,p in enumerate(yarray):
    if (xarray[i] > dTRANGE_L and xarray[i] < dTRANGE_R):
      tot_int += p
  return tot_int

def scale_y1_to_y2 (x1, y1, y2):
  """
  scales y1 to y2 by the number of counts in [SCALERANGE_L:SCALERANGE_R]
  and returns the result
  """
  ytmp = []
  refint = get_int(x1, y2)
  if refint < 1:
    print "ERROR: Bad reference histogram! Integral is smaller than 1. Check parameters!"
    exit(-1)
  tmpint = get_int(x1, y1)
  if tmpint < 1:
    return y1
  factor = refint/tmpint
  for i,y in enumerate(y1):
    ytmp.append(y*factor)
  return ytmp

def rebin_XY (x1, y1, fact):
  xout = []
  yout = []
  j = 0
  for i,x in enumerate(x1):
    if i>0 and i%fact == 0:
      xout[j] = x-fact/2.0
      yout[j] = ytot
      ytot = 0
      j+=1
  return xout, yout


############
## M A I N :
############


#initialising some things:

firstrun = 0
lastrun = 1

#Checking input parameters:
if( len(sys.argv) > 1 and len(sys.argv) >= 6):
  ofilename = sys.argv[1]
  rfilename = sys.argv[2]
  th2dname = sys.argv[3]
  refT = float(sys.argv[4])
  cutrange = int(sys.argv[5])
  if len(sys.argv) >= 7:
    firstrun = int(sys.argv[6])
  if len(sys.argv) >= 8:
    lastrun = int(sys.argv[7])
else:
  print "\n\nDetermines gain matching parameters from histograms. Needs a rootfile with TH2D"
  print "containing the data sorted by run or by timestamp in minutes."
  print "The TH2D name has to follow the naming convention:\n<somename><det#>, e.g. S480gm_det0 or gainmatch0"
  print "Usage:\n---------\n"
  print "Arguments in <> are mandatory, arguments in {} are optional"
  print "By time stamp range:"
  print "python2 ",sys.argv[0]," <outfilename.txt> <file.root> <TH2D name base> <value to shift to> <ts window in min> {lower time stamp limit} {upper time stamp limit}"
  exit(0)


#Set up the root environment, get the TH2D etc.
infile = ROOT.TFile.Open(rfilename, "read")
#th2dname = "{0}{1:02d}".format(th2dnamebase, detn)
gmTH2D = infile.Get(th2dname)
try:
    gmTH2D.GetName()
except ReferenceError:
    print "ERROR: Unable to extract histogram."
    print "       There is no histogram '{0}' in the given root file.\n".format(th2dname)
    exit(-1)
print "TH2D name is", gmTH2D.GetName()

#Setting up the plotting window:
plt.subplot (2,1,2)
plt.xlabel("dt (ps)")
plt.ylabel("counts")
plt.subplots_adjust(hspace=0.3)
plt.xlim(-4000,4000)
#Getting all the gainmatching parameters
runarray = []
##First check if the set ts range is valid:
if firstrun >= lastrun:
  print "Given time stamp range is invalid!"
  exit(-1)
maxts = int(math.ceil(gmTH2D.GetNbinsX()*gmTH2D.GetXaxis().GetBinWidth(0)))
maxsteps = int(maxts/float(cutrange))
if firstrun < gmTH2D.GetXaxis().GetBinLowEdge(1) or firstrun > maxts:
  firstrun = int(gmTH2D.GetXaxis().GetBinLowEdge(1)/float(cutrange))
else:
  firstrun = int(firstrun/float(cutrange))
if lastrun < gmTH2D.GetXaxis().GetBinLowEdge(2) or lastrun > maxts:
  lastrun = maxsteps
else:
  lastrun = int(lastrun/float(cutrange)+0.5) #+0.5 to catch potential last incomplete cut
for r in range(firstrun, lastrun):
  runarray.append(int( r*cutrange))

f = open(ofilename, 'w+')
X = []
Y = []
Y_ref = []
ts = []
Toffset = 0
offsetList = []

# Find range with best statistics to use for ref
# This is only used for plotting so that alignement can be checked easily.
maxint_i = 0
maxint = 0
for i, run in enumerate(runarray):
  X = []
  Y = []
  X,Y,Toffset = getoffset(gmTH2D, run, cutrange, refT)
  tmpint = get_int(X, Y)
  if maxint < tmpint:
    maxint = tmpint
    maxint_i = i
X,Y_ref,Toffset = getoffset(gmTH2D, runarray[maxint_i], cutrange, refT)

scaling = 1.0
for run in runarray:
  ts.append(run)
  X = []
  Y = []
  X,Y,Toffset = getoffset(gmTH2D, run, cutrange, refT)
  ts_range = "TSRANGE_{}_{}".format(run, run+cutrange)
  f.write("{} {} {:20}\n".format(run, run + cutrange, Toffset))
  if get_int(X,Y) < 100:
    offsetList.append(0)
    continue
  offsetList.append(Toffset)
  Y = scale_y1_to_y2(X,Y,Y_ref)
  plt.plot(X,Y)

f.close()

#Plotting the gm results vs respective timestamp
plt.subplot(2,1,1)
plt.xlabel("timestamp (min)")
plt.ylabel("dt offset (ps)")
plt.xlim([gmTH2D.GetXaxis().GetBinLowEdge(1)-cutrange, gmTH2D.GetXaxis().GetBinLowEdge(gmTH2D.GetNbinsX())+2*cutrange])
plt.ylim([min(offsetList)-100,max(offsetList)+200])
plt.plot(ts, offsetList, 'r*')
plt.savefig("sum.png")
plt.show()
