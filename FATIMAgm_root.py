#!/usr/bin/python2
"""
Creates a file containing gain shift parameters for FATIMA histograms.
Needs a root file as input which has a TH2D histogram with E on Y axis
and some parameter on X which is used to track (run number, timestamp)
SCRIPT
command line options:
  run:  It will interpret the x-axis as run numbers and determine
        gm factors for projection spectra from each bin.
  ts:   It will interpret the x-axis as time stamp value with unit minutes and determine
        gm factors for projections using the provided ts-range (in minutes).

Execute as follows ( <> are mandatory, {} are optional to set a run/ts
                     range, otherwise full range is used):
run:
  python2 FATIMAgm_root.py run <file.root> <TH2D name base> <detector number> <reference run number> {first run} {last run}
ts:
  python2 FATIMAgm_root.py ts  <file.root> <TH2D name base> <detector number> <ts window in min> <start of ref window in min> {first time in min}, {last time in min}


Output: Creates a file with gm factors per line. In case of run, first column
        is the run number, in case of ts, first column has information on the
        detector number as well as on the time stamp range where the factor is valid.

Note:   The script currently will only find a proportional factor for the gain shift
        (linear gain, no offset). Higher order polynomial shifts might be necessary,
        in this case the code should be extended with a more complicated
        multi-dimension minimisation algorithm.
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

SCALERANGE_L = 200   #low energy limit in keV
SCALERANGE_R = 3000  #upper energy limit in keV
BINNING = 1 #is adjusted to actual binning in the code. No need to change.


def gethist(inth2d, xarray, yarray, startbin, rangebins, boolrun, rebin=0):
  """
  Extracts 1D histogram from 2D histogram and changes xarray and yarray
  xarray -> bin unit, yarray -> bincontent

  Args: inth2d - TH2D root histogram
        xarray, yarray - empty arrays to be filled
        startbin - left bin of projection range on x in TH2D
        rangebins - width of the cut window for the projection
        rebin - TODO, add option for rebinning to account for high res digitisers

  Returns: nothing, changes mutable objects xarray and yarray
  """
  if(boolrun):
    maxbin = inth2d.GetNbinsX()
  else:
    maxbin = inth2d.GetXaxis().GetBinLowEdge(inth2d.GetNbinsX()) + inth2d.GetXaxis().GetBinWidth(0)
  print inth2d.GetName(), startbin, rangebins, maxbin
  if (startbin < 0 or startbin > (maxbin - 1)):
    print "Trying cut, but is out of range. Check input parameter of reference run/ts!"
    exit(0)
  if (startbin + rangebins > maxbin):
    endbin = maxbin
  else:
    endbin = startbin + rangebins

  inth2d.GetXaxis().SetRangeUser(startbin, endbin)
  histo = inth2d.ProjectionY()

  if (rebin < 2):
    for i in range (1, histo.GetNbinsX()):
      xarray.append(float(histo.GetBinLowEdge(i)))
      yarray.append(float(histo.GetBinContent(i)))
  else:
    #TODO: MR, add support for rebinning
    print "Rebinning not supported yet!"
    exit(0)
  return

def get_int (yarray, xmin=0, xmax=4000):
  tot_int = 0
  for i,p in enumerate(yarray):
    if (i>xmin and i < xmax):
      tot_int += p
  return tot_int

def get_max (yarray, xmin=0, xmax=4000):
  ymax = 0
  for i,p in enumerate(yarray):
    if (i>xmin and i < xmax):
      ymax = max(ymax, p)
  return ymax

def get_mean (yarray, xmin=0, xmax=4000):
  tot_mean = 0
  tot_int  = 0
  for i,p in enumerate(yarray):
    if (i>xmin and i<xmax):
      tot_mean += i*p
      tot_int  += p
  return tot_mean/tot_int

def scale_y1_to_y2 (y1, y2):
  """
  scales y1 to y2 by the number of counts in [SCALERANGE_L:SCALERANGE_R]
  and returns the result
  """
  ytmp = []
  refint = get_int(y2,SCALERANGE_L,SCALERANGE_R)
  if refint < 1:
    print "ERROR: Bad reference histogram! Integral is smaller than 1. Check parameters!"
    exit(-1)
  factor = refint/get_int(y1,SCALERANGE_L, SCALERANGE_R)
  for i,y in enumerate(y1):
    ytmp.append(y*factor)
  return ytmp


def apply_gm (factor, inx):
  outx = []
  for ix in inx:
    outx.append (factor*ix)
  return outx

def getFOM (xin, yin, yr):
  """
  Evaluates how close yin (input hist) and yr (reference hist)
  are by calculating a bin-wise square distance
  Returns: A number that is smaller if agreement is better. Also returns
    the interpolated histogram it used for the comparison.
  """
  tmpFOM = 100000
  tmpN = 0
  ycomp = 0
  yout = []
  xout = []
  a = 0
  b = 0
  for i, x in enumerate(xin):
    if i <SCALERANGE_L:
      continue
    if i > SCALERANGE_R:
      break
    iref = int(math.floor(x/BINNING))
    if iref > len(yr) - 1:
      break;
    a = (yin[i] - yin[i-1])/(xin[i] - xin[i-1])
    b = yin[i] - a*xin[i]
    ycomp = a*iref*BINNING + b
    tmpFOM += (yr[iref] - ycomp)**2
    tmpN += 1
    yout.append(ycomp)
    xout.append(iref*BINNING)
  FOM =np.sqrt(tmpFOM)
  #print("FOM:", FOM)
  return FOM, xout, yout

def find_gm_factor (thisth2d, leftbin, bins, inYref, boolrun):
  """
  Finds the best gm factor for a histogram produced from the input TH2D
  which is projected according to the input parameters.

  Returns: The best gain matching factor found, as well as the corresponding hist (as x,y arrays)
  """
  Xhist = []
  Yhist = []
  Yhistscaled = []
  Xgm = []
  Ygm = []
  gethist(thisth2d, Xhist, Yhist, leftbin, bins, boolrun)
  if (get_int(Yhist, SCALERANGE_L, SCALERANGE_R) < 0.2*get_int(inYref, SCALERANGE_L, SCALERANGE_R)):
    print "Not enough statistics in cut range/run. Will keep current gm factor."
    return 0, Xgm, Ygm
  Yhistscaled = scale_y1_to_y2(Yhist, inYref)
  if(boolrun):
    print "Gain-matching for run", leftbin
  else:
    print "Gain-matching for ts range", leftbin, leftbin + bins

  notfound = 1
  gmfactor = 1
  fom      = 10000000
  last_fom = 20000000
  best_fom = fom
  best_gmfactor = 1
  turns = 0
  up = 1
  bad_rounds=0
  iterations = 0
  while notfound == 1:
    fom, Xgm, Ygm = getFOM(apply_gm(gmfactor, Xhist), Yhistscaled,inYref)
    #print("it ", iterations, "  gm:", gmfactor, "  FOM:", fom, "turns: ", turns, "change is ", (2*up - 1)*0.002/2**turns)
    if (fom >= best_fom):
      if(fom >= last_fom):
        bad_rounds += 1
      if(bad_rounds>3):
        if(turns>2):
          break
        turns += 1
        bad_rounds=0
        #print("  turning ", turns)
        up = bool(not up)
        gmfactor = best_gmfactor
    else:
      bad_rounds = 0
      best_fom = fom
      best_gmfactor = gmfactor

    last_fom = fom
    iterations+=1
    gmfactor += (2*up - 1)*0.002/2**turns

  print "   found best factor", best_gmfactor, "after", iterations, "iterations"
  fom, Xgm, Ygm = getFOM(apply_gm(best_gmfactor, Xhist), Yhistscaled,inYref)
  return best_gmfactor, Xgm, Ygm



############
## M A I N :
############


#initialising some things:
Xref = []
Yref = []

Xhist = []
Yhist = []
Yhistscaled = []

Xgm1 = []
Ygm1 = []

ts = []
gmbyts = []

bestgm = 1

runs = 0

firstrun = 0
lastrun = 1

#Checking input parameters:
if( len(sys.argv) > 1 and sys.argv[1] == "run" and len(sys.argv) >= 6):
  runs = 1
  rfilename = sys.argv[2]
  th2dnamebase = sys.argv[3]
  detn = int(sys.argv[4])
  cutrange = 1
  refbin = int(sys.argv[5])
  if len(sys.argv) >= 7:
    firstrun = int(sys.argv[6])
  if len(sys.argv) >= 8:
    lastrun = int(sys.argv[7])

elif( len(sys.argv) > 1 and sys.argv[1] == "ts" and len(sys.argv) >= 7):
  runs = 0
  rfilename = sys.argv[2]
  th2dnamebase = sys.argv[3]
  detn = int(sys.argv[4])
  cutrange = int(sys.argv[5])
  refbin = int(sys.argv[6])
  if len(sys.argv) >= 8:
    firstrun = int(sys.argv[7])
  if len(sys.argv) >= 9:
    lastrun = int(sys.argv[8])
else:
  print "\n\nDetermines gain matching parameters from histograms. Needs a rootfile with TH2D"
  print "containing the data sorted by run or by timestamp in minutes."
  print "The TH2D name has to follow the naming convention:\n<somename><det#>, e.g. S480gm_det0 or gainmatch0"
  print "Usage:\n---------\n"
  print "Arguments in <> are mandatory, arguments in {} are optional"
  print "By run:"
  print "python2 ",sys.argv[0]," run <file.root> <TH2D name base> <detector number> <reference run number> {lower run limit} {upper run limit}"
  print "Example:"
  print "python2 code.py run gmfiles.root gainmatching 0 10"
  print "-------------------------------------------------------\n"
  print "By time stamp range:"
  print "python2 ",sys.argv[0]," ts <file.root> <TH2D name base> <detector number> <ts window in min> <start of ref window in min> {lower time stamp limit} {upper time stamp limit}"
  print "Example:"
  print "python2 code.py run gmfiles.root Histograms/FATIMA/gm/gainmatching 0 20 900"
  exit(0)


#Set up the root environment, get the TH2D etc.
infile = ROOT.TFile.Open(rfilename, "read")
th2dname = "{0}{1:02d}".format(th2dnamebase, detn)
#th2dname = "{0}{1:02d}_energy_vs_time".format(th2dnamebase, detn)
gmTH2D = infile.Get(th2dname)
try:
    gmTH2D.GetName()
except ReferenceError:
    print "ERROR: Unable to extract histogram."
    print "       There is no histogram '{0}' in the given root file.\n".format(th2dname)
    exit(-1)
print "TH2D name is", gmTH2D.GetName()

#Getting the reference histogram (and plotting it):
gethist(gmTH2D, Xref, Yref, refbin, cutrange, runs)
plt.subplot (2,1,2)
plt.xlabel("E (keV)")
plt.ylabel("counts")
plt.subplots_adjust(hspace=0.3)
plt.xlim(Xref[0], Xref[-1])
plt.plot(Xref, Yref, color='#E1B9A2', linestyle = '-', linewidth=7)
#Getting all the gainmatching parameters
runarray = []
##First check if the set run/ts range is valid:
if runs:
  if firstrun >= lastrun:
    print "Given run range is invalid!"
    exit(-1)
  if firstrun < 0 or firstrun > int(math.ceil(gmTH2D.GetNbinsX()/float(cutrange))):
    firstrun = 0
  else:
    firstrun = int(firstrun/float(cutrange))
  if lastrun < 1 or lastrun > int(math.ceil(gmTH2D.GetNbinsX()/float(cutrange))):
    lastrun = int(math.ceil(gmTH2D.GetNbinsX()/float(cutrange)))
  else:
    lasttrun = int(lasttrun/float(cutrange)+0.5) #+0.5 to catch potential last incomplete cut
  for r in range(firstrun, lastrun):
    runarray.append(r*cutrange)
else:
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

#Check histogram and binning:
if len(Xref) > 1:
  BINNING = Xref[1] - Xref[0]
else:
  print "Energy range is too small! Check the input histogram."
  exit(-1)
#Then check energy range value (set to bin numbers):
if BINNING > 1:
  SCALERANGE_L = int(SCALERANGE_L/BINNING)
  SCALERANGE_R = int(SCALERANGE_R/BINNING)

f = open("Det{0}.gmlist".format(detn), 'w')
prevgm = 1
detname = ""
for run in runarray:
  bestgm, Xgm1, Ygm1 = find_gm_factor(gmTH2D, run, cutrange, Yref, runs)
  ts.append(run)

  #if the fing_gm_factor doesn't return a valid array, keep the previous gm factor
  if (len(Xgm1) == 0):
    bestgm = prevgm
  else:
    prevgm =  bestgm

  gmbyts.append(bestgm)
  plt.plot(Xgm1, Ygm1)
  if(runs):
    f.write("{:25}  0  {:20}\n".format(run, bestgm))
  else:
    detname = "Det_{}_{}_{}".format(detn, run, run+cutrange)
    f.write("{:25}  0  {:20}\n".format(detname, bestgm))
f.close()

print "BINNING was {0}".format(BINNING)
#Plotting the gm results vs respective timestamp
plt.subplot(2,1,1)
if(runs):
  plt.xlabel("run number")
else:
  plt.xlabel("timestamp (min)")
plt.ylabel("gainmatching factor")
#if runs:
#  plt.xlim([firstrun -20, lastrun+20])
#else:
#  plt.xlim([gmTH2D.GetXaxis().GetBinLowEdge(1)-cutrange, gmTH2D.GetXaxis().GetBinLowEdge(gmTH2D.GetNbinsX())+2*cutrange])
plt.xlim(firstrun*cutrange - cutrange, lastrun*cutrange + cutrange)
plt.ylim([min(gmbyts)-0.02,max(gmbyts)+0.02])
plt.plot(ts, gmbyts, 'r*')
plt.savefig("gm_Det{0}.png".format(detn))
plt.show()
