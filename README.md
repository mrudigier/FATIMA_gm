# FATIMA_gm
Python programme to track gain drifts which occur in long measurements using LaBr3 or other scintillator detectors. Used during the DESPEC FATIMA campaign 2020/2021.

Works with python 2.7.

# Description:
Creates a file containing gain shift parameters for FATIMA histograms.
Needs a root file as input which has a TH2D histogram with E on Y axis
and some parameter on X which is used to track (run number, timestamp).
Apart from the output file with the gain drift correction factors, the
results are also displayed graphically to immediately check if they
are reasonable.

## Command line options:
### run:
It will interpret the x-axis as run numbers and determine
gm factors for projection spectra from each bin.
### ts:
It will interpret the x-axis as time stamp value with unit minutes and determine
gm factors for projections using the provided ts-range (in minutes).

# Execute as follows
( <> are mandatory, {} are optional to set a run/ts
                     range, otherwise full range is used):
## By run:
python2  FATIMAgm_root.py  run \<file.root\> \<TH2D name base\> \<detector number\> \<reference run number\> {lower run limit} {upper run limit}
### Example:
python2 code.py run gmfiles.root gainmatching 0 10

------------------------------------------------------

## By time stamp range:
python2  FATIMAgm_root.py  ts \<file.root\> \<TH2D name base\> \<detector number\> \<ts window in min\> \<start of ref window in min\> {lower time stamp limit} {upper time stamp limit}
### Example:
python2 code.py run gmfiles.root Histograms/FATIMA/gm/gainmatching 0 20 900


## Output:
Creates a file with gm factors per line. In case of run, first column
is the run number, in case of ts, first column has information on the
detector number as well as on the time stamp range where the factor is valid.

## Working principle:
The code determines a reference histogram based on the user input parameters. It will then produce histograms from the 
input 2D-histogram by cutting with the given window width. These histograms are then compared to the reference histogram
bin by bin. By applying a linear caliration this bin-by-bin difference is minimised. The best found linear coefficient is
kept as the gain drift correction factor and eventually written to the output file.

This way of matching the gain to a reference has some advantages over peak-tracking algorithms when the spectra to be
aligned have no clear peaks. In this case peak tracking algorithms have a hard time to produce consistent results.

## Note:
The script currently will only find a proportional factor for the gain shift
(linear gain, no offset). Higher order polynomial shifts might be necessary,
in this case the code should be extended with a more complicated
multi-dimension minimisation algorithm.

