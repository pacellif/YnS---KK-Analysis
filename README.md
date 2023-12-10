# $\Upsilon(2S)$+$\Phi$ Spectrum Analysis

In this framework I am presenting the script I wrote for an analysis about the searching of possible resonances in the spectrum of the $\Upsilon(2S) \rightarrow \mu\mu + \Phi \rightarrow KK$. 
All data was collected at CMS during the Run2 is from the _MuOnia_ dataset. 
The latter contains all the event collections which are useful to detect the production of quarkonium states (such as $J/\Psi$, $\Upsilon$, etc.) through the reconstruction of dimuons. 
In particular, dimuons in the range of $\Upsilon(2S)$ are selected, then ditracks were attached to form the Candidate. 
The CMSSW producer is responsible for this process and store all the objects involved, then CMSSW rootupler unpacks the kinematic quantities and vertex information. 

## The program

The .root files are read with two scripts: `Y2SKK.py` and `CandAnalysis.py`, which exploit PyROOT and are executable with the version 3.10 of Python.
```
python3 <script_name>.py
```
Each script opens one of the two trees which are present in the rootuple and builds a **RDataFrame**:

- The `UpsTree`, which stores the kinematics of the $\Upsilon(2S)$ and of each of its daughter muons.
- The `CandidateTree`, which stores the kinematics of the single muon, their combination (dimuon), the single tracks, their combination (ditrack), and the total candidate. 

Other `.py` files are modules where some definitions of functions are stored:
- In `declarations.py`, some algebric functions are defined in order to compute physically relevant quantities for the analysis, i. e. `ComputePseudoRapidity()`. They are defined directly using the `ROOT.gInterpreter()`
- In `definitions.py`, new columns are defined in the RDataFrame through the method `Define()`, in order to enrich the analysis.

`Y2SPhiRun2List.txt` contains the full paths of all the rootuples: since I downloaded the whole dataset on a hard disk, I filled the file by command-line using
```
ls path/to/hard/disk/*.root > Y2SPhiRun2List.txt
```

## Structure of the code

Each of the two macro can be divided in three parts:
1. The first one includes the importation of libraries and modules, opens the .root files and builds the RDataFrame.
2. The second part contains the definition of functions which will draw and save all the plots for the analysis.
3. In the last part, it is implemented an interface that allows to choose which plots to print.

### 1. Starting

The libraries and modules imported by the program there are
```
import ROOT 
import os
from time import time	#to get the computing time
from sys import argv	#to pass arguments by command line
from definitions import UpsTreeDefinitions #or CandTreeDefinitions for the other script
import declarations
```

I imported the `time` function because I was interested in the computing time. In addition, the `argv` list was very useful to run the script by inserting two different arguments from command line:
- an integer, corresponding to the number of .root files to open, in order to execute only a sample of the whole dataset, and speed up the debugging of the code or the correction of the plots
- a string, corresponding to the name of the folder of storage of the plots

If the argument is a number, the plots are automatically stored in a folder named "test", while if it is a string, the whole dataset is passed in the RDataFrame constructor.


















