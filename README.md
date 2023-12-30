# Analysis of $\Upsilon(2S)$ + $\phi$ spectrum

In this framework I am presenting the script I wrote for an analysis about the searching of possible resonances in the spectrum of the $\Upsilon(2S) \rightarrow \mu\mu + \phi \rightarrow KK$. 
All data was collected at CMS during the Run2, elaborated and stored in the _MuOnia_ dataset. 
The latter contains all the event collections which are useful to detect the production of quarkonium states ($J/\Psi$, $\Upsilon$, etc.) through the reconstruction of dimuons. 
In particular, dimuons in the range of $\Upsilon(2S)$ are selected, as it is done for ditracks in range of $\phi$. Then, all these objects were attached together to form the *candidate*. 
The CMSSW producer is responsible for this process and save all the objects involved, then CMSSW rootupler unpacks the kinematic quantities and vertex information in a .root file. 

This tool reads the rootuples and provides a large variety of plots that simply show distributions of specific quantities (such as p~T~, pseudorapidity and more) or allow to work on the mass spectra of the $\Upsilon$, the $\phi$ and of the whole candidate.

The spectrum in *Figure 1* represents the first three excited states of $\Upsilon$, for which I will focus only on the second peak. 
The $\phi$ spectrum is extracted from the ditrack, which are assumed to be pairs of $K^+K^-$, since CMS can not do PID for this kind of objects. 
In *Figure 2* there is an example of plot.

![YMass](https://upload.wikimedia.org/wikipedia/commons/e/e0/Upsilon_mesons_CMS.svg)
|---|

Figure 1: $\Upsilon$ spectrum from the dimuon invariant mass
|--|


![PhiMass](https://www.science20.com/files/images/phicms.png)
|---|

Figure 2: $\phi$ spectrum from the ditrack invariant mass. In the framework, ditracks are assumed to be Kaons.
|-----|


### Requirements

In order to be able to run the scripts it is necessary to have installed ROOT, with the PyROOT module and Python v3.

It is possible to download the repository from terminal with the command
```
cd $HOME
git clone https://github.com/pacellif/YnS-KK-Analysis.git
```
In the repository there are the README.md, the .py scripts and modules, a directory to store the data, which already contains a .txt file where to write the entire paths of the rootuples.


### Download the data

Since the .root files overcome the maximum size for the upload on GitHub, I stored a sample of 10 files on CERNBox folder (about 190 MB each). They must be downloaded by accessing with LXplus username and password and using the command

```
cd YnS-KK-Analysis
scp <username>@lxplus.cern.ch:/eos/user/f/fpacelli/data/*.root /data
```

Before running the scripts, fill the file `Y2SPhiRun2List.txt` by typing from command-line
```
ls $PWD/data/*.root > data/Y2SPhiRun2List.txt
```

Now you are ready to proceed.

## The framework

The .root files are read with two scripts: `Y2SKK.py` and `CandAnalysis.py`, which exploit PyROOT and are executable with Python 3.

```
python3 Y2SKK.py
```
or 
```
python3 CandAnalysis.py
```

Each script opens one of the two trees which are present in the rootuple and builds a **RDataFrame**:

- The `UpsTree`, which stores the kinematics of the $\Upsilon(2S)$ and of each of its daughter muons.
- The `CandidateTree`, which stores the kinematics of the single muon, their combination (dimuon), the single tracks, their combination (ditrack), and the total candidate. 

Other `.py` files are modules where some definitions of functions are stored:
- In `declarations.py`, some algebric functions are defined in order to compute physically relevant quantities for the analysis, i. e. `ComputePseudoRapidity()`. They are defined directly using the `ROOT.gInterpreter()`
- In `definitions.py`, new columns are defined in the RDataFrame through the method `Define()`, in order to enrich the analysis.


## Structure of the code

Each of the two macro can be divided in three parts:
1. The first one includes the importation of libraries and modules, opens the .root files and builds the RDataFrame.
2. The second part contains the definition of functions which will draw and save all the plots for the analysis (_plot functions_).
3. In the last part, it is implemented an interface that allows to choose which plots to print.

### 1. Starting

The libraries and modules imported by the program there are
```py
import ROOT 
import os
from time import time						#to get the computing time
from sys import argv						#to pass arguments by command line
from definitions import UpsTreeDefinitions	#or CandTreeDefinitions for the other script
import declarations
```
Right after, the `Y2SPhiRun2List.txt` is read, and all its lines are trimmed and stored in a list.
In order to execute only a sample of the whole dataset and speed up the debugging of the code or the correction of the plots, it is possible to pass an integer _N_ as command line argument, selecting the first _N_ .root files to open. At the end of the program, all the saved plots are stored in a directory named "test".

Alternatively, a string can be passed as command line argument, as the name of a different storing folder for the plots. This option is meant to be chosen when the code is correclty working and ready to analyse the whole dataset.

At this point, the sample of rootuples is passed as argument of RDataFrame constructor. Together with the dataframe, a TFile is opened to save the most elaborated plots, so that it is possible to access and edit them through the TBrowser or by the line command `rootbrowse <fileroot>.root`.

In the end, I imported the `time` function because I was interested in the computing time.

### 2. The core

In the middle part of the code, the plot functions are defined. 

I implemented a standard function `cprint` to simply draw and save a histogram which plots one of the columns of the dataframe.

In general, histograms were made using the RooFit package: 
`Histo1D` or `Histo2D` of the RDataFrame class are used as standard methods to print simple histograms, while to perform the fit of a histogram, RDataFrame is pythonized into a RooDataSet through the method `from_numpy`.

[//] # insert examples plots 


### 3. The interface

The interface is designed to allow the user to select which plots to work on. A menu is printed on the terminal at the beginning of the execution, displaying the keys to digit to select the plots. It is possible to choose multiple plots by separating the keys with space.

In case there is a non-valid key inserted, it is possible to correct it thanks to a while loop.

At the end of the last part the computing time is displayed in the terminal.

## Summary and conclusions

The purpose of this framework is to optimize work from local host, from the storage of the plots that are produced to the debugging and test of the code in order to avoid waste of time. 
In addition, thanks to a repetitive structure it should be easy to extend or improve the analysis.

Actually, the initial prospect was to work from remote on LXplus without downloading many GB of data, but I had some issues in the production of some plots. This was enough to test the difference in the remote computing time, which was not too much longer than the local.

For this reason, it could be useful to try solving that issue in order to work entirely from remote. 







