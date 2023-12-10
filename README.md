# Y(2S)+Phi Spectrum Analysis

In this framework I am presenting the script I wrote for an analysis about the searching of possible resonances in the spectrum of the $\Upsilon(2S) → \mu\mu + \Phi → KK$. 
All data was collected at CMS during the Run2 is from the _MuOnia_ dataset. 
The latter contains all the event collections which are useful to detect the production of quarkonium states (such as J/Psi, Y, etc.) through the reconstruction of dimuons. 
In particular, dimuons in the range of Y(2S) are selected, then ditracks were attached to form the Candidate. 
The CMSSW producer is responsible for this process and store all the objects involved, then CMSSW rootupler unpacks the kinematic quantities and vertex information. 

## The program

The .root files are read with two scripts: `Y2SKK.py` and `CandAnalysis.py`, which exploit PyROOT and are executable with the version 3.10 of Python.
```
python3 <script_name>.py
```
Each script opens one of the two trees which are present in the rootuple and builds a **RDataFrame**:

- The `UpsTree`, which stores the kinematics of the Y(2S) and of each of its daughter muons.
- The `CandidateTree`, which stores the kinematics of the single muon, their combination (dimuon), the single tracks, their combination (ditrack), and the total candidate. 

Other `.py` files are modules where some definitions of functions are stored:
- In `declarations.py`, some algebric functions are defined in order to compute physically relevant quantities for the analysis, i. e. `ComputePseudoRapidity()`. They are defined directly using the `ROOT.gInterpreter()`
- In `definitions.py`, new columns are defined in the RDataFrame through the method `Define()`, in order to enrich the analysis.

`Y2SPhiRun2List.txt` contains the full paths of all the rootuples: since I downloaded the whole dataset on a hard disk, I filled the file by command-line using
```
ls path/to/hard/disk/*.root > Y2SPhiRun2List.txt
```

## Structure of the code


