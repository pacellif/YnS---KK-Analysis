# Y(2S)+Phi Spectrum Analysis
In this framework I am presenting the script I wrote for an analysis about the searching of possible resonances 
in the spectrum of the Y(2S) → µµ + Phi → KK. 
All data is from the "MuOnia" dataset. The data of interest has been previously processed by a CMSSW producer, then stored in a rootuple through a CMSSW rootupler.

The rootuples contain two trees:
- The `UpsTree`, which stores the kinematics of the Y(2S) and of each of its daughter muons.
- The `CandidateTree`, which stores the kinematics of the single muon, their combination (dimuon), the single tracks, their combination (ditrack), and the total candidate. 

##The program
The scripts for the analysis are `Y2SKK.py` and `CandAnalysis.py`, which open respectively the UpsTree and the CandTree.

Other `.py` files are modules where some definitions of functions are stored.
`Y2SPhiRun2List.txt` contains the full paths of all the rootuples: since I downloaded the whole dataset on a hard disk, I filled the file by command-line using

```
ls path/to/hard/disk/*.root > Y2SPhiRun2List.txt
```
