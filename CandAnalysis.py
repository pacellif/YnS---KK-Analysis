"""
Software & Computing project
by Francesco Pacelli
	a. y. 22/23


'CandAnalysis.py' is devoted to the data reading and elaboration of the CandidateTree stored in the rootuples.
It exploit the PyROOT package and its structure follows 3 sections:
	1. Setup - Import of modules, custom function definitions, input management
 	2. Core  - Definition of 'Plot functions': each plot is associated to a function that will be executed
  	3. Menu  - Implementation of a menu to select plots by command line

Execution: simply type `python3 CandAnalysis.py` from command line, plus arguments (see TREES READING section)

"""


#--------------------------------------------------------------------------------------
#	SETUP
#--------------------------------------------------------------------------------------
print("\n^^^^^^^^\tY(2S)->µµ + Phi->KK SPECTRUM ANALYSER\t^^^^^^^^\n\nImporting modules...")

import ROOT
from time import time
import os
from sys import argv
from definitions import CandTreeDefinitions
import declarations
import pandas as pd
import matplotlib.pyplot as plt

#	FILES PATH READER

with open('data/Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()
    f.close()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")


#	TREES READING - using command line arguments
#	- no command line arguments: analysis of all data and store plots in a "test" directory
# 	- argv is an integer: analysis of first n root files, in a "test" directory
#	- argv is a string: analysis of all data and store plots in directory named <argv_value>

d = "test"
try: 
	sample = allFiles[:int(argv[1])] 
except IndexError:
	sample = allFiles[:]
except ValueError:
	sample = allFiles[:]
	d = argv[1]
											
if not os.path.isdir(f"./{d}/"):			
	os.system(f"mkdir {d}")		
	
# create a root file to store and edit plots, if necessary			
p = ROOT.TFile.Open(d+"/cand_plots.root","RECREATE")		
												

#	PLOT TEMPLATE
def cprint (hist, path, name, opt="", stats=False, root = p):
	c = ROOT.TCanvas(name, name)
	
	if stats==False: hist.SetStats(0)
	hist.Draw(opt)
		
	c.SaveAs(path+name+".pdf")
	c.Write(name)
#_____________________________________________ END OF DEF
					

#	OPEN THE SAMPLE
		
print("Creating Dataset...")

dataY2SKK = CandTreeDefinitions\
(ROOT.RDataFrame("rootuple/CandidateTree",sample))


#	KK filters to compare the distribution of interest (K+K-) with 
#	a random distribution (K+K+ or K-K-)

dataKKrs = dataY2SKK.Filter("track1_charge * track2_charge < 0")
dataKKws = dataY2SKK.Filter("track1_charge * track2_charge > 0")


#------------------------------------------------------------------	
#	CORE - DEFINITIONS OF THE PLOT FUNCTIONS
#------------------------------------------------------------------

#	LEADING KAON PT PLOT

def kaonL_pt ():
	hist = dataKKrs.Histo1D(("Leading Kaon", "#phi #rightarrow K^{+}K^{-};p_{T}(K^{L}) [GeV];Counts", 200, 0., 3.), "trackL_pT")
	cprint(hist, d+"/","LeadingK", stats=True)
#_______________________________________________________ END OF DEF

# 	SOFT KAON PT PLOT

def kaonS_pt():
	hist = dataKKrs.Histo1D(("Soft Kaon", "#phi #rightarrow K^{+}K^{-};p_{T}(K^{soft}) [GeV];Counts", 200, 0., 2.), "trackS_pT")
	cprint(hist, d+"/", "SoftK", stats=True)
#_______________________________________________________ END OF DEF

#	KK INVARIANT MASS PLOT, WITH CUTS

def m_kk():

#	The purpose of this plot is to find and show which are the best KK pT cuts
#	to set boundaries to the phi mass

		# cuts for leading K and soft K pTs
	LValues = [0., 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0]
	SValues = [0., 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8]
		
		#create array of histograms with all cuts  
	hists = []
	
	for i,j in zip(LValues, SValues):
		hists.append(dataKKrs\
		.Filter(f"trackL_pT > {i} & trackS_pT > {j}")\
		.Histo1D(("KK invariant mass", "#phi #rightarrow K^{+}K^{-};m(KK) [GeV];Counts", 100, 0.99, 1.06), "ditrack_mass"))
		

	c0 = ROOT.TCanvas()
	legend = ROOT.TLegend(0.46, 0.15, 0.89, 0.50)

	hists[0].SetStats(0)	#no stats

		#overlap plots with cuts
	for i, hist in enumerate(hists):
		
		hist.Scale(1./hist.Integral())		# normalization of plots
		hist.GetYaxis().SetRangeUser(0, 0.015)	# range of y axis
		hist.SetLineColor(i+1)					
		legend.AddEntry(hist.GetPtr(), "p_{T} lead > "+str(LValues[i])+"; p_{T} soft > "+str(SValues[i]), "l")	
		hist.Draw("same")
		
	legend.Draw("")

	c0.Draw("")
	c0.SaveAs(d+"/MassKK.pdf")
	p.WriteObject(c0,"MassKK")
	
		 # display only last cut
	hists[-1].SetTitle("#phi #rightarrow K^{+}K^{-}: p_{T}(K^{L}) > 1.0, p_{T}(K^{S}) > 0.8")
	cprint(hists[-1], d+"/", "MassKK_lastcut")
#___________________________________________________________ END OF DEF	


#	FIT OF KK_PT CUT	#####################
#	take the cut pt_kL > 1.2 and pt_kS > 1.0

#	Signal -> Voigtian
#	Background -> Polynomial

def quadrature (a,b):					# to determine boundaries of the phi mass ()
	return pow( (pow(a,2) + pow(b,2)), 0.5 )	# µ ± 2 Sigma with [ Sigma = ( (Gamma/2)^2 + sigma^2 )^(0.5) ]

def m_kk_fit(ptL = 1.2, ptS = 1.0):

		# create the filtered dataset AsNumpy
	dfkk = dataKKrs.Filter(f"trackL_pT > {ptL} & trackS_pT > {ptS}")\
	.Filter("ditrack_mass > 1.0 & ditrack_mass < 1.04")\
	.Filter("track1_pvAssocQ + track2_pvAssocQ > 11")\
	.AsNumpy(columns=["ditrack_mass"])
		
		# variable
	kkmass = ROOT.RooRealVar("ditrack_mass", "m(KK) [GeV]", 1.00, 1.04)

		# create the RooDataSet
	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
	kkroohist = kkroodata.binnedClone()
	
		#	Number of entries
	entries = kkroohist.sumEntries()
	
		#	Create frame
	phiframe = kkmass.frame(Title="Dikaon Candidate Mass")

	
		#	Signal parameters
	mean = ROOT.RooRealVar("#mu_{#phi}", "mean of gaussian", 1.019, 1.018, 1.020) #mean value, min value, max value
	sigma = ROOT.RooRealVar("#sigma_{#phi}", "resolution", 0.00125, 0.001, 0.002) 
	width = ROOT.RooRealVar("#Gamma_{#phi}", "width", 0.00439, 0.003, 0.006)
	
	sigma.setConstant(1)		# fixed res from MC????????

		#	Chebyshev coefficients
	f0 = ROOT.RooRealVar("f0", "f0", 0.2, 0., 2.)
	f1 = ROOT.RooRealVar("f1", "f1", -0.05, -2., 0.)

		#	Number of events
	Nbkg = ROOT.RooRealVar("N_{bkg}", "N bkg events", 50, 0., entries)
	Nsig = ROOT.RooRealVar("N_{sig}", "N sig events", 50, 0., entries)

		#	Model Functions
	sig = ROOT.RooVoigtian("signal", "signal", kkmass, mean, width, sigma)
	bkg = ROOT.RooChebychev("bkg", "Background", kkmass, [f0, f1]) 

		#	Total model
	model = ROOT.RooAddPdf("model", "voigt+cheb", [bkg, sig], [Nbkg, Nsig])
		
	model.fitTo(kkroohist)
		
		# print
	c0 = ROOT.TCanvas()

	kkroohist.plotOn(phiframe)
	model.plotOn(phiframe) # By default only fitted range is shown
	model.plotOn(phiframe, Components={sig}, LineStyle=":", LineColor="r")
	model.plotOn(phiframe, Components={bkg}, LineStyle=":", LineColor="g")
	model.paramOn(phiframe, ROOT.RooFit.Parameters([mean, width, Nbkg, Nsig, f0, f1]), ROOT.RooFit.Layout(0.65, 0.9, 0.9))

		# create boundaries and highlight them with lines
	xmin = mean.getVal() - 2*quadrature(width.getVal()/2, sigma.getVal())
	xmax = mean.getVal() + 2*quadrature(width.getVal()/2, sigma.getVal())

#		to highlight the boundaries of the signal uncomment and set an adequate height

#	line0 = ROOT.TLine(xmin, 0., xmin, 10000) 	
#	line1 = ROOT.TLine(xmax, 0., xmax, 10000)
#	line0.SetLineStyle(2)
#	line0.SetLineColor(7)
#	line0.SetLineWidth(4)
#	line1.SetLineStyle(2)
#	line1.SetLineColor(7)
#	line1.SetLineWidth(4)

	phiframe.Draw()
#	line0.Draw("same")
#	line1.Draw("same")
	
		# print and save

	p.WriteObject(phiframe,"phi_mass_fit")
	c0.SaveAs(d+"/PhiMassPlot.pdf")
#___________________________________________________________ END OF DEF

	
# PHI CANDIDATE PLOT	
# ws = wrong sign (take K+K+ and K-K-)
# compare the selection distribution with one random distribution


	#Γ φ = 0.00446 ± 0.00018
	#μφ  = 1.019445 ± 0.000036
	
	#	Y mass cuts: µ ± 2 sigma 
	#	φ mass cuts: µ ± 2 Sigma

# filters
ymumu_filter= "dimuon_mass > 9.881 & dimuon_mass < 10.147 & dimuon_pT > 18 & "
phiKKSelection = '''candidate_vProb > 0.1 &
ditrack_mass > 1.0150 &
ditrack_mass < 1.0239 &
trackL_pT > 1.0 & 
trackS_pT > 0.8'''
quality_filter = " & track1_pvAssocQ + track2_pvAssocQ > 11"
candProb = " & candidate_vProb > 0.05"

# ditrack mass: interval of phi rest mass
# vProb: vertex probability

binning = 250
edge = [10.8, 13.3]

	# cuts written on file
tagli = ymumu_filter + phiKKSelection + quality_filter + candProb
tagli = tagli.replace(" & ", "\n")
os.system(f"echo \"{tagli}\nbinning: {binning}\" > {d}/tagli.txt")



#	CANDIDATE INVARIANT MASS DISTRIBUTION PLOT

# now find the spectrum of Y + phi, comparing a spectrum where a real phi is considered (opposite sign kaons)
# with respect to a random distribution given by two kaons with same sign

def mumukk(zoom = False):

		# opposite sign kaons
	hist0 = dataKKrs\
	.Filter(ymumu_filter + phiKKSelection)\
	.Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", 
		  binning, edge[0], edge[1]), "candidate_vMass")

		# same sign kaons
	hist1 = dataKKws\
	.Filter(ymumu_filter + phiKKSelection)\
	.Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", 
		  binning, edge[0], edge[1]), "candidate_vMass")

	c0 = ROOT.TCanvas()

	hist0.SetStats(0)

	hist0.SetLineColor(1)
	hist1.SetLineColor(2)
	
	legend = ROOT.TLegend(0.7, 0.1, 0.89, 0.3) #(xmin,ymin,xmax,ymax)

	legend.AddEntry(hist0.GetPtr(), "RS kaons", "l")
	legend.AddEntry(hist1.GetPtr(), "WS kaons (norm)", "l")

	hist1.Scale(hist0.Integral()/hist1.Integral())
	        
	hist0.Draw("")
	hist1.Draw("same")
	legend.Draw("")
	c0.Draw("")
	c0.SaveAs(d+"/PhiCandidate.pdf")
	
	p.WriteObject(c0,"phi_candidate")

# eventually, it is possible to zoom in a particular region, in case
# abundancies of countings are seen in the spectrum
	
	if zoom:

		hist0.GetXaxis().SetRangeUser(11., 11.6)
		hist1.GetXaxis().SetRangeUser(11., 11.6)
		
		hist2 = hist0
		hist3 = hist1
		
		
		hist3.Scale(hist2.Integral()/hist3.Integral())
		
		hist2.SetTitle("Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-}) (Zoom)")
		
		hist2.Draw("")
		hist3.Draw("same")
		
		legend.Draw("")
		c0.Draw("")
		c0.SaveAs(d+"/PhiCandidateZoom.pdf")
		
		p.WriteObject(c0,"phi_candidate(zoom)")
#___________________________________________________________ END OF DEF

#		NUMBER OF CANDIDATES PER EVENT
#	from the combination of ditracks and dimuons, it is possible to have several many candidates from the same dimuon
#	by looking at the number of candidates it is possible to act on the multiplicity of candidates per event to 
#	clean the spectrum from further background countings. 
#	Moreover, it is computed the percentage of multiple candidates per event.

def Ncand():
	
	#	the selected candidates
	candSelectionRS = dataKKrs.Filter(ymumu_filter+phiKKSelection+ candProb + quality_filter)


	dfcands = candSelectionRS.AsNumpy(columns=["event", "run", "candidate_pT",
									 "candidate_vProb", "candidate_vMass"])
	
	candDF = pd.DataFrame(dfcands)
	candxEvent = candDF.groupby("event")["candidate_vProb"].count()
	#hist = candxEvent.plot(kind="hist", logy=True)
	
	plt.hist(candxEvent, bins=range(1,9), align="left", edgecolor='black')
	plt.title("Candidate multiplicity")
	plt.xlabel("Multiplicity")
	plt.ylabel("Counts")
	plt.yscale('log')

	
	plt.show()
	
	multiplicity_ratio = (candxEvent > 1).sum()/(candxEvent > 0).sum()
	print(f"\nmultiplicity ratio is {round(multiplicity_ratio*100, 2)}%")
	
#___________________________________________________________ END OF DEF
########################################################################

#----------------------------------------------------------------------
#	MENU
#----------------------------------------------------------------------

# The menu allows to select which plot(s) to print:
# from command line, insert the key(s) of the plot(s), 
# separated by space for multiple choice.

lang = input("\nSelect plots (Separate by spacing):\n" + 
	     "1. Leading Kaon pt\n" + 
	     "2. Soft Kaon pt\n" + 
	     "3. KK invariant mass (with K_pt cuts)\n" + 
	     "4. Fit KK invariant mass\n" + 
	     "5. Phi Candidate plot\n" +
	     "6. Phi candidate plot with Zoom\n" + 
	     "7. Candidate multiplicity\n" +
	     "ENTER: 1-5 Plots\n" + 
	     "Press \"q\" to EXIT.\n").split()

if "q" not in lang:
	print("Processing...") 

	# Start timer
start = time()

if not lang:	# print all plots

	kaonL_pt()
	kaonS_pt()
	m_kk()					
	m_kk_fit()			#save to .root
	mumukk()			#save to .root
else:
	for i in lang:
	
		#	avoid misdigit
		while i not in "1 2 3 4 5 6 7 q":
			i = input(f"\"{i}\" is not valid. Please insert a valid key: ")

		#	execution
		if i == "1":	kaonL_pt()
		if i == "2":	kaonS_pt()
		if i == "3":	m_kk()
		if i == "4":	m_kk_fit()		#save to .root
		if i == "5":	mumukk()		#save to .root
		if i == "6":	mumukk(zoom=True)
		if i == "7":	Ncand()
		if i == "q":
			print ("Bye Bye")
			exit()
	#_____________________________________ END OF LOOP

"""
################ for python version >= 3.10
	match i:
			case "1":
				kaonL_pt()
			case "2":
				kaonS_pt()
			case "3":
				m_kk()
			case "4":
				m_kk_fit()		#save to .root
			case "5":
				mumukk()		#save to .root
			case "6":		
				mumukk(zoom=True)
    			case "7": 
       				Ncand()
			case "q":
				print ("Bye Bye")
				exit()
"""

	# End timer
end = time()

	# Calculate elapsed time
elapsed = end - start
print("\nComputing time: ", elapsed, "\n") 


#	open the file with the TBrowser
if "7" not in lang:
	b = ROOT.TBrowser()
	p.Browse(b)


#	exit command
if input("Press \"q\" to exit from the framework \n") == "q":
	p.Close()
	exit()



