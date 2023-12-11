print("\n^^^^^^^^\tY(2S)->µµ + Phi->KK SPECTRUM ANALYSER\t^^^^^^^^\n\nImporting modules...")

import ROOT 
from time import time
import os
from sys import argv
from definitions import CandTreeDefinitions
import declarations

ROOT.gROOT.SetBatch(True)


#	PLOT TEMPLATE
def cprint (hist, name, opt="", stats=False):
	title= "Y2S "+name
	c = ROOT.TCanvas(title, title)
	
	if stats==False: hist.SetStats(0)
	hist.Draw(opt)
		
	c.SaveAs(name+".pdf")
	os.system("xdg-open "+name+".pdf")



with open('Y2SPhiRun2List.txt') as f:
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
p = ROOT.TFile.Open(d+"/cand.root","RECREATE")		
												

# creating datasets for the analysis		
print("Creating Dataset...")

dataY2SKK = CandTreeDefinitions\
(ROOT.RDataFrame("rootuple/CandidateTree",sample))#\
#.Define("mumupipi_mcmass", "ComputeMuMuPiPiInvMassMC(muonp_p4_fX, muonp_p4_fY, muonp_p4_fZ, muonn_p4_fX, muonn_p4_fY, muonn_p4_fZ, track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")	#Mass constraint


#KK filters
dataKKrs = dataY2SKK.Filter("track1_charge * track2_charge < 0")
dataKKws = dataY2SKK.Filter("track1_charge * track2_charge > 0")

#############################################################
#					CANDIDATES KK DECAY						#
#############################################################


#	DEFINING DIFFERENT PLOTS
"""
#	MASS PLOT OF CANDIDATE	##############################
	# Start timer
start = time()

YKKMass = ROOT.RooRealVar("YKKMass", "YKKMass", 11, 13.5)
YKKMass.setBins(1000)
rooDataSetYKK = dataY2SKK.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(YKKMass))), ["candidate_vMass"])
Cand_frame = YKKMass.frame(Title="Y(2S)KK Mass")
rooDataSetYKK.plotOn(Cand_frame,LineColor="b",MarkerSize=0.3)

c_Cand = ROOT.TCanvas("MassPlotY2SKK", "MassPlotY2SKK", 800, 400)
Cand_frame.Draw()
c_Cand.SaveAs("phikk_plots/MassPlotY2SKK.pdf")

	#candidate mass
hist = dataY2SKK.Histo1D(("Candidate Mass", "Candidate Mass; m(K^{L})+m(K^{S})+ m(Y(2S))[GeV/c^{2}];Counts", 500, 11., 13.6), "candidate_vMass")
cprint(hist, "phikk_plots/candidateMass")
	#dimuon mass
hist = dataY2SKK.Histo1D(("Dimuon Mass", "Candidate Mass; m[GeV/c^{2}];Counts", 500, 9.4, 10.5), "dimuon_vMass")
cprint(hist, "phikk_plots/dimuonMass", stats=True)
"""

#	LEADING KAON PT PLOT

def kaonL_pt ():
	hist = dataKKrs.Histo1D(("Leading Kaon", "#phi #rightarrow K^{+}K^{-};p_{T}(K^{L}) [GeV];Counts", 200, 0., 3.), "trackL_pT")
	cprint(hist, d+"/LeadingK", stats=True)

# 	SOFT KAON PT PLOT
def kaonS_pt():
	hist = dataKKrs.Histo1D(("Soft Kaon", "#phi #rightarrow K^{+}K^{-};p_{T}(K^{soft}) [GeV];Counts", 200, 0., 2.), "trackS_pT")
	cprint(hist, d+"/SoftK", stats=True)

#	KK INVARIANT MASS PLOT, WITH CUTS	#################
def m_kk():
		# Start timer
	start = time()

		#borders definitions
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
		
		hist.Scale(1./hist.Integral())	#normalization of plots
		hist.GetYaxis().SetRangeUser(0, 0.015)	#range of y axis
		hist.SetLineColor(i+1)					
		legend.AddEntry(hist.GetPtr(), "p_{T} lead > "+str(LValues[i])+"; p_{T} soft > "+str(SValues[i]), "l")	
		hist.Draw("same")	#non si distinguono molto bene i dati
		
	legend.Draw("")


	c0.Draw("")
	c0.SaveAs(d+"/MassKK.pdf")
	os.system(f"xdg-open {d}/MassKK.pdf")
		 #only last cut
	hists[-1].SetTitle("#phi #rightarrow K^{+}K^{-}: p_{T}(K^{L}) > 1.0, p_{T}(K^{S}) > 0.8")
	cprint(hists[-1], d+"/MassKK_lastcut")

#	FIT OF KK_PT CUT	#####################
#	take the cut pt_kL > 1.2 and pt_kS > 1.0

def m_kk_fit(ptL = 1.2, ptS = 1.0):
	dfkk = dataKKrs.Filter(f"trackL_pT > {ptL} & trackS_pT > {ptS}")\
	.Filter("ditrack_mass > 1.0 & ditrack_mass < 1.04")\
	.AsNumpy(columns=["ditrack_mass"])
		
		#variable
	kkmass = ROOT.RooRealVar("ditrack_mass", "m(KK) [GeV]", 1.00, 1.04)

	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
	kkroohist = kkroodata.binnedClone()
	phiframe = kkmass.frame(Title="Dikaon Candidate Mass")

	
		#model
	mean = ROOT.RooRealVar("#mu_{#phi}", "mean of gaussian", 1.019, 1.018, 1.020) #mean value, min value, max value
	sigma = ROOT.RooRealVar("#sigma_{#phi}", "resolution", 0.00125, 0.001, 0.002) 
	sigma.setConstant(1)		# fixed res from MC????????
	width = ROOT.RooRealVar("#Gamma_{#phi}", "width", 0.00439, 0.003, 0.006)

	f0 = ROOT.RooRealVar("f0", "f0", +0.2, -1, +1)
	f1 = ROOT.RooRealVar("f1", "f1", -0.05, -1, +1)

	bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.92, 0.8, 1.0)

	sig = ROOT.RooVoigtian("signal", "signal", kkmass, mean, width, sigma)
	bkg = ROOT.RooChebychev("bkg", "Background", kkmass, [f0, f1]) ## to increase the degree, just increase the coefficients

	model = ROOT.RooAddPdf("model", "voigt+cheb", [bkg, sig], [bkgfrac])
		
	model.fitTo(kkroohist)
		
		#print
	c0 = ROOT.TCanvas("canvas0", "canvas0", 1200, 600)

	kkroohist.plotOn(phiframe)
	model.plotOn(phiframe) # By default only fitted range is shown
	model.plotOn(phiframe, Components={sig}, LineStyle=":", LineColor="r")
	model.plotOn(phiframe, Components={bkg}, LineStyle=":", LineColor="g")
	model.paramOn(phiframe, ROOT.RooFit.Parameters([mean, width, bkgfrac]), ROOT.RooFit.Layout(0.65, 0.9, 0.9))

	xmin = mean.getVal() - width.getVal()
	xmax = mean.getVal() + width.getVal()

	line0 = ROOT.TLine(xmin, 0., xmin, 10000)
	line1 = ROOT.TLine(xmax, 0., xmax, 10000)
	line0.SetLineStyle(2)
	line0.SetLineColor(7)
	line0.SetLineWidth(4)
	line1.SetLineStyle(2)
	line1.SetLineColor(7)
	line1.SetLineWidth(4)

	phiframe.Draw()
	line0.Draw("same")
	line1.Draw("same")
	c0.Draw()
	p.WriteObject(c0,"phi_mass_fit")
	c0.SaveAs(d+"/PhiMassPlot.pdf")
	os.system(f"xdg-open {d}/PhiMassPlot.pdf")
		
		# End timer
	end = time()

		# Calculate elapsed time
	elapsed = end - start
	print("\nTime for Phi mass plot: ", elapsed, "\n") 


# PHI CANDIDATE PLOT	
#ws=wrong sign (take K+K+ and K-K-)
#compare the selection distribution with one random distribution


	#Γ φ = 0.00446 ± 0.00018
	#μφ  = 1.019445 ± 0.000036

ymumu_filter="dimuon_mass > 9.8 & dimuon_mass < 10.15 & dimuon_pT > 18 & "
phiKKSelection = '''candidate_vProb > 0.1 &
ditrack_mass > 1.0150 &
ditrack_mass < 1.0239 &
trackL_pT > 1.0 & 
trackS_pT > 0.8'''


# ditrack mass: interval of phi rest mass
# vProb: vertex probability

binning = 250
edge = [10.8, 13.3]

# new dir with cuts written on file

tagli = ymumu_filter+phiKKSelection
tagli = tagli.replace(" & ", "\n")
os.system(f"echo \"{tagli}\nbinning: {binning}\" > {d}/tagli.txt")

#histograms
def mumukk(zoom = False):
	hist0 = dataKKrs.Filter(ymumu_filter + phiKKSelection).Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", binning, edge[0], edge[1]), "candidate_vMass")
	hist1 = dataKKws.Filter(ymumu_filter + phiKKSelection).Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", binning, edge[0], edge[1]), "candidate_vMass")

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
	os.system(f"xdg-open {d}/PhiCandidate.pdf")
	
	p.WriteObject(c0,"phi_candidate")
	
	if zoom:
		hist0.GetXaxis().SetRangeUser(11.3, 11.6)
		hist1.GetXaxis().SetRangeUser(11.3, 11.6)
		
		hist2 = hist0
		hist3 = hist1
		
		hist3.Scale(hist2.Integral()/hist3.Integral())
		
		hist2.SetTitle("Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-}) (Zoom)")
		hist2.Draw("")
		hist3.Draw("same")
		legend.Draw("")
		c0.Draw("")
		c0.SaveAs(d+"/PhiCandidateZoom.pdf")
		os.system(f"xdg-open {d}/PhiCandidateZoom.pdf")

#	MENU

compute = {	"1" : kaonL_pt,
			"2" : kaonS_pt,
			"3" : m_kk,
			"4" : m_kk_fit,		
			"5" : mumukk,		
			"6" : mumukk,
			"q" : exit
		  }


lang = input("\nSelect plots (Separate by spacing):\n1. Leading Kaon pt\n2. Soft Kaon pt\n3. KK invariant mass (with K_pt cuts)\n4. Fit KK invariant mass\n5. Phi Candidate plot\n6. Phi candidate plot with Zoom\nENTER: 1-5 Plots\nPress \"q\" to EXIT.\n").split()

if "q" not in lang:
	print("Processing...") 

	# Start timer
start = time()

if not lang:	#print all plots		
	for func in compute.values() : func()

else:			#print only selected-by-key plots
	for i in lang: 
	
		while i not in compute.keys(): 	#in case of multiple misdigit
			i = input(f"\"{i}\" is not valid. Please insert a valid key:\n")
		
		if i != "6" : compute[i]()
		else : compute[i](zoom = True)

"""
if not lang:	#print all plots

		kaonL_pt()
		kaonS_pt()
		m_kk()					
		m_kk_fit()			#save to .root
		mumukk()			#save to .root
else:
	for i in lang:
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
			case "q":
				print ("Bye Bye")
				exit()
			case _:
				print("Not valid")
				exit()

"""
	# End timer
end = time()

	# Calculate elapsed time
elapsed = end - start
print("\nComputing time: ", elapsed, "\n") 
#os.system(f"echo {elapsed} > time.txt")
p.Close()
