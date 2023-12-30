print("\n^^^^^^^^\tY(2S)->µµ ANALYSER\t^^^^^^^^\n\nImporting modules...")


import ROOT 
import os
from time import time
from sys import argv
from definitions import UpsTreeDefinitions
import declarations

ROOT.gROOT.SetBatch(True)
#ROOT.ROOT.EnableImplicitMT(4)


#	FILE MANAGEMENT:	---------------------------------------------| 
#	- open .root dataset: it is possible to select a subset
#	- make custom directories in which store plots and .root files

with open('Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()
    f.close()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")


#	TREES READING - USING COMMAND LINE ARGUMENTS
#		it is possible to pass a number of rootuples to open 
#		or the name of the directory where to store the plots

d = "test"
try: 
	sample = allFiles[:int(argv[1])] 
	
except IndexError:			# in case argv[1] is empty: open all files
 	sample = allFiles[:]

except ValueError:			# in case argv[1] is a string: open all files
	sample = allFiles[:]	# and save them in a specific folder
	d = argv[1]
											
if not os.path.isdir(f"./{d}/"):	# check if the directory exists,	
	os.system(f"mkdir {d}")			# otherwise create one



#	plot template	-------------------------------------------------|
def cprint (hist, name, opt="", stats=False):
	title= "Y2S "+name
	c = ROOT.TCanvas(title, title)
	
	if stats==False: hist.SetStats(0)
	hist.Draw(opt)
		
	c.SaveAs(d+"/"+name+".pdf")
	os.system(f"xdg-open {d}/"+name+".pdf")
#------------------------------------------------------------------
def binCount (var, rangeName):
	varBinning = var.getBinning()
	a = var.getRange(rangeName)[0]
	b = var.getRange(rangeName)[1]
	nbins = varBinning.binNumber(b) - varBinning.binNumber(a)
	return nbins
#------------------------------------------------------------------	
#	OPENING THE SAMPLE
#------------------------------------------------------------------

print("Creating Dataset...")
dataY2S = UpsTreeDefinitions(ROOT.RDataFrame("rootuple/UpsTree",sample))


fileroot = ROOT.TFile.Open("ups_mass.root","RECREATE")

#------------------------------------------------------------------	
#	DEFINITIONS OF THE PLOT FUNCTIONS
#------------------------------------------------------------------

#	MUON TRANSVERSE MOMENTUM

def mu_pt():
	hist1 = dataY2S.Histo1D(("mu1_pt Distribution", "mu1_pt Distribution;p_{T}(#mu^{+});Counts", 500, 0., 40.), "mu1_pt")

	hist2 = dataY2S.Histo1D(("mu2_pt Distribution", "mu2_pt Distribution;p_{T}(#mu^{-});Counts", 500, 0., 40.), "mu2_pt")

#		#other two variables - same results
#	hist3 = dataY2S.Histo1D(("muonP_pT Distribution", "muonP_pT Distribution;p_{T}(#mu^{+});Counts", 500, 0., 40.), "muonP_pT")

#	hist4 = dataY2S.Histo1D(("muonN_pT Distribution", "muonN_pT Distribution;p_{T}(#mu^{-});Counts", 500, 0., 40.), "muonN_pT")

	c = ROOT.TCanvas("Muon pT distribution", "Muon pT distribution")

	c.Divide(1,2)

	c.cd(1)
	hist1.Draw()


	c.cd(2)
	hist2.Draw()


	c.SaveAs(d+"/muonPTdistribution.pdf")
	fileroot.WriteObject(c,"MuonPT")
	os.system(f"xdg-open {d}/muonPTdistribution.pdf")
#___________________________________________________________ END OF DEF

#	MASS PLOT

def m_Y2S():

	hist = dataY2S\
	.Filter("ups_vMass > 9.6 & ups_vMass < 10.3")\	#apply some filters
	.Filter("ups_pT > 15")\
	.Histo1D(("dimuon invariant mass", "Y(2S) #rightarrow #mu^{+}#mu^{-};m(#mu^{+}#mu^{-}) [GeV];Counts", 500, 9.6, 10.3), "ups_vMass")
	
	cprint(hist, "YinvMass")
#___________________________________________________________ END OF DEF

#	FIT	OF Y2S MASS PLOT
	
def fit_Y2S():	
	
		#	create the frame
	
	#	the variable
	upsmass = ROOT.RooRealVar("ups_vMass", "m(#mu#mu) [GeV]", 9.6, 10.3)
	upsmass.setBins(500)
	
	#alternatively
	#mumuroohist = dataY2S.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(YMass))), ["ups_vMass"])
	
	#	the histogram from the filtered dataset
	massDF = dataY2S\
	.Filter("ups_vMass > 9.6 & ups_vMass < 10.3")\
	.Filter("ups_pT > 15")\
	.AsNumpy(columns=["ups_vMass"])
	
	mumuroodata = ROOT.RooDataSet\
	.from_numpy({"ups_vMass": massDF["ups_vMass"]},ROOT.RooArgSet(upsmass))
	mumuroohist = mumuroodata.binnedClone()
	
	xframe = upsmass.frame(Title="Y(2S) Mass")
	
# mass Y2S PDG = 10.02326
# mass Y3S PDG = 10.35520

		#signal mean
	mean2s = ROOT.RooRealVar("#mu_{Y(2S)}", "mean of gaussians", 
							  10., 9.8, 10.2)

		#sigmas
	sigma2s = ROOT.RooRealVar("#sigma_{Y(2S)}", "width", 0.063, 0.01, 5.) 
		
		#Crystal Ball parameters
	alpha = ROOT.RooRealVar ("alpha","alpha", 1.62, 0., 5.)
	n= ROOT.RooRealVar ("n","n", 0.1, -1., 1.);

		#chebychev coefficients
	f0 = ROOT.RooRealVar("f0", "f0", 5., 0., 10.)
	f1 = ROOT.RooRealVar("f1", "f1", 0., -20., 1.)
	f2 = ROOT.RooRealVar("f2", "f2", 2., 0., 8.)
		
		#fractions
	bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.5, 0.001, 1.)

	# 	MODELS FOR MASS PLOT

		#signals
	sig2s1 = ROOT.RooGaussian("signal2s_1", "signal2s_1", upsmass, mean2s, sigma2s)
	cb = ROOT.RooCBShape("Double CB", "#upsilon(2s) Pdf", upsmass, mean2s, sigma2s, alpha, n);
		
		#backgrounds
	bkg2 = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1))
	bkg3 = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1,f2))

		#MODELS
	model1 = ROOT.RooAddPdf("model1", "Cheb2+Gaus", [bkg2,sig2s1], [bkgfrac]) 
	model2 = ROOT.RooAddPdf("model2", "Cheb2+CB", [bkg2,cb], [bkgfrac]) 
	model3 = ROOT.RooAddPdf("model3", "Cheb3+Gaus", [bkg3,sig2s1], [bkgfrac]) 
	model4 = ROOT.RooAddPdf("model4", "Cheb3+CB", [bkg3,cb], [bkgfrac]) 

	allModels = [model1, model2, model3, model4]

		#choose model to fit
	model = allModels[2]

	########	FIT RESULT AND CHI SQUARED COMPUTATION  ·················
	upsmass.setRange("range", 9.8 ,10.25)	 #set range before fitting
		
	fitResult = model.fitTo(mumuroohist, Range="range", Save=True)

#		#chi squared
#	chiSquared = int(-2 * fitResult.minNll())
#	ndof = binCount(upsmass,"range") - model.getParameters(mumuroodata).getSize()

#	reducedChiSquared = round(chiSquared/ndof)

	########	PLOTTING	···········································
	mumuroohist.plotOn(xframe,LineColor="b",MarkerSize=0.3)
	model.plotOn(xframe,LineColor="r")

	component = model.pdfList()


	model.plotOn(xframe, Components={component[0]}, LineStyle=":", LineColor="g")
	model.plotOn(xframe, Components={component[1]}, LineStyle=":", LineColor="b")
	model.plotOn(xframe,LineColor="r")	#total
	
	model.paramOn(xframe, ROOT.RooFit.Layout(0.1, 0.9, 0.9)) #print all parameters

#		#print chisquare
#	text_box = ROOT.TPaveText(0.65, 0.75, 0.9, 0.9, "NDC")
#	text_box.AddText( str(model.getTitle()) )	#type of fit
#	text_box.AddText(f"#chi^{2}/ndof = {reducedChiSquared}")
#	text_box.SetFillColor(0)
#	text_box.SetBorderSize(1)
	
		#print fit and pullplot
	c = ROOT.TCanvas("MassPlotY2S", "MassPlotY2S")
	c.Divide(1,2)
	c.cd(1)
	xframe.Draw()
	text_box.Draw()
	c.cd(2)
	xframe.pullHist().Draw()	#pull histogram

	fileroot.WriteObject(c,"UpsInvMass")
	c.SaveAs(d+"/MassPlotY2S.pdf")
	os.system(f"xdg-open {d}/MassPlotY2S.pdf")
#______________________________________________________________ END OF DEF

#######		TRANSVERSAL MOMENTUM PLOTS	#############


def Y_pt():

		#	plot the cumulative data from 2016 to 2018
		
	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")	#	2016 
	
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")	#	2017

	hist18 = dataY2S.Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum (stacked);p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50.), "ups_pT")	#	2018 

	c_pt = ROOT.TCanvas("Y(2S) pT", "Y(2S) pT")


	hist17.SetFillColor(5)		#	change color
	hist16.SetFillColor(3)


	hist18.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	c_pt.SaveAs(d+"/pt.pdf")
	os.system(f"xdg-open {d}/pt.pdf")
#_____________________________________________________________ END OF DEF

###########		PROBABILITY PLOT	###############

def Y_vProb():

	p_hist = dataY2S.Histo1D(("Y(2S) Probability", "Y(2S) Probability ;p;Counts", 500, 0., 1.), "ups_vProb")

	c_p = ROOT.TCanvas("Y2S prob", "Y2S prob")
	
	p_hist.Draw()
	c_p.SaveAs(d+"/prob.pdf")
	os.system(f"xdg-open {d}/prob.pdf")
#_____________________________________________________________ END OF DEF

##########		RAPIDITY PLOT	###########

def Y_rap():

		#	plot the cumulative data from 2016 to 2018

	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) Rapidity", "Y(2S) Rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")	# 2016
	
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) rapidity", "Y(2S) rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")	# 2017
	
	hist18 = dataY2S.Histo1D(("Y(2S) rapidity", "Y(2S) rapidity (stacked);y(#mu#mu);Counts", 500, -2.5, 2.5), "ups_rap")	# 2018

	hist17.SetFillColor(5)
	hist16.SetFillColor(3)


	c_rap = ROOT.TCanvas("Y2S Rapidity", "Y2S Rapidity")

	hist18.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	c_rap.SaveAs(d+"/rapidity.pdf")
	os.system(f"xdg-open {d}/rapidity.pdf")
#_____________________________________________________________ END OF DEF

#	PSEUDRAPIDITY PLOT	##############################

def Y_pseudorap():

		#	plot the cumulative data from 2016 to 2018

	hist18 = dataY2S.Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")	# 2018
	
	hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")	# 2016
	
	hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")	# 2017

	hist17.SetFillColor(5)	# change color to highlight
	hist16.SetFillColor(3)


	c = ROOT.TCanvas("Y2S Pseudorapidity", "Y2S Pseudorapidity")

	hist18.Draw("")
	hist17.Draw("same")
	hist16.Draw("same")

	c.SaveAs(d+"/Pseudorapidity.pdf")
	os.system(f"xdg-open {d}/Pseudorapidity.pdf")
#_____________________________________________________________ END OF DEF


#	DIMUON DECAY GEOMETRY		##############################
def dimuon_decay():

		#	angular phase plot
		
	hist1 = dataY2S.Histo2D(("dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#Delta#eta(#mu^{+}#mu^{-});#Delta#phi(#mu^{+}#mu^{-})", 50, 0., 2.5, 50, 0., 3.), "ups_deltaEta", "ups_deltaPhi")
	
	cprint(hist1, "Dimuon", "colz" )

		#	phase plot

	hist2 = dataY2S.Histo2D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#DeltaR(#mu^{+}#mu^{-});p_{T}(#mu^{+}#mu^{-}) [GeV]", 50, 0, 3, 50, 0, 150), "ups_deltaR", "ups_pT")
	
	cprint(hist2, "Dimuon2", "colz" )

		#	phase profile

	hist3 = dataY2S.Profile1D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};p_{T}(#mu^{+}#mu^{-}) [GeV];#DeltaR(#mu^{+}#mu^{-})", 50, 0, 150), "ups_pT", "ups_deltaR")
	
	cprint(hist3,"profile",stats=True)
#_____________________________________________________________ END OF DEF	


#	MENU

#	call functions by inserting the key from command line
#	couple the key with the plot function, then execute in a for loop
#	(since plot functions have no arguments, it reduces the code)
compute = {	"1" : mu_pt,
			"2" : m_Y2S,
			"3" : fit_Y2S,
			"4" : Y_pt,		
			"5" : Y_vProb,		
			"6" : Y_rap,
			"7" : Y_pseudorap,
			"8" : dimuon_decay,
			"q" : exit
		  }

	
lang = input("\nSelect plots (Separate by spacing):\n1. Single Muon pT Plot\n2. Y(2S) Mass Plot\n3. Fit of Y(2S) Mass Plot\n4. Y(2S) pT\n5. Y(2S) Vertex Probability\n6. Y(2S) Rapidity\n7. Y(2S) Pseudorapidity\n8. Plots for Geometry of Di-Muon Decay\nENTER to print all plots.\nPress \"q\" to EXIT.\n").split()

if "q" not in lang:
	print("Processing...") 
	

	# Start timer
start = time()

if not lang:	#print all plots		
	for func in compute.values() : func()

else:			#print only selected-by-key plots
	for i in lang: 
		
		#	avoid misdigit
		while i not in compute.keys(): 	
			i = input(f"\"{i}\" is not valid. Please insert a valid key:\n")
		
		#	execution
		compute[i]()
	
#	****** 	end of loop

	# End timer
end = time()

	# Calculate elapsed time
elapsed = end - start
print("\nComputing time: ", elapsed, "\n") 


fileroot.Close()

