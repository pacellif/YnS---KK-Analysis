import ROOT 
import pandas as pd
import os
import time
from sys import argv
from definitions import UpsTreeDefinitions
import declarations

ROOT.gROOT.SetBatch(True)


#directory: /lustre/cms/store/user/vmastrap/MuMuKKRun2/Y2SKKfromSametV4
# comando scp per copiare macro in account remoto
# cambiare il path per aprire i file


#	PLOT TEMPLATE
def cprint (hist, name, opt="", stats=False, x=300,y=200):
	title= "Y2S "+name
	c = ROOT.TCanvas(title, title, x, y)
	
	if stats==False: hist.SetStats(0)
	hist.Draw(opt)
		
	c.SaveAs(name+".pdf")
	os.system("xdg-open "+name+".pdf")

def binCount (var, rangeName):
	varBinning = var.getBinning()
	a = var.getRange(rangeName)[0]
	b = var.getRange(rangeName)[1]
	nbins = varBinning.binNumber(b) - varBinning.binNumber(a)
	return nbins
	

#ROOT.ROOT.EnableImplicitMT(4)

with open('Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")


#	TREES READING - using command line arguments
try: 
	sample = allFiles[:int(argv[1])] 
except IndexError:
	sample = allFiles[:]


dataY2S = UpsTreeDefinitions(ROOT.RDataFrame("rootuple/UpsTree",sample))
print("Data opened")

fileroot = ROOT.TFile.Open("ups_plots/ups_mass.root","RECREATE")

	# Start timer
start = time.time()


#	VARIABLES	###################################################

upsmass = ROOT.RooRealVar("ups_vMass", "m(#mu#mu) [GeV]", 9.6, 10.3)
upsmass.setBins(500)

VertexProb = ROOT.RooRealVar("ups_vProb", "ups_vProb", 0, 1)
VertexProb.setBins(1000)

YpT = ROOT.RooRealVar("ups_pT", "p_{T}(#mu#mu) [GeV/c]",0.0, 66.4)
rap = ROOT.RooRealVar("ups_rap", "rap(#mu#mu) ",-2.5, 2.5)



#	FILTERS + DICTIONARY OF FILTERS

###################################################################

data_filt = dataY2S\
.Filter("ups_vMass > 9.6 & ups_vMass < 10.3")\
.Filter("ups_pT > 25")

###################################################################

m_filt = data_filt.AsNumpy(columns=["ups_vMass"])






pt_filt = dataY2S.Filter("ups_pT < 97.4").AsNumpy(columns=["ups_pT"])
rap_filt= dataY2S.Filter("ups_rap < 2.45 & ups_rap > -2.5").AsNumpy(columns=["ups_rap"])
p_filt= dataY2S.Filter("ups_rap < 1 & ups_rap > 0").AsNumpy(columns=["ups_vProb"])

#rooDataSet = dataY2S.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(YMass,VertexProb, YpT))), ("ups_vMass","ups_vProb", "ups_pT"))
#rooDataSet = dataY2S.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(YMass))), ["ups_vMass"])

#------> DICTIONARY
d = {
"ups_vMass": m_filt["ups_vMass"],
"ups_pT": pt_filt["ups_pT"],
"ups_rap": rap_filt["ups_rap"],
"ups_vProb": p_filt["ups_vProb"]
}


#	DATASETS

mumuroodata = ROOT.RooDataSet.from_numpy(d, ROOT.RooArgSet(upsmass))
mumuroohist = mumuroodata.binnedClone()

#pt_ds = ROOT.RooDataSet.from_numpy(d, ROOT.RooArgSet(YpT))
#rap_ds = ROOT.RooDataSet.from_numpy(d, ROOT.RooArgSet(rap))
#prob_ds = ROOT.RooDataSet.from_numpy(d, ROOT.RooArgSet(VertexProb))


######	   MASS PLOT	#########

xframe = upsmass.frame(Title="Y(2S) Mass")
hist = data_filt\
.Histo1D(("dimuon invariant mass", "Y(2S) #rightarrow #mu^{+}#mu^{-};m(#mu^{+}#mu^{-}) [GeV];Counts", 500, 9.6, 10.3), "ups_vMass")



	#	FIT
	#define variables and parameters

# mass Y2S PDG = 10.02326

	#gaussian means
#pdgups3s= ROOT.RooRealVar ("pdgups3s","mass Y3S PDG",10.35520,9.,11.)#try fix and not to fix

mean2s = ROOT.RooRealVar("#mu_{Y(2S)}", "mean of gaussians", 10., 9.44, 10.3)
mean_reso = ROOT.RooRealVar("#mu_{res}", "mean of resolution", 10., 9.44, 10.3)
	#sigmas
sigma2s = ROOT.RooRealVar("#sigma^{Y(2S)}", "width1", 0.063, 0.01, 5.) 
sigma2s2 = ROOT.RooRealVar("#sigma_{res}^{Y(2S)}", "width2", 0.144, 0.01, 5.) 
#sigma3s = ROOT.RooRealVar("#sigma_{2}^{Y(3S)}", "width3", 0.144, 0.01, 5.)
	
	#chebychev coefficients
f0 = ROOT.RooRealVar("f0", "f0", 0.1, -8., 8.)
f1 = ROOT.RooRealVar("f1", "f1", 0.1, -8., 8.)
f2 = ROOT.RooRealVar("f2", "f2", 0.1, -8., 8.)
	
	#fractions
bkgfrac = ROOT.RooRealVar("f_{bkg}", "fraction of background", 0.02, 0.001, 1.)
resfrac = ROOT.RooRealVar("f_{reso}", "fraction of resolution", 0.45, 0., 1.)
#sig2frac = ROOT.RooRealVar("f_{signal}", "fraction of second signal of 2s", 0.45, 0., 1.)

alpha = ROOT.RooRealVar ("alpha","alpha", 1.62, 0., 5.)
n= ROOT.RooRealVar ("n","n", 0.1, -1., 1.);

  

# 		MODELS FOR MASS PLOT

	#signals
sig2s = ROOT.RooGaussian("signal2s", "signal2s", upsmass, mean2s, sigma2s)
cb = ROOT.RooCBShape("Double CB", "#upsilon(2s) Pdf", upsmass, mean2s, sigma2s, alpha, n)
	
	#backgrounds
bkg = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1,f2)) 

	#total
#model = ROOT.RooAddPdf("model", "bkg+sig", [bkg,sig2s], [bkgfrac]) 
model = ROOT.RooAddPdf("model", "bkg+sig_cb", [bkg,cb], [bkgfrac]) 



########	FIT RESULT AND CHI SQUARED COMPUTATION  ·················
upsmass.setRange("range", 9.8 ,10.15) #set range before fitting
	
fitResult = model.fitTo(mumuroohist, Range="range", Save=True)

	#chi squared
chiSquared = int(-2 * fitResult.minNll())
ndof = binCount(upsmass,"range") - model.getParameters(mumuroodata).getSize()

reducedChiSquared = round(chiSquared/ndof)

########	PLOTTING	···········································
mumuroohist.plotOn(xframe,LineColor="b",MarkerSize=0.3)
model.plotOn(xframe,LineColor="r")


model.plotOn(xframe,LineColor="r")
model.plotOn(xframe, Components={bkg}, LineStyle=":", LineColor="g")
model.plotOn(xframe, Components={sig2s}, LineStyle=":", LineColor="b")
model.paramOn(xframe, ROOT.RooFit.Layout(0.1, 0.9, 0.9)) #print all parameters

	#print chisquare
text_box = ROOT.TPaveText(0.65, 0.75, 0.9, 0.9, "NDC")
text_box.AddText("CB+Cheb3")	#type of fit
text_box.AddText(f"#chi^{2}/ndof = {reducedChiSquared}")
text_box.SetFillColor(0)
text_box.SetBorderSize(1)


c = ROOT.TCanvas("MassPlotY2S", "MassPlotY2S", 800, 800)
c.Divide(1,2)
c.cd(1)
xframe.Draw()
text_box.Draw()
c.cd(2)
hist.Draw()

fileroot.WriteObject(c,"UpsInvMass")
c.SaveAs("ups_plots/MassPlotY2S.pdf")
os.system("xdg-open ups_plots/MassPlotY2S.pdf")

	# End timer
end = time.time()

	# Calculate elapsed time
elapsed = end - start
print("\nTime for Mass Plot: ", elapsed,"\n") 

fileroot.Close()


"""
	# Start timer
start = time.time()

#######		TRANSVERSAL MOMENTUM PLOTS	#############

#pt_frame = YpT.frame(Title="Y(2S) Transverse Momentum")
#pt_ds.plotOn(pt_frame,LineColor="b",MarkerSize=0.3)

hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) transverse momentum", "Y(1S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")
hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) transverse momentum", "Y(1S) transverse momentum;p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50), "ups_pT")

hist18 = dataY2S.Histo1D(("Y(2S) transverse momentum", "Y(2S) transverse momentum (stacked);p_{T}(#mu#mu) [GeV];Counts", 500, 0., 50.), "ups_pT")

c_pt = ROOT.TCanvas("Y(2S) pT", "Y(2S) pT", 800, 800)


hist17.SetFillColor(5)
hist16.SetFillColor(3)


hist18.Draw("")
hist17.Draw("same")
hist16.Draw("same")

c_pt.SaveAs("ups_plots/pt.pdf")
os.system("xdg-open ups_plots/pt.pdf")

###########		PROBABILITY PLOT	###############
#pframe = VertexProb.frame(Title="Y(2S) Probability")
#prob_ds.plotOn(pframe, LineColor="b", MarkerSize= 0.3)

p_hist = dataY2S.Histo1D(("Y(2S) Probability", "Y(2S) Probability ;p;Counts", 500, 0., 1.), "ups_vProb")

c_p = ROOT.TCanvas("Y2S prob", "Y2S prob", 800, 800)
#c_p.Divide(1,2)
#c_p.cd(1)
#pframe.Draw()
#c_p.cd(2)
p_hist.Draw()
c_p.SaveAs("ups_plots/prob.pdf")
os.system("xdg-open ups_plots/prob.pdf")


##########		RAPIDITY PLOT	###########
#rap_frame = rap.frame(Title="Y(2S) Rapidity")
#rap_ds.plotOn(rap_frame,LineColor="b",MarkerSize=0.3)

hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) Rapidity", "Y(2S) Rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")

#cprint(hist16, "ups_plots/Rapidity")


hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) rapidity", "Y(2S) rapidity;y(#mu#mu);Counts", 500,-2.5,2.5), "ups_rap")
hist19 = dataY2S.Histo1D(("Y(2S) rapidity", "Y(2S) rapidity (stacked);y(#mu#mu);Counts", 500, -2.5, 2.5), "ups_rap")

hist17.SetFillColor(5)
hist16.SetFillColor(3)


c_rap = ROOT.TCanvas("Y2S Rapidity", "Y2S Rapidity", 800, 800)

hist19.Draw("")
hist17.Draw("same")
hist16.Draw("same")

c_rap.SaveAs("ups_plots/rapidity.pdf")
os.system("xdg-open ups_plots/rapidity.pdf")


#	PSEUDRAPIDITY PLOT	##############################
hist18 = dataY2S.Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")
hist16 = dataY2S.Filter("run < 290000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")
hist17 = dataY2S.Filter("run < 310000").Histo1D(("Y(2S) pseudorapidity", "Y(2S) Pseudorapidity (stacked);#eta(#mu#mu);Counts", 500, -2.0, 2.0), "ups_eta")

hist17.SetFillColor(5)
hist16.SetFillColor(3)


c_psrap = ROOT.TCanvas("Y2S Pseudorapidity", "Y2S Pseudorapidity", 800, 800)

hist18.Draw("")
hist17.Draw("same")
hist16.Draw("same")

c_psrap.SaveAs("ups_plots/Pseudorapidity.pdf")
os.system("xdg-open ups_plots/Pseudorapidity.pdf")

#cprint(hist18, "ups_plots/Pseudorapidity")

"""
"""
#	DIMUON DECAY - ANGULAR PHASE PLOT	##############################
hist = data_filt.Histo2D(("dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#Delta#eta(#mu^{+}#mu^{-});#Delta#phi(#mu^{+}#mu^{-})", 50, 0., 2.5, 50, 0., 3.), "ups_deltaEta", "ups_deltaPhi")
cprint(hist, "ups_plots/Dimuon", "colz" )


#	DIMUON DECAY - PHASE PLOT	 ###########################################
hist = data_filt.Histo2D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};#DeltaR(#mu^{+}#mu^{-});p_{T}(#mu^{+}#mu^{-}) [GeV]", 50, 0, 3, 50, 0, 150), "ups_deltaR", "ups_pT")
cprint(hist, "ups_plots/Dimuon2", "colz" )

hist = data_filt.Profile1D(("Dimuon decay", "Y(2S) #rightarrow #mu^{+}#mu^{-};p_{T}(#mu^{+}#mu^{-}) [GeV];#DeltaR(#mu^{+}#mu^{-})", 50, 0, 150), "ups_pT", "ups_deltaR")
cprint(hist,"ups_plots/profile",stats=True)
	
	# End timer
end = time.time()

	# Calculate elapsed time
elapsed = end - start
print("\nTime for other Y Plots: ", elapsed,"\n") 
"""


# c.Divide(2)
# c.cd(1)
# ROOT.gPad.SetLeftMargin(0.15)
# xframe.GetYaxis().SetTitleOffset(1.6)

# c.cd(2)
# ROOT.gPad.SetLeftMargin(0.15)
# xframe2.GetYaxis().SetTitleOffset(1.6)
# xframe2.Draw()





#From root to Pandas
# Y2S_tree = myfileY2SKK["rootuple/UpsTree;1"]
# print(Y2S_tree.keys())

# Y2SKK_tree = myfileY2SKK["rootuple/CandidateTree;1"]
# print(Y2SKK_tree.keys())

# data_Y2S = Y2S_tree.arrays(library="pd")
# print(data_Y2S.head())











#Roofit ToolBox

# x = ROOT.RooRealVar("x", "x", -10, 10)
# mean = ROOT.RooRealVar("mean", "mean of gaussian", 1, -10, 10)
# sigma = ROOT.RooRealVar("sigma", "width of gaussian", 1, 0.1, 10)
 
# # Build gaussian pdf in terms of x,mean and sigma
# gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)
 
# # Construct plot frame in 'x'
# xframe = x.frame(Title="Gaussian pdf")  # RooPlot
 
# # Plot model and change parameter values
# # ---------------------------------------------------------------------------
# # Plot gauss in frame (i.e. in x)
# gauss.plotOn(xframe)
 
# # Change the value of sigma to 3
# sigma.setVal(3)
 
# # Plot gauss in frame (i.e. in x) and draw frame on canvas
# gauss.plotOn(xframe, LineColor="r")

# data = gauss.generate({x}, 10000)  # ROOT.RooDataSet
 
# # Make a second plot frame in x and draw both the
# # data and the pdf in the frame
# xframe2 = x.frame(Title="Gaussian pdf with data")  # RooPlot
# data.plotOn(xframe2)
# gauss.plotOn(xframe2)
 
# # Fit model to data
# # -----------------------------
# # Fit pdf to data
# gauss.fitTo(data)
 
# # Print values of mean and sigma (that now reflect fitted values and
# # errors)
# mean.Print()
# sigma.Print()
 
# # Draw all frames on a canvas
# c = ROOT.TCanvas("rf101_basics", "rf101_basics", 800, 400)
# c.Divide(2)
# c.cd(1)
# ROOT.gPad.SetLeftMargin(0.15)
# xframe.GetYaxis().SetTitleOffset(1.6)
# xframe.Draw()
# c.cd(2)
# ROOT.gPad.SetLeftMargin(0.15)
# xframe2.GetYaxis().SetTitleOffset(1.6)
# xframe2.Draw()
 
# c.SaveAs("rf101_basics.png")

#Upsilon:
['run', 'event', 'numPrimaryVertices', 'trigger', 'ups_p4', 'muonP_p4', 'muonN_p4', 'iPVwithmuons_ups', 'ups_vertexWeight', 'ups_vProb', 'ups_vMass', 'ups_vNChi2', 'ups_DCA', 'ups_ctauPV', 'ups_ctauErrPV', 'ups_lxyPV', 'ups_lxyErrPV', 'ups_cosAlpha', 'ups_ctauBS', 'ups_ctauErrBS', 'ups_lxyBS', 'ups_lxyErrBS', 'mu1_pt', 'mu1_ptErr', 'mu1_d0', 'mu1_d0Err', 'mu1_dz', 'mu1_dzErr', 'mu1_dxy', 'mu1_dxyErr', 'mu1_nvsh', 'mu1_nvph', 'mu1_charge', 'mu2_pt', 'mu2_ptErr', 'mu2_d0', 'mu2_d0Err', 'mu2_dz', 'mu2_dzErr', 'mu2_dxy', 'mu2_dxyErr', 'mu2_nvsh', 'mu2_nvph', 'mu2_charge']
#Upsilon KKcandidate_vMass
['run', 'event', 'nCandPerEvent', 'numPrimaryVertices', 'trigger', 'candidate_p4', 'track1_p4', 'track2_p4', 'ditrack_p4', 'dimuon_p4', 'muonp_p4', 
'muonn_p4', 'iPVwithmuons', 'dimuon_diMuIndx', 'dimuon_vertexWeight', 'dimuon_vProb', 'dimuon_vMass', 'dimuon_vNChi2', 'dimuon_DCA', 'dimuon_ctauPV', 
'dimuon_ctauErrPV', 'dimuon_lxyPV', 'dimuon_lxyErrPV', 'dimuon_cosAlpha', 'dimuon_ctauBS', 'dimuon_ctauErrBS', 'dimuon_lxyBS', 'dimuon_lxyErrBS', 
'candidate_vMass', 'candidate_vProb', 'candidate_vChi2', 'candidate_cosAlpha', 'candidate_ctauPV', 'candidate_ctauErrPV', 'candidate_charge', 
'candidate_lxy', 'candidate_lxyErr', 'candidate_lxyz', 'candidate_lxyzErr', 'thePrimaryV_X', 'thePrimaryV_Y', 'thePrimaryV_Z', 
'TheDecayVertex_X', 'TheDecayVertex_Y', 'TheDecayVertex_Z', 'thePrimaryV_2D_position', 'thePrimaryV_3D_position', 'TheDecayVertex_2D_position', 
'TheDecayVertex_3D_position', 'TheVertexDistance_2D', 'TheVertexDistance_3D', 'track1_d0', 'track1_d0Err', 'track1_dz', 'track1_dzErr', 'track1_dxy',
'track1_dxyErr', 'track1_nvsh', 'track1_nvph', 'track1_dRdimuon', 'track1_charge', 'track1_PV', 'track1_refVtx', 'track1_pvAssocQ', 'track1_dzAssocPV',
'track2_d0', 'track2_d0Err', 'track2_dz', 'track2_dzErr', 'track2_dxy', 'track2_dxyErr', 'track2_nvsh', 'track2_nvph', 'track2_dRdimuon', 
'track2_charge', 'track2_PV', 'track2_refVtx', 'track2_pvAssocQ', 'track2_dzAssocPV', 'ditrack_dRdimuon']
