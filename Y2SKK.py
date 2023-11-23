import ROOT 
import uproot
import pandas as pd
import os
import time
from filereader import test
ROOT.gROOT.SetBatch(True)

#vertex probability of four tracks
#tagli al pt ()



#1. eseguire i plot per tutte le variabili del jupyter notebook per la Y2s
#2. esplorare i tagli per ogni plot (provare a plottare tutti i tagli sullo stesso frame)
#3. Capire se ci sono dei tagli particolari per cui non si vede la 2s nel plot della 1s

#4. APRIRE PIU' FILE INSIEME PER AUMENTARE LA STATISTICA - GUARDARE PRIME RIGHE DI JUPYTER

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



#ROOT.ROOT.EnableImplicitMT(4)

	# Start timer
start = time.time()


#Ipotesi di aprire i file direttamente con ssh e non da disco esterno
#vedere chatgpt -> pacchetto paramiko
"""
with open('Y2SPhiRun2List.txt') as f:
    allFiles = f.readlines()

for i in range(len(allFiles)):			
    allFiles[i] = allFiles[i].replace("\n", "")
"""


#	COMPUTING FUNCTIONS
ROOT.gInterpreter.Declare('''
Double_t ComputePseudoRapidity(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    const ROOT::Math::PxPyPzE4D p(px, py, pz, ene);
    return p.Eta();
}
Double_t ComputePhi(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    const ROOT::Math::PxPyPzE4D p(px, py, pz, ene);
    return p.Phi();
}
Double_t ComputeRapidity(Double_t px, Double_t py, Double_t pz, Double_t ene) {
    return 0.5 * log( (ene+pz)/(ene-pz) );
}
''')

ROOT.gInterpreter.Declare('''
Double_t ComputeMuMuPiPiInvMass(Double_t px1, Double_t py1, Double_t pz1, 
                                    Double_t px2, Double_t py2, Double_t pz2,
                                    Double_t px3, Double_t py3, Double_t pz3,
                                    Double_t px4, Double_t py4, Double_t pz4) {
    const Double_t muonMass = 0.105658;
    const Double_t pionMass = 0.139570;
    const ROOT::Math::PxPyPzMVector mu1(px1, py1, pz1, muonMass);
    const ROOT::Math::PxPyPzMVector mu2(px2, py2, pz2, muonMass);
    const ROOT::Math::PxPyPzMVector pi1(px3, py3, pz3, pionMass);
    const ROOT::Math::PxPyPzMVector pi2(px4, py4, pz4, pionMass);
    return (mu1+mu2+pi1+pi2).mass();
}
Double_t ComputePiPiInvMass(Double_t px3, Double_t py3, Double_t pz3,
                     Double_t px4, Double_t py4, Double_t pz4) {
    const Double_t pionMass = 0.139570;
    const ROOT::Math::PxPyPzMVector pi1(px3, py3, pz3, pionMass);
    const ROOT::Math::PxPyPzMVector pi2(px4, py4, pz4, pionMass);
    return (pi1+pi2).mass();
}
''')


	# End timer
end = time.time()

	# Calculate elapsed time
elapsed = end - start
print("Files opened in", elapsed,"s\n")

#	TREES READING
sample = test(4)


dataY2SKK = ROOT.RDataFrame("rootuple/CandidateTree", sample)
if not dataY2SKK:
	print("Connect hard disk")
	exit()
dataY2S = ROOT.RDataFrame("rootuple/UpsTree", sample)

"""
#https://root.cern/doc/master/classROOT_1_1RDataFrame.html

	# Start timer
start = time.time()
	
#	Definition of NEW COLUMNS
dataY2S = dataY2S.Define("ups_pT", "sqrt(ups_p4_fX*ups_p4_fX + ups_p4_fY*ups_p4_fY)")\
.Define("ups_eta", "ComputePseudoRapidity(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
.Define("ups_phi", "ComputePhi(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
.Define("ups_rap", "ComputeRapidity(ups_p4_fX, ups_p4_fY, ups_p4_fZ, ups_p4_fE)")\
.Define("muonP_pT", "sqrt(muonP_p4_fX*muonP_p4_fX + muonP_p4_fY*muonP_p4_fY)")\
.Define("muonP_eta", "ComputePseudoRapidity(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
.Define("muonP_phi", "ComputePhi(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
.Define("muonP_rap", "ComputeRapidity(muonP_p4_fX, muonP_p4_fY, muonP_p4_fZ, muonP_p4_fE)")\
.Define("muonN_pT", "sqrt(muonN_p4_fX*muonN_p4_fX + muonN_p4_fY*muonN_p4_fY)")\
.Define("muonN_eta", "ComputePseudoRapidity(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
.Define("muonN_phi", "ComputePhi(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
.Define("muonN_rap", "ComputeRapidity(muonN_p4_fX, muonN_p4_fY, muonN_p4_fZ, muonN_p4_fE)")\
.Define("ups_deltaEta", "abs(muonP_eta - muonN_eta)")\
.Define("ups_deltaPhi", "Double_t angle = abs(muonP_phi - muonN_phi); if (angle > ROOT::Math::Pi()) angle = 2*ROOT::Math::Pi() - angle; return angle;")\
.Define("ups_deltaR", "sqrt(ups_deltaEta*ups_deltaEta + ups_deltaPhi*ups_deltaPhi)")
	
	# End timer
end = time.time()

	# Calculate elapsed time
elapsed = end - start
print("\nTime for definition of columns: ", elapsed,"\n") 


fileroot = ROOT.TFile.Open("ups_mass.root","RECREATE")


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
.Filter("ups_rap < 0.7 & ups_rap > -0.7")

###################################################################

m_filt = data_filt.AsNumpy(columns=["ups_vMass"])



"""
"""

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
#### COME COSTRUAMO IL MODELLO PER IL FIT? SEGUIAMO IL 
#### GAUS+CHEBISHEV O ALTRO? 
#### AGGIUNGERE UNA GAUSSIANA CHE RAPPRESENTI LA Y3S, CERCARE VALORE NOMINALE DELLA 3S E IMPOSTARLO COME MEDIA FISSA DELLA SECONDA GAUSSIANA
	#define variables and parameters

# mass Y2S PDG = 10.02326)

	#gaussian means
#pdgups3s= ROOT.RooRealVar ("pdgups3s","mass Y3S PDG",10.35520,9.,11.)#try fix and not to fix

mean2s = ROOT.RooRealVar("#mu_{Y(2S)}", "mean of gaussians", 10., 9.44, 10.3)
mean_reso = ROOT.RooRealVar("#mu_{res}", "mean of resolution", 10., 9.44, 10.3)
	#sigmas
sigma2s1 = ROOT.RooRealVar("#sigma_{1}^{Y(2S)}", "width1", 0.063, 0.01, 5.) 
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

#alpha = ROOT.RooRealVar ("alpha","alpha", 1.62,0.,2.)
#n= ROOT.RooRealVar ("n","n", 2.5,0.,5.);
#ups2sOne = ROOT.RooCBShape("Double CB", "#upsilon(2s) Pdf",YMass,mean,sigma,alpha,n);
  

# 		MODELS FOR MASS PLOT

	#signals
sig2s1 = ROOT.RooGaussian("signal2s_1", "signal2s_1", upsmass, mean2s, sigma2s1)
sig2s2 = ROOT.RooGaussian("signal2s_2", "signal2s_2", upsmass, mean2s, sigma2s2)
	
	#backgrounds
bkg = ROOT.RooChebychev("bkg", "Background", upsmass, ROOT.RooArgList(f0,f1,f2)) 

	#total
#model = ROOT.RooAddPdf("model", "bkg+reso+sig", [bkg,sig2s2,sig2s1], [bkgfrac,resfrac]) #with bkg
model = ROOT.RooAddPdf("model", "reso+sig", [sig2s2,sig2s1], [resfrac]) #no bkg

upsmass.setRange("range", 9.8 ,10.15)
model.fitTo(mumuroohist, Range="range")

#		PLOTTING
mumuroohist.plotOn(xframe,LineColor="b",MarkerSize=0.3)
model.plotOn(xframe,LineColor="r")
#model.plotOn(xframe, Components={bkg}, LineStyle=":", LineColor="r")
model.plotOn(xframe, Components={sig2s2}, LineStyle=":", LineColor="g")
model.plotOn(xframe, Components={sig2s1}, LineStyle=":", LineColor="b")
model.paramOn(xframe, ROOT.RooFit.Parameters([mean2s, sigma2s1, sigma2s2, resfrac]), ROOT.RooFit.Layout(0.1, 0.9, 0.9))



c = ROOT.TCanvas("MassPlotY2S", "MassPlotY2S", 800, 800)
c.Divide(1,2)
c.cd(1)
xframe.Draw()
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


################################################################
#	CANDIDATES KK DECAY	########################################
################################################################


#	DEFINITIONS
dataY2SKK = dataY2SKK\
.Define("ditrack_mass", "sqrt(ditrack_p4_fE*ditrack_p4_fE - ditrack_p4_fX*ditrack_p4_fX - ditrack_p4_fY*ditrack_p4_fY - ditrack_p4_fZ*ditrack_p4_fZ)")\
.Define("dimuon_mass", "sqrt(dimuon_p4_fE*dimuon_p4_fE - dimuon_p4_fX*dimuon_p4_fX - dimuon_p4_fY*dimuon_p4_fY - dimuon_p4_fZ*dimuon_p4_fZ)")\
.Define("track1_pT", "sqrt(track1_p4_fX*track1_p4_fX + track1_p4_fY*track1_p4_fY)")\
.Define("track1_eta", "ComputePseudoRapidity(track1_p4_fX, track1_p4_fY, track1_p4_fZ, track1_p4_fE)")\
.Define("track2_pT", "sqrt(track2_p4_fX*track2_p4_fX + track2_p4_fY*track2_p4_fY)")\
.Define("track2_eta", "ComputePseudoRapidity(track2_p4_fX, track2_p4_fY, track2_p4_fZ, track2_p4_fE)")\
.Define("ditrack_pT", "sqrt(ditrack_p4_fX*ditrack_p4_fX + ditrack_p4_fY*ditrack_p4_fY)")\
.Define("ditrack_eta", "ComputePseudoRapidity(ditrack_p4_fX, ditrack_p4_fY, ditrack_p4_fZ, ditrack_p4_fE)")\
.Define("dimuon_pT", "sqrt(dimuon_p4_fX*dimuon_p4_fX + dimuon_p4_fY*dimuon_p4_fY)")\
.Define("dimuon_eta", "ComputePseudoRapidity(dimuon_p4_fX, dimuon_p4_fY, dimuon_p4_fZ, dimuon_p4_fE)")\
.Define("candidate_pT", "sqrt(candidate_p4_fX*candidate_p4_fX + candidate_p4_fY*candidate_p4_fY)")\
.Define("candidate_eta", "ComputePseudoRapidity(candidate_p4_fX, candidate_p4_fY, candidate_p4_fZ, candidate_p4_fE)")\
.Define("trackL_pT", "max(track1_pT, track2_pT)")\
.Define("trackS_pT", "min(track1_pT, track2_pT)")\
.Define("mumupipi_mass", "ComputeMuMuPiPiInvMass(muonp_p4_fX, muonp_p4_fY, muonp_p4_fZ, muonn_p4_fX, muonn_p4_fY, muonn_p4_fZ, track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")\
.Define("pipi_mass", "ComputePiPiInvMass(track1_p4_fX, track1_p4_fY, track1_p4_fZ, track2_p4_fX, track2_p4_fY, track2_p4_fZ)")


#############################################
											#
d="time_test"								#
if not os.path.isdir(f"./{d}/"):			#
	os.system(f"mkdir {d}")					#
# ROOT FILE									#
p = ROOT.TFile.Open(d+"/h.root","RECREATE")	#
											#
#############################################


"""
#	MASS PLOT OF CANDIDATE	##############################
	# Start timer
start = time.time()

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

#KK filters
dataKKrs = dataY2SKK.Filter("track1_charge * track2_charge < 0")
dataKKws = dataY2SKK.Filter("track1_charge * track2_charge > 0")

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
	start = time.time()

		#borders definitions
	maxValues = [0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0]
	minValues = [0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8]
		
		#cuts definitions
	cuts = ["track1_charge * track2_charge < 0"]

	for i in range(len(maxValues)):
		    cuts.append(cuts[0] +" & trackL_pT > " + str(maxValues[i]) + " & trackS_pT > " + str(minValues[i]))
		  
		#create array of histograms with all cuts  
	hists = []

	for cut in cuts:
		hists.append(dataY2SKK.Filter(cut).Histo1D(("KK invariant mass", "#phi #rightarrow K^{+}K^{-};m(KK) [GeV];Counts", 100, 0.99, 1.06), "ditrack_mass"))


	a = [0., 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0]
	b = [0., 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8]

	c0 = ROOT.TCanvas()

	legend = ROOT.TLegend(0.46, 0.15, 0.89, 0.50)

	hists[0].SetStats(0)	#no stats
	print("type of hists: ", type(hists[0]))

		#overlap plots with cuts
	for i, hist in enumerate(hists):
		
		hist.Scale(1./hist.Integral())	#normalization of plots
		hist.GetYaxis().SetRangeUser(0, 0.015)	#range of y axis
		hist.SetLineColor(i+1)					
		legend.AddEntry(hist.GetPtr(), "p_{T} lead > " + str(a[i]) + "; p_{T} soft > " + str(b[i]), "l")	
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
def m_kk_fit():
	dfkk = dataKKrs.Filter("trackL_pT > 1.2 & trackS_pT > 1.0")\
	.Filter("ditrack_mass > 1.0 & ditrack_mass < 1.04").AsNumpy(columns=["ditrack_mass"])
		
		#variable
	kkmass = ROOT.RooRealVar("ditrack_mass", "m(KK) [GeV]", 1.00, 1.04)

	kkroodata = ROOT.RooDataSet.from_numpy({"ditrack_mass": dfkk["ditrack_mass"]}, [kkmass])
	kkroohist = kkroodata.binnedClone()
	phiframe = kkmass.frame(Title="Dikaon Candidate Mass")

	
		#model
	mean = ROOT.RooRealVar("#mu_{#phi}", "mean of gaussian", 1.019, 1.018, 1.020) #mean value, min value, max value
	sigma = ROOT.RooRealVar("#sigma_{#phi}", "resolution", 0.00125, 0.001, 0.002) # fixed res from MC
	sigma.setConstant(1)
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
	end = time.time()

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
phiKKwsSelection = '''candidate_vProb > 0.1 &
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
def phi_cand(zoom = False):
	hist0 = dataKKrs.Filter(ymumu_filter + phiKKSelection).Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", binning, edge[0], edge[1]), "candidate_vMass")
	hist1 = dataKKws.Filter(ymumu_filter + phiKKwsSelection).Histo1D(("MuMuKK cands", "Y(2S)(#rightarrow #mu^{+}#mu^{-})#phi(#rightarrow K^{+}K^{-});m(#mu#muKK) - m(#mu#mu) + m^{PDG}(Y) [GeV];Counts", binning, edge[0], edge[1]), "candidate_vMass")

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
		os.system(f"xdg-open {d}/PhiCandidateZoom.pdf")


#	MENU
lang = input("Select plots (Separate by spacing):\n1. Leading Kaon pt\n2. Soft Kaon pt\n3. KK invariant mass (with K_pt cuts)\n4. Fit KK invariant mass\n5. Phi Candidate plot\n6. Phi candidate plot with Zoom\nENTER: 1-5 Plots\n").split()
#print(lang)
#exit()
print("Processing...")

	# Start timer
start = time.time()
if not lang:
#print all plots
		kaonL_pt()
		kaonS_pt()
		m_kk()					
		m_kk_fit()				#save to .root
		phi_cand()				#save to .root
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
				phi_cand()		#save to .root
			case "6":			
				phi_cand(zoom=True)
			case _:
				print("Not valid")
				exit()


	# End timer
end = time.time()

	# Calculate elapsed time
elapsed = end - start
print("\nComputing time: ", elapsed, "\n") 
os.system(f"echo {elapsed} > time.txt")
p.Close()

#blablabla
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
