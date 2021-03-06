
import ROOT
import sys,os
import array
from math import *
from ROOT import *
from array import *
import runlist
from runlist import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

#Elev_cut_lower = 55
#Elev_cut_upper = 85
#Elev_cut_lower = 65
#Elev_cut_upper = 70
#Elev_cut_lower = 70
#Elev_cut_upper = 75
#Elev_cut_lower = 75
#Elev_cut_upper = 85

UseLooseControlRegions = True
#UseLooseControlRegions = False

#target = '2ndCrab'
#target = 'PKS1424'
#target = '3C264'
#target = 'H1426'
#target = 'Ton599'
#target = 'IC443'

#method = 'Depth'
#method = 'MSCW'
#method = 'MSCL'

#region = 'SR'
#region = 'VR'

target = sys.argv[1]
method = sys.argv[2]
region = sys.argv[3]
Elev_cut_lower = float(sys.argv[4])
Elev_cut_upper = Elev_cut_lower+5
Elev_cut_lower = 55
Elev_cut_upper = 85

Hist_MSCW_Crab_On_Alpha = ROOT.TH1D("Hist_MSCW_Crab_On_Alpha","",110,-1,10)
Hist_MSCW_Crab_Off_Alpha = ROOT.TH1D("Hist_MSCW_Crab_Off_Alpha","",110,-1,10)
Hist_Loose_MSCW_Crab_On = ROOT.TH1D("Hist_Loose_MSCW_Crab_On","",110,-1,10)
Hist_Loose_MSCW_Crab_Off = ROOT.TH1D("Hist_Loose_MSCW_Crab_Off","",110,-1,10)
Hist_Loose_MSCW_Target_On = ROOT.TH1D("Hist_Loose_MSCW_Target_On","",110,-1,10)
Hist_Loose_MSCL_Crab_On = ROOT.TH1D("Hist_Loose_MSCL_Crab_On","",110,-1,10)
Hist_Loose_MSCL_Crab_Off = ROOT.TH1D("Hist_Loose_MSCL_Crab_Off","",110,-1,10)
Hist_Loose_MSCL_Target_On = ROOT.TH1D("Hist_Loose_MSCL_Target_On","",110,-1,10)
Hist_Loose_MSCW_Crab_On_TeV = ROOT.TH1D("Hist_Loose_MSCW_Crab_On_TeV","",110,-1,10)
Hist_Loose_MSCW_Crab_Off_TeV = ROOT.TH1D("Hist_Loose_MSCW_Crab_Off_TeV","",110,-1,10)
Hist_Loose_MSCW_Target_On_TeV = ROOT.TH1D("Hist_Loose_MSCW_Target_On_TeV","",110,-1,10)
Hist_Loose_MSCL_Crab_On_TeV = ROOT.TH1D("Hist_Loose_MSCL_Crab_On_TeV","",110,-1,10)
Hist_Loose_MSCL_Crab_Off_TeV = ROOT.TH1D("Hist_Loose_MSCL_Crab_Off_TeV","",110,-1,10)
Hist_Loose_MSCL_Target_On_TeV = ROOT.TH1D("Hist_Loose_MSCL_Target_On_TeV","",110,-1,10)
Hist_Loose_SlantDepth_Crab_On = ROOT.TH1D("Hist_Loose_SlantDepth_Crab_On","",100,0,30)
Hist_Loose_SlantDepth_Crab_Off = ROOT.TH1D("Hist_Loose_SlantDepth_Crab_Off","",100,0,30)
Hist_Loose_SlantDepth_Target_On = ROOT.TH1D("Hist_Loose_SlantDepth_Target_On","",100,0,30)
Hist_Loose_SlantDepth_Crab_On_TeV = ROOT.TH1D("Hist_Loose_SlantDepth_Crab_On_TeV","",100,0,30)
Hist_Loose_SlantDepth_Crab_Off_TeV = ROOT.TH1D("Hist_Loose_SlantDepth_Crab_Off_TeV","",100,0,30)
Hist_Loose_SlantDepth_Target_On_TeV = ROOT.TH1D("Hist_Loose_SlantDepth_Target_On_TeV","",100,0,30)
Hist_Tight_MSCW_Crab_On = ROOT.TH1D("Hist_Tight_MSCW_Crab_On","",110,-1,10)
Hist_Tight_MSCW_Crab_Off = ROOT.TH1D("Hist_Tight_MSCW_Crab_Off","",110,-1,10)
Hist_Tight_MSCW_Target_On = ROOT.TH1D("Hist_Tight_MSCW_Target_On","",110,-1,10)
Hist_Tight_MSCL_Crab_On = ROOT.TH1D("Hist_Tight_MSCL_Crab_On","",110,-1,10)
Hist_Tight_MSCL_Crab_Off = ROOT.TH1D("Hist_Tight_MSCL_Crab_Off","",110,-1,10)
Hist_Tight_MSCL_Target_On = ROOT.TH1D("Hist_Tight_MSCL_Target_On","",110,-1,10)
Hist_Tight_MSCW_Crab_On_TeV = ROOT.TH1D("Hist_Tight_MSCW_Crab_On_TeV","",110,-1,10)
Hist_Tight_MSCW_Crab_Off_TeV = ROOT.TH1D("Hist_Tight_MSCW_Crab_Off_TeV","",110,-1,10)
Hist_Tight_MSCW_Target_On_TeV = ROOT.TH1D("Hist_Tight_MSCW_Target_On_TeV","",110,-1,10)
Hist_Tight_MSCL_Crab_On_TeV = ROOT.TH1D("Hist_Tight_MSCL_Crab_On_TeV","",110,-1,10)
Hist_Tight_MSCL_Crab_Off_TeV = ROOT.TH1D("Hist_Tight_MSCL_Crab_Off_TeV","",110,-1,10)
Hist_Tight_MSCL_Target_On_TeV = ROOT.TH1D("Hist_Tight_MSCL_Target_On_TeV","",110,-1,10)
Hist_Tight_SlantDepth_Crab_On = ROOT.TH1D("Hist_Tight_SlantDepth_Crab_On","",100,0,30)
Hist_Tight_SlantDepth_Crab_Off = ROOT.TH1D("Hist_Tight_SlantDepth_Crab_Off","",100,0,30)
Hist_Tight_SlantDepth_Target_On = ROOT.TH1D("Hist_Tight_SlantDepth_Target_On","",100,0,30)
Hist_Tight_SlantDepth_Crab_On_TeV = ROOT.TH1D("Hist_Tight_SlantDepth_Crab_On_TeV","",100,0,30)
Hist_Tight_SlantDepth_Crab_Off_TeV = ROOT.TH1D("Hist_Tight_SlantDepth_Crab_Off_TeV","",100,0,30)
Hist_Tight_SlantDepth_Target_On_TeV = ROOT.TH1D("Hist_Tight_SlantDepth_Target_On_TeV","",100,0,30)

Hist2D_TelElev_vs_SlantDepth_Target = TH2D("Hist2D_TelElev_vs_SlantDepth_Target","",20,30,85,48,0,36)
Hist2D_TelElev_vs_EmissionHeight_Target = TH2D("Hist2D_TelElev_vs_EmissionHeight_Target","",20,30,85,40,0,40)

Hist_TelElevAzim_Counts_Target = TH2D("Hist_TelElevAzim_Counts_Target","",100,30,85,100,0,360)
Hist_TelElevAzim_Counts_Crab = TH2D("Hist_TelElevAzim_Counts_Crab","",100,30,85,100,0,360)

energy_bins = []
#energy_bins += [pow(10,2.2),pow(10,2.4),pow(10,2.6),pow(10,2.8)]
#energy_bins += [pow(10,3.0),pow(10,3.2),pow(10,3.4),pow(10,3.6),pow(10,3.8)]
#energy_bins += [pow(10,4.0),pow(10,4.2),pow(10,4.4)]
energy_bins += [pow(10,2.1),pow(10,2.2),pow(10,2.3),pow(10,2.4),pow(10,2.5),pow(10,2.6),pow(10,2.7),pow(10,2.8),pow(10,2.9)]
energy_bins += [pow(10,3.0),pow(10,3.1),pow(10,3.2),pow(10,3.3),pow(10,3.4),pow(10,3.5),pow(10,3.6),pow(10,3.7),pow(10,3.8)]
energy_bins += [pow(10,3.9),pow(10,4.0),pow(10,4.1),pow(10,4.2),pow(10,4.3),pow(10,4.4)]
Hist_ErecS_Target_SR = []
Hist_ErecS_Target_CR = []
Hist_ErecS_Target_Bkg = []
for energy in range(0,len(energy_bins)-1):
    Hist_ErecS_Target_SR += [ROOT.TH1D("Hist_ErecS_Target_SR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",len(energy_bins)-1,array('d',energy_bins))]
    Hist_ErecS_Target_CR += [ROOT.TH1D("Hist_ErecS_Target_CR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",len(energy_bins)-1,array('d',energy_bins))]
    Hist_ErecS_Target_Bkg += [ROOT.TH1D("Hist_ErecS_Target_Bkg_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",len(energy_bins)-1,array('d',energy_bins))]

Hist2D_Energy_vs_SlantDepth_Target = TH2D("Hist2D_Energy_vs_SlantDepth_Target","",len(energy_bins)-1,array('d',energy_bins),60,0,36)
Hist2D_Energy_vs_SlantDepth_Crab_On = TH2D("Hist2D_Energy_vs_SlantDepth_Crab_On","",len(energy_bins)-1,array('d',energy_bins),60,0,36)
Hist2D_Energy_vs_SlantDepth_Crab_Off = TH2D("Hist2D_Energy_vs_SlantDepth_Crab_Off","",len(energy_bins)-1,array('d',energy_bins),60,0,36)

Hist_On_SR_Count = []
Hist_On_CR_Count = []
Hist_SigAtt = []
Hist_Off_SR_Count = []
Hist_Off_CR_Count = []
Hist_BkgAtt = []
for energy in range(0,len(energy_bins)-1):
    Hist_On_SR_Count += [ROOT.TH1D("Hist_On_SR_Count_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]
    Hist_On_CR_Count += [ROOT.TH1D("Hist_On_CR_Count_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]
    Hist_SigAtt += [ROOT.TH1D("Hist_SigAtt_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]
    Hist_Off_SR_Count += [ROOT.TH1D("Hist_Off_SR_Count_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]
    Hist_Off_CR_Count += [ROOT.TH1D("Hist_Off_CR_Count_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]
    Hist_BkgAtt += [ROOT.TH1D("Hist_BkgAtt_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",1,0,1)]

Hist_Theta2_Target_SR = []
Hist_Theta2_Target_CR = []
Hist_Theta2_Target_Bkg = []
for energy in range(0,len(energy_bins)-1):
    Hist_Theta2_Target_SR += [ROOT.TH1D("Hist_Theta2_Target_SR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,4)]
    Hist_Theta2_Target_CR += [ROOT.TH1D("Hist_Theta2_Target_CR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,4)]
    Hist_Theta2_Target_Bkg += [ROOT.TH1D("Hist_Theta2_Target_Bkg_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,4)]

Hist_Theta2ZoomIn_Target_SR = []
Hist_Theta2ZoomIn_Target_CR = []
Hist_Theta2ZoomIn_Target_Bkg = []
for energy in range(0,len(energy_bins)-1):
    Hist_Theta2ZoomIn_Target_SR += [ROOT.TH1D("Hist_Theta2ZoomIn_Target_SR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,0.5)]
    Hist_Theta2ZoomIn_Target_CR += [ROOT.TH1D("Hist_Theta2ZoomIn_Target_CR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,0.5)]
    Hist_Theta2ZoomIn_Target_Bkg += [ROOT.TH1D("Hist_Theta2ZoomIn_Target_Bkg_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",20,0,0.5)]

Xcenter = 0
Ycenter = 0
Window = 2.
Nbins = 100
Hist_Xoff_vs_Yoff_SR = []
Hist_Xoff_vs_Yoff_CR = []
Hist_Xoff_vs_Yoff_Bkg = []
for energy in range(0,len(energy_bins)-1):
    Hist_Xoff_vs_Yoff_SR += [ROOT.TH2D("Hist_Xoff_vs_Yoff_SR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)]
    Hist_Xoff_vs_Yoff_CR += [ROOT.TH2D("Hist_Xoff_vs_Yoff_CR_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)]
    Hist_Xoff_vs_Yoff_Bkg += [ROOT.TH2D("Hist_Xoff_vs_Yoff_Bkg_ErecS%dto%d"%(energy_bins[energy],energy_bins[energy+1]),"",Nbins,Xcenter-Window,Xcenter+Window,Nbins,Ycenter-Window,Ycenter+Window)]


MSCW_cut_lower = -0.5
MSCW_cut_upper = 0.6
MSCL_cut_lower = -0.5
MSCL_cut_upper = 0.6
Depth_cut_lower = 7
Depth_cut_upper = 13

if method=='DepthUpper':
    Depth_cut_control = 13
if method=='DepthLower':
    Depth_cut_control = 7
MSCW_cut_control = 2.0
MSCL_cut_control = 2.0

def QualitySelection(tree):
    if tree.EmissionHeightChi2/tree.EmissionHeight>0.2: return False
    if tree.EChi2S/tree.ErecS>5: return False
    return True

def SignalSelection(tree):
    if not region=='VR':
        if SignalRegion(tree): return True
        return False
    else:
        if SignalRegion(tree): return False
        if ControlSelection(tree): return False
        if method=='MSCW':
            if tree.MSCW>MSCW_cut_upper*1.0: return False
            if tree.MSCL>MSCL_cut_upper*4.0: return False
        if method=='DepthUpper':
            if tree.MSCW>MSCW_cut_upper*4.0: return False
            if tree.MSCW<MSCW_cut_lower: return False
            if tree.MSCL>MSCL_cut_upper*1.0: return False
            if tree.MSCL<MSCL_cut_lower: return False
            if tree.SlantDepth*100./37.>Depth_cut_lower: return False
        return True

def SignalRegion(tree):
    if tree.MSCW>MSCW_cut_upper: return False
    if tree.MSCW<MSCW_cut_lower: return False
    if tree.MSCL>MSCL_cut_upper: return False
    if tree.MSCL<MSCL_cut_lower: return False
    if tree.SlantDepth*100./37.<Depth_cut_lower: return False
    if tree.SlantDepth*100./37.>Depth_cut_upper: return False
    return True

def ControlSelection(tree):
    if tree.MSCW>MSCW_cut_upper*3.0: return False
    if tree.MSCL>MSCL_cut_upper*3.0: return False
    if method=='DepthUpper':
        if tree.MSCW>MSCW_cut_upper*3.0: return False
        if tree.MSCL>MSCL_cut_upper*1.: return False
        if tree.MSCW<MSCW_cut_lower: return False
        if tree.MSCL<MSCL_cut_lower: return False
        if tree.SlantDepth*100./37.<Depth_cut_upper: return False
        #if tree.SlantDepth*100./37.>Depth_cut_lower and tree.SlantDepth*100./37.<Depth_cut_upper: return False
    if method == 'MSCW':
        if tree.MSCL>MSCL_cut_upper*1.25: return False
        if tree.MSCL<MSCL_cut_lower: return False
        if tree.SlantDepth*100./37.<Depth_cut_lower: return False
        if tree.SlantDepth*100./37.>Depth_cut_upper: return False
        if tree.MSCW<MSCW_cut_upper*1.5: return False
    if method == 'MSCL':
        if tree.MSCW>MSCW_cut_upper*1.5: return False
        if tree.MSCW<MSCW_cut_lower: return False
        if tree.SlantDepth*100./37.<Depth_cut_lower: return False
        if tree.SlantDepth*100./37.>Depth_cut_upper: return False
        if tree.MSCL<MSCL_cut_upper*1.5: return False
    if method=='DepthLower':
        if not UseLooseControlRegions:
            if tree.MSCW>MSCW_cut_upper: return False
            if tree.MSCW<MSCW_cut_lower: return False
            if tree.MSCL>MSCL_cut_upper: return False
            if tree.MSCL<MSCL_cut_lower: return False
        if tree.SlantDepth*100./37.>Depth_cut_control: return False
    return True

def DiagnosticTightSelection(tree):
    if method=='DepthUpper':
        if tree.MSCW>MSCW_cut_upper: return False
        if tree.MSCW<MSCW_cut_lower: return False
        if tree.MSCL>MSCL_cut_upper: return False
        if tree.MSCL<MSCL_cut_lower: return False
    if method=='DepthLower':
        if tree.MSCW>MSCW_cut_upper: return False
        if tree.MSCW<MSCW_cut_lower: return False
        if tree.MSCL>MSCL_cut_upper: return False
        if tree.MSCL<MSCL_cut_lower: return False
    if method == 'MSCW':
        if tree.MSCL>MSCL_cut_upper: return False
        if tree.MSCL<MSCL_cut_lower: return False
        if tree.SlantDepth*100./37.<Depth_cut_lower: return False
        if tree.SlantDepth*100./37.>Depth_cut_upper: return False
    if method == 'MSCL':
        if tree.MSCW>MSCW_cut_upper: return False
        if tree.MSCW<MSCW_cut_lower: return False
        if tree.SlantDepth*100./37.<Depth_cut_lower: return False
        if tree.SlantDepth*100./37.>Depth_cut_upper: return False
    return True

def DiagnosticLooseSelection(tree):
    return True

def FOVSelection(tree):
    if tree.theta2<0.3: return False
    #if tree.theta2>1.0: return False
    return True

def CalculateBkgMSCWMethod(Hist_CR,Hist_Data,Hist_CR_Raw,CR_attenu,err_CR_attenu):

    Hist_CR.Reset()
    if not CR_attenu==0:
        Hist_CR.Add(Hist_CR_Raw,1./CR_attenu)
        for b in range(0,Hist_CR.GetNbinsX()):
            old_err = Hist_CR.GetBinError(b+1)
            add_err = Hist_CR.GetBinContent(b+1)*(err_CR_attenu/CR_attenu)
            Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
        for b in range(0,Hist_CR.GetNbinsX()):
            Hist_CR.SetBinContent(b+1,max(0,Hist_CR.GetBinContent(b+1)))

def CalculateBkgHeightMethod(Hist_CR,Hist_Data,Hist_CR_Raw,CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu):

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = a*N_egamma + b*N_CR, a = gamma attenuation, b = CR attenuation
    #N_control - b*N_signal = (a-b)*N_egamma
    #N_egamma = (N_control-b*N_signal)/(a-b)
    #N_CR = N_signal-N_egamma

    #if NoSignalAssumption: gamma_attenu = 0

    print 'sig_att = %f +/- %f'%(gamma_attenu,err_gamma_attenu)
    print 'bkg_att = %f +/- %f'%(CR_attenu,err_CR_attenu)

    Hist_CR.Reset()
    if not CR_attenu==0 and not gamma_attenu/CR_attenu==1:
        Hist_CR.Add(Hist_Data)
        Hist_CR.Add(Hist_CR_Raw,-1./CR_attenu)
        for b in range(0,Hist_CR.GetNbinsX()):
            old_err = Hist_CR.GetBinError(b+1)
            add_err = Hist_CR_Raw.GetBinContent(b+1)/CR_attenu*(err_CR_attenu/CR_attenu)
            Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
        Hist_CR.Scale(1./(1.-gamma_attenu/CR_attenu))
        if not gamma_attenu==0:
            for b in range(0,Hist_CR.GetNbinsX()):
                old_err = Hist_CR.GetBinError(b+1)
                numerator = pow(err_gamma_attenu/gamma_attenu,2)+pow(err_CR_attenu/CR_attenu,2)
                numerator = pow(numerator,0.5)*gamma_attenu/CR_attenu
                add_err = numerator/(1.-gamma_attenu/CR_attenu)*(1./(1.-gamma_attenu/CR_attenu))
                Hist_CR.SetBinError(b+1,pow(old_err*old_err+add_err*add_err,0.5))
        Hist_CR.Scale(-1.)
        Hist_CR.Add(Hist_Data)
        for b in range(0,Hist_CR.GetNbinsX()):
            Hist_CR.SetBinContent(b+1,max(0,Hist_CR.GetBinContent(b+1)))

def AnalyzeTarget(target):
    runlist_Target = []
    if target=='3C264': runlist_Target = runlist_3C264
    if target=='H1426': runlist_Target = runlist_H1426
    if target=='Ton599': runlist_Target = runlist_Ton599
    if target=='IC443': runlist_Target = runlist_IC443
    if target=='PKS1424': runlist_Target = PKS1424_runlist
    if target=='2ndCrab': runlist_Target = SecondCrab_runlist
    for run in runlist_Target:
        print 'getting Target run %s'%(run)
        filename = target
        if target=='2ndCrab': filename = 'Crab'
        Target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/%s_V6_Moderate-TMVA-BDT.RB.%s.root"%(filename,run))
        pointing_tree = Target_file.Get("run_%s/stereo/pointingDataReduced"%(run))
        pointing_tree.GetEntry(0)
        TelElev_begin = pointing_tree.TelElevation
        pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
        TelElev_end = pointing_tree.TelElevation
        if (TelElev_begin+TelElev_end)/2.<Elev_cut_lower: continue
        if (TelElev_begin+TelElev_end)/2.>Elev_cut_upper: continue

        pointing_tree.GetEntry(0)
        time_begin = pointing_tree.Time
        pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
        time_end = pointing_tree.Time
        hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
        for i in range(0,pointing_tree.GetEntries()):
            pointing_tree.GetEntry(i)
            hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
            hist_Azi_vs_time.Fill(pointing_tree.Time,pointing_tree.TelAzimuth)

        Target_tree_On = Target_file.Get("run_%s/stereo/data_on"%(run))
        for entry in range(0,Target_tree_On.GetEntries()):
            Target_tree_On.GetEntry(entry)
            if not QualitySelection(Target_tree_On): continue
            if Target_tree_On.ErecS*1000.<energy_bins[0]: continue
            if Target_tree_On.ErecS*1000.>energy_bins[len(energy_bins)-1]: continue
            energy_bin = Hist_ErecS_Target_SR[0].FindBin(Target_tree_On.ErecS*1000.)-1
            if energy_bin>len(energy_bins)-2: continue
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(Target_tree_On.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(Target_tree_On.Time))
            Hist2D_TelElev_vs_SlantDepth_Target.Fill(elevation,Target_tree_On.SlantDepth*100./37.)
            Hist2D_TelElev_vs_EmissionHeight_Target.Fill(elevation,Target_tree_On.EmissionHeight)
            Hist_TelElevAzim_Counts_Target.Fill(elevation,azimuth)
            if SignalSelection(Target_tree_On):
                if FOVSelection(Target_tree_On):
                    Hist_ErecS_Target_SR[energy_bin].Fill(Target_tree_On.ErecS*1000.)
                Hist_Theta2_Target_SR[energy_bin].Fill(Target_tree_On.theta2)
                Hist_Theta2ZoomIn_Target_SR[energy_bin].Fill(Target_tree_On.theta2)
                Hist_Xoff_vs_Yoff_SR[energy_bin].Fill(Target_tree_On.ra,Target_tree_On.dec)
            if ControlSelection(Target_tree_On):
                if FOVSelection(Target_tree_On):
                    Hist_ErecS_Target_CR[energy_bin].Fill(Target_tree_On.ErecS*1000.)
                Hist_Theta2_Target_CR[energy_bin].Fill(Target_tree_On.theta2)
                Hist_Theta2ZoomIn_Target_CR[energy_bin].Fill(Target_tree_On.theta2)
                Hist_Xoff_vs_Yoff_CR[energy_bin].Fill(Target_tree_On.ra,Target_tree_On.dec)
            if DiagnosticLooseSelection(Target_tree_On):
                if FOVSelection(Target_tree_On):
                    Hist_Loose_SlantDepth_Target_On.Fill(Target_tree_On.SlantDepth*100./37.)
                    Hist_Loose_MSCW_Target_On.Fill(Target_tree_On.MSCW)
                    Hist_Loose_MSCL_Target_On.Fill(Target_tree_On.MSCL)
                    if Target_tree_On.ErecS*1000.>1000.:
                        Hist_Loose_SlantDepth_Target_On_TeV.Fill(Target_tree_On.SlantDepth*100./37.)
                        Hist_Loose_MSCW_Target_On_TeV.Fill(Target_tree_On.MSCW)
                        Hist_Loose_MSCL_Target_On_TeV.Fill(Target_tree_On.MSCL)
            if DiagnosticTightSelection(Target_tree_On):
                Hist2D_Energy_vs_SlantDepth_Target.Fill(Target_tree_On.ErecS*1000.,Target_tree_On.SlantDepth*100./37.)
                if FOVSelection(Target_tree_On):
                    Hist_Tight_SlantDepth_Target_On.Fill(Target_tree_On.SlantDepth*100./37.)
                    Hist_Tight_MSCW_Target_On.Fill(Target_tree_On.MSCW)
                    Hist_Tight_MSCL_Target_On.Fill(Target_tree_On.MSCL)
                    if Target_tree_On.ErecS*1000.>1000.:
                        Hist_Tight_SlantDepth_Target_On_TeV.Fill(Target_tree_On.SlantDepth*100./37.)
                        Hist_Tight_MSCW_Target_On_TeV.Fill(Target_tree_On.MSCW)
                        Hist_Tight_MSCL_Target_On_TeV.Fill(Target_tree_On.MSCL)
        Target_file.Close()


for run in Crab_runlist:
    print 'getting Crab run %s'%(run)
    Crab_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    pointing_tree = Crab_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    pointing_tree.GetEntry(0)
    TelElev_begin = pointing_tree.TelElevation
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    TelElev_end = pointing_tree.TelElevation
    if (TelElev_begin+TelElev_end)/2.<Elev_cut_lower: continue
    if (TelElev_begin+TelElev_end)/2.>Elev_cut_upper: continue

    pointing_tree.GetEntry(0)
    time_begin = pointing_tree.Time
    pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
    time_end = pointing_tree.Time
    hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,pointing_tree.TelAzimuth)

    Crab_tree_On = Crab_file.Get("run_%s/stereo/data_on"%(run))
    for entry in range(0,Crab_tree_On.GetEntries()):
        Crab_tree_On.GetEntry(entry)
        if not QualitySelection(Crab_tree_On): continue
        if Crab_tree_On.ErecS*1000.<energy_bins[0]: continue
        if Crab_tree_On.ErecS*1000.>energy_bins[len(energy_bins)-1]: continue
        energy_bin = Hist_ErecS_Target_SR[0].FindBin(Crab_tree_On.ErecS*1000.)-1
        if energy_bin>len(energy_bins)-2: continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(Crab_tree_On.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(Crab_tree_On.Time))
        if SignalSelection(Crab_tree_On):
            Hist_On_SR_Count[energy_bin].Fill(0)
        if ControlSelection(Crab_tree_On):
            Hist_On_CR_Count[energy_bin].Fill(0)
        Hist_MSCW_Crab_On_Alpha.Fill(Crab_tree_On.MSCW)
        if DiagnosticLooseSelection(Crab_tree_On):
            Hist_Loose_SlantDepth_Crab_On.Fill(Crab_tree_On.SlantDepth*100./37.)
            Hist_Loose_MSCW_Crab_On.Fill(Crab_tree_On.MSCW)
            Hist_Loose_MSCL_Crab_On.Fill(Crab_tree_On.MSCL)
            if Crab_tree_On.ErecS*1000.>1000.:
                Hist_Loose_SlantDepth_Crab_On_TeV.Fill(Crab_tree_On.SlantDepth*100./37.)
                Hist_Loose_MSCW_Crab_On_TeV.Fill(Crab_tree_On.MSCW)
                Hist_Loose_MSCL_Crab_On_TeV.Fill(Crab_tree_On.MSCL)
        if DiagnosticTightSelection(Crab_tree_On):
            Hist2D_Energy_vs_SlantDepth_Crab_On.Fill(Crab_tree_On.ErecS*1000.,Crab_tree_On.SlantDepth*100./37.)
            Hist_Tight_SlantDepth_Crab_On.Fill(Crab_tree_On.SlantDepth*100./37.)
            Hist_Tight_MSCW_Crab_On.Fill(Crab_tree_On.MSCW)
            Hist_Tight_MSCL_Crab_On.Fill(Crab_tree_On.MSCL)
            if Crab_tree_On.ErecS*1000.>1000.:
                Hist_Tight_SlantDepth_Crab_On_TeV.Fill(Crab_tree_On.SlantDepth*100./37.)
                Hist_Tight_MSCW_Crab_On_TeV.Fill(Crab_tree_On.MSCW)
                Hist_Tight_MSCL_Crab_On_TeV.Fill(Crab_tree_On.MSCL)

    Crab_tree_Off = Crab_file.Get("run_%s/stereo/data_off"%(run))
    for entry in range(0,Crab_tree_Off.GetEntries()):
        Crab_tree_Off.GetEntry(entry)
        if not QualitySelection(Crab_tree_Off): continue
        if Crab_tree_Off.ErecS*1000.<energy_bins[0]: continue
        if Crab_tree_Off.ErecS*1000.>energy_bins[len(energy_bins)-1]: continue
        energy_bin = Hist_ErecS_Target_SR[0].FindBin(Crab_tree_Off.ErecS*1000.)-1
        if energy_bin>len(energy_bins)-2: continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(Crab_tree_Off.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(Crab_tree_Off.Time))
        Hist_TelElevAzim_Counts_Crab.Fill(elevation,azimuth)
        if SignalSelection(Crab_tree_Off):
            Hist_Off_SR_Count[energy_bin].Fill(0)
        if ControlSelection(Crab_tree_Off):
            Hist_Off_CR_Count[energy_bin].Fill(0)
        Hist_MSCW_Crab_Off_Alpha.Fill(Crab_tree_Off.MSCW)
        if DiagnosticLooseSelection(Crab_tree_Off):
            Hist_Loose_SlantDepth_Crab_Off.Fill(Crab_tree_Off.SlantDepth*100./37.)
            Hist_Loose_MSCW_Crab_Off.Fill(Crab_tree_Off.MSCW)
            Hist_Loose_MSCL_Crab_Off.Fill(Crab_tree_Off.MSCL)
            if Crab_tree_Off.ErecS*1000.>1000.:
                Hist_Loose_SlantDepth_Crab_Off_TeV.Fill(Crab_tree_Off.SlantDepth*100./37.)
                Hist_Loose_MSCW_Crab_Off_TeV.Fill(Crab_tree_Off.MSCW)
                Hist_Loose_MSCL_Crab_Off_TeV.Fill(Crab_tree_Off.MSCL)
        if DiagnosticTightSelection(Crab_tree_Off):
            Hist2D_Energy_vs_SlantDepth_Crab_Off.Fill(Crab_tree_Off.ErecS*1000.,Crab_tree_Off.SlantDepth*100./37.)
            Hist_Tight_SlantDepth_Crab_Off.Fill(Crab_tree_Off.SlantDepth*100./37.)
            Hist_Tight_MSCW_Crab_Off.Fill(Crab_tree_Off.MSCW)
            Hist_Tight_MSCL_Crab_Off.Fill(Crab_tree_Off.MSCL)
            if Crab_tree_Off.ErecS*1000.>1000.:
                Hist_Tight_SlantDepth_Crab_Off_TeV.Fill(Crab_tree_Off.SlantDepth*100./37.)
                Hist_Tight_MSCW_Crab_Off_TeV.Fill(Crab_tree_Off.MSCW)
                Hist_Tight_MSCL_Crab_Off_TeV.Fill(Crab_tree_Off.MSCL)

    Crab_file.Close()

AnalyzeTarget(target)

bin_begin = Hist_Loose_MSCW_Crab_On.FindBin(4)
bin_end = Hist_Loose_MSCW_Crab_On.FindBin(10)
alpha = 0
if not Hist_MSCW_Crab_Off_Alpha.Integral(bin_begin,bin_end)==0:
    alpha = Hist_MSCW_Crab_On_Alpha.Integral(bin_begin,bin_end)/Hist_MSCW_Crab_Off_Alpha.Integral(bin_begin,bin_end)
Hist_Loose_SlantDepth_Crab_Off.Scale(alpha)
Hist_Loose_MSCW_Crab_Off.Scale(alpha)
Hist_Loose_MSCL_Crab_Off.Scale(alpha)
Hist_Loose_SlantDepth_Crab_Off_TeV.Scale(alpha)
Hist_Loose_MSCW_Crab_Off_TeV.Scale(alpha)
Hist_Loose_MSCL_Crab_Off_TeV.Scale(alpha)
Hist_Tight_SlantDepth_Crab_Off.Scale(alpha)
Hist_Tight_MSCW_Crab_Off.Scale(alpha)
Hist_Tight_MSCL_Crab_Off.Scale(alpha)
Hist_Tight_SlantDepth_Crab_Off_TeV.Scale(alpha)
Hist_Tight_MSCW_Crab_Off_TeV.Scale(alpha)
Hist_Tight_MSCL_Crab_Off_TeV.Scale(alpha)

Hist2D_Energy_vs_SlantDepth_Crab_Off.Scale(alpha)
Hist2D_Energy_vs_SlantDepth_Crab_On.Add(Hist2D_Energy_vs_SlantDepth_Crab_Off,-1.)

for energy in range(0,len(energy_bins)-1):

    Hist_Off_SR_Count[energy].Scale(alpha)
    Hist_Off_CR_Count[energy].Scale(alpha)
    Hist_On_SR_Count[energy].Add(Hist_Off_SR_Count[energy],-1)
    Hist_On_CR_Count[energy].Add(Hist_Off_CR_Count[energy],-1)
    Hist_SigAtt[energy].Reset()
    Hist_SigAtt[energy].Add(Hist_On_CR_Count[energy])
    Hist_SigAtt[energy].Divide(Hist_On_SR_Count[energy])
    Hist_BkgAtt[energy].Reset()
    Hist_BkgAtt[energy].Add(Hist_Off_CR_Count[energy])
    Hist_BkgAtt[energy].Divide(Hist_Off_SR_Count[energy])
    sig_att = Hist_SigAtt[energy].GetBinContent(1)
    err_sig_att = Hist_SigAtt[energy].GetBinError(1)
    bkg_att = Hist_BkgAtt[energy].GetBinContent(1)
    err_bkg_att = Hist_BkgAtt[energy].GetBinError(1)

    if method == 'DepthUpper' or method == 'DepthLower':
        CalculateBkgHeightMethod(Hist_ErecS_Target_Bkg[energy],Hist_ErecS_Target_SR[energy],Hist_ErecS_Target_CR[energy],bkg_att,err_bkg_att,sig_att,err_sig_att)
        CalculateBkgHeightMethod(Hist_Theta2_Target_Bkg[energy],Hist_Theta2_Target_SR[energy],Hist_Theta2_Target_CR[energy],bkg_att,err_bkg_att,sig_att,err_sig_att)
        CalculateBkgHeightMethod(Hist_Theta2ZoomIn_Target_Bkg[energy],Hist_Theta2ZoomIn_Target_SR[energy],Hist_Theta2ZoomIn_Target_CR[energy],bkg_att,err_bkg_att,sig_att,err_sig_att)
        CalculateBkgHeightMethod(Hist_Xoff_vs_Yoff_Bkg[energy],Hist_Xoff_vs_Yoff_SR[energy],Hist_Xoff_vs_Yoff_CR[energy],bkg_att,err_bkg_att,sig_att,err_sig_att)
    if method == 'MSCW' or method == 'MSCL':
        CalculateBkgMSCWMethod(Hist_ErecS_Target_Bkg[energy],Hist_ErecS_Target_SR[energy],Hist_ErecS_Target_CR[energy],bkg_att,err_bkg_att)
        CalculateBkgMSCWMethod(Hist_Theta2_Target_Bkg[energy],Hist_Theta2_Target_SR[energy],Hist_Theta2_Target_CR[energy],bkg_att,err_bkg_att)
        CalculateBkgMSCWMethod(Hist_Theta2ZoomIn_Target_Bkg[energy],Hist_Theta2ZoomIn_Target_SR[energy],Hist_Theta2ZoomIn_Target_CR[energy],bkg_att,err_bkg_att)
        CalculateBkgMSCWMethod(Hist_Xoff_vs_Yoff_Bkg[energy],Hist_Xoff_vs_Yoff_SR[energy],Hist_Xoff_vs_Yoff_CR[energy],bkg_att,err_bkg_att)

if method=='DepthUpper':
    bin_begin = Hist_Loose_SlantDepth_Crab_On.FindBin(Depth_cut_control)
    bin_end = Hist_Loose_SlantDepth_Crab_On.FindBin(1000.)
if method=='DepthLower':
    bin_begin = Hist_Loose_SlantDepth_Crab_On.FindBin(0)
    bin_end = Hist_Loose_SlantDepth_Crab_On.FindBin(Depth_cut_control)
scale = 0
if not Hist_Loose_SlantDepth_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_SlantDepth_Crab_Off.Integral(bin_begin,bin_end)/Hist_Loose_SlantDepth_Target_On.Integral(bin_begin,bin_end)
Hist_Loose_SlantDepth_Target_On.Scale(scale)

if method=='DepthUpper':
    bin_begin = Hist_Tight_SlantDepth_Crab_On.FindBin(Depth_cut_control)
    bin_end = Hist_Tight_SlantDepth_Crab_On.FindBin(1000.)
if method=='DepthLower':
    bin_begin = Hist_Tight_SlantDepth_Crab_On.FindBin(0)
    bin_end = Hist_Tight_SlantDepth_Crab_On.FindBin(Depth_cut_control)
scale = 0
if not Hist_Tight_SlantDepth_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_SlantDepth_Crab_Off.Integral(bin_begin,bin_end)/Hist_Tight_SlantDepth_Target_On.Integral(bin_begin,bin_end)
Hist_Tight_SlantDepth_Target_On.Scale(scale)

bin_begin = Hist_Loose_MSCW_Crab_On.FindBin(MSCW_cut_control)
bin_end = Hist_Loose_MSCW_Crab_On.FindBin(100)
scale = 0
if not Hist_Loose_MSCW_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_MSCW_Crab_Off.Integral(bin_begin,bin_end)/Hist_Loose_MSCW_Target_On.Integral(bin_begin,bin_end)
Hist_Loose_MSCW_Target_On.Scale(scale)

bin_begin = Hist_Tight_MSCW_Crab_On.FindBin(MSCW_cut_control)
bin_end = Hist_Tight_MSCW_Crab_On.FindBin(100)
scale = 0
if not Hist_Tight_MSCW_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_MSCW_Crab_Off.Integral(bin_begin,bin_end)/Hist_Tight_MSCW_Target_On.Integral(bin_begin,bin_end)
Hist_Tight_MSCW_Target_On.Scale(scale)

bin_begin = Hist_Loose_MSCL_Crab_On.FindBin(MSCL_cut_control)
bin_end = Hist_Loose_MSCL_Crab_On.FindBin(100)
scale = 0
if not Hist_Loose_MSCL_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_MSCL_Crab_Off.Integral(bin_begin,bin_end)/Hist_Loose_MSCL_Target_On.Integral(bin_begin,bin_end)
Hist_Loose_MSCL_Target_On.Scale(scale)

bin_begin = Hist_Tight_MSCL_Crab_On.FindBin(MSCL_cut_control)
bin_end = Hist_Tight_MSCL_Crab_On.FindBin(100)
scale = 0
if not Hist_Tight_MSCL_Target_On.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_MSCL_Crab_Off.Integral(bin_begin,bin_end)/Hist_Tight_MSCL_Target_On.Integral(bin_begin,bin_end)
Hist_Tight_MSCL_Target_On.Scale(scale)

if method=='DepthUpper':
    bin_begin = Hist_Loose_SlantDepth_Crab_On_TeV.FindBin(Depth_cut_control)
    bin_end = Hist_Loose_SlantDepth_Crab_On_TeV.FindBin(1000.)
if method=='DepthLower':
    bin_begin = Hist_Loose_SlantDepth_Crab_On_TeV.FindBin(0)
    bin_end = Hist_Loose_SlantDepth_Crab_On_TeV.FindBin(Depth_cut_control)
scale = 0
if not Hist_Loose_SlantDepth_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_SlantDepth_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Loose_SlantDepth_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Loose_SlantDepth_Target_On_TeV.Scale(scale)

if method=='DepthUpper':
    bin_begin = Hist_Tight_SlantDepth_Crab_On_TeV.FindBin(Depth_cut_control)
    bin_end = Hist_Tight_SlantDepth_Crab_On_TeV.FindBin(1000.)
if method=='DepthLower':
    bin_begin = Hist_Tight_SlantDepth_Crab_On_TeV.FindBin(0)
    bin_end = Hist_Tight_SlantDepth_Crab_On_TeV.FindBin(Depth_cut_control)
scale = 0
if not Hist_Tight_SlantDepth_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_SlantDepth_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Tight_SlantDepth_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Tight_SlantDepth_Target_On_TeV.Scale(scale)

bin_begin = Hist_Loose_MSCW_Crab_On_TeV.FindBin(MSCW_cut_control)
bin_end = Hist_Loose_MSCW_Crab_On_TeV.FindBin(100)
scale = 0
if not Hist_Loose_MSCW_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_MSCW_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Loose_MSCW_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Loose_MSCW_Target_On_TeV.Scale(scale)

bin_begin = Hist_Tight_MSCW_Crab_On_TeV.FindBin(MSCW_cut_control)
bin_end = Hist_Tight_MSCW_Crab_On_TeV.FindBin(100)
scale = 0
if not Hist_Tight_MSCW_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_MSCW_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Tight_MSCW_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Tight_MSCW_Target_On_TeV.Scale(scale)

bin_begin = Hist_Loose_MSCL_Crab_On_TeV.FindBin(MSCL_cut_control)
bin_end = Hist_Loose_MSCL_Crab_On_TeV.FindBin(100)
scale = 0
if not Hist_Loose_MSCL_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Loose_MSCL_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Loose_MSCL_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Loose_MSCL_Target_On_TeV.Scale(scale)

bin_begin = Hist_Tight_MSCL_Crab_On_TeV.FindBin(MSCL_cut_control)
bin_end = Hist_Tight_MSCL_Crab_On_TeV.FindBin(100)
scale = 0
if not Hist_Tight_MSCL_Target_On_TeV.Integral(bin_begin,bin_end)==0:
    scale = Hist_Tight_MSCL_Crab_Off_TeV.Integral(bin_begin,bin_end)/Hist_Tight_MSCL_Target_On_TeV.Integral(bin_begin,bin_end)
Hist_Tight_MSCL_Target_On_TeV.Scale(scale)


if not UseLooseControlRegions:
    OutputFile = ROOT.TFile("output/Histograms_%s_%s_%s_Elev%dto%d_Tight.root"%(target,method,region,Elev_cut_lower,Elev_cut_upper),"recreate")
else:
    OutputFile = ROOT.TFile("output/Histograms_%s_%s_%s_Elev%dto%d.root"%(target,method,region,Elev_cut_lower,Elev_cut_upper),"recreate")
for energy in range(0,len(energy_bins)-1):
    Hist_ErecS_Target_SR[energy].Write()
    Hist_ErecS_Target_Bkg[energy].Write()
    Hist_Theta2_Target_SR[energy].Write()
    Hist_Theta2_Target_Bkg[energy].Write()
    Hist_Theta2ZoomIn_Target_SR[energy].Write()
    Hist_Theta2ZoomIn_Target_Bkg[energy].Write()
    Hist_Xoff_vs_Yoff_SR[energy].Write()
    Hist_Xoff_vs_Yoff_Bkg[energy].Write()
Hist_Loose_MSCW_Crab_On.Write()
Hist_Loose_MSCW_Crab_Off.Write()
Hist_Loose_MSCW_Target_On.Write()
Hist_Loose_MSCL_Crab_On.Write()
Hist_Loose_MSCL_Crab_Off.Write()
Hist_Loose_MSCL_Target_On.Write()
Hist_Loose_MSCW_Crab_On_TeV.Write()
Hist_Loose_MSCW_Crab_Off_TeV.Write()
Hist_Loose_MSCW_Target_On_TeV.Write()
Hist_Loose_MSCL_Crab_On_TeV.Write()
Hist_Loose_MSCL_Crab_Off_TeV.Write()
Hist_Loose_MSCL_Target_On_TeV.Write()
Hist_Loose_SlantDepth_Crab_On.Write()
Hist_Loose_SlantDepth_Crab_Off.Write()
Hist_Loose_SlantDepth_Target_On.Write()
Hist_Loose_SlantDepth_Crab_On_TeV.Write()
Hist_Loose_SlantDepth_Crab_Off_TeV.Write()
Hist_Loose_SlantDepth_Target_On_TeV.Write()
Hist_Tight_MSCW_Crab_On.Write()
Hist_Tight_MSCW_Crab_Off.Write()
Hist_Tight_MSCW_Target_On.Write()
Hist_Tight_MSCL_Crab_On.Write()
Hist_Tight_MSCL_Crab_Off.Write()
Hist_Tight_MSCL_Target_On.Write()
Hist_Tight_MSCW_Crab_On_TeV.Write()
Hist_Tight_MSCW_Crab_Off_TeV.Write()
Hist_Tight_MSCW_Target_On_TeV.Write()
Hist_Tight_MSCL_Crab_On_TeV.Write()
Hist_Tight_MSCL_Crab_Off_TeV.Write()
Hist_Tight_MSCL_Target_On_TeV.Write()
Hist_Tight_SlantDepth_Crab_On.Write()
Hist_Tight_SlantDepth_Crab_Off.Write()
Hist_Tight_SlantDepth_Target_On.Write()
Hist_Tight_SlantDepth_Crab_On_TeV.Write()
Hist_Tight_SlantDepth_Crab_Off_TeV.Write()
Hist_Tight_SlantDepth_Target_On_TeV.Write()
Hist2D_TelElev_vs_SlantDepth_Target.Write()
Hist2D_TelElev_vs_EmissionHeight_Target.Write()
Hist2D_Energy_vs_SlantDepth_Target.Write()
Hist2D_Energy_vs_SlantDepth_Crab_On.Write()
Hist2D_Energy_vs_SlantDepth_Crab_Off.Write()
Hist_TelElevAzim_Counts_Target.Write()
Hist_TelElevAzim_Counts_Crab.Write()
OutputFile.Close()

