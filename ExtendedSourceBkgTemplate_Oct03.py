
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
ROOT.gStyle.SetPaintTextFormat("0.3f")

Hist_TelElev_Counts = TH1D("Hist_TelElev_Counts","",45,0,90)
Hist_TelElev_Counts_target = TH1D("Hist_TelElev_Counts_target","",18,0,90)
Hist_TelElev_Counts_source = TH1D("Hist_TelElev_Counts_source","",18,0,90)
Hist_TelElev_Counts_scale = TH1D("Hist_TelElev_Counts_scale","",18,0,90)

def ControlRegionHeightSelection(tree,elev):
    
    cut_upper = 10.+control_width
    cut_lower = 10.-control_width

    if (tree.EmissionHeight<cut_upper and tree.EmissionHeight>cut_lower): return False
    return True

def SignalRegionHeightSelection(tree,elev):
    
    cut_upper = 10.+signal_width
    cut_lower = 10.-signal_width

    if (tree.EmissionHeight>cut_upper or tree.EmissionHeight<cut_lower): return False
    return True

def CommonSelection(tree):
    if (tree.MSCL<MSCL_cut_lower): return False
    if (tree.MSCL>MSCL_cut_upper): return False
    if (tree.Erec<Erec_cut_lower): return False
    if (tree.Erec>Erec_cut_upper): return False
    return True

def MakeATag():
    global tag
    tag = ''
    tag += '_Elev%sto%s'%(Elev_cut_lower,Elev_cut_upper)
    tag += '_MSCL%sto%s'%(MSCL_cut_lower,MSCL_cut_upper)
    tag += '_Erec%sto%s'%(Erec_cut_lower,Erec_cut_upper)
    tag += '_'+Height_cut
    tag += '_'+RW_method
    tag += '_'+Fake_Removal
    tag += '_'+source+'_'+field


def GetAnttenuationHistogram(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if region == 'on':
            anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
            ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        else:
            anasum_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("data")

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

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        if not CommonSelection(anasum_tree): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        scale = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        if RW_method=='NoRW': 
            scale = 1
        if ControlRegionHeightSelection(anasum_tree,elevation):
            if region == 'on':
                Hist_OnData_Control.Fill(anasum_tree.MSCW,scale)
            else:
                Hist_OffData_Control.Fill(anasum_tree.MSCW,scale)
        if SignalRegionHeightSelection(anasum_tree,elevation):
            if region == 'on':
                Hist_OnData_Signal.Fill(anasum_tree.MSCW,scale)
            else:
                Hist_OffData_Signal.Fill(anasum_tree.MSCW,scale)

def GetSourceElevation(run):

    if source == 'Crab':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == '2ndCrab':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == 'PKS1424':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == '3C264':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))

    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        Hist_TelElev_Counts_source.Fill(pointing_tree.TelElevation)

def GetTargetElevation(run):

    if target == 'Crab':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == '2ndCrab':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'PKS1424':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == '3C264':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    target_pointing_tree = target_file.Get("run_%s/stereo/pointingDataReduced"%(run))

    for i in range(0,target_pointing_tree.GetEntries()):
        target_pointing_tree.GetEntry(i)
        Hist_TelElev_Counts_target.Fill(target_pointing_tree.TelElevation)


def MeasureAttenuation():

    Rate_CR_Attenuation = 0
    Err_Rate_CR_Attenuation = 0
    Rate_Gamma_Attenuation = 0
    Err_Rate_Gamma_Attenuation = 0
    
    if source == 'Crab':
        for run in Crab_runlist:
            GetSourceElevation(run)
    if target == 'Crab':
        for run in Crab_runlist:
            GetTargetElevation(run)
    if target == '2ndCrab':
        for run in SecondCrab_runlist:
            GetTargetElevation(run)
    if target == 'PKS1424':
        for run in PKS1424_runlist:
            GetTargetElevation(run)
    if target == '3C264':
        for run in runlist_3C264:
            GetTargetElevation(run)
    #print Hist_TelElev_Counts_target.Print('All')
    for b in range(0,Hist_TelElev_Counts_scale.GetNbinsX()):
        Hist_TelElev_Counts_scale.SetBinContent(b+1,Hist_TelElev_Counts_target.GetBinContent(b+1))
    Hist_TelElev_Counts_scale.Divide(Hist_TelElev_Counts_source)
    #print Hist_TelElev_Counts_scale.Print('All')
    if source == 'Crab':
        for run in Crab_runlist:
            GetAnttenuationHistogram('Crab','on',run,'anasum')
            GetAnttenuationHistogram('Crab','off',run,'anasum')
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            GetAnttenuationHistogram('Crab','on',run,'anasum')
            GetAnttenuationHistogram('Crab','off',run,'anasum')
    Hist_TelElev_Counts_scale.Reset()
    Hist_TelElev_Counts_source.Reset()
    Hist_TelElev_Counts_target.Reset()
    
    print 'Hist_OnData_Signal.GetEntries() = %s'%(Hist_OnData_Signal.GetEntries())
    print 'Hist_OnData_Control.GetEntries() = %s'%(Hist_OnData_Control.GetEntries())
    print 'Hist_OffData_Signal.GetEntries() = %s'%(Hist_OffData_Signal.GetEntries())
    print 'Hist_OffData_Control.GetEntries() = %s'%(Hist_OffData_Control.GetEntries())
    print 'Hist_OnData_Signal.Integral() = %s'%(Hist_OnData_Signal.GetBinContent(Hist_OnData_Signal.FindBin(0)))
    print 'Hist_OnData_Control.Integral() = %s'%(Hist_OnData_Control.GetBinContent(Hist_OnData_Signal.FindBin(0)))
    print 'Hist_OffData_Signal.Integral() = %s'%(Hist_OffData_Signal.GetBinContent(Hist_OnData_Signal.FindBin(0)))

    if Hist_OffData_Signal.GetBinContent(len(MSCW_Bin)-1)<10:
        return 0,0,0,0,0,0
    if Hist_OffData_Control.GetBinContent(len(MSCW_Bin)-1)<10:
        return 0,0,0,0,0,0
    Normalization_Factor_Signal = Hist_OnData_Signal.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Signal = Normalization_Factor_Signal/Hist_OffData_Signal.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Control = Hist_OnData_Control.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Control = Normalization_Factor_Control/Hist_OffData_Control.GetBinContent(len(MSCW_Bin)-1)
    Hist_OffData_Signal.Scale(Normalization_Factor_Signal)
    Hist_OffData_Control.Scale(Normalization_Factor_Signal)
    
    for b in range(0,Hist_CR_Attenuation.GetNbinsX()):
        Hist_CR_Attenuation.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_CR_Attenuation.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_CR_Attenuation.Divide(Hist_OffData_Signal)
    Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinContent(Hist_CR_Attenuation.FindBin(0))
    Err_Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinError(Hist_CR_Attenuation.FindBin(0))
    print 'Rate_CR_Attenuation = %s +/- %s'%(Rate_CR_Attenuation,Err_Rate_CR_Attenuation)
    
    for b in range(0,Hist_Gamma_Signal.GetNbinsX()):
        Hist_Gamma_Signal.SetBinContent(b+1,Hist_OnData_Signal.GetBinContent(b+1))
        Hist_Gamma_Signal.SetBinError(b+1,Hist_OnData_Signal.GetBinError(b+1))
    Hist_Gamma_Signal.Add(Hist_OffData_Signal,-1)
    for b in range(0,Hist_Gamma_Control.GetNbinsX()):
        Hist_Gamma_Control.SetBinContent(b+1,Hist_OnData_Control.GetBinContent(b+1))
        Hist_Gamma_Control.SetBinError(b+1,Hist_OnData_Control.GetBinError(b+1))
    Hist_Gamma_Control.Add(Hist_OffData_Control,-1)
    for b in range(0,Hist_Gamma_Attenuation.GetNbinsX()):
        Hist_Gamma_Attenuation.SetBinContent(b+1,Hist_Gamma_Control.GetBinContent(b+1))
        Hist_Gamma_Attenuation.SetBinError(b+1,Hist_Gamma_Control.GetBinError(b+1))
    Hist_Gamma_Attenuation.Divide(Hist_Gamma_Signal)
    Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinContent(Hist_Gamma_Attenuation.FindBin(0))
    Err_Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinError(Hist_Gamma_Attenuation.FindBin(0))
    print 'Rate_Gamma_Attenuation = %s +/- %s'%(Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation)

    Hist_OnData_Signal.Reset()
    Hist_OnData_Control.Reset()
    Hist_OffData_Signal.Reset()
    Hist_OffData_Control.Reset()
    Hist_Gamma_Signal.Reset()
    Hist_Gamma_Control.Reset()
    Hist_CR_Attenuation.Reset()
    Hist_Gamma_Attenuation.Reset()

    return Rate_CR_Attenuation,Err_Rate_CR_Attenuation,Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation

def SingleRunAnalysis(source,region,run,data_type):

    global N_GammaMC_Fake
    global N_GammaData_Fake
    global N_GammaMC_ES
    global N_GammaData_ES

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if region == 'on':
            anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
            ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        else:
            anasum_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("data")

    if not source=="MC":
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

        for i in range(0,anasum_tree.GetEntries()):
            anasum_tree.GetEntry(i)
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
            Hist_TelElev_Counts.Fill(elevation)
            if not CommonSelection(anasum_tree): continue
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            scale = 1
            if ControlRegionHeightSelection(anasum_tree,elevation):
                if (anasum_tree.MSCW<1.):
                    Hist_Erec_CR.Fill(anasum_tree.Erec,scale)
                    Hist_Erec_CR_Raw.Fill(anasum_tree.Erec,1.)
            if SignalRegionHeightSelection(anasum_tree,elevation):
                if (anasum_tree.MSCW<1.):
                    Hist_Erec_Data.Fill(anasum_tree.Erec)
                if (anasum_tree.MSCW>3.):
                    Hist_Norm_Data.Fill(0)

    if region == 'on':
        for i in range(0,ring_tree.GetEntries()):
            ring_tree.GetEntry(i)
            if not CommonSelection(ring_tree): continue
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(ring_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(ring_tree.Time))
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            if SignalRegionHeightSelection(ring_tree,elevation):
                if (ring_tree.MSCW<1.):
                    Hist_Erec_Ring.Fill(ring_tree.Erec)
                if (ring_tree.MSCW>3.):
                    Hist_Norm_Ring.Fill(0)

def RunExtendedSourceAnalysis(CR_attenu,gamma_attenu):

    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print 'Erec_cut_lower = %s, Erec_cut_upper = %s'%(Erec_cut_lower,Erec_cut_upper)
    print 'CR_attenu = %s, gamma_attenu = %s'%(CR_attenu,gamma_attenu)

    if source == 'Crab':
        for run in Crab_runlist:
            SingleRunAnalysis('Crab',field,run,'anasum')
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            SingleRunAnalysis('Crab',field,run,'anasum')
    if source == 'PKS1424':
        for run in PKS1424_runlist:
            SingleRunAnalysis('PKS1424',field,run,'anasum')
    if source == '3C264':
        for run in runlist_3C264:
            SingleRunAnalysis('3C264',field,run,'anasum')
    
    if field == 'on':
        Normalization_Factor_Ring = Hist_Norm_Data.Integral()
        Normalization_Factor_Ring = Normalization_Factor_Ring/Hist_Norm_Ring.Integral()
        Hist_Erec_Ring.Scale(Normalization_Factor_Ring)

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = a*N_egamma + b*N_CR, a = gamma attenuation, b = CR attenuation
    #N_control - b*N_signal = (a-b)*N_egamma
    #N_egamma = (N_control-b*N_signal)/(a-b)
    #N_CR = N_signal-N_egamma

    print 'Hist_Erec_CR_Raw.Integral() = %s'%(Hist_Erec_CR_Raw.Integral())
    print 'Hist_Erec_Data.Integral() = %s'%(Hist_Erec_Data.Integral())
    if not Fake_Removal=='NoFakeRemoved':
        Hist_Erec_CR.Reset()
        Hist_Erec_CR.Add(Hist_Erec_Data)
        Hist_Erec_CR.Add(Hist_Erec_CR_Raw,-1./CR_attenu)
        Hist_Erec_CR.Scale(1./(1.-gamma_attenu/CR_attenu))
        Hist_Erec_CR.Scale(-1.)
        Hist_Erec_CR.Add(Hist_Erec_Data)
        #Hist_Erec_CR.Reset()
        #Hist_Erec_CR.Add(Hist_Erec_CR_Raw)
        #Hist_Erec_CR.Add(Hist_Erec_Data,-CR_attenu)
        #Hist_Erec_CR.Scale(1./(CR_attenu-gamma_attenu))
        #Hist_Erec_CR.Add(Hist_Erec_Data)
    Err_N_CR = ROOT.Double(0.)
    N_CR = Hist_Erec_CR.IntegralAndError(1,Hist_Erec_CR.GetNbinsX(),Err_N_CR)
    print 'Hist_Erec_CR.Integral() = (CR_raw - gamma_att*data)/(CR_att-gamma_att) = %s +/- %s'%(N_CR,Err_N_CR)
    print 'Hist_Erec_Ring.Integral() = %s'%(Hist_Erec_Ring.Integral())
    
    
    Hist_Erec_Data_Sum.Add(Hist_Erec_Data)
    Hist_Erec_Ring_Sum.Add(Hist_Erec_Ring)
    Hist_Erec_CR_Sum.Add(Hist_Erec_CR)

    Hist_Norm_Data.Reset()
    Hist_Norm_Ring.Reset()
    Hist_Erec_CR.Reset()
    Hist_Erec_CR_Raw.Reset()
    Hist_Erec_Data.Reset()
    Hist_Erec_Ring.Reset()

Hist_Norm_Data = TH1D("Hist_Norm_Data","",1,0,1)
Hist_Norm_Ring = TH1D("Hist_Norm_Ring","",1,0,1)

Elev_Bin = [55,85]
#Elev_Bin = [30,40,50,60,70,80,90]
Hist_Elev_Bins = TH1D("Hist_Elev_Bins","",len(Elev_Bin)-1,array('d',Elev_Bin))

MSCW_Bin = [-100,1,3,100]  # the last bin has to be the CR normalization region!!
Hist_OnData_Signal = TH1D("Hist_OnData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OnData_Control = TH1D("Hist_OnData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Signal = TH1D("Hist_OffData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Control = TH1D("Hist_OffData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Signal = TH1D("Hist_Gamma_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Control = TH1D("Hist_Gamma_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Attenuation = TH1D("Hist_CR_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Attenuation = TH1D("Hist_Gamma_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))

energy_bins = [0.1,0.2,0.3,1000]

Hist_Erec_CR = TH1D("Hist_Erec_CR","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Raw = TH1D("Hist_Erec_CR_Raw","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data = TH1D("Hist_Erec_Data","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring = TH1D("Hist_Erec_Ring","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data_Sum = TH1D("Hist_Erec_Data_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring_Sum = TH1D("Hist_Erec_Ring_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Sum = TH1D("Hist_Erec_CR_Sum","",len(energy_bins)-1,array('d',energy_bins))


source = 'Crab'
field = 'on'

#target = 'Crab'
#target = '2ndCrab'
target = '3C264'
#target = 'PKS1424'
#target_field = 'on'
target_field = 'off'

Height_cut = ''
RW_method = ''
Fake_Removal = ''
Height_cut = ''
#RW_method = 'NoRW'
#Fake_Removal = 'NoFakeRemoved'


Elev_cut_lower = 60
Elev_cut_upper = 90
MSCL_cut_lower = -2
MSCL_cut_upper = 1
Erec_cut_lower = 0.1
Erec_cut_upper = 1000
signal_width = 1000.0
control_width = 4.0

tag = ''
Hist_e2p_Ratio = []
Hist_Gamma_Attenu = []
Hist_CR_Attenu = []
for energy in range(0,len(energy_bins)-1):
    Hist_Gamma_Attenu += [Hist_Elev_Bins.Clone()]
    Hist_CR_Attenu += [Hist_Elev_Bins.Clone()]
    for elev in range(1,Hist_Elev_Bins.GetNbinsX()+1):
        Elev_cut_lower = Hist_Elev_Bins.GetBinLowEdge(elev)
        Elev_cut_upper = Hist_Elev_Bins.GetBinLowEdge(elev+1)
        Erec_cut_lower = energy_bins[energy]
        Erec_cut_upper = energy_bins[energy+1]
        #Erec_cut_upper = 1000.
        MakeATag()
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print 'selection = %s'%(tag)
        CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu = MeasureAttenuation()
        Hist_Gamma_Attenu[energy].SetBinContent(Hist_Elev_Bins.FindBin(Elev_cut_lower),gamma_attenu)
        Hist_Gamma_Attenu[energy].SetBinError(Hist_Elev_Bins.FindBin(Elev_cut_lower),err_gamma_attenu)
        Hist_CR_Attenu[energy].SetBinContent(Hist_Elev_Bins.FindBin(Elev_cut_lower),CR_attenu)
        Hist_CR_Attenu[energy].SetBinError(Hist_Elev_Bins.FindBin(Elev_cut_lower),err_CR_attenu)

source = target
field = target_field

Elev_cut_lower = Elev_Bin[0]
Elev_cut_upper = Elev_Bin[len(Elev_Bin)-1]

#Do_analysis = False
Do_analysis = True
if Do_analysis:
    for energy in range(0,len(energy_bins)-1):
        Erec_cut_lower = energy_bins[energy]
        Erec_cut_upper = energy_bins[energy+1]
        MakeATag()
        cr_att = 1
        gamma_att = 0
        if not Fake_Removal=='NoFakeRemoved':
            cr_att = Hist_CR_Attenu[energy].GetBinContent(Hist_CR_Attenu[energy].FindBin(70.))
            gamma_att = Hist_Gamma_Attenu[energy].GetBinContent(Hist_Gamma_Attenu[energy].FindBin(70.))
        RunExtendedSourceAnalysis(cr_att,gamma_att)
