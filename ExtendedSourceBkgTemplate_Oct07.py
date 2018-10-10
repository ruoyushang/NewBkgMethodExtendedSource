
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

def ControlRegionMSCWSelection(tree,elev):
    
    if (tree.MSCW<control_width): return False
    return True

def SignalRegionMSCWSelection(tree,elev):
    
    if (tree.MSCW>MSCW_cut_upper): return False
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
    #tag += '_Erec%sto%s'%(Erec_cut_lower,Erec_cut_upper)
    tag += '_'+Height_cut
    tag += '_'+RW_method
    tag += '_'+Fake_Removal
    tag += '_'+source+'_'+field
    tag += '_MethodMSCW'


def GetEfficiencyHistogram(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '2ndCrab':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'PKS1424':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == '3C264':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'Ton599':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
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
    hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
    for i in range(0,pointing_tree.GetEntries()):
        pointing_tree.GetEntry(i)
        hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
        hist_Azi_vs_time.Fill(pointing_tree.Time,pointing_tree.TelAzimuth)
        hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
        hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

    for i in range(0,anasum_tree.GetEntries()):
        anasum_tree.GetEntry(i)
        if not CommonSelection(anasum_tree): continue
        elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
        azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
        Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        scale = Hist_TelElev_Counts_scale.GetBinContent(Hist_TelElev_Counts_scale.FindBin(elevation))
        if RW_method=='NoRW': 
            scale = 1
        if ControlRegionMSCWSelection(anasum_tree,elevation):
            Hist_OffData_Control.Fill(0,scale)
        if SignalRegionMSCWSelection(anasum_tree,elevation):
            Hist_OffData_Signal.Fill(0,scale)

def GetSourceElevation(run):

    if source == 'Crab':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == '2ndCrab':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Crab_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == 'PKS1424':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/PKS1424_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == '3C264':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/3C264_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == 'H1426':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if source == 'Ton599':
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
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
    if target == 'H1426':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    if target == 'Ton599':
        target_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
    target_pointing_tree = target_file.Get("run_%s/stereo/pointingDataReduced"%(run))

    for i in range(0,target_pointing_tree.GetEntries()):
        target_pointing_tree.GetEntry(i)
        Hist_TelElev_Counts_target.Fill(target_pointing_tree.TelElevation)


def MeasureEfficiency():

    Rate_CR_Efficiency = 0
    Err_Rate_CR_Efficiency = 0
    
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
    if target == 'H1426':
        for run in runlist_H1426:
            GetTargetElevation(run)
    if target == 'Ton599':
        for run in runlist_Ton599:
            GetTargetElevation(run)
    #print Hist_TelElev_Counts_target.Print('All')
    for b in range(0,Hist_TelElev_Counts_scale.GetNbinsX()):
        Hist_TelElev_Counts_scale.SetBinContent(b+1,Hist_TelElev_Counts_target.GetBinContent(b+1))
    Hist_TelElev_Counts_scale.Divide(Hist_TelElev_Counts_source)
    #print Hist_TelElev_Counts_scale.Print('All')
    if source == 'Crab':
        for run in Crab_runlist:
            GetEfficiencyHistogram('Crab','off',run,'anasum')
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            GetEfficiencyHistogram('Crab','off',run,'anasum')
    Hist_TelElev_Counts_scale.Reset()
    Hist_TelElev_Counts_source.Reset()
    Hist_TelElev_Counts_target.Reset()
    
    print 'Hist_OffData_Signal.GetEntries() = %s'%(Hist_OffData_Signal.GetEntries())
    print 'Hist_OffData_Control.GetEntries() = %s'%(Hist_OffData_Control.GetEntries())
    print 'Hist_OffData_Signal.Integral() = %s'%(Hist_OffData_Signal.GetBinContent(Hist_OffData_Signal.FindBin(0)))

    
    for b in range(0,Hist_CR_Efficiency.GetNbinsX()):
        Hist_CR_Efficiency.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_CR_Efficiency.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_CR_Efficiency.Divide(Hist_OffData_Signal)
    Rate_CR_Efficiency = Hist_CR_Efficiency.GetBinContent(Hist_CR_Efficiency.FindBin(0))
    Err_Rate_CR_Efficiency = Hist_CR_Efficiency.GetBinError(Hist_CR_Efficiency.FindBin(0))
    print 'Rate_CR_Efficiency = %s +/- %s'%(Rate_CR_Efficiency,Err_Rate_CR_Efficiency)
    
    Hist_OffData_Signal.Reset()
    Hist_OffData_Control.Reset()
    Hist_CR_Efficiency.Reset()

    return Rate_CR_Efficiency,Err_Rate_CR_Efficiency

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
        if source == 'H1426':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/H1426_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if source == 'Ton599':
            anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/Ton599_V6_Moderate-TMVA-BDT.RB.%s.root"%(run))
        if region == 'on':
            anasum_tree = anasum_file.Get("run_%s/stereo/data_on"%(run))
            ring_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        else:
            anasum_tree = anasum_file.Get("run_%s/stereo/data_off"%(run))
        pointing_tree = anasum_file.Get("run_%s/stereo/pointingDataReduced"%(run))
    else:
        anasum_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
        anasum_tree = anasum_file.Get("data")
        pointing_tree = anasum_file.Get("pointingDataReduced")

    if not source=="MC":
        pointing_tree.GetEntry(0)
        time_begin = pointing_tree.Time
        pointing_tree.GetEntry(pointing_tree.GetEntries()-1)
        time_end = pointing_tree.Time
        hist_Ele_vs_time = ROOT.TProfile("hist_Ele_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Azi_vs_time = ROOT.TProfile("hist_Azi_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),-400,400)
        hist_RA_vs_time = ROOT.TProfile("hist_RA_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        hist_Dec_vs_time = ROOT.TProfile("hist_Dec_vs_time","",100,time_begin,time_end+0.01*(time_end-time_begin),0,100)
        for i in range(0,pointing_tree.GetEntries()):
            pointing_tree.GetEntry(i)
            hist_Ele_vs_time.Fill(pointing_tree.Time,pointing_tree.TelElevation)
            hist_Azi_vs_time.Fill(pointing_tree.Time,pointing_tree.TelAzimuth)
            hist_RA_vs_time.Fill(pointing_tree.Time,pointing_tree.TelRAJ2000)
            hist_Dec_vs_time.Fill(pointing_tree.Time,pointing_tree.TelDecJ2000)

        for i in range(0,anasum_tree.GetEntries()):
            anasum_tree.GetEntry(i)
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
            RA = hist_RA_vs_time.GetBinContent(hist_RA_vs_time.FindBin(anasum_tree.Time))
            Dec = hist_Dec_vs_time.GetBinContent(hist_Dec_vs_time.FindBin(anasum_tree.Time))
            Hist_TelElev_Counts.Fill(elevation)
            if not CommonSelection(anasum_tree): continue
            if (anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_TelElev_vs_EmissionHeight_Data_All.Fill(elevation,anasum_tree.EmissionHeight)
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            if (anasum_tree.MSCW<MSCW_cut_upper):
                Hist2D_Erec_vs_EmissionHeight_Data_All.Fill(anasum_tree.Erec,anasum_tree.EmissionHeight)
            scale = 1
            Xoff_derot = anasum_tree.Xoff_derot
            Yoff_derot = anasum_tree.Yoff_derot
            Xoff = anasum_tree.Xoff
            Yoff = anasum_tree.Yoff
            Xcore = anasum_tree.Xcore
            Ycore = anasum_tree.Ycore
            if ControlRegionMSCWSelection(anasum_tree,elevation):
                if (anasum_tree.theta2<0.02):
                    Hist_Erec_CR_Raw.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec),1.)
		Hist2D_Xoff_vs_Yoff_CR_Raw.Fill(anasum_tree.ra,anasum_tree.dec)
		Hist_theta2_CR_Raw.Fill(anasum_tree.theta2)
            if SignalRegionMSCWSelection(anasum_tree,elevation):
                if (anasum_tree.theta2<0.02):
                    Hist_Erec_Data.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
                elif (anasum_tree.theta2>0.1 and anasum_tree.theta2<0.5):
                    Hist_Erec_Ring.Fill(min(energy_bins[len(energy_bins)-2],anasum_tree.Erec))
		Hist2D_Xoff_vs_Yoff_Data.Fill(anasum_tree.ra,anasum_tree.dec)
		Hist_theta2_Data.Fill(anasum_tree.theta2)
            if (anasum_tree.MSCW>3.):
                if (anasum_tree.theta2<0.02):
                    Hist_Norm_Data.Fill(0)
                elif (anasum_tree.theta2>0.1 and anasum_tree.theta2<0.5):
                    Hist_Norm_Ring.Fill(0)


def RunExtendedSourceAnalysis(data_type,CR_Efficiency):

    print 'Erec_cut_lower = %s, Erec_cut_upper = %s'%(Erec_cut_lower,Erec_cut_upper)
    print 'CR_Efficiency = %s'%(CR_Efficiency)

    if source == 'Crab':
        for run in Crab_runlist:
            SingleRunAnalysis('Crab',field,run,data_type)
    if source == '2ndCrab':
        for run in SecondCrab_runlist:
            SingleRunAnalysis('Crab',field,run,data_type)
    if source == 'PKS1424':
        for run in PKS1424_runlist:
            SingleRunAnalysis('PKS1424',field,run,data_type)
    if source == '3C264':
        for run in runlist_3C264:
            SingleRunAnalysis('3C264',field,run,data_type)
    if source == 'H1426':
        for run in runlist_H1426:
            SingleRunAnalysis('H1426',field,run,data_type)
    if source == 'Ton599':
        for run in runlist_Ton599:
            SingleRunAnalysis('Ton599',field,run,data_type)
    
    if field == 'on' and data_type=='anasum':
        Normalization_Factor_Ring = 0
        if not Hist_Norm_Ring.Integral()==0:
            Normalization_Factor_Ring = Hist_Norm_Data.Integral()
            Normalization_Factor_Ring = Normalization_Factor_Ring/Hist_Norm_Ring.Integral()
        Hist_Erec_Ring.Scale(Normalization_Factor_Ring)

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = 0*N_egamma + b*N_CR, b = CR efficiency
    #N_CR = N_control/b

    print 'Hist_Erec_CR_Raw.Integral() = %s'%(Hist_Erec_CR_Raw.Integral())
    print 'Hist_Erec_Data.Integral() = %s'%(Hist_Erec_Data.Integral())
    if not Fake_Removal=='NoFakeRemoved':
        Hist_Erec_CR.Reset()
        Hist_Erec_CR.Add(Hist_Erec_CR_Raw,1./CR_Efficiency)
        Hist2D_Xoff_vs_Yoff_CR.Reset()
        Hist2D_Xoff_vs_Yoff_CR.Add(Hist2D_Xoff_vs_Yoff_CR_Raw,1./CR_Efficiency)
        Hist_theta2_CR.Reset()
        Hist_theta2_CR.Add(Hist_theta2_CR_Raw,1./CR_Efficiency)
    Err_N_CR = ROOT.Double(0.)
    N_CR = Hist_Erec_CR.IntegralAndError(1,Hist_Erec_CR.GetNbinsX(),Err_N_CR)
    print 'Hist_Erec_CR.Integral() = %s +/- %s'%(N_CR,Err_N_CR)
    print 'Hist_Erec_Ring.Integral() = %s'%(Hist_Erec_Ring.Integral())
    #print Hist_theta2_CR_Raw.Print("All")
    #print Hist_theta2_CR.Print("All")
    
    
    Hist_Erec_Data_Sum.Add(Hist_Erec_Data)
    Hist_Erec_Ring_Sum.Add(Hist_Erec_Ring)
    Hist_Erec_CR_Sum.Add(Hist_Erec_CR)

    Hist2D_Xoff_vs_Yoff_Data_Sum.Add(Hist2D_Xoff_vs_Yoff_Data)
    Hist2D_Xoff_vs_Yoff_CR_Sum.Add(Hist2D_Xoff_vs_Yoff_CR)
    Hist_theta2_Data_Sum.Add(Hist_theta2_Data)
    Hist_theta2_CR_Sum.Add(Hist_theta2_CR)

    Hist_Norm_Data.Reset()
    Hist_Norm_Ring.Reset()
    Hist_Erec_CR.Reset()
    Hist_Erec_CR_Raw.Reset()
    Hist_Erec_Data.Reset()
    Hist_Erec_Ring.Reset()
    Hist2D_Xoff_vs_Yoff_Data.Reset()
    Hist2D_Xoff_vs_Yoff_CR.Reset()
    Hist2D_Xoff_vs_Yoff_CR_Raw.Reset()
    Hist_theta2_Data.Reset()
    Hist_theta2_CR.Reset()
    Hist_theta2_CR_Raw.Reset()

def EfficiencyRateAtDifferentHeight():

    global Erec_cut_lower
    global Erec_cut_upper
    global MSCL_cut_lower
    global MSCL_cut_upper
    global control_width
    global Hist_CR_Efficiency_list
    global Hist_Inflation

    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Efficiency_list += [Hist_ControlWidth_Bins.Clone()]
        Hist_Inflation += [Hist_ControlWidth_Bins.Clone()]
        for width in range(1,Hist_ControlWidth_Bins.GetNbinsX()+1):
            control_width = Hist_ControlWidth_Bins.GetBinLowEdge(width)
            Erec_cut_lower = energy_bins[energy]
            Erec_cut_upper = energy_bins[energy+1]
            if energy==len(energy_bins)-2: Erec_cut_upper = 1e10
            MakeATag()
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print 'selection = %s, control region = %s'%(tag, control_width)
            CR_Efficiency,err_CR_Efficiency = MeasureEfficiency()
            Hist_CR_Efficiency_list[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),CR_Efficiency)
            Hist_CR_Efficiency_list[energy].SetBinError(Hist_ControlWidth_Bins.FindBin(control_width),err_CR_Efficiency)
            inflation = 1./(CR_Efficiency)
            print 'inflation = %s'%(inflation)
            Hist_Inflation[energy].SetBinContent(Hist_ControlWidth_Bins.FindBin(control_width),inflation)

    canvas2 = ROOT.TCanvas("canvas2","canvas2", 200, 10, 600, 600)
    pad2 = ROOT.TPad("pad2","pad2",0.1,0.8,1,1)
    pad2.SetBottomMargin(0.0)
    pad2.SetTopMargin(0.03)
    pad2.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0.1,0,1,0.8)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    legend = ROOT.TLegend(0.2,0.1,0.7,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)

    pad1.cd()
    Hist_CR_Efficiency_list[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_CR_Efficiency_list[0].GetYaxis().SetTitle('Bkg Efficiency')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_CR_Efficiency_list[energy].GetMaximum(),maxContent)
    Hist_CR_Efficiency_list[0].SetMaximum(1.2*maxContent)
    Hist_CR_Efficiency_list[0].SetMinimum(0)
    Hist_CR_Efficiency_list[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Efficiency_list[energy].SetLineColor(energy+2)
        Hist_CR_Efficiency_list[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_CR_Efficiency_list[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/CR_Efficiency_%s.pdf'%(tag))

    pad1.cd()
    Hist_Inflation[0].GetXaxis().SetTitle('Emmission Height of Control Region')
    Hist_Inflation[0].GetYaxis().SetTitle('Unc. inflation factor')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_Inflation[energy].GetMaximum(),maxContent)
    Hist_Inflation[0].SetMaximum(1.2*maxContent)
    Hist_Inflation[0].SetMinimum(0)
    Hist_Inflation[0].Draw("E")
    for energy in range(0,len(energy_bins)-1):
        Hist_Inflation[energy].SetLineColor(energy+2)
        Hist_Inflation[energy].Draw("E same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_Inflation[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/Inflation_%s.pdf'%(tag))

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeStackPlot(Hist_Data,Hist_CR,Hist_Ring,title):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,0.8)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
    pad2.SetBottomMargin(0.39)
    pad2.SetTopMargin(0.0)
    pad2.SetBorderMode(0)
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    Hist_Data.SetLineColor(1)
    Hist_Data.SetLineWidth(3)
    set_histStyle( Hist_CR , 30)
    stack = ROOT.THStack("stack", "")
    stack.Add( Hist_CR )
    if Hist_Data.GetMaximum()>stack.GetMaximum():
        if Hist_Data.GetMinimum()>stack.GetMinimum():
            Hist_Data.SetMinimum(max(1.,stack.GetMinimum()/2.))
        Hist_Data.Draw("E")
    else:
        if Hist_Data.GetMinimum()<stack.GetMinimum():
            stack.SetMinimum(max(1.,Hist_Data.GetMinimum()/2.))
        stack.Draw("hist")
    stack.Draw("hist same")
    Hist_Data.Draw("E same")
    if not Hist_Ring.Integral()==0:
        Hist_Ring.SetLineColor(2)
        Hist_Ring.SetLineWidth(2)
        Hist_Ring.Draw("E same")
    Hist_All = Hist_CR.Clone()
    Hist_All.SetFillColor(kBlack)
    Hist_All.SetFillStyle(3004)
    Hist_All.SetMarkerSize(0)
    Hist_All.Draw("e2 same")
    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_Data,'Data',"pl")
    legend.AddEntry(Hist_CR,'Background',"f")
    if not Hist_Ring.Integral()==0:
        legend.AddEntry(Hist_Ring,'Ring model',"pl")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.6,'%0.2f < Tel elev. < %0.2f'%(Elev_cut_lower,Elev_cut_upper) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.45,'%0.2f < MSCL < %0.2f'%(MSCL_cut_lower,MSCL_cut_upper) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.3,'%0.2f < MSCW < %0.2f'%(MSCL_cut_lower,MSCL_cut_upper) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    pad2.cd()
    Hist_Band = Hist_All.Clone()
    Hist_Band.Divide(Hist_All)
    Hist_Band.SetFillColor(kBlack)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    if not Hist_Ring.Integral()==0:
        Hist_Band.GetYaxis().SetTitle("Ring/Bkg")
    else:
        Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.SetMaximum(2.5)
    Hist_Band.SetMinimum(0)
    Hist_Band.Draw("e2")
    if not Hist_Ring.Integral()==0:
        Hist_Ratio = Hist_Ring.Clone()
    else:
        Hist_Ratio = Hist_Data.Clone()
    Hist_Ratio.Divide(Hist_All)
    Hist_Ratio.Draw("E same")
    if title=='Erec':
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()
    c_both.SaveAs('output/%s_%s.pdf'%(title,tag))
    pad1.SetLogy(0)
    pad1.SetLogx(0)
    pad2.SetLogx(0)

Hist_Norm_Data = TH1D("Hist_Norm_Data","",1,0,1)
Hist_Norm_Ring = TH1D("Hist_Norm_Ring","",1,0,1)

Elev_Bin = [55,85]
#Elev_Bin = [30,40,50,60,70,80,90]
Hist_Elev_Bins = TH1D("Hist_Elev_Bins","",len(Elev_Bin)-1,array('d',Elev_Bin))

ControlWidth_Bin = [3,4]
Hist_ControlWidth_Bins = TH1D("Hist_ControlWidth_Bins","",len(ControlWidth_Bin)-1,array('d',ControlWidth_Bin))

MSCW_Bin = [0,1]
Hist_OffData_Signal = TH1D("Hist_OffData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Control = TH1D("Hist_OffData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Efficiency = TH1D("Hist_CR_Efficiency","",len(MSCW_Bin)-1,array('d',MSCW_Bin))

energy_bins = [pow(10,0),pow(10,0.2),pow(10,0.4),pow(10,0.6),pow(10,0.8),pow(10,1.0)]
ControlWidthAtThisEnergy = [3,3,3,3,3,3]

Hist_Erec_CR = TH1D("Hist_Erec_CR","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Raw = TH1D("Hist_Erec_CR_Raw","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data = TH1D("Hist_Erec_Data","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring = TH1D("Hist_Erec_Ring","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Data_Sum = TH1D("Hist_Erec_Data_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_Ring_Sum = TH1D("Hist_Erec_Ring_Sum","",len(energy_bins)-1,array('d',energy_bins))
Hist_Erec_CR_Sum = TH1D("Hist_Erec_CR_Sum","",len(energy_bins)-1,array('d',energy_bins))

Xcenter = 0
Ycenter = 0
Window = 2.
Hist2D_Xoff_vs_Yoff_Data = TH2D("Hist2D_Xoff_vs_Yoff_Data","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR = TH2D("Hist2D_Xoff_vs_Yoff_CR","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Raw = TH2D("Hist2D_Xoff_vs_Yoff_CR_Raw","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Data_Sum = TH2D("Hist2D_Xoff_vs_Yoff_Data_Sum","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_CR_Sum = TH2D("Hist2D_Xoff_vs_Yoff_CR_Sum","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist2D_Xoff_vs_Yoff_Sig = TH2D("Hist2D_Xoff_vs_Yoff_Sig","",20,Xcenter-Window,Xcenter+Window,20,Ycenter-Window,Ycenter+Window)
Hist_theta2_Data = TH1D("Hist_theta2_Data","",20,0,0.5)
Hist_theta2_CR = TH1D("Hist_theta2_CR","",20,0,0.5)
Hist_theta2_CR_Raw = TH1D("Hist_theta2_CR_Raw","",20,0,0.5)
Hist_theta2_Data_Sum = TH1D("Hist_theta2_Data_Sum","",20,0,0.5)
Hist_theta2_CR_Sum = TH1D("Hist_theta2_CR_Sum","",20,0,0.5)
Hist_theta2_Ring_Sum = TH1D("Hist_theta2_Ring_Sum","",20,0,0.5)

variable_bins = []
b = -1
while (b<2):
    variable_bins += [pow(10,b)]
    b += 0.2
Hist2D_Erec_vs_EmissionHeight_Data_All = TH2D("Hist2D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins),40,0,40)
Hist1D_Erec_vs_EmissionHeight_Data_All = TH1D("Hist1D_Erec_vs_EmissionHeight_Data_All","",len(variable_bins)-1,array('d',variable_bins))
Hist2D_TelElev_vs_EmissionHeight_Data_All = TH2D("Hist2D_TelElev_vs_EmissionHeight_Data_All","",18,0,90,40,0,40)
Hist1D_TelElev_vs_EmissionHeight_Data_All = TH1D("Hist1D_TelElev_vs_EmissionHeight_Data_All","",18,0,90)

source = 'Crab'
field = 'on'

#target = 'Crab'
#target = '2ndCrab'
target = '3C264'
#target = 'PKS1424'
#target = 'H1426'
#target = 'Ton599'
target_field = 'on'
#target_field = 'off'

Height_cut = ''
RW_method = ''
Fake_Removal = ''
Height_cut = ''
#RW_method = 'NoRW'
#Fake_Removal = 'NoFakeRemoved'


Elev_cut_lower = Elev_Bin[0]
Elev_cut_upper = Elev_Bin[len(Elev_Bin)-1]
MSCL_cut_lower = -2
MSCL_cut_upper = 0.7
MSCW_cut_lower = -2
MSCW_cut_upper = 0.7
Erec_cut_lower = 0.1
Erec_cut_upper = 1000

Hist_CR_Efficiency_list = []
Hist_Inflation = []

signal_width = 1000.0
control_width = 4.0
EfficiencyRateAtDifferentHeight()

source = target
field = target_field


#Do_analysis = False
Do_analysis = True
if Do_analysis:
    for energy in range(0,len(energy_bins)-1):
        Erec_cut_lower = energy_bins[energy]
        Erec_cut_upper = energy_bins[energy+1]
        if energy==len(energy_bins)-2: Erec_cut_upper = 1e10
        MakeATag()
        cr_eff = 1
        control_width = ControlWidthAtThisEnergy[energy]
        cr_eff = Hist_CR_Efficiency_list[energy].GetBinContent(Hist_CR_Efficiency_list[energy].FindBin(control_width))
        print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print 'control region = %s'%(control_width)
        inflation = 1./cr_eff
        print 'inflation = %s'%(inflation)
        RunExtendedSourceAnalysis('anasum',cr_eff)
    MakeStackPlot(Hist_Erec_Data_Sum,Hist_Erec_CR_Sum,Hist_Erec_Ring_Sum,'Erec')
    MakeStackPlot(Hist_theta2_Data_Sum,Hist_theta2_CR_Sum,Hist_theta2_Ring_Sum,'theta2')
    #print Hist_theta2_CR_Sum.Print("All")

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    Hist2D_TelElev_vs_EmissionHeight_Data_All.GetYaxis().SetTitle('EmissionHeight')
    Hist2D_TelElev_vs_EmissionHeight_Data_All.GetXaxis().SetTitle('Tel. elevation')
    Hist2D_TelElev_vs_EmissionHeight_Data_All.Draw("COL4Z")
    for b in range(0,Hist1D_TelElev_vs_EmissionHeight_Data_All.GetNbinsX()):
        hist_temp = Hist2D_TelElev_vs_EmissionHeight_Data_All.ProjectionY("hist_temp",b+1,b+1)
        Hist1D_TelElev_vs_EmissionHeight_Data_All.SetBinContent(b+1,hist_temp.GetMean())
        Hist1D_TelElev_vs_EmissionHeight_Data_All.SetBinError(b+1,hist_temp.GetRMS())
    Hist1D_TelElev_vs_EmissionHeight_Data_All.SetLineColor(2)
    Hist1D_TelElev_vs_EmissionHeight_Data_All.Draw("E same")
    #f_cut_mean = ROOT.TF1("f_cut_mean","10.+10.*pow((x-90)/50.,4)",0,90)
    #f_cut_mean.SetLineColor(2)
    #f_cut_mean.SetLineStyle(2)
    #f_cut_mean.Draw("same")
    #f_cut_upper = ROOT.TF1("f_cut_upper","3.0+10.+10.*pow((x-90)/47.,4)",0,90)
    #f_cut_upper.SetLineColor(2)
    #f_cut_upper.Draw("same")
    #f_cut_lower = ROOT.TF1("f_cut_lower","-3.0+10.+10.*pow((x-90)/53.,4)",0,90)
    #f_cut_lower.SetLineColor(2)
    #f_cut_lower.Draw("same")
    canvas.SetLogz()
    canvas.SaveAs('output/TelElev_vs_EmissionHeight_All_%s.pdf'%(tag))

    Hist2D_Erec_vs_EmissionHeight_Data_All.GetYaxis().SetTitle('EmissionHeight')
    Hist2D_Erec_vs_EmissionHeight_Data_All.GetXaxis().SetTitle('Erec')
    Hist2D_Erec_vs_EmissionHeight_Data_All.Draw("COL4Z")
    for b in range(0,Hist1D_Erec_vs_EmissionHeight_Data_All.GetNbinsX()):
        hist_temp = Hist2D_Erec_vs_EmissionHeight_Data_All.ProjectionY("hist_temp",b+1,b+1)
        Hist1D_Erec_vs_EmissionHeight_Data_All.SetBinContent(b+1,hist_temp.GetMean())
        Hist1D_Erec_vs_EmissionHeight_Data_All.SetBinError(b+1,hist_temp.GetRMS())
    Hist1D_Erec_vs_EmissionHeight_Data_All.SetLineColor(2)
    Hist1D_Erec_vs_EmissionHeight_Data_All.Draw("E same")
    canvas.SetLogz()
    canvas.SetLogx()
    canvas.SaveAs('output/Erec_vs_EmissionHeight_All_%s.pdf'%(tag))

    canvas.SetLogz(0)
    canvas.SetLogx(0)
    canvas.SetLogy(0)
    for bx in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsX()):
        for by in range(0,Hist2D_Xoff_vs_Yoff_Sig.GetNbinsY()):
            this_data = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinContent(bx,by)
            this_data_err = Hist2D_Xoff_vs_Yoff_Data_Sum.GetBinError(bx,by)
            this_bkg = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinContent(bx,by)
            this_bkg_err = Hist2D_Xoff_vs_Yoff_CR_Sum.GetBinError(bx,by)
            value = 0
            if not this_bkg_err==0:
                value = (this_data-this_bkg)/pow(this_bkg_err*this_bkg_err+this_data_err*this_data_err,0.5)
            Hist2D_Xoff_vs_Yoff_Sig.SetBinContent(bx,by,value)
    Hist2D_Xoff_vs_Yoff_Sig.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_Sig.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_Sig.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_Sig_%s.pdf'%(tag))
    Hist2D_Xoff_vs_Yoff_Data_Sum.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_Data_Sum.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_Data_Sum.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_Data_%s.pdf'%(tag))
    Hist2D_Xoff_vs_Yoff_CR_Sum.GetYaxis().SetTitle('Yoff')
    Hist2D_Xoff_vs_Yoff_CR_Sum.GetXaxis().SetTitle('Xoff')
    Hist2D_Xoff_vs_Yoff_CR_Sum.Draw("COL4Z")
    canvas.SaveAs('output/Xoff_vs_Yoff_CR_%s.pdf'%(tag))
