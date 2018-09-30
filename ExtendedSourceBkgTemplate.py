
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

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeStackPlot(Hist_Data,Hist_CR,Hist_ES,Hist_Fake,Hist_Ring,title):
    
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
    set_histStyle( Hist_ES , 38)
    stack = ROOT.THStack("stack", "")
    stack.Add( Hist_CR )
    stack.Add( Hist_ES )
    if Hist_Data.GetMaximum()>stack.GetMaximum():
        Hist_Data.Draw("E")
    else:
        stack.Draw("hist")
    stack.Draw("hist same")
    Hist_Data.Draw("E same")
    #Hist_Fake.SetLineColor(9)
    #Hist_Fake.Draw("E same")
    if not Hist_Ring.Integral()==0:
        Hist_Ring.SetLineColor(2)
        Hist_Ring.SetLineWidth(2)
        Hist_Ring.SetLineStyle(2)
        Hist_Ring.Draw("E same")
    Hist_All = Hist_CR.Clone()
    Hist_All.Add(Hist_ES)
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
    legend.AddEntry(Hist_Data,'Data at low alt.',"pl")
    legend.AddEntry(Hist_CR,'Hadron',"f")
    legend.AddEntry(Hist_ES,'e + others',"f")
    #legend.AddEntry(Hist_Fake,'e/#gamma contamin.',"pl")
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
    lumilab4 = ROOT.TLatex(0.15,0.3,'%0.2f < Erec < %0.2f'%(Erec_cut_lower,Erec_cut_upper) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    #lumilab5 = ROOT.TLatex(0.15,0.15,'e/p = %0.2f +/- %0.2f'%(Ratio_ES2CR,Err_Ratio_ES2CR) )
    #lumilab5.SetNDC()
    #lumilab5.SetTextSize(0.15)
    #lumilab5.Draw()
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
    Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.SetMaximum(2.5)
    Hist_Band.SetMinimum(0)
    Hist_Band.Draw("e2")
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

def CommonSelection(tree):
    #if (tree.theta2>3.0): return False
    if (tree.MSCL<MSCL_cut_lower): return False
    if (tree.MSCL>MSCL_cut_upper): return False
    if (tree.Erec<Erec_cut_lower): return False
    if (tree.Erec>Erec_cut_upper): return False
    return True

def ControlRegionHeightSelection(tree,elev):
    # best
    #width = 3.
    #if elev>60:
    #    mean_height = 10.
    #    cut_upper = mean_height+width
    #    cut_lower = mean_height-width
    #else:
    #    mean_height = -0.5*(elev-20.)+30
    #    cut_upper = mean_height+(width+(mean_height-10.)*0.15)
    #    cut_lower = mean_height-(width+(mean_height-10.)*0.15)

    # best
    cut_upper = 3.0+10.+10.*pow((elev-90)/47.,4)-(1+log10(tree.Erec))
    cut_lower = -3.0+10.+10.*pow((elev-90)/53.,4)-(1+log10(tree.Erec))

    if Height_cut=='simple':
        cut_upper = 13.
        cut_lower = 7.

    if (tree.EmissionHeight<cut_upper and tree.EmissionHeight>cut_lower): return False
    return True

def FillReweightingHistogram(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
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
        if (anasum_tree.MSCW<3): continue
	elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
	azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
        if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
        if ControlRegionHeightSelection(anasum_tree,elevation):
            Hist_TelElev_Control.Fill(elevation)
            Hist_MSCL_Control.Fill(anasum_tree.MSCL)
            Hist_Erec_Control.Fill(anasum_tree.Erec)
        else:
            Hist_TelElev_Signal.Fill(elevation)
            Hist_MSCL_Signal.Fill(anasum_tree.MSCL)
            Hist_Erec_Signal.Fill(anasum_tree.Erec)
    #print Hist_TelElev_Control.Print('All')
        
def FillHistBkgElectron(run,scale_data_MC):

    input_file = ROOT.TFile.Open("$VERITAS_USER_DATA_DIR/analysis/Results/%s.mscw.root"%(run))
    tree = input_file.Get("data")

    total_entries = tree.GetEntries()
    total_entries = int(float(tree.GetEntries())/10.)
    for i in range(0,total_entries):
    	tree.GetEntry(i)
        if not CommonSelection(tree): continue
        if not ControlRegionHeightSelection(tree,90.):
            scale = scale_data_MC
            Hist_MSCW_ES.Fill((tree.MSCW+1.0),scale)
            Hist_Erec_ES.Fill(tree.Erec,scale)

def GetAnttenuation(source,region,run,data_type):

    if data_type=='anasum':
        if source == 'Crab':
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
        if ControlRegionHeightSelection(anasum_tree,elevation):
            if region == 'on':
                Hist_OnData_Control.Fill(anasum_tree.MSCW)
            else:
                Hist_OffData_Control.Fill(anasum_tree.MSCW)
        else:
            if region == 'on':
                Hist_OnData_Signal.Fill(anasum_tree.MSCW)
            else:
                Hist_OffData_Signal.Fill(anasum_tree.MSCW)
        

def SingleRunAnalysis(source,region,run,data_type):

    global N_GammaMC_Fake
    global N_GammaData_Fake
    global N_GammaMC_ES
    global N_GammaData_ES

    if data_type=='anasum':
        if source == 'Crab':
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

        #print Hist_TelElev_RW_Factor.Print("All")
        for i in range(0,anasum_tree.GetEntries()):
            anasum_tree.GetEntry(i)
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(anasum_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(anasum_tree.Time))
            Hist_TelElev_Counts.Fill(elevation)
            if not CommonSelection(anasum_tree): continue
            #if (anasum_tree.MSCW<3): continue
            if (anasum_tree.MSCW<1 and anasum_tree.MSCL<1):
                TProf_TelElev_vs_EmissionHeight_Data.Fill(elevation,anasum_tree.EmissionHeight)
                Hist2D_TelElev_vs_EmissionHeight_Data.Fill(elevation,anasum_tree.EmissionHeight)
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            Hist2D_Erec_vs_EmissionHeight_Data.Fill(anasum_tree.Erec,anasum_tree.EmissionHeight)
            scale = Hist_TelElev_RW_Factor.GetBinContent(Hist_TelElev_RW_Factor.FindBin(elevation))
            scale = scale*Hist_MSCL_RW_Factor.GetBinContent(Hist_MSCL_RW_Factor.FindBin(anasum_tree.MSCL))
            scale = scale*Hist_Erec_RW_Factor.GetBinContent(Hist_Erec_RW_Factor.FindBin(anasum_tree.Erec))
            if RW_method=='NoRW': 
                scale = 1
            if ControlRegionHeightSelection(anasum_tree,elevation):
                Hist_MSCW_CR.Fill(anasum_tree.MSCW,scale)
                Hist_MSCW_CR_Raw.Fill(anasum_tree.MSCW,1.)
                Hist_Erec_CR.Fill(anasum_tree.Erec,scale)
                Hist_Erec_CR_Raw.Fill(anasum_tree.Erec,1.)
                Hist2D_MSCW_vs_EmissionHeight_Data.Fill(anasum_tree.MSCW,anasum_tree.EmissionHeight,scale)
                if (anasum_tree.MSCW<Fake_cut_MSCW and anasum_tree.MSCL<Fake_cut_MSCL):
                    N_GammaData_Fake += scale 
            else:
                Hist_MSCW_Data.Fill(anasum_tree.MSCW)
                Hist_Erec_Data.Fill(anasum_tree.Erec)
                Hist2D_MSCW_vs_EmissionHeight_Data.Fill(anasum_tree.MSCW,anasum_tree.EmissionHeight)
                N_GammaData_ES += 1.0 
    else:
        total_entries = anasum_tree.GetEntries()
        total_entries = int(float(anasum_tree.GetEntries())/10.)
        for i in range(0,total_entries):
	    anasum_tree.GetEntry(i)
            if not CommonSelection(anasum_tree): continue
            if ControlRegionHeightSelection(anasum_tree,90.):
                if ((anasum_tree.MSCW+1.0)<Fake_cut_MSCW and anasum_tree.MSCL<Fake_cut_MSCL):
                    N_GammaMC_Fake += 1 
            else:
                Hist_MSCW_ES_Raw.Fill(anasum_tree.MSCW+1.0)
                Hist_Erec_ES_Raw.Fill(anasum_tree.Erec)

    if region == 'on':
        for i in range(0,ring_tree.GetEntries()):
            ring_tree.GetEntry(i)
            if not CommonSelection(ring_tree): continue
            #if (ring_tree.MSCW<3): continue
            elevation = hist_Ele_vs_time.GetBinContent(hist_Ele_vs_time.FindBin(ring_tree.Time))
            azimuth = hist_Azi_vs_time.GetBinContent(hist_Azi_vs_time.FindBin(ring_tree.Time))
            if elevation<Elev_cut_lower or elevation>Elev_cut_upper: continue
            if not ControlRegionHeightSelection(ring_tree,elevation):
                Hist_MSCW_Ring.Fill(ring_tree.MSCW)
                Hist_Erec_Ring.Fill(ring_tree.Erec)
        
def RunExtendedSourceAnalysis(e2p_ratio,CR_attenu,gamma_attenu):

    if source == 'Crab':
        for run in Crab_runlist:
            FillReweightingHistogram('Crab',field,run,'anasum')
    if source == 'PKS1424':
        for run in PKS1424_runlist:
            FillReweightingHistogram('PKS1424',field,run,'anasum')
    if source == '3C264':
        for run in runlist_3C264:
            FillReweightingHistogram('3C264',field,run,'anasum')
    for b in range(0,Hist_TelElev_RW_Factor.GetNbinsX()):
        Hist_TelElev_RW_Factor.SetBinContent(b+1,Hist_TelElev_Signal.GetBinContent(b+1))
        Hist_TelElev_RW_Factor.SetBinError(b+1,Hist_TelElev_Signal.GetBinError(b+1))
    Hist_TelElev_RW_Factor.Divide(Hist_TelElev_Control)
    for b in range(0,Hist_MSCL_RW_Factor.GetNbinsX()):
        Hist_MSCL_RW_Factor.SetBinContent(b+1,Hist_MSCL_Signal.GetBinContent(b+1))
        Hist_MSCL_RW_Factor.SetBinError(b+1,Hist_MSCL_Signal.GetBinError(b+1))
    Hist_MSCL_RW_Factor.Divide(Hist_MSCL_Control)
    for b in range(0,Hist_Erec_RW_Factor.GetNbinsX()):
        Hist_Erec_RW_Factor.SetBinContent(b+1,Hist_Erec_Signal.GetBinContent(b+1))
        Hist_Erec_RW_Factor.SetBinError(b+1,Hist_Erec_Signal.GetBinError(b+1))
    Hist_Erec_RW_Factor.Divide(Hist_Erec_Control)
    if source == 'Crab':
        for run in Crab_runlist:
            SingleRunAnalysis('Crab',field,run,'anasum')
    if source == 'PKS1424':
        for run in PKS1424_runlist:
            SingleRunAnalysis('PKS1424',field,run,'anasum')
    if source == '3C264':
        for run in runlist_3C264:
            SingleRunAnalysis('3C264',field,run,'anasum')
    
    Normalization_Factor = Hist_MSCW_Data.Integral(Hist_MSCW_Data.FindBin(3),Hist_MSCW_Data.FindBin(100))
    Normalization_Factor = Normalization_Factor/Hist_MSCW_CR.Integral(Hist_MSCW_CR.FindBin(3),Hist_MSCW_CR.FindBin(100))
    print 'Normalization_Factor = %s'%(Normalization_Factor)
    Hist_MSCW_CR.Scale(Normalization_Factor)
    Hist_Erec_CR.Scale(Normalization_Factor)
    if field == 'on':
        Normalization_Factor_Ring = Hist_MSCW_Data.Integral(Hist_MSCW_Data.FindBin(3),Hist_MSCW_Data.FindBin(100))
        Normalization_Factor_Ring = Normalization_Factor_Ring/Hist_MSCW_Ring.Integral(Hist_MSCW_Ring.FindBin(3),Hist_MSCW_Ring.FindBin(100))
        Hist_MSCW_Ring.Scale(Normalization_Factor_Ring)
        Hist_Erec_Ring.Scale(Normalization_Factor_Ring)

    #SingleRunAnalysis('MC','off',"20130",'mscw')
    #N_Data_All = Hist_MSCW_Data.Integral()
    #N_CR_All = Hist_MSCW_CR.Integral()
    #N_ES_All = e2p_ratio*N_CR_All
    #Scale_ES = N_ES_All/Hist_MSCW_ES_Raw.Integral()
    #if Add_Electron=='NoEle':
    #    Scale_ES = 0
    #FillHistBkgElectron("20130",Scale_ES)

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = a*N_egamma + b*N_CR, a = gamma attenuation, b = CR attenuation
    #N_control - b*N_signal = (a-b)*N_egamma
    #N_egamma = (N_control-b*N_signal)/(a-b)
    #N_CR = N_signal-N_egamma

    if not Fake_Removal=='NoFakeRemoved':
        #CR_attenu = Hist_MSCW_CR_Raw.Integral(Hist_MSCW_CR_Raw.FindBin(3),Hist_MSCW_CR_Raw.FindBin(100))
        #CR_attenu = CR_attenu/Hist_MSCW_Data.Integral(Hist_MSCW_Data.FindBin(3),Hist_MSCW_Data.FindBin(100))
        #for b in range(0,Hist_MSCW_CR_Raw.GetNbinsX()):
        for b in range(Hist_MSCW_CR_Raw.FindBin(-2.5),Hist_MSCW_CR_Raw.FindBin(2.5)):
            N_egamma = (Hist_MSCW_CR_Raw.GetBinContent(b+1)-CR_attenu*Hist_MSCW_Data.GetBinContent(b+1))
            N_egamma = N_egamma/(gamma_attenu-CR_attenu)
            if N_egamma>0:
                #Hist_MSCW_CR.SetBinContent(b+1,(Hist_MSCW_Data.GetBinContent(b+1)-N_egamma)*Normalization_Factor)
                Hist_MSCW_CR.SetBinContent(b+1,max(0,Hist_MSCW_Data.GetBinContent(b+1)-N_egamma))
                Hist_MSCW_Fake.SetBinContent(b+1,N_egamma*gamma_attenu*Normalization_Factor)
        #Hist_MSCW_CR.Add(Hist_MSCW_Fake,-1.)
        for b in range(0,Hist_Erec_CR_Raw.GetNbinsX()):
            N_egamma = (Hist_Erec_CR_Raw.GetBinContent(b+1)-CR_attenu*Hist_Erec_Data.GetBinContent(b+1))
            N_egamma = N_egamma/(gamma_attenu-CR_attenu)
            if N_egamma>0:
                #Hist_Erec_CR.SetBinContent(b+1,(Hist_Erec_Data.GetBinContent(b+1)-N_egamma)*Normalization_Factor)
                Hist_Erec_CR.SetBinContent(b+1,max(0,Hist_Erec_Data.GetBinContent(b+1)-N_egamma))
                Hist_Erec_Fake.SetBinContent(b+1,N_egamma*gamma_attenu*Normalization_Factor)
        #Hist_Erec_CR.Add(Hist_Erec_Fake,-1.)
    
    SingleRunAnalysis('MC','off',"20130",'mscw')
    N_Data_All = Hist_MSCW_Data.Integral()
    N_CR_All = Hist_MSCW_CR.Integral()
    N_ES_All = e2p_ratio*N_CR_All
    Scale_ES = N_ES_All/Hist_MSCW_ES_Raw.Integral()
    if Add_Electron=='NoEle':
        Scale_ES = 0
    FillHistBkgElectron("20130",Scale_ES)

    MakeStackPlot(Hist_MSCW_Data,Hist_MSCW_CR,Hist_MSCW_ES,Hist_MSCW_Fake,Hist_MSCW_Ring,'MSCW')
    #MakeStackPlot(Hist_Erec_Data,Hist_Erec_CR,Hist_Erec_ES,Hist_Erec_Fake,Hist_Erec_Ring,'Erec')
    
    Hist_MSCW_Data_Sum.Add(Hist_MSCW_Data)
    Hist_MSCW_Ring_Sum.Add(Hist_MSCW_Ring)
    Hist_MSCW_Fake_Sum.Add(Hist_MSCW_Fake)
    Hist_MSCW_ES_Sum.Add(Hist_MSCW_ES)
    Hist_MSCW_CR_Sum.Add(Hist_MSCW_CR)
    Hist_Erec_Data_Sum.Add(Hist_Erec_Data)
    Hist_Erec_Ring_Sum.Add(Hist_Erec_Ring)
    Hist_Erec_Fake_Sum.Add(Hist_Erec_Fake)
    Hist_Erec_ES_Sum.Add(Hist_Erec_ES)
    Hist_Erec_CR_Sum.Add(Hist_Erec_CR)


    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    #Hist2D_MSCW_vs_EmissionHeight_Data.GetYaxis().SetTitle('EmissionHeight')
    #Hist2D_MSCW_vs_EmissionHeight_Data.GetXaxis().SetTitle('MSCW')
    #Hist2D_MSCW_vs_EmissionHeight_Data.Draw("COL4Z")
    #Hist2D_MSCW_vs_EmissionHeight_Data.SetLineColor(1)
    #Hist2D_MSCW_vs_EmissionHeight_Data.Draw("CONT3 same")
    #canvas.SetLogz()
    #canvas.SaveAs('output/MSCW_vs_EmissionHeight_%s.pdf'%(tag))
    Hist2D_TelElev_vs_EmissionHeight_Data.GetYaxis().SetTitle('EmissionHeight')
    Hist2D_TelElev_vs_EmissionHeight_Data.GetXaxis().SetTitle('Tel. elevation')
    Hist2D_TelElev_vs_EmissionHeight_Data.Draw("COL4Z")
    for b in range(0,Hist1D_TelElev_vs_EmissionHeight_Data.GetNbinsX()):
        hist_temp = Hist2D_TelElev_vs_EmissionHeight_Data.ProjectionY("hist_temp",b+1,b+1)
        Hist1D_TelElev_vs_EmissionHeight_Data.SetBinContent(b+1,hist_temp.GetMean())
        Hist1D_TelElev_vs_EmissionHeight_Data.SetBinError(b+1,hist_temp.GetRMS())
    Hist1D_TelElev_vs_EmissionHeight_Data.Draw("E same")
    f_cut_mean = ROOT.TF1("f_cut_mean","10.+10.*pow((x-90)/50.,4)",0,90)
    f_cut_mean.SetLineColor(2)
    f_cut_mean.SetLineStyle(2)
    f_cut_mean.Draw("same")
    f_cut_upper = ROOT.TF1("f_cut_upper","3.0+10.+10.*pow((x-90)/47.,4)",0,90)
    f_cut_upper.SetLineColor(2)
    f_cut_upper.Draw("same")
    f_cut_lower = ROOT.TF1("f_cut_lower","-3.0+10.+10.*pow((x-90)/53.,4)",0,90)
    f_cut_lower.SetLineColor(2)
    f_cut_lower.Draw("same")
    canvas.SetLogz()
    canvas.SaveAs('output/TelElev_vs_EmissionHeight_%s.pdf'%(tag))
    #Hist2D_Erec_vs_EmissionHeight_Data.GetYaxis().SetTitle('EmissionHeight')
    #Hist2D_Erec_vs_EmissionHeight_Data.GetXaxis().SetTitle('Erec')
    #Hist2D_Erec_vs_EmissionHeight_Data.Draw("COL4Z")
    #Hist2D_Erec_vs_EmissionHeight_Data.SetLineColor(1)
    #Hist2D_Erec_vs_EmissionHeight_Data.Draw("CONT3 same")
    #canvas.SetLogx()
    #canvas.SetLogz()
    #canvas.SaveAs('output/Erec_vs_EmissionHeight_%s.pdf'%(tag))

    Hist_MSCW_CR.Reset()
    Hist_MSCW_CR_Raw.Reset()
    Hist_MSCW_Fake.Reset()
    Hist_MSCW_ES.Reset()
    Hist_MSCW_ES_Raw.Reset()
    Hist_MSCW_Data.Reset()
    Hist_MSCW_Ring.Reset()
    Hist_Erec_CR.Reset()
    Hist_Erec_CR_Raw.Reset()
    Hist_Erec_Fake.Reset()
    Hist_Erec_ES.Reset()
    Hist_Erec_ES_Raw.Reset()
    Hist_Erec_Data.Reset()
    Hist_Erec_Ring.Reset()
    Hist2D_MSCW_vs_EmissionHeight_Data.Reset()
    TProf_TelElev_vs_EmissionHeight_Data.Reset()
    Hist1D_TelElev_vs_EmissionHeight_Data.Reset()
    Hist2D_TelElev_vs_EmissionHeight_Data.Reset()
    Hist2D_Erec_vs_EmissionHeight_Data.Reset()

def MakeATag():
    global tag
    tag = ''
    tag += '_Elev%sto%s'%(Elev_cut_lower,Elev_cut_upper)
    tag += '_MSCL%sto%s'%(MSCL_cut_lower,MSCL_cut_upper)
    tag += '_Erec%sto%s'%(Erec_cut_lower,Erec_cut_upper)
    tag += '_'+Height_cut
    tag += '_'+RW_method
    tag += '_'+Add_Electron
    tag += '_'+Fake_Removal
    tag += '_'+source+'_'+field


def MeasureElectronComponent():

    Rate_CR_Attenuation = 0
    Err_Rate_CR_Attenuation = 0
    Rate_Gamma_Attenuation = 0
    Err_Rate_Gamma_Attenuation = 0
    
    if source == 'Crab':
        for run in Crab_runlist:
            GetAnttenuation('Crab','on',run,'anasum')
            GetAnttenuation('Crab','off',run,'anasum')
    if source == 'PKS1424':
        for run in PKS1424_runlist:
            GetAnttenuation('PKS1424','on',run,'anasum')
            GetAnttenuation('PKS1424','off',run,'anasum')
    if source == '3C264':
        for run in runlist_3C264:
            GetAnttenuation('3C264','on',run,'anasum')
            GetAnttenuation('3C264','off',run,'anasum')
    
    print 'Hist_OnData_Signal.GetEntries() = %s'%(Hist_OnData_Signal.GetEntries())
    print 'Hist_OnData_Control.GetEntries() = %s'%(Hist_OnData_Control.GetEntries())
    print 'Hist_OffData_Signal.GetEntries() = %s'%(Hist_OffData_Signal.GetEntries())
    print 'Hist_OffData_Control.GetEntries() = %s'%(Hist_OffData_Control.GetEntries())

    Normalization_Factor_Signal = Hist_OnData_Signal.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Signal = Normalization_Factor_Signal/Hist_OffData_Signal.GetBinContent(len(MSCW_Bin)-1)
    #print 'Normalization_Factor_Signal = %s'%(Normalization_Factor_Signal)
    Normalization_Factor_Control = Hist_OnData_Control.GetBinContent(len(MSCW_Bin)-1)
    Normalization_Factor_Control = Normalization_Factor_Control/Hist_OffData_Control.GetBinContent(len(MSCW_Bin)-1)
    #print 'Normalization_Factor_Control = %s'%(Normalization_Factor_Control)
    Hist_OffData_Signal.Scale(Normalization_Factor_Signal)
    Hist_OffData_Control.Scale(Normalization_Factor_Signal)
    
    for b in range(0,Hist_CR_Attenuation.GetNbinsX()):
        Hist_CR_Attenuation.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_CR_Attenuation.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_CR_Attenuation.Divide(Hist_OffData_Signal)
    Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinContent(len(MSCW_Bin)-1)
    Err_Rate_CR_Attenuation = Hist_CR_Attenuation.GetBinError(len(MSCW_Bin)-1)
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
    Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinContent(Hist_Gamma_Attenuation.FindBin(0.5))
    Err_Rate_Gamma_Attenuation = Hist_Gamma_Attenuation.GetBinError(Hist_Gamma_Attenuation.FindBin(0.5))
    print 'Rate_Gamma_Attenuation = %s +/- %s'%(Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation)
    

    #N_signal (number of events in signal region) = N_egamma + N_CR
    #N_control (number of events in control region) = a*N_egamma + b*N_CR, a = gamma attenuation, b = CR attenuation
    #N_control - b*N_signal = (a-b)*N_egamma
    #N_egamma = (N_control-b*N_signal)/(a-b)
    #N_CR = N_signal-N_egamma
    for b in range(0,Hist_eGamma_Signal.GetNbinsX()):
        Hist_eGamma_Signal.SetBinContent(b+1,Hist_OnData_Control.GetBinContent(b+1))
        Hist_eGamma_Signal.SetBinError(b+1,Hist_OnData_Control.GetBinError(b+1))
    Hist_eGamma_Signal.Add(Hist_OnData_Signal,-1.*Rate_CR_Attenuation)
    Hist_eGamma_Signal.Scale(1./(Rate_Gamma_Attenuation-Rate_CR_Attenuation))
    for b in range(0,Hist_eGamma_Signal.GetNbinsX()):
        Hist_eGamma_Signal.SetBinContent(b+1,max(0.,Hist_eGamma_Signal.GetBinContent(b+1)))
    print Hist_eGamma_Signal.Print("All")

    for b in range(0,Hist_El_Signal.GetNbinsX()):
        Hist_El_Signal.SetBinContent(b+1,Hist_OffData_Control.GetBinContent(b+1))
        Hist_El_Signal.SetBinError(b+1,Hist_OffData_Control.GetBinError(b+1))
    Hist_El_Signal.Add(Hist_OffData_Signal,-1.*Rate_CR_Attenuation)
    Hist_El_Signal.Scale(1./(Rate_Gamma_Attenuation-Rate_CR_Attenuation))
    for b in range(0,Hist_El_Signal.GetNbinsX()):
        Hist_El_Signal.SetBinContent(b+1,max(0.,Hist_El_Signal.GetBinContent(b+1)))
    print Hist_El_Signal.Print("All")

    for b in range(0,Hist_CR_Signal.GetNbinsX()):
        Hist_CR_Signal.SetBinContent(b+1,Hist_OffData_Signal.GetBinContent(b+1))
        Hist_CR_Signal.SetBinError(b+1,Hist_OffData_Signal.GetBinError(b+1))
    Hist_CR_Signal.Add(Hist_El_Signal,-1.)
    for b in range(0,Hist_CR_Signal.GetNbinsX()):
        Hist_CR_Signal.SetBinContent(b+1,max(0.,Hist_CR_Signal.GetBinContent(b+1)))
    print Hist_CR_Signal.Print("All")

    Err_N_CR = ROOT.Double(0.)
    Err_N_El = ROOT.Double(0.)
    N_CR = Hist_CR_Signal.IntegralAndError(1,Hist_CR_Signal.GetNbinsX(),Err_N_CR)
    N_El = Hist_El_Signal.IntegralAndError(1,Hist_El_Signal.GetNbinsX(),Err_N_El)
    #N_CR = Hist_CR_Signal.IntegralAndError(Hist_CR_Signal.FindBin(-0.5),Hist_CR_Signal.FindBin(1.5),Err_N_CR)
    #N_El = Hist_El_Signal.IntegralAndError(Hist_El_Signal.FindBin(-0.5),Hist_El_Signal.FindBin(1.5),Err_N_El)
    Hist_Nel_temp = TH1D("Hist_Nel_temp","",1,0,1)
    Hist_Ncr_temp = TH1D("Hist_Ncr_temp","",1,0,1)
    Hist_e2p_temp = TH1D("Hist_e2p_temp","",1,0,1)
    Hist_Nel_temp.SetBinContent(1,N_El)
    Hist_Nel_temp.SetBinError(1,Err_N_El)
    Hist_Ncr_temp.SetBinContent(1,N_CR)
    Hist_Ncr_temp.SetBinError(1,Err_N_CR)
    Hist_e2p_temp = Hist_Nel_temp.Clone()
    Hist_e2p_temp.Divide(Hist_Ncr_temp)
    e2p_ratio_temp = Hist_e2p_temp.GetBinContent(1)
    err_e2p_ratio_temp = Hist_e2p_temp.GetBinError(1)
    #print 'OnData = %s'%(Hist_OnData_Signal.Integral(Hist_OnData_Signal.FindBin(-0.5),Hist_OnData_Signal.FindBin(1.5)))
    #print 'OffData = %s'%(Hist_OffData_Signal.Integral(Hist_OffData_Signal.FindBin(-0.5),Hist_OffData_Signal.FindBin(1.5)))
    print 'N_CR = %s +/- %s, N_EL = %s +/- %s'%(N_CR,Err_N_CR,N_El,Err_N_El)
    print 'N_EL/N_CR = %s +/- %s'%(e2p_ratio_temp,err_e2p_ratio_temp)

    Hist_OnData_Signal.Reset()
    Hist_OnData_Control.Reset()
    Hist_OffData_Signal.Reset()
    Hist_OffData_Control.Reset()
    Hist_Gamma_Signal.Reset()
    Hist_Gamma_Control.Reset()
    Hist_CR_Attenuation.Reset()
    Hist_Gamma_Attenuation.Reset()
    Hist_eGamma_Signal.Reset()
    Hist_CR_Signal.Reset()
    Hist_El_Signal.Reset()

    return Rate_CR_Attenuation,Err_Rate_CR_Attenuation,Rate_Gamma_Attenuation,Err_Rate_Gamma_Attenuation,e2p_ratio_temp,err_e2p_ratio_temp


#TelElev_Bin = [0,10,20,30,40,50,90]
TelElev_Bin = [0,10,20,30,40,50,60,70,80,90]
Hist_TelElev_RW_Factor = TH1D("Hist_TelElev_RW_Factor","",len(TelElev_Bin)-1,array('d',TelElev_Bin))
Hist_TelElev_Signal = TH1D("Hist_TelElev_Signal","",len(TelElev_Bin)-1,array('d',TelElev_Bin))
Hist_TelElev_Control = TH1D("Hist_TelElev_Control","",len(TelElev_Bin)-1,array('d',TelElev_Bin))
MSCL_Bin = [-10,0,1,2,3,4,5,100]
Hist_MSCL_RW_Factor = TH1D("Hist_MSCL_RW_Factor","",len(MSCL_Bin)-1,array('d',MSCL_Bin))
Hist_MSCL_Signal = TH1D("Hist_MSCL_Signal","",len(MSCL_Bin)-1,array('d',MSCL_Bin))
Hist_MSCL_Control = TH1D("Hist_MSCL_Control","",len(MSCL_Bin)-1,array('d',MSCL_Bin))
Erec_Bin = [0,0.5,1,2,3,5,10,100]
Hist_Erec_RW_Factor = TH1D("Hist_Erec_RW_Factor","",len(Erec_Bin)-1,array('d',Erec_Bin))
Hist_Erec_Signal = TH1D("Hist_Erec_Signal","",len(Erec_Bin)-1,array('d',Erec_Bin))
Hist_Erec_Control = TH1D("Hist_Erec_Control","",len(Erec_Bin)-1,array('d',Erec_Bin))

Hist2D_MSCW_vs_EmissionHeight_Data = TH2D("Hist2D_MSCW_vs_EmissionHeight_Data","",24,-4,20,40,0,40)
Hist2D_TelElev_vs_EmissionHeight_Data = TH2D("Hist2D_TelElev_vs_EmissionHeight_Data","",18,0,90,40,0,40)
Hist1D_TelElev_vs_EmissionHeight_Data = TH1D("Hist1D_TelElev_vs_EmissionHeight_Data","",18,0,90)
TProf_TelElev_vs_EmissionHeight_Data = TProfile("TProf_TelElev_vs_EmissionHeight_Data","",18,0,90,0,40)

nbins = 12
Hist_MSCW_CR = TH1D("Hist_MSCW_CR","",nbins,-5,19)
Hist_MSCW_CR_Raw = TH1D("Hist_MSCW_CR_Raw","",nbins,-5,19)
Hist_MSCW_Fake = TH1D("Hist_MSCW_Fake","",nbins,-5,19)
Hist_MSCW_ES = TH1D("Hist_MSCW_ES","",nbins,-5,19)
Hist_MSCW_ES_Raw = TH1D("Hist_MSCW_ES_Raw","",nbins,-5,19)
Hist_MSCW_Data = TH1D("Hist_MSCW_Data","",nbins,-5,19)
Hist_MSCW_Ring = TH1D("Hist_MSCW_Ring","",nbins,-5,19)
Hist_MSCW_Data_Sum = TH1D("Hist_MSCW_Data_Sum","",nbins,-5,19)
Hist_MSCW_Ring_Sum = TH1D("Hist_MSCW_Ring_Sum","",nbins,-5,19)
Hist_MSCW_Fake_Sum = TH1D("Hist_MSCW_Fake_Sum","",nbins,-5,19)
Hist_MSCW_ES_Sum = TH1D("Hist_MSCW_ES_Sum","",nbins,-5,19)
Hist_MSCW_CR_Sum = TH1D("Hist_MSCW_CR_Sum","",nbins,-5,19)

variable_bins = []
b = -1
while (b<2):
    variable_bins += [pow(10,b)]
    b += 0.5
Hist_Erec_CR = TH1D("Hist_Erec_CR","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_CR_Raw = TH1D("Hist_Erec_CR_Raw","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Fake = TH1D("Hist_Erec_Fake","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_ES = TH1D("Hist_Erec_ES","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_ES_Raw = TH1D("Hist_Erec_ES_Raw","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Data = TH1D("Hist_Erec_Data","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Ring = TH1D("Hist_Erec_Ring","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Data_Sum = TH1D("Hist_Erec_Data_Sum","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Ring_Sum = TH1D("Hist_Erec_Ring_Sum","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_Fake_Sum = TH1D("Hist_Erec_Fake_Sum","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_ES_Sum = TH1D("Hist_Erec_ES_Sum","",len(variable_bins)-1,array('d',variable_bins))
Hist_Erec_CR_Sum = TH1D("Hist_Erec_CR_Sum","",len(variable_bins)-1,array('d',variable_bins))

N_GammaData_Fake = 0
N_GammaMC_Fake = 0
N_GammaData_ES = 0
N_GammaMC_ES = 0

MSCW_Bin = [-100,3,100]  # the last bin has to be the CR normalization region!!
Hist_OnData_Signal = TH1D("Hist_OnData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OnData_Control = TH1D("Hist_OnData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Signal = TH1D("Hist_OffData_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_OffData_Control = TH1D("Hist_OffData_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Signal = TH1D("Hist_Gamma_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Control = TH1D("Hist_Gamma_Control","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Attenuation = TH1D("Hist_CR_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_Gamma_Attenuation = TH1D("Hist_Gamma_Attenuation","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_eGamma_Signal = TH1D("Hist_eGamma_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_CR_Signal = TH1D("Hist_CR_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))
Hist_El_Signal = TH1D("Hist_El_Signal","",len(MSCW_Bin)-1,array('d',MSCW_Bin))

#Elev_Bin = [75,90]
Elev_Bin = [50,55,60,65,70,75,80,85]
Hist_e2p_Ratio_LowE = TH1D("Hist_e2p_Ratio_LowE","",len(Elev_Bin)-1,array('d',Elev_Bin))
Hist_Gamma_Attenu_LowE = TH1D("Hist_Gamma_Attenu_LowE","",len(Elev_Bin)-1,array('d',Elev_Bin))
Hist_CR_Attenu_LowE = TH1D("Hist_CR_Attenu_LowE","",len(Elev_Bin)-1,array('d',Elev_Bin))

Hist_TelElev_Counts = TH1D("Hist_TelElev_Counts","",45,0,90)

variable_bins = []
b = -1
while (b<2):
    variable_bins += [pow(10,b)]
    b += 0.2
Hist2D_Erec_vs_EmissionHeight_Data = TH2D("Hist2D_Erec_vs_EmissionHeight_Data","",len(variable_bins)-1,array('d',variable_bins),40,0,40)

source = 'Crab'
#source = 'PKS1424'
#source = '3C264'
#field = 'off'
field = 'on'

RW_method = ''
Fake_Removal = ''
Add_Electron = ''
Height_cut = ''
Elev_cut_lower = 0
Elev_cut_upper = 90
MSCL_cut_lower = -2
MSCL_cut_upper = 1
Erec_cut_lower = 0.1
Erec_cut_upper = 100

energy_bins = [0.1,0.3,1.0,1000]
#gamma_att = [0.094,0.084,0.092]
#cr_att = [0.81,0.52,0.36]
#e2p_ratio = [0.68,0.36,0.07]

Fake_cut_MSCL = 0.
Fake_cut_MSCW = 0.

tag = ''

#Do_e2p_Measurement = False
Do_e2p_Measurement = True
Hist_e2p_Ratio = []
Hist_Gamma_Attenu = []
Hist_CR_Attenu = []
if Do_e2p_Measurement:
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
    for energy in range(0,len(energy_bins)-1):
        Hist_e2p_Ratio += [Hist_e2p_Ratio_LowE.Clone()]
        Hist_Gamma_Attenu += [Hist_Gamma_Attenu_LowE.Clone()]
        Hist_CR_Attenu += [Hist_CR_Attenu_LowE.Clone()]
        for elev in range(1,Hist_e2p_Ratio_LowE.GetNbinsX()+1):
            Elev_cut_lower = Hist_e2p_Ratio_LowE.GetBinLowEdge(elev)
            Elev_cut_upper = Hist_e2p_Ratio_LowE.GetBinLowEdge(elev+1)
            Erec_cut_lower = energy_bins[energy]
            Erec_cut_upper = energy_bins[energy+1]
            #Erec_cut_upper = 1000.
            MakeATag()
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print 'selection = %s'%(tag)
            CR_attenu,err_CR_attenu,gamma_attenu,err_gamma_attenu,e2p_ratio,err_e2p_ratio = MeasureElectronComponent()
            Hist_e2p_Ratio[energy].SetBinContent(Hist_e2p_Ratio_LowE.FindBin(Elev_cut_lower),e2p_ratio)
            Hist_e2p_Ratio[energy].SetBinError(Hist_e2p_Ratio_LowE.FindBin(Elev_cut_lower),err_e2p_ratio)
            Hist_Gamma_Attenu[energy].SetBinContent(Hist_Gamma_Attenu_LowE.FindBin(Elev_cut_lower),gamma_attenu)
            Hist_Gamma_Attenu[energy].SetBinError(Hist_Gamma_Attenu_LowE.FindBin(Elev_cut_lower),err_gamma_attenu)
            Hist_CR_Attenu[energy].SetBinContent(Hist_CR_Attenu_LowE.FindBin(Elev_cut_lower),CR_attenu)
            Hist_CR_Attenu[energy].SetBinError(Hist_CR_Attenu_LowE.FindBin(Elev_cut_lower),err_CR_attenu)
    pad1.cd()
    #f_LowE = ROOT.TF1("f_LowE","[0]",30,90)
    #Hist_e2p_Ratio_LowE.Fit("f_LowE","R")
    #f_LowE.SetLineColor(4)
    Hist_e2p_Ratio[0].GetXaxis().SetTitle('Tel Elevation')
    Hist_e2p_Ratio[0].GetYaxis().SetTitle('e/p ratio')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_e2p_Ratio[energy].GetMaximum(),maxContent)
    Hist_e2p_Ratio[0].SetMaximum(1.2*maxContent)
    Hist_e2p_Ratio[0].SetMinimum(0)
    Hist_e2p_Ratio[0].Draw("E TEXT")
    for energy in range(0,len(energy_bins)-1):
        Hist_e2p_Ratio[energy].SetLineColor(energy+2)
        Hist_e2p_Ratio[energy].Draw("E TEXT same")
    #f_LowE.Draw("same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_e2p_Ratio[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/e2p_Ratio_%s.pdf'%(source))
    pad1.cd()
    #f_LowE = ROOT.TF1("f_LowE","[0]",30,90)
    #Hist_Gamma_Attenu_LowE.Fit("f_LowE","R")
    #f_LowE.SetLineColor(4)
    Hist_Gamma_Attenu[0].GetXaxis().SetTitle('Tel Elevation')
    Hist_Gamma_Attenu[0].GetYaxis().SetTitle('#gamma attenuation')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_Gamma_Attenu[energy].GetMaximum(),maxContent)
    Hist_Gamma_Attenu[0].SetMaximum(1.2*maxContent)
    Hist_Gamma_Attenu[0].SetMinimum(0)
    Hist_Gamma_Attenu[0].Draw("E TEXT")
    for energy in range(0,len(energy_bins)-1):
        Hist_Gamma_Attenu[energy].SetLineColor(energy+2)
        Hist_Gamma_Attenu[energy].Draw("E TEXT same")
    #f_LowE.Draw("same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_Gamma_Attenu[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/Gamma_Attenu_%s.pdf'%(source))
    pad1.cd()
    #f_LowE = ROOT.TF1("f_LowE","[0]",30,90)
    #Hist_CR_Attenu_LowE.Fit("f_LowE","R")
    #f_LowE.SetLineColor(4)
    Hist_CR_Attenu[0].GetXaxis().SetTitle('Tel Elevation')
    Hist_CR_Attenu[0].GetYaxis().SetTitle('CR attenuation')
    maxContent = 0.
    for energy in range(0,len(energy_bins)-1):
        maxContent = max(Hist_CR_Attenu[energy].GetMaximum(),maxContent)
    Hist_CR_Attenu[0].SetMaximum(1.2*maxContent)
    Hist_CR_Attenu[0].SetMinimum(0)
    Hist_CR_Attenu[0].Draw("E TEXT")
    for energy in range(0,len(energy_bins)-1):
        Hist_CR_Attenu[energy].SetLineColor(energy+2)
        Hist_CR_Attenu[energy].Draw("E TEXT same")
    #f_LowE.Draw("same")
    pad2.cd()
    legend.Clear()
    for energy in range(0,len(energy_bins)-1):
        legend.AddEntry(Hist_CR_Attenu[energy],'E %s-%s TeV'%(energy_bins[energy],energy_bins[energy+1]),"pl")
    legend.Draw("SAME")
    canvas2.SaveAs('output/CR_Attenu_%s.pdf'%(source))

RW_method = ''
Fake_Removal = ''
Add_Electron = ''
Height_cut = ''
#RW_method = 'NoRW'
#Fake_Removal = 'NoFakeRemoved'
#Add_Electron = 'NoEle'
#Height_cut = 'simple'

source = 'Crab'
#source = 'PKS1424'
#source = '3C264'
#field = 'off'
field = 'on'

Elev_cut_lower = 60
Elev_cut_upper = 90
#MSCL_cut_lower = -2
#MSCL_cut_upper = 1

#Do_analysis = False
Do_analysis = True
if Do_analysis:
    for energy in range(0,len(energy_bins)-1):
        Erec_cut_lower = energy_bins[energy]
        Erec_cut_upper = energy_bins[energy+1]
        MakeATag()
        #RunExtendedSourceAnalysis(e2p_ratio[energy],cr_att[energy],gamma_att[energy])
        e2p_ratio = Hist_e2p_Ratio[energy].GetBinContent(Hist_e2p_Ratio[energy].FindBin(70.))
        cr_att = Hist_CR_Attenu[energy].GetBinContent(Hist_CR_Attenu[energy].FindBin(70.))
        gamma_att = Hist_Gamma_Attenu[energy].GetBinContent(Hist_Gamma_Attenu[energy].FindBin(70.))
        RunExtendedSourceAnalysis(e2p_ratio,cr_att,gamma_att)
    
    Erec_cut_lower = energy_bins[0]
    Erec_cut_upper = energy_bins[len(energy_bins)-1]
    MakeATag()
    MakeStackPlot(Hist_MSCW_Data_Sum,Hist_MSCW_CR_Sum,Hist_MSCW_ES_Sum,Hist_MSCW_Fake_Sum,Hist_MSCW_Ring_Sum,'MSCW')
    MakeStackPlot(Hist_Erec_Data_Sum,Hist_Erec_CR_Sum,Hist_Erec_ES_Sum,Hist_Erec_Fake_Sum,Hist_Erec_Ring_Sum,'Erec')

canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
Hist_TelElev_Counts.GetXaxis().SetTitle('Tel. elevation')
Hist_TelElev_Counts.Draw("E")
canvas.SaveAs('output/TelElev_Counts_%s.pdf'%(tag))

