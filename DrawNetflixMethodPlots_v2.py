
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")



FileFolder = []
FileTag = []
#FileFolder += ['output_unblind_4x4_dNSBp0_tight']
#FileTag += ['unblind_4x4_dNSBp0_tight']
#FileFolder += ['output_unblind_4x4_dNSBp0_medium']
#FileTag += ['unblind_4x4_dNSBp0_medium']
FileFolder += ['output_unblind_4x4_dNSBp0_loose']
FileTag += ['unblind_4x4_dNSBp0_loose']
#FileFolder += ['output_unblind_4x4_dNSBp3_loose']
#FileTag += ['unblind_4x4_dNSBp3_loose']

elev_range = []
elev_range += [[75,85]]
elev_range += [[65,75]]
elev_range += [[55,65]]
elev_range += [[45,55]]
#elev_range += [[35,45]]
#elev_range += [[25,35]]

source_list = []
sky_coord = []
source_list += ['Crab']
sky_coord += ['05 34 31.97 +22 00 52.1']
source_list += ['Mrk421']
sky_coord += ['11 04 19 +38 11 41']
source_list += ['H1426']
sky_coord += ['14 28 32.609 +42 40 21.05']
source_list += ['PKS1424']
sky_coord += ['14 27 00 +23 47 00']
source_list += ['3C264']
sky_coord += ['11 45 5.009 +19 36 22.74']
source_list += ['OJ287V6']
sky_coord += ['08 54 49.1 +20 05 58.89']
source_list += ['1ES0229']
sky_coord += ['02 32 53.2 +20 16 21']
source_list += ['S3_1227_V6']
sky_coord += ['12 30 14.1 +25 18 07']
source_list += ['MS1221V6']
sky_coord += ['12 24 24.2 +24 36 24']
source_list += ['PKS1441V6']
sky_coord += ['14 43 56.9 +25 01 44']
source_list += ['RBS0413V6']
sky_coord += ['03 19 47 +18 45 42']
source_list += ['PG1553V6']
sky_coord += ['15 55 44.7 +11 11 41']
source_list += ['Segue1V6']
sky_coord += ['10 07 04 +16 04 55']
source_list += ['1ES1011V6']
sky_coord += ['10 15 4.139 +49 26 0.71']
source_list += ['NGC1275V6']
sky_coord += ['03 19 48.1 +41 30 42']
source_list += ['1ES0647V6']
sky_coord += ['06 50 46.490 +25 02 59.62']
source_list += ['1ES1440V6']
sky_coord += ['14 42 48.277 +12 00 40.37']
source_list += ['1ES1741V6']
sky_coord += ['17 44 01.2 +19 32 47']
source_list += ['IC443HotSpot']
sky_coord += ['06 18 2.700 +22 39 36.00']
source_list += ['RGBJ0710']
sky_coord += ['07 10 26.4 +59 09 00']
source_list += ['CasA']
sky_coord += ['23 23 13.8 +58 48 26']
source_list += ['WComaeV6']
sky_coord += ['12 21 31.7 +28 13 59']
source_list += ['M82']
sky_coord += ['09 55 52.7 +69 40 46']
#source_list += ['G079']
#sky_coord += ['20 32 28.56 +40 19 41.52']
#source_list += ['1ES1218V6']
#sky_coord += ['12 21 26.3 +30 11 29']
#source_list += ['MGRO_J1908_V6']
#sky_coord += ['19 07 54 +06 16 07']
#source_list += ['MGRO_J1908_V5']
#sky_coord += ['19 07 54 +06 16 07']
#source_list += ['Segue1V5']
#sky_coord += ['10 07 04 +16 04 55']
#source_list += ['IC443HotSpotV5']
#sky_coord += ['06 18 2.700 +22 39 36.00']

energy_list = []
energy_list += [200]
energy_list += [237]
energy_list += [282]
energy_list += [335]
energy_list += [398] # works from here
energy_list += [473]
energy_list += [562]
energy_list += [794]
energy_list += [1122]
energy_list += [2239]
energy_list += [4467]
#energy_list += [8913]

exposure_hours = 0.
exposure_hours_dark = 0.
NSB_avg = 0.
NSB_avg_dark = 0.
ErecS_lower_cut = 0
ErecS_upper_cut = 0

N_bins_for_deconv = 40
MSCW_plot_lower = -1.
MSCW_plot_upper = 3.
MSCL_plot_lower = -1.
MSCL_plot_upper = 3.

MSCW_lower_cut = -1.0
MSCW_upper_cut = 1.0
MSCW_blind_cut = 1.0
MSCL_lower_cut = -1.0
MSCL_upper_cut = 1.0
MSCL_blind_cut = 1.0

Hist2D_Data = ROOT.TH2D("Hist2D_Data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Dark = ROOT.TH2D("Hist2D_Dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Bkgd = ROOT.TH2D("Hist2D_Bkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_GammaDark = ROOT.TH2D("Hist2D_GammaDark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_GammaRDBM = ROOT.TH2D("Hist2D_GammaRDBM","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_TrueBkgd = ROOT.TH2D("Hist2D_TrueBkgd","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaDark_MSCW = ROOT.TH1D("Hist_GammaDark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaRDBM_MSCW = ROOT.TH1D("Hist_GammaRDBM_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_MSCW = ROOT.TH1D("Hist_Data_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Bkgd_MSCW = ROOT.TH1D("Hist_Bkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_TrueBkgd_MSCW = ROOT.TH1D("Hist_TrueBkgd_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Dark_MSCW = ROOT.TH1D("Hist_Dark_MSCW","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaDark_MSCL = ROOT.TH1D("Hist_GammaDark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_GammaRDBM_MSCL = ROOT.TH1D("Hist_GammaRDBM_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Data_MSCL = ROOT.TH1D("Hist_Data_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Bkgd_MSCL = ROOT.TH1D("Hist_Bkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_TrueBkgd_MSCL = ROOT.TH1D("Hist_TrueBkgd_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Dark_MSCL = ROOT.TH1D("Hist_Dark_MSCL","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Data_Theta2 = ROOT.TH1D("Hist_Data_Theta2","",64,0,10)
Hist_TrueBkgd_Theta2 = ROOT.TH1D("Hist_TrueBkgd_Theta2","",64,0,10)
Hist_Dark_Theta2 = ROOT.TH1D("Hist_Dark_Theta2","",64,0,10)
Hist_Bkgd_Theta2 = ROOT.TH1D("Hist_Bkgd_Theta2","",64,0,10)
Hist_Bkgd_Theta2_Raw = ROOT.TH1D("Hist_Bkgd_Theta2_Raw","",64,0,10)

Hist2D_Data_SumE = ROOT.TH2D("Hist2D_Data_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Dark_SumE = ROOT.TH2D("Hist2D_Dark_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_Bkgd_SumE = ROOT.TH2D("Hist2D_Bkgd_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_GammaDark_SumE = ROOT.TH2D("Hist2D_GammaDark_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_GammaRDBM_SumE = ROOT.TH2D("Hist2D_GammaRDBM_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist2D_TrueBkgd_SumE = ROOT.TH2D("Hist2D_TrueBkgd_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaDark_MSCW_SumE = ROOT.TH1D("Hist_GammaDark_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaRDBM_MSCW_SumE = ROOT.TH1D("Hist_GammaRDBM_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_MSCW_SumE = ROOT.TH1D("Hist_Data_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Bkgd_MSCW_SumE = ROOT.TH1D("Hist_Bkgd_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_TrueBkgd_MSCW_SumE = ROOT.TH1D("Hist_TrueBkgd_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Dark_MSCW_SumE = ROOT.TH1D("Hist_Dark_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_GammaDark_MSCL_SumE = ROOT.TH1D("Hist_GammaDark_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_GammaRDBM_MSCL_SumE = ROOT.TH1D("Hist_GammaRDBM_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Data_MSCL_SumE = ROOT.TH1D("Hist_Data_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Bkgd_MSCL_SumE = ROOT.TH1D("Hist_Bkgd_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_TrueBkgd_MSCL_SumE = ROOT.TH1D("Hist_TrueBkgd_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Dark_MSCL_SumE = ROOT.TH1D("Hist_Dark_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
Hist_Data_Theta2_SumE = ROOT.TH1D("Hist_Data_Theta2_SumE","",64,0,10)
Hist_Dark_Theta2_SumE = ROOT.TH1D("Hist_Dark_Theta2_SumE","",64,0,10)
Hist_TrueBkgd_Theta2_SumE = ROOT.TH1D("Hist_TrueBkgd_Theta2_SumE","",64,0,10)
Hist_Bkgd_Theta2_SumE = ROOT.TH1D("Hist_Bkgd_Theta2_SumE","",64,0,10)
Hist_Bkgd_Theta2_Raw_SumE = ROOT.TH1D("Hist_Bkgd_Theta2_Raw_SumE","",64,0,10)

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def CalculateSignificance(s,b,err):
    if (b*b+(s+b)*err*err)==0.: return 0.
    if (s+b)*(b+err*err)==0.: return 0.
    first_term = (s+b)*math.log((s+b)*(b+err*err)/(b*b+(s+b)*err*err))
    if err>0. and b>0:
        second_term = b*b/(err*err)*math.log(1.+err*err*s/(b*(b+err*err)))
    else: 
        second_term = 0.
    result = 0.
    if first_term>second_term: result = pow(2*(first_term-second_term),0.5)
    else: result = pow(2*(-first_term+second_term),0.5)
    if s>0: return result
    else: return -1.*result

def IntegralAndError(Hist,bin1,bin2):
    
    integral = 0
    error = 0
    for b in range(1,Hist.GetNbinsX()+1):
        if b<bin1: continue
        if b>bin2: continue
        integral += Hist.GetBinContent(b)
        error += pow(Hist.GetBinError(b),2)
    error = pow(error,0.5)
    if math.isnan(integral) or math.isnan(error):
        integral = 0
        error = 0
    return integral, error

def IntegralAndSystError(Hist,bin1,bin2,syst):
    
    integral = 0
    error = 0
    sum_syst_err = 0.
    sum_stat_err = 0.
    for b in range(1,Hist.GetNbinsX()+1):
        if b<bin1: continue
        if b>bin2: continue
        integral += Hist.GetBinContent(b)
        total_err = Hist.GetBinError(b)
        syst_err = Hist.GetBinContent(b)*syst
        stat_err = 0.
        if syst<0:  # 100% linear syst
            syst_err = total_err
        if not total_err*total_err-syst_err*syst_err<0:
            stat_err = pow(total_err*total_err-syst_err*syst_err,0.5)
        else:
            stat_err = 0.
        sum_syst_err += syst_err
        sum_stat_err += stat_err*stat_err
    if syst<0:
        sum_stat_err = 0
        error = sum_syst_err
    else:
        sum_stat_err = pow(sum_stat_err,0.5)
        error = pow(sum_syst_err*sum_syst_err+sum_stat_err*sum_stat_err,0.5)
    return integral, error

def Variation_ratio(Hist_SR, Hist_Bkg,range_lower,range_upper):

    #range_lower = 0.4

    norm_bin_low_target = Hist_SR.FindBin(range_lower)
    norm_bin_up_target = Hist_Bkg.FindBin(range_upper)-1
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hist_SR,norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target)
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0 and not data_SR==0:
        sbratio = (data_SR)/(predict_bkg)
        sbratio_err = (data_SR)/(predict_bkg)*pow(pow(err_SR/data_SR,2)+pow(err_bkg/predict_bkg,2),0.5)
    return sbratio, sbratio_err

def MakeChi2Plot(Hists,legends,colors,title,name,doSum,doNorm,range_lower,range_upper,syst):
    
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
    pad1.SetGrid()
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    norm_bin_low_target = Hists[0].FindBin(range_lower)
    norm_bin_up_target = Hists[0].FindBin(range_upper)-1
    if 'fine_theta2' in name:
        norm_bin_low_target = Hists[0].FindBin(0)
        norm_bin_up_target = Hists[0].FindBin(0.5)-1
    if norm_bin_up_target<norm_bin_low_target: norm_bin_up_target+=1

    if doNorm:
        for h in range(1,len(Hists)):
            bin_low = Hists[0].FindBin(1.)
            bin_up = Hists[0].FindBin(3.)-1
            scale = Hists[0].Integral(bin_low,bin_up)/Hists[h].Integral(bin_low,bin_up)
            Hists[h].Scale(scale)

    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    bkg_idx = 1
    for h in range(0,len(Hists)):
        mean += [Hists[h].GetMean()]
        rms += [Hists[h].GetRMS()]
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            if legends[h]=='predict. bkg.':
                bkg_idx = h
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #if not 'Energy' in name:
    #    Hists[max_hist].SetMinimum(0)
    #if not 'ZoomIn' in name:
    #    low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
    #    high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
    #    low_end = -1.
    #    high_end = 3.
    #    Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    #else:
    #    low_end = range_lower
    #    high_end = range_upper
    #    Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    #Hists[0].SetMinimum(0)
    #Hists[0].Draw("E")
    Hists[max_hist].Draw("E")

    fill_color = [0,0,0,30,38,46]
    if doSum:
        stack = ROOT.THStack("stack", "")
        Hist_Sum = Hists[1].Clone()
        Hist_Sum.Reset()
        Hist_Sys = Hists[1].Clone()
        Hist_Sys.Reset()
        for h in range(0,len(Hists)):
            if colors[h]==1: continue
            if colors[h]==2: continue
            if not legends[h]=='syst.':
                set_histStyle( Hists[h] , fill_color[colors[h]])
                stack.Add( Hists[h] )
                Hist_Sum.Add( Hists[h] )
            else:
                Hist_Sys.Add( Hists[h] )
        stack.Draw("hist same")
        Hist_Err = Hist_Sum.Clone()
        for binx in range(0,Hist_Err.GetNbinsX()):
            old_err = Hist_Err.GetBinError(binx+1)
            new_err = Hist_Sys.GetBinContent(binx+1)
            Hist_Err.SetBinError(binx+1,pow(old_err*old_err+new_err*new_err,0.5))
        Hist_Err.SetFillColor(1)
        Hist_Err.SetFillStyle(3004)
        Hist_Err.SetMarkerSize(0)
        Hist_Err.Draw("e2 same")

    #Hist_Excess = Hists[0].Clone()
    #Hist_Excess.Add(Hists[1],-1.)
    #Hist_Excess.SetFillColor(0)
    #Hist_Excess.SetFillStyle(0)
    #Hist_Excess.SetLineColor(2)
    #Hist_Excess.Draw("e2 same")

    line1 = ROOT.TLine(Hists[max_hist].GetBinLowEdge(norm_bin_up_target+1),Hists[max_hist].GetMaximum()*2,Hists[max_hist].GetBinLowEdge(norm_bin_up_target+1),0)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()


    for h in range(0,len(Hists)):
        if colors[h]==1 or colors[h]==2:
            Hists[h].SetLineWidth(3)
            Hists[h].Draw("E same")
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    #legend.SetTextSize(0.2)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hists[0],legends[0],"pl")
    for h in range(1,len(Hists)):
        if Hists[h]!=0:
            if 'NoElec' in name and legends[h]=='electron': continue
            legend.AddEntry(Hists[h],legends[h],"f")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.80,'E >%0.1f GeV (%.1f hrs)'%(ErecS_lower_cut,exposure_hours) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hists[0],norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hists[bkg_idx],norm_bin_low_target,norm_bin_up_target)
    predict_sys, err_sys = IntegralAndSystError(Hist_Err,norm_bin_low_target,norm_bin_up_target,-1)
    err_bkg = pow(err_bkg*err_bkg+err_sys*err_sys,0.5)
    Sig = 1.*CalculateSignificance(data_SR-predict_bkg,predict_bkg,err_bkg)
    lumilab2 = ROOT.TLatex(0.15,0.60,'Excess = %0.1f#pm%0.1f (%0.1f#sigma)'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5),Sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.40,'Bkg = %0.1f#pm%0.1f'%(predict_bkg,err_bkg) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0: 
        sbratio = (data_SR-predict_bkg)/(predict_bkg)
    if not data_SR-predict_bkg==0 and not predict_bkg==0:
        sbratio_err = sbratio*pow(pow(pow(err_SR*err_SR+err_bkg*err_bkg,0.5)/(data_SR-predict_bkg),2)+pow(err_bkg/predict_bkg,2),0.5)
    else: sbratio_err = 0
    lumilab4 = ROOT.TLatex(0.15,0.20,'S/B = %0.3f#pm%0.3f'%(sbratio,sbratio_err) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()

    pad2.cd()
    Hist_Band = Hist_Sum.Clone()
    #Hist_Band.Divide(Hists[0])
    for b in range(0,Hist_Band.GetNbinsX()):
        #Hist_Band.SetBinContent(b+1,0)
        #Hist_Band.SetBinError(b+1,3)
        Hist_Band.SetBinContent(b+1,1)
        Hist_Band.SetBinError(b+1,0.1)
    Hist_Band.SetFillColor(1)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    #Hist_Band.GetYaxis().SetTitle("#sqrt{#chi^{2}} (per bin)")
    Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    Hist_Band.GetYaxis().SetTitleSize(0.10)
    Hist_Band.GetYaxis().SetNdivisions(505)
    #Hist_Band.SetMaximum(5)
    #Hist_Band.SetMinimum(-5)
    Hist_Band.SetMaximum(1.2)
    Hist_Band.SetMinimum(0.8)
    #Hist_Band.GetXaxis().SetRangeUser(low_end,high_end)
    Hist_Band.Draw("e2")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    Hist_Ratio = []
    Hist_Bkg = Hists[1].Clone()
    for h in range(0,len(Hists)):
        if colors[h]==4:
            Hist_Bkg = Hists[h].Clone()
    for h in range(0,len(Hists)):
        if colors[h]==1:
            Hist_Ratio += [Hists[h].Clone()]
            Hist_Ratio[h].Divide(Hist_Bkg)
            Hist_Ratio[h].SetLineWidth(2)
            Hist_Ratio[h].Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeComparisonPlotTwoColumn(Hists,legends,colors,title,name,minheight,maxheight,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.2)
    pad1.SetTopMargin(0.0)
    pad1.SetLeftMargin(0.2)
    pad1.SetBorderMode(0)
    if logy: pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,len(Hists)):
        mean += [0]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].GetXaxis().SetTitle(title)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    if 'MSCW' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        for h in range(0,len(Hists)):
            Hists[h].GetXaxis().SetRangeUser(low_end,high_end)

    if not logy: Hists[max_hist].SetMinimum(0)
    if not maxheight==0: Hists[max_hist].SetMaximum(maxheight)
    Hists[max_hist].SetMinimum(minheight)
    Hists[max_hist].GetXaxis().SetTitle(title)
    Hists[max_hist].GetXaxis().SetTitleOffset(0.6)
    Hists[max_hist].GetXaxis().SetTitleSize(0.06)
    Hists[max_hist].GetXaxis().SetLabelSize(0.06)
    Hists[max_hist].GetYaxis().SetLabelSize(0.06)
    Hists[max_hist].GetYaxis().SetTitleOffset(1.0)
    Hists[max_hist].GetYaxis().SetTitle("R = true / predict")
    Hists[max_hist].GetYaxis().SetTitleSize(0.06)
    Hists[max_hist].Draw("E")

    line1 = ROOT.TLine(Hists[max_hist].GetBinLowEdge(1),1,Hists[max_hist].GetBinLowEdge(Hists[max_hist].GetNbinsX()+1),1)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw("same")

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].SetMarkerStyle(22)
            Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.9,0.9)
    legend.SetNColumns(2)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    pad1.cd()
    syst = []
    syst_err = []
    for binx in range(0,Hists[0].GetNbinsX()):
        weight_this_bin = 0.
        raw_weight_this_bin = 0.
        syst_this_bin = 0.
        syst_err_this_bin = 0.
        for h in range(0,len(Hists)):
            #weight = 1./pow(Hists[h].GetBinError(binx+1),2)
            weight = 1.
            syst_this_bin += weight*pow(1.-Hists[h].GetBinContent(binx+1),2)
            weight_this_bin += weight
            raw_weight_this_bin += 1.
        syst_this_bin = syst_this_bin/weight_this_bin
        syst_this_bin = pow(syst_this_bin,0.5)
        syst += [syst_this_bin]
        syst_err += [syst_this_bin/pow(raw_weight_this_bin,0.5)]

    lumilab1 = ROOT.TLatex(0.25,0.9,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[0],100.*syst_err[0]) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.05)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.60,0.9,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[1],100.*syst_err[1]) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.05)
    lumilab2.Draw()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeComparisonPlotSigDist(Hists,legends,colors,title,name,minheight,maxheight,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.4,1,0.8)
    pad1.SetBottomMargin(0.1)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.4)
    pad2.SetBottomMargin(0.2)
    pad2.SetTopMargin(0.1)
    pad2.SetBorderMode(0)
    if logy: pad1.SetGrid()
    pad2.Draw()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()

    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,len(Hists)):
        mean += [0]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].GetXaxis().SetTitle(title)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    if 'MSCW' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        for h in range(0,len(Hists)):
            Hists[h].GetXaxis().SetRangeUser(low_end,high_end)

    if not logy: Hists[max_hist].SetMinimum(0)
    if not maxheight==0: Hists[max_hist].SetMaximum(maxheight)
    Hists[max_hist].SetMinimum(minheight)
    Hists[max_hist].GetXaxis().SetTitle(title)
    Hists[max_hist].GetXaxis().SetTitleOffset(0.6)
    Hists[max_hist].GetXaxis().SetTitleSize(0.06)
    Hists[max_hist].GetXaxis().SetLabelSize(0.06)
    Hists[max_hist].GetYaxis().SetLabelSize(0.06)
    Hists[max_hist].GetYaxis().SetTitleOffset(0.5)
    Hists[max_hist].GetYaxis().SetTitle("R = true / predict")
    Hists[max_hist].GetYaxis().SetTitleSize(0.08)
    Hists[max_hist].Draw("E")

    line1 = ROOT.TLine(Hists[max_hist].GetBinLowEdge(1),1,Hists[max_hist].GetBinLowEdge(Hists[max_hist].GetNbinsX()+1),1)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw("same")

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].SetMarkerStyle(22)
            Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.9,0.9)
    legend.SetNColumns(2)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    pad2.cd()
    shift = []
    syst = []
    syst_err = []
    for binx in range(0,Hists[0].GetNbinsX()):
        weight_this_bin = 0.
        raw_weight_this_bin = 0.
        shift_this_bin = 0.
        syst_this_bin = 0.
        syst_err_this_bin = 0.
        for h in range(0,len(Hists)):
            if Hists[h].GetBinContent(binx+1)==0.: continue
            #weight = 1./pow(Hists[h].GetBinError(binx+1),2)
            weight = 1.
            shift_this_bin += weight*(1.-Hists[h].GetBinContent(binx+1))
            syst_this_bin += weight*pow(1.-Hists[h].GetBinContent(binx+1),2)
            weight_this_bin += weight
            raw_weight_this_bin += 1.
        if weight_this_bin==0.:
            shift_this_bin = 1.
            syst_this_bin = 1.
            syst_err_this_bin = 0.
        else:
            shift_this_bin = shift_this_bin/weight_this_bin
            syst_this_bin = syst_this_bin/weight_this_bin
            syst_err_this_bin = syst_this_bin/pow(raw_weight_this_bin,0.5)
        syst_this_bin = pow(syst_this_bin,0.5)
        shift += [shift_this_bin]
        syst += [syst_this_bin]
        syst_err += [syst_err_this_bin]

    for binx in range(0,Hists[0].GetNbinsX()):
        print 'MDM_correction += [%0.3f]'%(shift[binx])

    Hist_Syst = Hists[0].Clone()
    for b in range(0,Hists[0].GetNbinsX()):
        Hist_Syst.SetBinContent(b+1,100.*syst[b])
        Hist_Syst.SetBinError(b+1,100.*syst_err[b])
    Hist_Syst.SetMaximum(15.)
    Hist_Syst.SetMinimum(0.)
    Hist_Syst.GetXaxis().SetLabelSize(0.06)
    Hist_Syst.GetYaxis().SetLabelSize(0.06)
    Hist_Syst.GetXaxis().SetTitleSize(0.06)
    Hist_Syst.GetXaxis().SetTitle(title)
    Hist_Syst.GetYaxis().SetTitleOffset(0.5)
    Hist_Syst.GetYaxis().SetTitleSize(0.08)
    Hist_Syst.GetYaxis().SetTitle('Syst. Unc. %')
    Hist_Syst.GetXaxis().SetTitleOffset(1.0)
    Hist_Syst.Draw("E")

    if logx: 
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def GetSourceInfo(file_list):

    global MSCW_blind_cut
    global MSCL_blind_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark
    exposure_hours = 0.
    exposure_hours_dark = 0.
    NSB_avg = 0.
    NSB_avg_dark = 0.
    for path in range(1,len(file_list)):
        if not os.path.isfile(file_list[path]):continue
        InputFile = ROOT.TFile(file_list[path])
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        exposure_hours += InfoTree.exposure_hours
        exposure_hours_dark += InfoTree.exposure_hours_dark
        NSB_avg += InfoTree.exposure_hours*InfoTree.NSB
        NSB_avg_dark += InfoTree.exposure_hours_dark*InfoTree.NSB_dark
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind
        InputFile.Close()
    NSB_avg = NSB_avg/exposure_hours
    NSB_avg_dark = NSB_avg_dark/exposure_hours_dark

    return exposure_hours, NSB_avg, exposure_hours_dark, NSB_avg_dark

def ResetStackedShowerHistograms():

    Hist2D_Data_SumE.Reset()
    Hist2D_Dark_SumE.Reset()
    Hist2D_Bkgd_SumE.Reset()
    Hist2D_TrueBkgd_SumE.Reset()
    Hist2D_GammaDark_SumE.Reset()
    Hist2D_GammaRDBM_SumE.Reset()

    Hist_Data_MSCL_SumE.Reset()
    Hist_Dark_MSCL_SumE.Reset()
    Hist_Bkgd_MSCL_SumE.Reset()
    Hist_TrueBkgd_MSCL_SumE.Reset()
    Hist_GammaDark_MSCL_SumE.Reset()
    Hist_GammaRDBM_MSCL_SumE.Reset()

    Hist_Data_MSCW_SumE.Reset()
    Hist_Dark_MSCW_SumE.Reset()
    Hist_Bkgd_MSCW_SumE.Reset()
    Hist_TrueBkgd_MSCW_SumE.Reset()
    Hist_GammaDark_MSCW_SumE.Reset()
    Hist_GammaRDBM_MSCW_SumE.Reset()

    Hist_Data_Theta2_SumE.Reset()
    Hist_TrueBkgd_Theta2_SumE.Reset()
    Hist_Dark_Theta2_SumE.Reset()
    Hist_Bkgd_Theta2_SumE.Reset()
    Hist_Bkgd_Theta2_Raw_SumE.Reset()

def StackShowerHistograms():

    Hist2D_Data_SumE.Add(Hist2D_Data)
    Hist2D_Dark_SumE.Add(Hist2D_Dark)
    Hist2D_Bkgd_SumE.Add(Hist2D_Bkgd)
    Hist2D_TrueBkgd_SumE.Add(Hist2D_TrueBkgd)
    Hist2D_GammaDark_SumE.Add(Hist2D_GammaDark)
    Hist2D_GammaRDBM_SumE.Add(Hist2D_GammaRDBM)

    Hist_Data_MSCL_SumE.Add(Hist_Data_MSCL)
    Hist_Dark_MSCL_SumE.Add(Hist_Dark_MSCL)
    Hist_Bkgd_MSCL_SumE.Add(Hist_Bkgd_MSCL)
    Hist_TrueBkgd_MSCL_SumE.Add(Hist_TrueBkgd_MSCL)
    Hist_GammaDark_MSCL_SumE.Add(Hist_GammaDark_MSCL)
    Hist_GammaRDBM_MSCL_SumE.Add(Hist_GammaRDBM_MSCL)

    Hist_Data_MSCW_SumE.Add(Hist_Data_MSCW)
    Hist_Dark_MSCW_SumE.Add(Hist_Dark_MSCW)
    Hist_Bkgd_MSCW_SumE.Add(Hist_Bkgd_MSCW)
    Hist_TrueBkgd_MSCW_SumE.Add(Hist_TrueBkgd_MSCW)
    Hist_GammaDark_MSCW_SumE.Add(Hist_GammaDark_MSCW)
    Hist_GammaRDBM_MSCW_SumE.Add(Hist_GammaRDBM_MSCW)

def Theta2HistScale(Hist,scale,scale_err):

    for b in range(1,Hist.GetNbinsX()+1):
        old_content = Hist.GetBinContent(b)
        old_error = Hist.GetBinError(b)
        new_content = old_content*scale
        new_error = 0
        if old_content>0 and scale>0:
            new_error = new_content*(old_error/old_content)
        Hist.SetBinContent(b,new_content)
        Hist.SetBinError(b,new_error)

def StackTheta2Histograms():

    Hist_Data_Theta2_SumE.Add(Hist_Data_Theta2)
    Hist_TrueBkgd_Theta2_SumE.Add(Hist_TrueBkgd_Theta2)
    Hist_Dark_Theta2_SumE.Add(Hist_Dark_Theta2)
    Hist_Bkgd_Theta2_SumE.Add(Hist_Bkgd_Theta2)
    Hist_Bkgd_Theta2_Raw_SumE.Add(Hist_Bkgd_Theta2_Raw)

def NormalizeTheta2Histograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1

    HistName = "Hist_Data_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Theta2.Reset()
    Hist_Data_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_TrueBkgd_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_TrueBkgd_Theta2.Reset()
    Hist_TrueBkgd_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Dark_Theta2.Reset()
    Hist_Dark_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_Theta2.Reset()
    Hist_Bkgd_Theta2.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_SelectFoV_Theta2_Raw_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_Theta2_Raw.Reset()
    Hist_Bkgd_Theta2_Raw.Add(InputFile.Get(HistName))

    bkg_total, bkg_err = IntegralAndError(Hist_Bkgd_MSCW,bin_lower,bin_upper)
    old_integral = Hist_Bkgd_Theta2.Integral()
    old_integral_raw = Hist_Bkgd_Theta2_Raw.Integral()
    bkgd_scale = 0
    bkgd_scale_err = 0
    bkgd_scale_raw = 0
    bkgd_scale_raw_err = 0
    if not bkg_total==0 and not old_integral==0:
        bkgd_scale = bkg_total/old_integral
        bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        bkgd_scale_raw = bkg_total/old_integral_raw
        bkgd_scale_raw_err = bkgd_scale_raw*(bkg_err/bkg_total)
    else:
        bkgd_scale = 0
        bkgd_scale_err = 0
        bkgd_scale_raw = 0
        bkgd_scale_raw_err = 0
    if not bkg_total==0:
        Theta2HistScale(Hist_Bkgd_Theta2,bkgd_scale,bkgd_scale_err)
        Theta2HistScale(Hist_Bkgd_Theta2_Raw,bkgd_scale_raw,bkgd_scale_raw_err)
    else:
        Hist_Bkgd_Theta2.Scale(0)
        Hist_Bkgd_Theta2_Raw.Scale(0)

    bkg_total, bkg_err = IntegralAndError(Hist_Dark_MSCW,bin_lower,bin_upper)
    old_integral = Hist_Dark_Theta2.Integral()
    bkgd_scale = 0
    bkgd_scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        bkgd_scale = bkg_total/old_integral
        bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
    else:
        bkgd_scale = 0
        bkgd_scale_err = 0
    if not bkg_total==0:
        Theta2HistScale(Hist_Dark_Theta2,bkgd_scale,bkgd_scale_err)
    else:
        Hist_Dark_Theta2.Scale(0)

def GetShowerHistogramsFromFile(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower_x = Hist2D_Data.GetXaxis().FindBin(MSCL_lower_cut)
    bin_upper_x = Hist2D_Data.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper_y = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1

    HistName = "Hist_Data_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_Data.Reset()
    Hist2D_Data.Add(InputFile.Get(HistName))
    HistName = "Hist_Dark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_Dark.Reset()
    Hist2D_Dark.Add(InputFile.Get(HistName))
    HistName = "Hist_TrueBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_TrueBkgd.Reset()
    Hist2D_TrueBkgd.Add(InputFile.Get(HistName))
    HistName = "Hist_GammaDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_GammaDark.Reset()
    Hist2D_GammaDark.Add(InputFile.Get(HistName))
    HistName = "Hist_GammaRDBM_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_GammaRDBM.Reset()
    Hist2D_GammaRDBM.Add(InputFile.Get(HistName))
    HistName = "Hist_Bkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist2D_Bkgd.Reset()
    Hist2D_Bkgd.Add(InputFile.Get(HistName))
    for bx in range(1,Hist2D_Bkgd.GetNbinsX()+1):
        for by in range(1,Hist2D_Bkgd.GetNbinsY()+1):
            data_err = pow(Hist2D_Bkgd.GetBinContent(bx,by),0.5)
            Hist2D_Bkgd.SetBinError(bx,by,data_err)

    Hist_Data_MSCL.Reset()
    Hist_Data_MSCL.Add(Hist2D_Data.ProjectionX("Hist1D_Data_MSCL",bin_lower_y,bin_upper_y))
    Hist_Data_MSCW.Reset()
    Hist_Data_MSCW.Add(Hist2D_Data.ProjectionY("Hist1D_Data_MSCW",bin_lower_x,bin_upper_x))
    Hist_Dark_MSCL.Reset()
    Hist_Dark_MSCL.Add(Hist2D_Dark.ProjectionX("Hist1D_Dark_MSCL",bin_lower_y,bin_upper_y))
    Hist_Dark_MSCW.Reset()
    Hist_Dark_MSCW.Add(Hist2D_Dark.ProjectionY("Hist1D_Dark_MSCW",bin_lower_x,bin_upper_x))
    Hist_TrueBkgd_MSCL.Reset()
    Hist_TrueBkgd_MSCL.Add(Hist2D_TrueBkgd.ProjectionX("Hist1D_TrueBkgd_MSCL",bin_lower_y,bin_upper_y))
    Hist_TrueBkgd_MSCW.Reset()
    Hist_TrueBkgd_MSCW.Add(Hist2D_TrueBkgd.ProjectionY("Hist1D_TrueBkgd_MSCW",bin_lower_x,bin_upper_x))
    Hist_GammaDark_MSCL.Reset()
    Hist_GammaDark_MSCL.Add(Hist2D_GammaDark.ProjectionX("Hist1D_GammaDark_MSCL",bin_lower_y,bin_upper_y))
    Hist_GammaDark_MSCW.Reset()
    Hist_GammaDark_MSCW.Add(Hist2D_GammaDark.ProjectionY("Hist1D_GammaDark_MSCW",bin_lower_x,bin_upper_x))
    Hist_GammaRDBM_MSCL.Reset()
    Hist_GammaRDBM_MSCL.Add(Hist2D_GammaRDBM.ProjectionX("Hist1D_GammaRDBM_MSCL",bin_lower_y,bin_upper_y))
    Hist_GammaRDBM_MSCW.Reset()
    Hist_GammaRDBM_MSCW.Add(Hist2D_GammaRDBM.ProjectionY("Hist1D_GammaRDBM_MSCW",bin_lower_x,bin_upper_x))
    Hist_Bkgd_MSCL.Reset()
    Hist_Bkgd_MSCL.Add(Hist2D_Bkgd.ProjectionX("Hist1D_Bkgd_MSCL",bin_lower_y,bin_upper_y))
    Hist_Bkgd_MSCW.Reset()
    Hist_Bkgd_MSCW.Add(Hist2D_Bkgd.ProjectionY("Hist1D_Bkgd_MSCW",bin_lower_x,bin_upper_x))

    InputFile.Close()

def PlotsStackedHistograms(tag):

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_MSCW_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_MSCW_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Bkgd_MSCW_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_GammaRDBM_MSCW_SumE]
    legends += ['predict. #gamma']
    colors += [3]
    plotname = 'Stack_MSCW_MDM_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_Theta2_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_Theta2_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Bkgd_Theta2_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    plotname = 'Stack_theta2_MDM_%s'%(tag)
    title = '#theta^{2}'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,0.,20.,-1)

def PlotsPerEnergyBin():

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_MSCW]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_MSCW]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Bkgd_MSCW]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_GammaRDBM_MSCW]
    legends += ['predict. #gamma']
    colors += [3]
    plotname = 'MSCW_MDM'
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)

def SystAsFunctionOfEnergy():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    PercentCrab = "_Crab0"
    ONOFF = "OFF"
    folder_tag = FileTag[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[source] += [FilePath]

    Hist_MDM_S2B_Energy = []
    Hist_Dark_S2B_Energy = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        exposure_hours, NSB_avg, exposure_hours_dark, NSB_avg_dark = GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        Hist_MDM_S2B_Energy += [ROOT.TH1D("Hist_MDM_S2B_Energy_%s"%(source_name),"",len(energy_list)-1,array('d',energy_list))]
        Hist_Dark_S2B_Energy += [ROOT.TH1D("Hist_Dark_S2B_Energy_%s"%(source_name),"",len(energy_list)-1,array('d',energy_list))]
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            ResetStackedShowerHistograms()
            for path in range(1,len(FilePath_Folder0[source])):
                if not os.path.isfile(FilePath_Folder0[source][path]):continue
                GetShowerHistogramsFromFile(FilePath_Folder0[source][path])
                StackShowerHistograms()
                NormalizeTheta2Histograms(FilePath_Folder0[source][path])
                StackTheta2Histograms()
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            Hist_MDM_S2B_Energy[len(Hist_MDM_S2B_Energy)-1].SetBinContent(e+1,s2b)
            Hist_MDM_S2B_Energy[len(Hist_MDM_S2B_Energy)-1].SetBinError(e+1,s2b_err)
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            Hist_Dark_S2B_Energy[len(Hist_Dark_S2B_Energy)-1].SetBinContent(e+1,s2b)
            Hist_Dark_S2B_Energy[len(Hist_Dark_S2B_Energy)-1].SetBinError(e+1,s2b_err)
    MakeComparisonPlotSigDist(Hist_MDM_S2B_Energy,legend_S2B,color_S2B,'E [GeV]','Mismodeling_MDM_systematics_%s'%(folder_tag),0.6,1.4,True,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_Energy,legend_S2B,color_S2B,'E [GeV]','Mismodeling_Dark_systematics_%s'%(folder_tag),0.6,1.4,True,False)

def SystDarkVsMDM():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    PercentCrab = "_Crab0"
    ONOFF = "OFF"
    folder_tag = FileTag[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[source] += [FilePath]

    Hist_S2B = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        exposure_hours, NSB_avg, exposure_hours_dark, NSB_avg_dark = GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        Hist_S2B += [ROOT.TH1D("Hist_S2B_%s"%(source_name),"",2,0,2)]
        Hist_S2B[len(Hist_S2B)-1].GetXaxis().SetBinLabel(1,'Initial')
        Hist_S2B[len(Hist_S2B)-1].GetXaxis().SetBinLabel(2,'MDM')
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        ResetStackedShowerHistograms()
        for path in range(1,len(FilePath_Folder0[source])):
            if not os.path.isfile(FilePath_Folder0[source][path]):continue
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                GetShowerHistogramsFromFile(FilePath_Folder0[source][path])
                StackShowerHistograms()
                NormalizeTheta2Histograms(FilePath_Folder0[source][path])
                StackTheta2Histograms()
            #PlotsStackedHistograms('%s_%s_%s'%(source_name,path,folder_tag))
        s2b = 0.
        s2b_err = 0.
        if not "_Crab0" in PercentCrab:
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
        Hist_S2B[len(Hist_S2B)-1].SetBinContent(2,s2b)
        Hist_S2B[len(Hist_S2B)-1].SetBinError(2,s2b_err)
        s2b = 0.
        s2b_err = 0.
        if not "_Crab0" in PercentCrab:
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        Hist_S2B[len(Hist_S2B)-1].SetBinContent(1,s2b)
        Hist_S2B[len(Hist_S2B)-1].SetBinError(1,s2b_err)
    MakeComparisonPlotTwoColumn(Hist_S2B,legend_S2B,color_S2B,'methods','Mismodeling_systematics_%s'%(folder_tag),0.6,1.4,False,False)

#SystAsFunctionOfEnergy()
SystDarkVsMDM()

