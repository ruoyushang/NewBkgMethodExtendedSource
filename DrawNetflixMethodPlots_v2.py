
import os
import sys,ROOT
import array
import math
from array import *
from ROOT import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from scipy import special

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")


#ONOFF = "ON"
ONOFF = "OFF"

PercentCrab = "_Crab0"
#PercentCrab = "_Crab10"
#PercentCrab = "_Crab20"
#PercentCrab = "_Crab50"
#PercentCrab = "_Crab100"
#PercentCrab = "_Crab200"

FileFolder = []
FileTag = []
FileLabel = []

Syst_MDM = 0.
Syst_Dark = 0.

#FileFolder += ['output_test']
#FileTag += ['test']
#FileLabel += ['test']

#FileFolder += ['output_unblind_4x4_nominal_tight']
#FileTag += ['tight']
#FileLabel += ['tight']
#FileFolder += ['output_unblind_4x4_nominal_medium']
#FileTag += ['medium']
#FileLabel += ['medium']
FileFolder += ['output_unblind_4x4_nominal_loose']
FileTag += ['loose']
FileLabel += ['loose']

#FileFolder += ['output_unblind_4x4_dNSBm3_loose']
#FileTag += ['dNSBm3']
#FileLabel += ['#Delta NSB = -3']
#FileFolder += ['output_unblind_4x4_nominal_loose']
#FileTag += ['dNSBp0']
#FileLabel += ['#Delta NSB = 0']
#FileFolder += ['output_unblind_4x4_dNSBp3_loose']
#FileTag += ['dNSBp3']
#FileLabel += ['#Delta NSB = 3']

#FileFolder += ['output_unblind_4x4_nominal_loose']
#FileTag += ['dElev0']
#FileLabel += ['#Delta elev. = 0']
#FileFolder += ['output_unblind_4x4_dElevm10_loose']
#FileTag += ['dElevm20']
#FileLabel += ['#Delta elev. = -10']
#FileFolder += ['output_unblind_4x4_dElevm20_loose']
#FileTag += ['dElevm20']
#FileLabel += ['#Delta elev. = -20']

elev_range = []
#elev_range += [[25,35]]
#elev_range += [[35,45]]
elev_range += [[45,55]]
elev_range += [[55,65]]
elev_range += [[65,75]]
elev_range += [[75,85]]

source_list = []
sky_coord = []
for_syst = []
source_list += ['Proton_NSB200']
sky_coord += ['10 07 04 +16 04 55']
for_syst += [False]
source_list += ['Proton_NSB750']
sky_coord += ['10 07 04 +16 04 55']
for_syst += [False]
source_list += ['Crab']
sky_coord += ['05 34 31.97 +22 00 52.1']
for_syst += [True]
source_list += ['Mrk421']
sky_coord += ['11 04 19 +38 11 41']
for_syst += [True]
source_list += ['H1426']
sky_coord += ['14 28 32.609 +42 40 21.05']
for_syst += [True]
source_list += ['PKS1424']
sky_coord += ['14 27 00 +23 47 00']
for_syst += [True]
source_list += ['3C264']
sky_coord += ['11 45 5.009 +19 36 22.74']
for_syst += [True]
source_list += ['OJ287V6']
sky_coord += ['08 54 49.1 +20 05 58.89']
for_syst += [True]
source_list += ['1ES0229']
sky_coord += ['02 32 53.2 +20 16 21']
for_syst += [True]
source_list += ['S3_1227_V6']
sky_coord += ['12 30 14.1 +25 18 07']
for_syst += [True]
source_list += ['MS1221V6']
sky_coord += ['12 24 24.2 +24 36 24']
for_syst += [True]
source_list += ['PKS1441V6']
sky_coord += ['14 43 56.9 +25 01 44']
for_syst += [True]
source_list += ['RBS0413V6']
sky_coord += ['03 19 47 +18 45 42']
for_syst += [True]
source_list += ['PG1553V6']
sky_coord += ['15 55 44.7 +11 11 41']
for_syst += [True]
source_list += ['Segue1V6']
sky_coord += ['10 07 04 +16 04 55']
for_syst += [True]
source_list += ['1ES1011V6']
sky_coord += ['10 15 4.139 +49 26 0.71']
for_syst += [True]
source_list += ['NGC1275V6']
sky_coord += ['03 19 48.1 +41 30 42']
for_syst += [True]
source_list += ['1ES0647V6']
sky_coord += ['06 50 46.490 +25 02 59.62']
for_syst += [True]
source_list += ['1ES1440V6']
sky_coord += ['14 42 48.277 +12 00 40.37']
for_syst += [True]
source_list += ['1ES1741V6']
sky_coord += ['17 44 01.2 +19 32 47']
for_syst += [True]
source_list += ['IC443HotSpot']
sky_coord += ['06 18 2.700 +22 39 36.00']
for_syst += [True]
source_list += ['RGBJ0710']
sky_coord += ['07 10 26.4 +59 09 00']
for_syst += [True]
source_list += ['CasA']
sky_coord += ['23 23 13.8 +58 48 26']
for_syst += [True]
source_list += ['WComaeV6']
sky_coord += ['12 21 31.7 +28 13 59']
for_syst += [True]
source_list += ['M82']
sky_coord += ['09 55 52.7 +69 40 46']
for_syst += [True]
source_list += ['G079']
sky_coord += ['20 32 28.56 +40 19 41.52']
for_syst += [False]
source_list += ['1ES1218V6']
sky_coord += ['12 21 26.3 +30 11 29']
for_syst += [False]
source_list += ['MGRO_J1908_V6']
sky_coord += ['19 07 54 +06 16 07']
for_syst += [False]
source_list += ['MGRO_J1908_V5']
sky_coord += ['19 07 54 +06 16 07']
for_syst += [False]
source_list += ['Segue1V5']
sky_coord += ['10 07 04 +16 04 55']
for_syst += [False]
source_list += ['IC443HotSpotV5']
sky_coord += ['06 18 2.700 +22 39 36.00']
for_syst += [False]
source_list += ['GemingaV6']
sky_coord += ['06 32 28 +17 22 00']
for_syst += [False]
source_list += ['GemingaV5']
sky_coord += ['06 32 28 +17 22 00']
for_syst += [False]

gal_coord = SkyCoord(sky_coord[:], unit=(u.hourangle, u.deg))

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
Hist_Data_Eigenvalues_real = ROOT.TH1D("Hist_Data_Eigenvalues_real","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist_Data_Eigenvector_0_real = ROOT.TH1D("Hist_Data_Eigenvector_0_real","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_Eigenvector_1_real = ROOT.TH1D("Hist_Data_Eigenvector_1_real","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_Eigenvector_2_real = ROOT.TH1D("Hist_Data_Eigenvector_2_real","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

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
Hist_Data_Eigenvalues_real_SumE = ROOT.TH1D("Hist_Data_Eigenvalues_real_SumE","",N_bins_for_deconv,0,N_bins_for_deconv)
Hist_Data_Eigenvector_0_real_SumE = ROOT.TH1D("Hist_Data_Eigenvector_0_real_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_Eigenvector_1_real_SumE = ROOT.TH1D("Hist_Data_Eigenvector_1_real_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
Hist_Data_Eigenvector_2_real_SumE = ROOT.TH1D("Hist_Data_Eigenvector_2_real_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)

NRGBs = 5
#NCont = 512
NCont = 104
stops = [0.00,0.45,0.50,0.55,1.00]
red =   [0.00,1.00,1.00,1.00,1.00]
green = [0.00,1.00,1.00,1.00,0.00]
blue =  [1.00,1.00,1.00,1.00,0.00]
ROOT.TColor.CreateGradientColorTable(NRGBs,array('d',stops),array('d',red),array('d',green),array('d',blue),NCont)
ROOT.gStyle.SetNumberContours(NCont)


def FindSourceIndex(source_name):
    for source in range(0,len(source_list)):
        if source_name==source_list[source]:
            return source
    return 0

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

def Signal_Significance(Hist_SR, Hist_Bkg,range_lower,range_upper,syst):

    #range_lower = 0.4

    norm_bin_low_target = Hist_SR.FindBin(range_lower)
    norm_bin_up_target = Hist_Bkg.FindBin(range_upper)-1
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hist_SR,norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target)
    err_bkg = pow(err_bkg,2)+pow(predict_bkg*syst,2)
    err_bkg = pow(err_bkg,0.5)
    Sig = 1.*CalculateSignificance(data_SR-predict_bkg,predict_bkg,err_bkg)
    return Sig

def MakeGaussComparisonPlot(Hists,legends,colors,title,name):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.15)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

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
            #Hists[h].GetXaxis().SetRangeUser(0,8)
            Hists[h].GetXaxis().SetRangeUser(-5,8)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #Hists[max_hist].SetMinimum(0)
    #Hists[max_hist].SetMaximum(1)
    Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            if h==0:
                Hists[h].Draw("hist same")
            else:
                Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.1,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            #legend.AddEntry(Hists[h],'%s, mean = %.5f, RMS = %.5f'%(legends[h],Hists[h].GetMean(),Hists[h].GetRMS()),"pl")
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")

    pad1.SetLogy()

    c_both.SaveAs('output_plots/%s.png'%(name))

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
        if 'MSCW' in name: Hists[h].GetXaxis().SetRangeUser(-1,MSCW_blind_cut+1)
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
    #err_bkg = pow(err_bkg*err_bkg+err_sys*err_sys,0.5)
    err_bkg = pow(err_bkg*err_bkg+(Syst_MDM*predict_bkg)*(Syst_MDM*predict_bkg),0.5)
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
    if 'MSCW' in name: line2 = ROOT.TLine(-1,1,MSCW_blind_cut+1,1)
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
        Hist_Ratio += [Hists[h].Clone()]
        Hist_Ratio[h].Divide(Hist_Bkg)
        Hist_Ratio[h].SetLineWidth(2)
    for h in range(0,len(Hists)):
        if colors[h]==2:
            Hist_Ratio[h].Draw("B same")
    for h in range(0,len(Hists)):
        if colors[h]==1:
            Hist_Ratio[h].Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s.png'%(name))

def Smooth2DMap(Hist_Old,smooth_size,addLinearly):

    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(4*smooth_size/bin_size) + 1
    for bx1 in range(1,Hist_Old.GetNbinsX()+1):
        for by1 in range(1,Hist_Old.GetNbinsY()+1):
            bin_content = 0
            bin_error = 0
            locationx1 = Hist_Old.GetXaxis().GetBinCenter(bx1)
            locationy1 = Hist_Old.GetYaxis().GetBinCenter(by1)
            for bx2 in range(bx1-nbin_smooth,bx1+nbin_smooth):
                for by2 in range(by1-nbin_smooth,by1+nbin_smooth):
                    if bx2>=1 and bx2<=Hist_Old.GetNbinsX():
                        if by2>=1 and by2<=Hist_Old.GetNbinsY():
                            locationx2 = Hist_Old.GetXaxis().GetBinCenter(bx2)
                            locationy2 = Hist_Old.GetYaxis().GetBinCenter(by2)
                            distance = pow(pow(locationx1-locationx2,2)+pow(locationy1-locationy2,2),0.5)
                            bin_content += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinContent(bx2,by2)
                            if not addLinearly:
                                bin_error += pow(ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2),2)
                            else:
                                bin_error += ROOT.TMath.Gaus(distance,0,smooth_size)*Hist_Old.GetBinError(bx2,by2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            if not addLinearly:
                Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
            else:
                Hist_Smooth.SetBinError(bx1,by1,bin_error)
    return Hist_Smooth

def Make2DSignificanceCameraFoVPlot(exposure,syst_method,Hist_SR,Hist_Bkg,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            theta2 = pow(Hist_SR.GetXaxis().GetBinCenter(bx+1),2)+pow(Hist_SR.GetYaxis().GetBinCenter(by+1),2)

    Hist_Skymap = Hist_SR.Clone()
    Hist_Excess = Hist_SR.Clone()
    syst_fov = 0.
    n_bins = 0
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            n_bins += 1
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+pow(syst_method*Hist_Bkg.GetBinContent(bx+1,by+1),2),0.5)
            syst_fov += pow(NSR/NBkg-1.,2)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
            Hist_Excess.SetBinContent(bx+1,by+1,NSR/NBkg-1.)
    syst_fov = syst_fov/n_bins
    syst_fov = pow(syst_fov,0.5)

    total_bins = 0.
    sum_bins = 0.
    sum_data = 0.
    sum_bkg = 0.
    sum_bkg_err = 0.
    bias_sum_bins = 0.
    bias_sum_data = 0.
    bias_sum_bkg = 0.
    bias_sum_bkg_err = 0.
    global_sum_data = 0.
    global_sum_bkg = 0.
    global_sum_bkg_err = 0.
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            global_sum_data += Hist_SR.GetBinContent(bx+1,by+1)
            global_sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
            global_sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
            bin_ra = Hist_SR.GetXaxis().GetBinCenter(bx+1)
            bin_dec = Hist_SR.GetYaxis().GetBinCenter(by+1)
            bin_theta2 = pow(pow(0.-bin_ra,2)+pow(0.-bin_dec,2),0.5)
            if bin_theta2<theta2_range: total_bins += 1.
            if bin_theta2<theta2_range and Hist_Skymap.GetBinContent(bx+1,by+1)>=highlight_threshold: 
                bias_sum_bins += 1.
                bias_sum_data += Hist_SR.GetBinContent(bx+1,by+1)
                bias_sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
                bias_sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
            if bin_theta2<theta2_range and abs(Hist_Skymap.GetBinContent(bx+1,by+1))>=highlight_threshold: 
                sum_bins += 1.
                sum_data += Hist_SR.GetBinContent(bx+1,by+1)
                sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
                sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
    sum_bkg_err = pow(sum_bkg_err+pow(syst_method*sum_bkg,2),0.5)
    sum_sig = 1.*CalculateSignificance(sum_data-sum_bkg,sum_bkg,sum_bkg_err)
    bias_sum_bkg_err = pow(bias_sum_bkg_err+pow(syst_method*bias_sum_bkg,2),0.5)
    bias_sum_sig = 1.*CalculateSignificance(bias_sum_data-bias_sum_bkg,bias_sum_bkg,bias_sum_bkg_err)
    global_sum_bkg_err = pow(global_sum_bkg_err+pow(syst_method*global_sum_bkg,2),0.5)
    global_sum_sig = 1.*CalculateSignificance(global_sum_data-global_sum_bkg,global_sum_bkg,global_sum_bkg_err)

    binx_size = Hist_SR.GetXaxis().GetBinCenter(2)-Hist_SR.GetXaxis().GetBinCenter(1)
    biny_size = Hist_SR.GetYaxis().GetBinCenter(2)-Hist_SR.GetYaxis().GetBinCenter(1)
    sum_bins = max(1.,sum_bins)
    source_size = sum_bins*binx_size*biny_size
    print 'sum_sig = %s'%(sum_sig)
    pvalue = (total_bins/sum_bins)*special.erfc(sum_sig/pow(2,0.5))
    print 'p-value = %s'%(pvalue)
    sum_sig_lee = max(0,special.erfinv(1.-pvalue)*pow(2,0.5))
    print 'sum_sig_lee = %s'%(sum_sig_lee)


    Hist_Skymap.GetYaxis().SetTitle(ytitle)
    Hist_Skymap.GetXaxis().SetTitle(xtitle)
    Hist_Skymap.SetMaximum(5)
    Hist_Skymap.SetMinimum(-5)
    Hist_Skymap.Draw("COL4Z")
    lumilab4 = ROOT.TLatex(0.2,0.85,'exposure = %0.1f hours'%(exposure) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.04)
    lumilab4.Draw()
    lumilab1 = ROOT.TLatex(0.2,0.8,'global %0.1f#sigma (syst = %0.1f%%)'%(global_sum_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.04)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.2,0.75,'local %0.1f#sigma'%(bias_sum_sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.04)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.2,0.70,'local %0.1f#sigma (unbiased)'%(sum_sig) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.04)
    lumilab3.Draw()
    if not math.isinf(sum_sig_lee):
        lumilab5 = ROOT.TLatex(0.2,0.65,'local %0.1f#sigma (LEE)'%(sum_sig_lee) )
        lumilab5.SetNDC()
        lumilab5.SetTextSize(0.04)
        lumilab5.Draw()
    canvas.SaveAs('output_plots/SkymapSig_%s.png'%(name))

    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.SetMaximum(0.02)
    Hist_Excess.SetMinimum(-0.02)
    Hist_Excess.Draw("COL4Z")
    lumilab4 = ROOT.TLatex(0.2,0.85,'exposure = %0.1f hours'%(exposure) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.04)
    lumilab4.Draw()
    lumilab3 = ROOT.TLatex(0.2,0.8,'rad. accept. syst. = %0.2f%%'%(syst_fov*100.) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.04)
    lumilab3.Draw()
    canvas.SaveAs('output_plots/SkymapExcess_%s.png'%(name))

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    Hist_Sig_cut = ROOT.TH1D("Hist_Sig_cut","",65,-5,8)
    for bx in range(0,Hist_Skymap.GetNbinsX()):
        for by in range(0,Hist_Skymap.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0 and not Hist_Bkg.GetBinContent(bx+1,by+1)==0:
                content = Hist_Skymap.GetBinContent(bx+1,by+1)
                Hist_Sig.Fill(content)
                locationx1 = Hist_Skymap.GetXaxis().GetBinCenter(bx)
                locationy1 = Hist_Skymap.GetYaxis().GetBinCenter(by)
                if (pow(locationx1*locationx1+locationy1*locationy1,0.5)>0.5):
                    Hist_Sig_cut.Fill(content)
    Hist_Model = ROOT.TH1D("Hist_Model","",65,-5,8)
    Hist_Model.FillRandom("func",10000*int(Hist_Sig.GetEntries()))
    Hist_Model.Scale(1./10000.)
    Hist_Model.SetMinimum(0.5)
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Perfect Gaussian']
    color_list += [2]
    Hist_list += [Hist_Sig]
    legend_list += ['Data']
    color_list += [4]
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))

def Make2DSignificancePlot(syst_method,Hist_SR,Hist_Bkg,Hist_SR_smooth,Hist_Bkg_smooth,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            theta2 = pow(Hist_SR.GetXaxis().GetBinCenter(bx+1),2)+pow(Hist_SR.GetYaxis().GetBinCenter(by+1),2)

    Hist_Skymap = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+pow(syst_method*Hist_Bkg.GetBinContent(bx+1,by+1),2),0.5)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap.SetBinContent(bx+1,by+1,Sig)
    Hist_Skymap_smooth = Hist_SR_smooth.Clone()
    for bx in range(0,Hist_SR_smooth.GetNbinsX()):
        for by in range(0,Hist_SR_smooth.GetNbinsY()):
            if Hist_Bkg_smooth.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR_smooth.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR_smooth.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg_smooth.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg_smooth.GetBinError(bx+1,by+1),2)+pow(syst_method*Hist_Bkg_smooth.GetBinContent(bx+1,by+1),2),0.5)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Skymap_smooth.SetBinContent(bx+1,by+1,Sig)

    total_bins = 0.
    sum_bins = 0.
    sum_data = 0.
    sum_bkg = 0.
    sum_bkg_err = 0.
    bias_sum_bins = 0.
    bias_sum_data = 0.
    bias_sum_bkg = 0.
    bias_sum_bkg_err = 0.
    global_sum_data = 0.
    global_sum_bkg = 0.
    global_sum_bkg_err = 0.
    Hist_Ratio = Hist_SR_smooth.Clone()
    Hist_Ratio.Add(Hist_Bkg_smooth,-1.)
    Hist_Ratio.Divide(Hist_Bkg_smooth)
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            Hist_Highlight_Skymap.SetBinContent(bx+1,by+1,0.)
            Hist_Highlight_Bias_Skymap.SetBinContent(bx+1,by+1,0.)
            global_sum_data += Hist_SR.GetBinContent(bx+1,by+1)
            global_sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
            global_sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
            bin_ra = Hist_SR.GetXaxis().GetBinCenter(bx+1)
            bin_dec = Hist_SR.GetYaxis().GetBinCenter(by+1)
            bin_theta2 = pow(pow(source_ra-bin_ra,2)+pow(source_dec-bin_dec,2),0.5)
            if bin_theta2<theta2_range: total_bins += 1.
            if bin_theta2<theta2_range and Hist_Skymap_smooth.GetBinContent(bx+1,by+1)>=highlight_threshold: 
                local_sig = Hist_Skymap_smooth.GetBinContent(bx+1,by+1);
                #Hist_Highlight_Bias_Skymap.SetBinContent(bx+1,by+1,int(local_sig))
                Hist_Highlight_Bias_Skymap.SetBinContent(bx+1,by+1,1.)
                bias_sum_bins += 1.
                bias_sum_data += Hist_SR.GetBinContent(bx+1,by+1)
                bias_sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
                bias_sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
            if bin_theta2<theta2_range and abs(Hist_Skymap_smooth.GetBinContent(bx+1,by+1))>=highlight_threshold: 
                local_sig = Hist_Skymap_smooth.GetBinContent(bx+1,by+1);
                #Hist_Highlight_Skymap.SetBinContent(bx+1,by+1,int(abs(local_sig)))
                Hist_Highlight_Skymap.SetBinContent(bx+1,by+1,1.)
                sum_bins += 1.
                sum_data += Hist_SR.GetBinContent(bx+1,by+1)
                sum_bkg += Hist_Bkg.GetBinContent(bx+1,by+1)
                sum_bkg_err += pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+Hist_Bkg.GetBinContent(bx+1,by+1)
            if Hist_Bkg_smooth.GetBinContent(bx+1,by+1)<100: 
                Hist_Ratio.SetBinContent(bx+1,by+1,0.)
    sum_bkg_err = pow(sum_bkg_err+pow(syst_method*sum_bkg,2),0.5)
    sum_sig = 1.*CalculateSignificance(sum_data-sum_bkg,sum_bkg,sum_bkg_err)
    bias_sum_bkg_err = pow(bias_sum_bkg_err+pow(syst_method*bias_sum_bkg,2),0.5)
    bias_sum_sig = 1.*CalculateSignificance(bias_sum_data-bias_sum_bkg,bias_sum_bkg,bias_sum_bkg_err)
    global_sum_bkg_err = pow(global_sum_bkg_err+pow(syst_method*global_sum_bkg,2),0.5)
    global_sum_sig = 1.*CalculateSignificance(global_sum_data-global_sum_bkg,global_sum_bkg,global_sum_bkg_err)

    binx_size = Hist_SR.GetXaxis().GetBinCenter(2)-Hist_SR.GetXaxis().GetBinCenter(1)
    biny_size = Hist_SR.GetYaxis().GetBinCenter(2)-Hist_SR.GetYaxis().GetBinCenter(1)
    sum_bins = max(1.,sum_bins)
    source_size = sum_bins*binx_size*biny_size
    print 'sum_sig = %s'%(sum_sig)
    pvalue = (total_bins/sum_bins)*special.erfc(sum_sig/pow(2,0.5))
    print 'p-value = %s'%(pvalue)
    sum_sig_lee = max(0,special.erfinv(1.-pvalue)*pow(2,0.5))
    print 'sum_sig_lee = %s'%(sum_sig_lee)

    max_sig = 0.
    max_x = 0
    max_y = 0
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            sig = Hist_Skymap_smooth.GetBinContent(bx+1,by+1)
            if max_sig < sig:
                max_sig = sig
                max_x = Hist_Skymap_smooth.GetXaxis().GetBinCenter(bx+1)
                max_y = Hist_Skymap_smooth.GetYaxis().GetBinCenter(by+2)

    Hist_Skymap.GetYaxis().SetTitle(ytitle)
    Hist_Skymap.GetXaxis().SetTitle(xtitle)
    Hist_Skymap.SetMaximum(5)
    Hist_Skymap.SetMinimum(-5)
    Hist_Skymap.Draw("COL4Z")
    Hist_Highlight_Bias_Skymap.Draw("CONT3 same")
    lumilab1 = ROOT.TLatex(0.2,0.85,'global %0.1f#sigma (syst = %0.1f%%)'%(global_sum_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.04)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.2,0.80,'local %0.1f#sigma'%(bias_sum_sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.04)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.2,0.75,'local %0.1f#sigma (unbiased)'%(sum_sig) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.04)
    lumilab3.Draw()
    if not math.isinf(sum_sig_lee):
        lumilab5 = ROOT.TLatex(0.2,0.70,'local %0.1f#sigma (LEE)'%(sum_sig_lee) )
        lumilab5.SetNDC()
        lumilab5.SetTextSize(0.04)
        lumilab5.Draw()
    canvas.SaveAs('output_plots/SkymapSig_%s.png'%(name))

    Hist_Skymap_smooth.GetYaxis().SetTitle(ytitle)
    Hist_Skymap_smooth.GetXaxis().SetTitle(xtitle)
    Hist_Skymap_smooth.SetMaximum(5)
    Hist_Skymap_smooth.SetMinimum(-5)
    Hist_Skymap_smooth.Draw("COL4Z")
    Hist_Highlight_Bias_Skymap.Draw("CONT3 same")
    lumilab1 = ROOT.TLatex(0.2,0.85,'global %0.1f#sigma (syst = %0.1f%%)'%(global_sum_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.04)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.2,0.80,'local %0.1f#sigma'%(bias_sum_sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.04)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.2,0.75,'local %0.1f#sigma (unbiased)'%(sum_sig) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.04)
    lumilab3.Draw()
    if not math.isinf(sum_sig_lee):
        lumilab5 = ROOT.TLatex(0.2,0.70,'local %0.1f#sigma (LEE)'%(sum_sig_lee) )
        lumilab5.SetNDC()
        lumilab5.SetTextSize(0.04)
        lumilab5.Draw()
    canvas.SaveAs('output_plots/SkymapSmoothSig_%s.png'%(name))

    #Hist_Ratio.GetYaxis().SetTitle(ytitle)
    #Hist_Ratio.GetXaxis().SetTitle(xtitle)
    ##Hist_Ratio.SetMaximum(0.1)
    ##Hist_Ratio.SetMinimum(-0.1)
    #Hist_Ratio.Draw("COL4Z")
    #Hist_Highlight_Bias_Skymap.Draw("CONT3 same")
    #canvas.SaveAs('output_plots/SkymapRatio_%s.png'%(name))

    #Hist_Excess = Hist_SR.Clone()
    #Hist_Excess.GetYaxis().SetTitle(ytitle)
    #Hist_Excess.GetXaxis().SetTitle(xtitle)
    #Hist_Excess.Draw("COL4Z")
    #lumilab1 = ROOT.TLatex(0.15,0.90,'Integral = %0.1f'%(Hist_Excess.Integral()) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.05)
    #lumilab1.Draw()
    #canvas.SaveAs('output_plots/SkymapTotal_%s_%s_%s.png'%(name,target,tag))

    #Hist_Excess = Hist_Bkg.Clone()
    #Hist_Excess.GetYaxis().SetTitle(ytitle)
    #Hist_Excess.GetXaxis().SetTitle(xtitle)
    #Hist_Excess.Draw("COL4Z")
    #lumilab2 = ROOT.TLatex(0.15,0.90,'Integral = %0.1f'%(Hist_Excess.Integral()) )
    #lumilab2.SetNDC()
    #lumilab2.SetTextSize(0.05)
    #lumilab2.Draw()
    #canvas.SaveAs('output_plots/SkymapBkg_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR_smooth.Clone()
    Hist_Excess.Add(Hist_Bkg_smooth,-1.)
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    Hist_Highlight_Bias_Skymap.Draw("CONT3 same")
    lumilab1 = ROOT.TLatex(0.2,0.85,'global %0.1f#sigma (syst = %0.1f%%)'%(global_sum_sig,syst_method*100.) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.04)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.2,0.80,'local %0.1f#sigma'%(bias_sum_sig) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.04)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.2,0.75,'local %0.1f#sigma (unbiased)'%(sum_sig) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.04)
    lumilab3.Draw()
    if not math.isinf(sum_sig_lee):
        lumilab5 = ROOT.TLatex(0.2,0.70,'local %0.1f#sigma (LEE)'%(sum_sig_lee) )
        lumilab5.SetNDC()
        lumilab5.SetTextSize(0.04)
        lumilab5.Draw()
    canvas.SaveAs('output_plots/SkymapSmoothExcess_%s.png'%(name))

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    Hist_Sig_cut = ROOT.TH1D("Hist_Sig_cut","",65,-5,8)
    for bx in range(0,Hist_Skymap.GetNbinsX()):
        for by in range(0,Hist_Skymap.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0 and not Hist_Bkg.GetBinContent(bx+1,by+1)==0:
                content = Hist_Skymap.GetBinContent(bx+1,by+1)
                Hist_Sig.Fill(content)
                locationx1 = Hist_Skymap.GetXaxis().GetBinCenter(bx)
                locationy1 = Hist_Skymap.GetYaxis().GetBinCenter(by)
                if (pow(locationx1*locationx1+locationy1*locationy1,0.5)>0.5):
                    Hist_Sig_cut.Fill(content)
    Hist_Model = ROOT.TH1D("Hist_Model","",65,-5,8)
    Hist_Model.FillRandom("func",10000*int(Hist_Sig.GetEntries()))
    Hist_Model.Scale(1./10000.)
    Hist_Model.SetMinimum(0.5)
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Perfect Gaussian']
    color_list += [2]
    Hist_list += [Hist_Sig]
    legend_list += ['Data']
    color_list += [4]
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))

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
            weight = 1./pow(Hists[h].GetBinError(binx+1),2)
            #weight = 1.
            syst_this_bin += weight*pow(1.-Hists[h].GetBinContent(binx+1),2)
            weight_this_bin += weight
            raw_weight_this_bin += 1.
        syst_this_bin = syst_this_bin/weight_this_bin
        syst_this_bin = pow(syst_this_bin,0.5)
        syst += [syst_this_bin]
        syst_err += [syst_this_bin/pow(raw_weight_this_bin,0.5)]

    #lumilab1 = ROOT.TLatex(0.25,0.9,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[0],100.*syst_err[0]) )
    #lumilab1.SetNDC()
    #lumilab1.SetTextSize(0.05)
    #lumilab1.Draw()
    #lumilab2 = ROOT.TLatex(0.60,0.9,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[1],100.*syst_err[1]) )
    #lumilab2.SetNDC()
    #lumilab2.SetTextSize(0.05)
    #lumilab2.Draw()
    lumilab = []
    delta = Hists[0].GetBinCenter(2)-Hists[0].GetBinCenter(1)
    for binx in range(0,Hists[0].GetNbinsX()):
        #lumilab += [ROOT.TLatex(float(binx+0.5*delta)/length,0.9,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[binx],100.*syst_err[binx]))]
        #lumilab[binx].SetNDC()
        lumilab += [ROOT.TLatex(Hists[0].GetBinCenter(binx+1)-0.25*delta,1.3,'Syst = %0.1f#pm%0.1f %%'%(100.*syst[binx],100.*syst_err[binx]))]
        lumilab[binx].SetTextSize(0.04)
        lumilab[binx].Draw()

    c_both.SaveAs('output_plots/%s.png'%(name))

    return syst[0], syst[1]

def MakeComparisonPlotSigDist(Hists,legends,colors,title,name,minheight,maxheight,logx,logy):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0.35,1,0.7)
    pad1.SetBottomMargin(0.0)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35)
    pad2.SetBottomMargin(0.25)
    pad2.SetTopMargin(0.0)
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
            weight = 1./pow(Hists[h].GetBinError(binx+1),2)
            #weight = 1.
            shift_this_bin += weight*(1.-Hists[h].GetBinContent(binx+1))
            syst_this_bin += weight*pow(1.-Hists[h].GetBinContent(binx+1),2)
            weight_this_bin += weight
            raw_weight_this_bin += 1.
        if weight_this_bin==0.:
            shift_this_bin = 0.
            syst_this_bin = 0.
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
    Hist_Syst.GetXaxis().SetLabelSize(0.1)
    Hist_Syst.GetYaxis().SetLabelSize(0.06)
    Hist_Syst.GetXaxis().SetTitleSize(0.1)
    Hist_Syst.GetXaxis().SetTitle(title)
    Hist_Syst.GetYaxis().SetTitleOffset(0.5)
    Hist_Syst.GetYaxis().SetTitleSize(0.08)
    Hist_Syst.GetYaxis().SetTitle('Syst. Unc. %')
    Hist_Syst.GetXaxis().SetTitleOffset(1.0)
    Hist_Syst.SetMarkerStyle(0)
    Hist_Syst.SetFillColor(1)
    Hist_Syst.SetFillStyle(3004)
    Hist_Syst.Draw("hist")

    if logx: 
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s.png'%(name))

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,y_min,y_max,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.7,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.7)
    pad1.SetBottomMargin(0.1)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    if logy: pad1.SetGrid()
    pad1.Draw()
    pad3.Draw()

    pad1.cd()
    if logy: pad1.SetLogy()


    min_heigh = 0
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
            Hists[h].GetXaxis().SetTitle(title_x)
            Hists[h].GetYaxis().SetTitle(title_y)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
            if min_heigh > Hists[h].GetMinimum(): 
                min_heigh = Hists[h].GetMinimum()

    if not y_max==0. and not y_min==0.:
        Hists[0].SetMaximum(y_max)
        Hists[0].SetMinimum(y_min)
        Hists[0].Draw("E")
    else:
        Hists[max_hist].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].SetLineWidth(2)
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
            legend.AddEntry(Hists[h],'%s'%(legends[h]),"pl")
    legend.Draw("SAME")


    if logx: 
        pad1.SetLogx()

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
    for path in range(0,len(file_list)):
        print 'Read file: %s'%(file_list[path])
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
    if not exposure_hours==0:
        NSB_avg = NSB_avg/exposure_hours
        NSB_avg_dark = NSB_avg_dark/exposure_hours_dark
    else:
        NSB_avg = 0
        NSB_avg_dark = 0

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

    Hist_Data_Eigenvalues_real_SumE.Reset()
    Hist_Data_Eigenvector_0_real_SumE.Reset()
    Hist_Data_Eigenvector_1_real_SumE.Reset()
    Hist_Data_Eigenvector_2_real_SumE.Reset()

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

    Hist_Data_Eigenvalues_real_SumE.Add(Hist_Data_Eigenvalues_real)
    Hist_Data_Eigenvector_0_real.Scale(1./Hist_Data_Eigenvector_0_real.Integral())
    Hist_Data_Eigenvector_0_real.Scale(Hist_Data_Eigenvalues_real.GetBinContent(1))
    Hist_Data_Eigenvector_0_real_SumE.Add(Hist_Data_Eigenvector_0_real)
    Hist_Data_Eigenvector_1_real.Scale(1./Hist_Data_Eigenvector_1_real.Integral())
    Hist_Data_Eigenvector_1_real.Scale(Hist_Data_Eigenvalues_real.GetBinContent(2))
    Hist_Data_Eigenvector_1_real_SumE.Add(Hist_Data_Eigenvector_1_real)
    Hist_Data_Eigenvector_2_real.Scale(1./Hist_Data_Eigenvector_2_real.Integral())
    Hist_Data_Eigenvector_2_real.Scale(Hist_Data_Eigenvalues_real.GetBinContent(3))
    Hist_Data_Eigenvector_2_real_SumE.Add(Hist_Data_Eigenvector_2_real)

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

def StackCameraFoVHistograms():

    Hist_Data_CameraFoV_SumE.Add(Hist_Data_CameraFoV)
    Hist_Bkgd_CameraFoV_SumE.Add(Hist_Bkgd_CameraFoV)
    Hist_Bkgd_CameraFoV_Raw_SumE.Add(Hist_Bkgd_CameraFoV_Raw)

def StackSkymapHistograms():

    Hist_Data_Skymap_SumE.Add(Hist_Data_Skymap)
    Hist_Dark_Skymap_SumE.Add(Hist_Dark_Skymap)
    Hist_Bkgd_Skymap_SumE.Add(Hist_Bkgd_Skymap)
    Hist_Bkgd_Skymap_Raw_SumE.Add(Hist_Bkgd_Skymap_Raw)

def StackTheta2Histograms():

    Hist_Data_Theta2_SumE.Add(Hist_Data_Theta2)
    Hist_TrueBkgd_Theta2_SumE.Add(Hist_TrueBkgd_Theta2)
    Hist_Dark_Theta2_SumE.Add(Hist_Dark_Theta2)
    Hist_Bkgd_Theta2_SumE.Add(Hist_Bkgd_Theta2)
    Hist_Bkgd_Theta2_Raw_SumE.Add(Hist_Bkgd_Theta2_Raw)

def RaDecHistScale(Hist,scale,scale_err):

    for bx in range(1,Hist.GetNbinsX()+1):
        for by in range(1,Hist.GetNbinsY()+1):
            old_content = Hist.GetBinContent(bx,by)
            old_error = Hist.GetBinError(bx,by)
            new_content = old_content*scale
            new_error = 0
            if old_content>0 and scale>0:
                new_error = new_content*(old_error/old_content)
            Hist.SetBinContent(bx,by,new_content)
            Hist.SetBinError(bx,by,new_error)
            new_error = 0
            if old_content>0 and scale>0:
                new_error = new_content*(scale_err/scale)

def NormalizeCameraFoVHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1

    HistName = "Hist_Data_SR_CameraFoV_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_CameraFoV.Reset()
    Hist_Data_CameraFoV.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_CameraFoV_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_CameraFoV.Reset()
    Hist_Bkgd_CameraFoV.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_CameraFoV_Raw_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_CameraFoV_Raw.Reset()
    Hist_Bkgd_CameraFoV_Raw.Add(InputFile.Get(HistName))

    #bkg_total, bkg_err = IntegralAndError(Hist_Data_MSCW,bin_lower,bin_upper)
    bkg_total = Hist_Data_CameraFoV.Integral()
    bkg_err = 0

    old_integral = Hist_Bkgd_CameraFoV.Integral()
    scale = 0
    scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        scale = bkg_total/old_integral
        scale_err = scale*(bkg_err/bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not bkg_total==0:
        RaDecHistScale(Hist_Bkgd_CameraFoV,scale,scale_err)
    else:
        Hist_Bkgd_CameraFoV.Scale(0)

    old_integral = Hist_Bkgd_CameraFoV_Raw.Integral()
    scale = 0
    scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        scale = bkg_total/old_integral
        scale_err = scale*(bkg_err/bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not bkg_total==0:
        RaDecHistScale(Hist_Bkgd_CameraFoV_Raw,scale,scale_err)
    else:
        Hist_Bkgd_CameraFoV_Raw.Scale(0)

def NormalizeSkyMapHistograms(FilePath):

    global MSCW_blind_cut
    global MSCL_blind_cut
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_blind_cut = InfoTree.MSCW_cut_blind
    MSCL_blind_cut = InfoTree.MSCL_cut_blind
    bin_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1

    HistName = "Hist_Data_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Skymap.Reset()
    Hist_Data_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Dark_Skymap.Reset()
    Hist_Dark_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_Skymap.Reset()
    Hist_Bkgd_Skymap.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_CR_Skymap_Raw_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Bkgd_Skymap_Raw.Reset()
    Hist_Bkgd_Skymap_Raw.Add(InputFile.Get(HistName))

    bkg_total, bkg_err = IntegralAndError(Hist_Bkgd_MSCW,bin_lower,bin_upper)
    dark_bkg_total, dark_bkg_err = IntegralAndError(Hist_Dark_MSCW,bin_lower,bin_upper)

    old_integral = Hist_Dark_Skymap.Integral()
    scale = 0
    scale_err = 0
    if not dark_bkg_total==0 and not old_integral==0:
        scale = dark_bkg_total/old_integral
        scale_err = scale*(dark_bkg_err/dark_bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not dark_bkg_total==0:
        RaDecHistScale(Hist_Dark_Skymap,scale,scale_err)
    else:
        Hist_Dark_Skymap.Scale(0)

    old_integral = Hist_Bkgd_Skymap.Integral()
    scale = 0
    scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        scale = bkg_total/old_integral
        scale_err = scale*(bkg_err/bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not bkg_total==0:
        RaDecHistScale(Hist_Bkgd_Skymap,scale,scale_err)
    else:
        Hist_Bkgd_Skymap.Scale(0)

    old_integral = Hist_Bkgd_Skymap_Raw.Integral()
    scale = 0
    scale_err = 0
    if not bkg_total==0 and not old_integral==0:
        scale = bkg_total/old_integral
        scale_err = scale*(bkg_err/bkg_total)
    else:
        scale = 0
        scale_err = 0
    if not bkg_total==0:
        RaDecHistScale(Hist_Bkgd_Skymap_Raw,scale,scale_err)
    else:
        Hist_Bkgd_Skymap_Raw.Scale(0)

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

    HistName = "Hist_Data_Eigenvalues_real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Eigenvalues_real.Reset()
    Hist_Data_Eigenvalues_real.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_InvEigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Eigenvector_0_real.Reset()
    Hist_Data_Eigenvector_0_real.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_InvEigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Eigenvector_1_real.Reset()
    Hist_Data_Eigenvector_1_real.Add(InputFile.Get(HistName))
    HistName = "Hist_Data_InvEigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
    Hist_Data_Eigenvector_2_real.Reset()
    Hist_Data_Eigenvector_2_real.Add(InputFile.Get(HistName))

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
    Hists += [Hist_Data_MSCL_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_MSCL_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Bkgd_MSCL_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_GammaRDBM_MSCL_SumE]
    legends += ['predict. #gamma']
    colors += [3]
    plotname = 'Stack_MSCL_MDM_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_MSCW_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_MSCW_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Dark_MSCW_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_GammaDark_MSCW_SumE]
    legends += ['predict. #gamma']
    colors += [3]
    plotname = 'Stack_MSCW_Dark_%s'%(tag)
    title = 'MSCW'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_MSCL_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_MSCL_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Dark_MSCL_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    Hists += [Hist_GammaDark_MSCL_SumE]
    legends += ['predict. #gamma']
    colors += [3]
    plotname = 'Stack_MSCL_Dark_%s'%(tag)
    title = 'MSCL'
    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)

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

    Hists = []
    legends = []
    colors = []
    Hists += [Hist_Data_Theta2_SumE]
    legends += ['obs. data']
    colors += [1]
    Hists += [Hist_TrueBkgd_Theta2_SumE]
    legends += ['true bkg.']
    colors += [2]
    Hists += [Hist_Bkgd_Theta2_Raw_SumE]
    legends += ['predict. bkg.']
    colors += [4]
    plotname = 'Stack_theta2_Raw_MDM_%s'%(tag)
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

def SystAsFunctionOfSignal():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    SignalFlux = []
    SignalFlux_label = []
    SignalFlux += ['_Crab0']
    SignalFlux_label += ['0']
    SignalFlux += ['_Crab20']
    SignalFlux_label += ['0.2']
    SignalFlux += ['_Crab50']
    SignalFlux_label += ['0.5']
    SignalFlux += ['_Crab100']
    SignalFlux_label += ['1.0']
    SignalFlux += ['_Crab200']
    SignalFlux_label += ['2.0']
    SignalFlux += ['_Crab400']
    SignalFlux_label += ['4.0']

    Hist_MDM_S2B_Flux = []
    Hist_Dark_S2B_Flux = []
    Hist_Eigenvalue_2ndRank = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        source_name = source_list[source]
        Hist_MDM_S2B_Flux += [ROOT.TH1D("Hist_MDM_S2B_Flux_%s"%(source_name),"",len(SignalFlux),0,len(SignalFlux))]
        Hist_Dark_S2B_Flux += [ROOT.TH1D("Hist_Dark_S2B_Flux_%s"%(source_name),"",len(SignalFlux),0,len(SignalFlux))]
        Hist_Eigenvalue_2ndRank += [ROOT.TH1D("Hist_Eigenvalue_2ndRank_%s"%(source_name),"",len(SignalFlux),0,len(SignalFlux))]
        legend_S2B += ['%s'%(source_name)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        for flux in range(0,len(SignalFlux)):
            ResetStackedShowerHistograms()
            for elev in range(0,len(elev_range)):
                file_elev_lower = elev_range[elev][0]
                file_elev_upper = elev_range[elev][1]
                FilePath_Folder0 = []
                exposure_hours = 0.
                exposure_hours_dark = 0.
                NSB_avg = 0.
                NSB_avg_dark = 0.
                FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+SignalFlux[flux]+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
                FilePath_Folder0 += [FilePath]
                if not os.path.isfile(FilePath_Folder0[0]):continue
                GetSourceInfo(FilePath_Folder0)
                for e in range(0,len(energy_list)-1):
                    ErecS_lower_cut = energy_list[e]
                    ErecS_upper_cut = energy_list[e+1]
                    GetShowerHistogramsFromFile(FilePath_Folder0[0])
                    StackShowerHistograms()
                    NormalizeTheta2Histograms(FilePath_Folder0[0])
                    StackTheta2Histograms()
            if Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)==0.: continue
            eigenvalue_err = pow(Hist_Data_Eigenvalues_real_SumE.GetBinContent(1),0.5)
            scale = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)
            Hist_Data_Eigenvalues_real_SumE.Scale(scale)
            eigenvalue_err = eigenvalue_err*scale
            Hist_Eigenvalue_2ndRank[len(Hist_Eigenvalue_2ndRank)-1].SetBinContent(flux+1,Hist_Data_Eigenvalues_real_SumE.GetBinContent(2))
            Hist_Eigenvalue_2ndRank[len(Hist_Eigenvalue_2ndRank)-1].SetBinError(flux+1,eigenvalue_err)
            Hist_Eigenvalue_2ndRank[len(Hist_Eigenvalue_2ndRank)-1].GetXaxis().SetBinLabel(flux+1,SignalFlux_label[flux])
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in SignalFlux[flux]:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            Hist_MDM_S2B_Flux[len(Hist_MDM_S2B_Flux)-1].SetBinContent(flux+1,s2b)
            Hist_MDM_S2B_Flux[len(Hist_MDM_S2B_Flux)-1].SetBinError(flux+1,s2b_err)
            Hist_MDM_S2B_Flux[len(Hist_MDM_S2B_Flux)-1].GetXaxis().SetBinLabel(flux+1,SignalFlux_label[flux])
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            Hist_Dark_S2B_Flux[len(Hist_Dark_S2B_Flux)-1].SetBinContent(flux+1,s2b)
            Hist_Dark_S2B_Flux[len(Hist_Dark_S2B_Flux)-1].SetBinError(flux+1,s2b_err)
            Hist_Dark_S2B_Flux[len(Hist_Dark_S2B_Flux)-1].GetXaxis().SetBinLabel(flux+1,SignalFlux_label[flux])
    MakeComparisonPlotSigDist(Hist_MDM_S2B_Flux,legend_S2B,color_S2B,'signal flux in Crab unit','Systematics_MDM_Flux_%s'%(folder_tag),0.6,1.4,False,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_Flux,legend_S2B,color_S2B,'signal flux in Crab unit','Systematics_Dark_Flux_%s'%(folder_tag),0.6,1.4,False,False)
    MakeComparisonPlot(Hist_Eigenvalue_2ndRank,legend_S2B,color_S2B,'signal flux in Crab unit','size of eigenvalues','Eigenvalue_2ndRank_Flux_%s'%(folder_tag),0.1,0.4,False,False)

    Hist_Eigenvalue_Rank = []
    Hist_Eigenvector_0 = []
    Hist_Eigenvector_1 = []
    Hist_Eigenvector_2 = []
    Hist_Shower_MSCW = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for flux in range(0,len(SignalFlux)):
        Hist_Eigenvalue_Rank += [ROOT.TH1D("Hist_Eigenvalue_Rank_%s"%(flux),"",N_bins_for_deconv,0,N_bins_for_deconv)]
        Hist_Eigenvector_0 += [ROOT.TH1D("Hist_Eigenvector_0_%s"%(flux),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Eigenvector_1 += [ROOT.TH1D("Hist_Eigenvector_1_%s"%(flux),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Eigenvector_2 += [ROOT.TH1D("Hist_Eigenvector_2_%s"%(flux),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Shower_MSCW += [ROOT.TH1D("Hist_Shower_MSCW_%s"%(flux),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        legend_S2B += ['Flux %s Crab unit'%(SignalFlux_label[flux])]
        color_S2B += [color_code[(flux % 3)]+int(flux/3.)]
        ResetStackedShowerHistograms()
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            for source in range(0,len(source_list)):
                if for_syst[source]==False: continue
                source_name = source_list[source]
                FilePath_Folder0 = []
                exposure_hours = 0.
                exposure_hours_dark = 0.
                NSB_avg = 0.
                NSB_avg_dark = 0.
                FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+SignalFlux[flux]+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
                FilePath_Folder0 += [FilePath]
                if not os.path.isfile(FilePath_Folder0[0]):continue
                GetSourceInfo(FilePath_Folder0)
                for e in range(0,len(energy_list)-1):
                    ErecS_lower_cut = energy_list[e]
                    ErecS_upper_cut = energy_list[e+1]
                    GetShowerHistogramsFromFile(FilePath_Folder0[0])
                    StackShowerHistograms()
                    NormalizeTheta2Histograms(FilePath_Folder0[0])
                    StackTheta2Histograms()

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_0_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_0_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_0[len(Hist_Eigenvector_0)-1].Add(Hist_Data_Eigenvector_0_real_SumE)

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(2)
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_1_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_1_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_1[len(Hist_Eigenvector_1)-1].Add(Hist_Data_Eigenvector_1_real_SumE)

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(3)
        Hist_Data_Eigenvector_2_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_2_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_2_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_2_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_2[len(Hist_Eigenvector_2)-1].Add(Hist_Data_Eigenvector_2_real_SumE)

        bin_lower = Hist2D_Data.GetYaxis().FindBin(1.)
        bin_upper = Hist2D_Data.GetYaxis().FindBin(3.)-1
        mscw_norm = Hist_Data_MSCW_SumE.Integral(bin_lower,bin_upper)
        scale_norm = 1./mscw_norm
        Hist_Data_MSCW_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Shower_MSCW[len(Hist_Shower_MSCW)-1].Add(Hist_Data_MSCW_SumE)

        Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Add(Hist_Data_Eigenvalues_real_SumE)
        if Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)>0:
            scale = 1./Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)
            Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Scale(scale)

    for hist in range(0,len(Hist_Eigenvalue_Rank)):
        Hist_Eigenvalue_Rank[hist].GetXaxis().SetRangeUser(0,16)
    MakeComparisonPlot(Hist_Eigenvalue_Rank,legend_S2B,color_S2B,'rank','size of eigenvalues','Eigenvalue_Rank_Flux_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_0,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_0_Flux_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_1,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_1_Flux_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_2,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_2_Flux_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Shower_MSCW,legend_S2B,color_S2B,'MSCW','normalized counts','ShowerMSCW_Flux_%s'%(folder_tag),0.,0.,False,False)

def RadialAcceptance():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    ResetStackedShowerHistograms()
    exposure_hours_total = 0.
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        source_name = source_list[source]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath_Folder0 = []
            exposure_hours = 0.
            exposure_hours_dark = 0.
            NSB_avg = 0.
            NSB_avg_dark = 0.
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0 += [FilePath]
            if not os.path.isfile(FilePath_Folder0[0]):continue
            GetSourceInfo(FilePath_Folder0)
            exposure_hours_total += exposure_hours
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                GetShowerHistogramsFromFile(FilePath_Folder0[0])
                StackShowerHistograms()
                NormalizeCameraFoVHistograms(FilePath_Folder0[0])
                StackCameraFoVHistograms()

    Hist_Data_CameraFoV_SumE.Rebin2D(n_rebin,n_rebin)
    Hist_Bkgd_CameraFoV_SumE.Rebin2D(n_rebin,n_rebin)
    Hist_Bkgd_CameraFoV_Raw_SumE.Rebin2D(n_rebin,n_rebin)

    Make2DSignificanceCameraFoVPlot(exposure_hours_total,Syst_MDM,Hist_Data_CameraFoV_SumE,Hist_Bkgd_CameraFoV_SumE,'X','Y','CameraFoV%s_%s_%s'%(PercentCrab,ONOFF,folder_tag))
    Make2DSignificanceCameraFoVPlot(exposure_hours_total,Syst_MDM,Hist_Data_CameraFoV_SumE,Hist_Bkgd_CameraFoV_Raw_SumE,'X','Y','CameraFoV_Raw%s_%s_%s'%(PercentCrab,ONOFF,folder_tag))

def SystAsFunctionOfElevation():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    Hist_MDM_S2B_Elev = []
    Hist_Dark_S2B_Elev = []
    Hist_Eigenvalue_2ndRank = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        source_name = source_list[source]
        Hist_MDM_S2B_Elev += [ROOT.TH1D("Hist_MDM_S2B_Elev_%s"%(source_name),"",len(elev_range),elev_range[0][0],elev_range[len(elev_range)-1][1])]
        Hist_Dark_S2B_Elev += [ROOT.TH1D("Hist_Dark_S2B_Elev_%s"%(source_name),"",len(elev_range),elev_range[0][0],elev_range[len(elev_range)-1][1])]
        Hist_Eigenvalue_2ndRank += [ROOT.TH1D("Hist_Eigenvalue_2ndRank_%s"%(source_name),"",len(elev_range),elev_range[0][0],elev_range[len(elev_range)-1][1])]
        legend_S2B += ['%s'%(source_name)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            ResetStackedShowerHistograms()
            FilePath_Folder0 = []
            exposure_hours = 0.
            exposure_hours_dark = 0.
            NSB_avg = 0.
            NSB_avg_dark = 0.
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0 += [FilePath]
            if not os.path.isfile(FilePath_Folder0[0]):continue
            GetSourceInfo(FilePath_Folder0)
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                GetShowerHistogramsFromFile(FilePath_Folder0[0])
                StackShowerHistograms()
                NormalizeTheta2Histograms(FilePath_Folder0[0])
                StackTheta2Histograms()
            eigenvalue_err = pow(Hist_Data_Eigenvalues_real_SumE.GetBinContent(1),0.5)
            scale = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)
            Hist_Data_Eigenvalues_real_SumE.Scale(scale)
            eigenvalue_err = eigenvalue_err*scale
            Hist_Eigenvalue_2ndRank[len(Hist_Eigenvalue_2ndRank)-1].SetBinContent(elev+1,Hist_Data_Eigenvalues_real_SumE.GetBinContent(2))
            Hist_Eigenvalue_2ndRank[len(Hist_Eigenvalue_2ndRank)-1].SetBinError(elev+1,eigenvalue_err)
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            Hist_MDM_S2B_Elev[len(Hist_MDM_S2B_Elev)-1].SetBinContent(elev+1,s2b)
            Hist_MDM_S2B_Elev[len(Hist_MDM_S2B_Elev)-1].SetBinError(elev+1,s2b_err)
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            Hist_Dark_S2B_Elev[len(Hist_Dark_S2B_Elev)-1].SetBinContent(elev+1,s2b)
            Hist_Dark_S2B_Elev[len(Hist_Dark_S2B_Elev)-1].SetBinError(elev+1,s2b_err)

    MakeComparisonPlotSigDist(Hist_MDM_S2B_Elev,legend_S2B,color_S2B,'Elev.','Systematics_MDM_Elev_%s'%(folder_tag),0.6,1.4,False,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_Elev,legend_S2B,color_S2B,'Elev.','Systematics_Dark_Elev_%s'%(folder_tag),0.6,1.4,False,False)
    MakeComparisonPlot(Hist_Eigenvalue_2ndRank,legend_S2B,color_S2B,'Elev.','size of eigenvalues','Eigenvalue_2ndRank_Elev_%s'%(folder_tag),0.1,0.4,False,False)

    Hist_Eigenvalue_Rank = []
    Hist_Eigenvector_0 = []
    Hist_Eigenvector_1 = []
    Hist_Eigenvector_2 = []
    Hist_Shower_MSCW = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for elev in range(0,len(elev_range)):
        file_elev_lower = elev_range[elev][0]
        file_elev_upper = elev_range[elev][1]
        Hist_Eigenvalue_Rank += [ROOT.TH1D("Hist_Eigenvalue_Rank_%s"%(elev),"",N_bins_for_deconv,0,N_bins_for_deconv)]
        Hist_Eigenvector_0 += [ROOT.TH1D("Hist_Eigenvector_0_%s"%(elev),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Eigenvector_1 += [ROOT.TH1D("Hist_Eigenvector_1_%s"%(elev),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Eigenvector_2 += [ROOT.TH1D("Hist_Eigenvector_2_%s"%(elev),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Shower_MSCW += [ROOT.TH1D("Hist_Shower_MSCW_%s"%(elev),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        legend_S2B += ['Elev. %s-%s'%(file_elev_lower,file_elev_upper)]
        color_S2B += [color_code[(elev % 3)]+int(elev/3.)]
        ResetStackedShowerHistograms()
        for source in range(0,len(source_list)):
            if for_syst[source]==False: continue
            source_name = source_list[source]
            FilePath_Folder0 = []
            exposure_hours = 0.
            exposure_hours_dark = 0.
            NSB_avg = 0.
            NSB_avg_dark = 0.
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0 += [FilePath]
            if not os.path.isfile(FilePath_Folder0[0]):continue
            GetSourceInfo(FilePath_Folder0)
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                GetShowerHistogramsFromFile(FilePath_Folder0[0])
                StackShowerHistograms()
                NormalizeTheta2Histograms(FilePath_Folder0[0])
                StackTheta2Histograms()

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_0_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_0_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_0[len(Hist_Eigenvector_0)-1].Add(Hist_Data_Eigenvector_0_real_SumE)

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(2)
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_1_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_1_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_1[len(Hist_Eigenvector_1)-1].Add(Hist_Data_Eigenvector_1_real_SumE)

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(3)
        Hist_Data_Eigenvector_2_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_2_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_2_real_SumE.GetBinContent(binx+1))
        scale_norm = 1./abs_integral
        Hist_Data_Eigenvector_2_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_2[len(Hist_Eigenvector_2)-1].Add(Hist_Data_Eigenvector_2_real_SumE)

        bin_lower = Hist2D_Data.GetYaxis().FindBin(1.)
        bin_upper = Hist2D_Data.GetYaxis().FindBin(3.)-1
        mscw_norm = Hist_Data_MSCW_SumE.Integral(bin_lower,bin_upper)
        scale_norm = 1./mscw_norm
        Hist_Data_MSCW_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Shower_MSCW[len(Hist_Shower_MSCW)-1].Add(Hist_Data_MSCW_SumE)

        Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Add(Hist_Data_Eigenvalues_real_SumE)
        if Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)>0:
            scale = 1./Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)
            Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Scale(scale)

    for hist in range(0,len(Hist_Eigenvalue_Rank)):
        Hist_Eigenvalue_Rank[hist].GetXaxis().SetRangeUser(0,16)
    MakeComparisonPlot(Hist_Eigenvalue_Rank,legend_S2B,color_S2B,'rank','size of eigenvalues','Eigenvalue_Rank_Elev_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_0,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_0_Elev_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_1,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_1_Elev_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_2,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_2_Elev_%s'%(folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Shower_MSCW,legend_S2B,color_S2B,'MSCW','normalized counts','ShowerMSCW_Elev_%s'%(folder_tag),0.,0.,False,False)

def SystAsFunctionOfGalLat():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

    Hist_MDM_S2B_GalLat = []
    Hist_Dark_S2B_GalLat = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        gal_x = gal_coord.galactic.l.degree[source]
        gal_y = gal_coord.galactic.b.degree[source]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        Hist_MDM_S2B_GalLat += [ROOT.TH1D("Hist_MDM_S2B_GalLat_%s"%(source_name),"",9,0,90)]
        Hist_Dark_S2B_GalLat += [ROOT.TH1D("Hist_Dark_S2B_GalLat_%s"%(source_name),"",9,0,90)]
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        ResetStackedShowerHistograms()
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
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
        Hist_MDM_S2B_GalLat[len(Hist_MDM_S2B_GalLat)-1].SetBinContent(Hist_MDM_S2B_GalLat[len(Hist_MDM_S2B_GalLat)-1].FindBin(abs(gal_y)),s2b)
        Hist_MDM_S2B_GalLat[len(Hist_MDM_S2B_GalLat)-1].SetBinError(Hist_MDM_S2B_GalLat[len(Hist_MDM_S2B_GalLat)-1].FindBin(abs(gal_y)),s2b_err)
        s2b = 0.
        s2b_err = 0.
        if not "_Crab0" in PercentCrab:
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        Hist_Dark_S2B_GalLat[len(Hist_Dark_S2B_GalLat)-1].SetBinContent(Hist_Dark_S2B_GalLat[len(Hist_Dark_S2B_GalLat)-1].FindBin(abs(gal_y)),s2b)
        Hist_Dark_S2B_GalLat[len(Hist_Dark_S2B_GalLat)-1].SetBinError(Hist_Dark_S2B_GalLat[len(Hist_Dark_S2B_GalLat)-1].FindBin(abs(gal_y)),s2b_err)
    MakeComparisonPlotSigDist(Hist_MDM_S2B_GalLat,legend_S2B,color_S2B,'gal. lat.','Systematics_MDM_GalLat%s_%s'%(PercentCrab,folder_tag),0.6,1.4,False,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_GalLat,legend_S2B,color_S2B,'gal. lat.','Systematics_Dark_GalLat%s_%s'%(PercentCrab,folder_tag),0.6,1.4,False,False)

def SensitivityAsFunctionOfExposure():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

    Hist_MDM_FOM_time = []
    Hist_Dark_FOM_time = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        Hist_MDM_FOM_time += [ROOT.TH1D("Hist_MDM_FOM_time_%s"%(source_name),"",10,0,100)]
        Hist_Dark_FOM_time += [ROOT.TH1D("Hist_Dark_FOM_time_%s"%(source_name),"",10,0,100)]
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        ResetStackedShowerHistograms()
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            for path in range(1,len(FilePath_Folder0[source])):
                if not os.path.isfile(FilePath_Folder0[source][path]):continue
                GetShowerHistogramsFromFile(FilePath_Folder0[source][path])
                StackShowerHistograms()
                NormalizeTheta2Histograms(FilePath_Folder0[source][path])
                StackTheta2Histograms()
        significance = Signal_Significance(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.,10.,Syst_MDM)
        Hist_MDM_FOM_time[len(Hist_MDM_FOM_time)-1].SetBinContent(Hist_MDM_FOM_time[len(Hist_MDM_FOM_time)-1].FindBin(exposure_hours),significance)
        significance = Signal_Significance(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.,10.,Syst_Dark)
        Hist_Dark_FOM_time[len(Hist_Dark_FOM_time)-1].SetBinContent(Hist_Dark_FOM_time[len(Hist_Dark_FOM_time)-1].FindBin(exposure_hours),significance)
    MakeComparisonPlot(Hist_MDM_FOM_time,legend_S2B,color_S2B,'exposure hours','detection significance','Sensitivity_MDM_time%s_%s'%(PercentCrab,folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Dark_FOM_time,legend_S2B,color_S2B,'exposure hours','detection significance','Sensitivity_Dark_time%s_%s'%(PercentCrab,folder_tag),0.,0.,False,False)

def SystAsFunctionOfNSB():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

    Hist_MDM_S2B_NSB = []
    Hist_Dark_S2B_NSB = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        Hist_MDM_S2B_NSB += [ROOT.TH1D("Hist_MDM_S2B_NSB_%s"%(source_name),"",4,4.5,8.5)]
        Hist_Dark_S2B_NSB += [ROOT.TH1D("Hist_Dark_S2B_NSB_%s"%(source_name),"",4,4.5,8.5)]
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        ResetStackedShowerHistograms()
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
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
        Hist_MDM_S2B_NSB[len(Hist_MDM_S2B_NSB)-1].SetBinContent(Hist_MDM_S2B_NSB[len(Hist_MDM_S2B_NSB)-1].FindBin(NSB_avg),s2b)
        Hist_MDM_S2B_NSB[len(Hist_MDM_S2B_NSB)-1].SetBinError(Hist_MDM_S2B_NSB[len(Hist_MDM_S2B_NSB)-1].FindBin(NSB_avg),s2b_err)
        s2b = 0.
        s2b_err = 0.
        if not "_Crab0" in PercentCrab:
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
        Hist_Dark_S2B_NSB[len(Hist_Dark_S2B_NSB)-1].SetBinContent(Hist_Dark_S2B_NSB[len(Hist_Dark_S2B_NSB)-1].FindBin(NSB_avg),s2b)
        Hist_Dark_S2B_NSB[len(Hist_Dark_S2B_NSB)-1].SetBinError(Hist_Dark_S2B_NSB[len(Hist_Dark_S2B_NSB)-1].FindBin(NSB_avg),s2b_err)
    MakeComparisonPlotSigDist(Hist_MDM_S2B_NSB,legend_S2B,color_S2B,'NSB','Systematics_MDM_NSB%s_%s'%(PercentCrab,folder_tag),0.6,1.4,False,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_NSB,legend_S2B,color_S2B,'NSB','Systematics_Dark_NSB%s_%s'%(PercentCrab,folder_tag),0.6,1.4,False,False)

def SystAsFunctionOfEnergy():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

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
        GetSourceInfo(FilePath_Folder0[source])
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
    MakeComparisonPlotSigDist(Hist_MDM_S2B_Energy,legend_S2B,color_S2B,'E [GeV]','Systematics_MDM_Energy%s_%s'%(PercentCrab,folder_tag),0.6,1.4,True,False)
    MakeComparisonPlotSigDist(Hist_Dark_S2B_Energy,legend_S2B,color_S2B,'E [GeV]','Systematics_Dark_Energy%s_%s'%(PercentCrab,folder_tag),0.6,1.4,True,False)

def SystDarkVsMDM():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark
    global Syst_MDM
    global Syst_Dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if for_syst[source]==False: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

    Hist_Eigenvalue_Rank = []
    Hist_Eigenvector_0 = []
    Hist_Eigenvector_1 = []
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
        GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        if exposure_hours==0: continue
        Hist_S2B += [ROOT.TH1D("Hist_S2B_%s"%(source_name),"",2,0,2)]
        Hist_S2B[len(Hist_S2B)-1].GetXaxis().SetBinLabel(1,'Initial')
        Hist_S2B[len(Hist_S2B)-1].GetXaxis().SetBinLabel(2,'MDM')
        Hist_Eigenvalue_Rank += [ROOT.TH1D("Hist_Eigenvalue_Rank_%s"%(source_name),"",N_bins_for_deconv,0,N_bins_for_deconv)]
        Hist_Eigenvector_0 += [ROOT.TH1D("Hist_Eigenvector_0_%s"%(source_name),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
        Hist_Eigenvector_1 += [ROOT.TH1D("Hist_Eigenvector_1_%s"%(source_name),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]
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

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(1)
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_0_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_0_real_SumE.GetBinContent(binx+1))
        scale_norm += 1./abs_integral
        Hist_Data_Eigenvector_0_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_0[len(Hist_Eigenvector_0)-1].Add(Hist_Data_Eigenvector_0_real_SumE)

        scale_norm = 1./Hist_Data_Eigenvalues_real_SumE.GetBinContent(2)
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same sign
        abs_integral = 0.
        for binx in range(0,Hist_Data_Eigenvector_1_real_SumE.GetNbinsX()):
            abs_integral += abs(Hist_Data_Eigenvector_1_real_SumE.GetBinContent(binx+1))
        scale_norm += 1./abs_integral
        Hist_Data_Eigenvector_1_real_SumE.Scale(scale_norm) # this makes sure they all have the same normalization
        Hist_Eigenvector_1[len(Hist_Eigenvector_1)-1].Add(Hist_Data_Eigenvector_1_real_SumE)

        Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Add(Hist_Data_Eigenvalues_real_SumE)
        if Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)>0:
            scale = 1./Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)
            Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Scale(scale)

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

    Syst_Dark, Syst_MDM = MakeComparisonPlotTwoColumn(Hist_S2B,legend_S2B,color_S2B,'methods','Systematics_DarkVsMDM%s_%s'%(PercentCrab,folder_tag),0.6,1.4,False,False)
    for hist in range(0,len(Hist_Eigenvalue_Rank)):
        Hist_Eigenvalue_Rank[hist].GetXaxis().SetRangeUser(0,16)
    MakeComparisonPlot(Hist_Eigenvalue_Rank,legend_S2B,color_S2B,'rank','size of eigenvalues','Eigenvalue_Rank%s_%s'%(PercentCrab,folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_0,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_0%s_%s'%(PercentCrab,folder_tag),0.,0.,False,False)
    MakeComparisonPlot(Hist_Eigenvector_1,legend_S2B,color_S2B,'columns','size of eigenvector','Eigenvector_1%s_%s'%(PercentCrab,folder_tag),0.,0.,False,False)

def CompareDiffFolders():

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    FilePath_AllFolders = []
    for folder in range(0,len(FileFolder)):
        folder_path = FileFolder[folder]
        folder_tag = FileTag[folder]
        folder_label = FileLabel[folder]
        FilePath_Folder0 = []
        for source in range(0,len(source_list)):
            if for_syst[source]==False: continue
            FilePath_Folder0 += [[source_list[source]]]
            for elev in range(0,len(elev_range)):
                file_elev_lower = elev_range[elev][0]
                file_elev_upper = elev_range[elev][1]
                FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
                FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]
        FilePath_AllFolders += [FilePath_Folder0]

    Hist_MDM_S2B = []
    Hist_Dark_S2B = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]
    for source in range(0,len(FilePath_AllFolders[0])):
        source_name = FilePath_AllFolders[0][source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        GetSourceInfo(FilePath_AllFolders[0][source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        if exposure_hours==0: continue
        Hist_MDM_S2B += [ROOT.TH1D("Hist_MDM_S2B_%s"%(source_name),"",len(FilePath_AllFolders),0,len(FilePath_AllFolders))]
        Hist_Dark_S2B += [ROOT.TH1D("Hist_Dark_S2B_%s"%(source_name),"",len(FilePath_AllFolders),0,len(FilePath_AllFolders))]
        legend_S2B += ['%s (%0.1f hrs)'%(source_name,exposure_hours)]
        color_S2B += [color_code[(source % 3)]+int(source/3.)]
        for folder in range(0,len(FileFolder)):
            folder_path = FileFolder[folder]
            folder_tag = FileTag[folder]
            folder_label = FileLabel[folder]
            ResetStackedShowerHistograms()
            for path in range(1,len(FilePath_AllFolders[folder][source])):
                if not os.path.isfile(FilePath_AllFolders[folder][source][path]):continue
                for e in range(0,len(energy_list)-1):
                    ErecS_lower_cut = energy_list[e]
                    ErecS_upper_cut = energy_list[e+1]
                    GetShowerHistogramsFromFile(FilePath_AllFolders[folder][source][path])
                    StackShowerHistograms()
                    NormalizeTheta2Histograms(FilePath_AllFolders[folder][source][path])
                    StackTheta2Histograms()
                #PlotsStackedHistograms('%s_%s_%s'%(source_name,path,folder_tag))
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Bkgd_Theta2_SumE,0.2,10.)
            Hist_MDM_S2B[len(Hist_MDM_S2B)-1].SetBinContent(folder+1,s2b)
            Hist_MDM_S2B[len(Hist_MDM_S2B)-1].SetBinError(folder+1,s2b_err)
            Hist_MDM_S2B[len(Hist_MDM_S2B)-1].GetXaxis().SetBinLabel(folder+1,folder_label)
            s2b = 0.
            s2b_err = 0.
            if not "_Crab0" in PercentCrab:
                s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            else:
                s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumE,Hist_Dark_Theta2_SumE,0.2,10.)
            Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinContent(folder+1,s2b)
            Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinError(folder+1,s2b_err)
            Hist_Dark_S2B[len(Hist_Dark_S2B)-1].GetXaxis().SetBinLabel(folder+1,folder_label)
            Hist_Dark_S2B[len(Hist_Dark_S2B)-1].GetXaxis().SetLabelSize(0.1)

    plot_name = 'Systematics_MDM_DiffCuts'
    for ntag in range(0,len(FileTag)):
        plot_name += '_%s'%(FileTag[ntag])
    MakeComparisonPlotSigDist(Hist_MDM_S2B,legend_S2B,color_S2B,'',plot_name,0.6,1.4,False,False)
    plot_name = 'Systematics_Dark_DiffCuts'
    for ntag in range(0,len(FileTag)):
        plot_name += '_%s'%(FileTag[ntag])
    MakeComparisonPlotSigDist(Hist_Dark_S2B,legend_S2B,color_S2B,'',plot_name,0.6,1.4,False,False)

def SkymapHighlightIntegral(hist_data,hist_highlight):

    integral = 0.
    for bx in range(0,hist_data.GetNbinsX()):
        for by in range(0,hist_data.GetNbinsY()):
            if not hist_highlight.GetBinContent(bx+1,by+1)==0.:
                integral += hist_data.GetBinContent(bx+1,by+1)

    return integral

def Event_rate(count,time):

    err = pow(count,0.5)
    return count/time, err/time

def SingleSourceSpectrum(source_name_input):

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    Hist_MDM_NormSig_Rate = []
    Hist_MDM_InclBkg_Rate = []
    Hist_MDM_Sig_Rate = []
    Hist_MDM_Bkg_Rate = []
    legend_S2B = []
    color_S2B = []
    color_code = [ROOT.kBlue,ROOT.kGreen,ROOT.kRed]

    source_name = source_name_input
    binx_lower = Hist2D_Data.GetYaxis().FindBin(MSCL_lower_cut)
    binx_upper = Hist2D_Data.GetYaxis().FindBin(MSCL_blind_cut)-1
    biny_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
    biny_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1

    for elev in range(0,len(elev_range)):
        file_elev_lower = elev_range[elev][0]
        file_elev_upper = elev_range[elev][1]
        FilePath_Folder0 = []
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        FilePath = "%s/Netflix_"%(folder_path)+source_name+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
        FilePath_Folder0 += [FilePath]
        if not os.path.isfile(FilePath_Folder0[0]): continue
        GetSourceInfo(FilePath_Folder0)
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        if exposure_hours==0: continue
        Hist_MDM_InclBkg_Rate += [ROOT.TH1D("Hist_MDM_InclBkg_Rate_%s"%(file_elev_upper),"",len(energy_list)-1,array('d',energy_list))]
        Hist_MDM_Bkg_Rate += [ROOT.TH1D("Hist_MDM_Bkg_Rate_%s"%(file_elev_upper),"",len(energy_list)-1,array('d',energy_list))]
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            GetShowerHistogramsFromFile(FilePath_Folder0[0])
            NormalizeSkyMapHistograms(FilePath_Folder0[0])
            incl_bkg_integral = Hist2D_Bkgd.Integral(binx_lower,binx_upper,biny_lower,biny_upper)
            incl_bkg_rate, incl_bkg_rate_err = Event_rate(incl_bkg_integral,exposure_hours)
            Hist_MDM_InclBkg_Rate[len(Hist_MDM_InclBkg_Rate)-1].SetBinContent(e+1,incl_bkg_rate)
            Hist_MDM_InclBkg_Rate[len(Hist_MDM_InclBkg_Rate)-1].SetBinError(e+1,incl_bkg_rate_err)
            Hist_temp = Hist_Bkgd_Skymap.Clone()
            Hist_temp.Rebin2D(n_rebin,n_rebin)
            bkg_integral = SkymapHighlightIntegral(Hist_temp,Hist_Highlight_Bias_Skymap)
            bkg_rate, bkg_rate_err = Event_rate(bkg_integral,exposure_hours)
            Hist_MDM_Bkg_Rate[len(Hist_MDM_Bkg_Rate)-1].SetBinContent(e+1,bkg_rate)
            Hist_MDM_Bkg_Rate[len(Hist_MDM_Bkg_Rate)-1].SetBinError(e+1,bkg_rate_err)

    for elev in range(0,len(elev_range)):
        file_elev_lower = elev_range[elev][0]
        file_elev_upper = elev_range[elev][1]
        FilePath_Folder0 = []
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        FilePath = "%s/Netflix_"%(folder_path)+source_name+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
        FilePath_Folder0 += [FilePath]
        if not os.path.isfile(FilePath_Folder0[0]): continue
        GetSourceInfo(FilePath_Folder0)
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        if exposure_hours==0: continue
        Hist_MDM_NormSig_Rate += [ROOT.TH1D("Hist_MDM_NormSig_Rate_%s"%(file_elev_upper),"",len(energy_list)-1,array('d',energy_list))]
        Hist_MDM_Sig_Rate += [ROOT.TH1D("Hist_MDM_Sig_Rate_%s"%(file_elev_upper),"",len(energy_list)-1,array('d',energy_list))]
        legend_S2B += ['%s elev. %s-%s'%(source_name,file_elev_lower,file_elev_upper)]
        color_S2B += [color_code[(elev % 3)]-2*int(elev/3.)]
        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            GetShowerHistogramsFromFile(FilePath_Folder0[0])
            NormalizeSkyMapHistograms(FilePath_Folder0[0])
            Hist_temp = Hist_Data_Skymap.Clone()
            Hist_temp.Rebin2D(n_rebin,n_rebin)
            data_integral = SkymapHighlightIntegral(Hist_temp,Hist_Highlight_Bias_Skymap)
            Hist_temp = Hist_Bkgd_Skymap.Clone()
            Hist_temp.Rebin2D(n_rebin,n_rebin)
            bkg_integral = SkymapHighlightIntegral(Hist_temp,Hist_Highlight_Bias_Skymap)
            incl_bkg_integral = Hist2D_Bkgd.Integral(binx_lower,binx_upper,biny_lower,biny_upper)
            data_rate, data_rate_err = Event_rate(data_integral,exposure_hours)
            bkg_rate, bkg_rate_err = Event_rate(bkg_integral,exposure_hours)
            incl_bkg_rate, incl_bkg_rate_err = Event_rate(incl_bkg_integral,exposure_hours)
            sig_rate = data_rate-bkg_rate
            sig_rate_err = pow(data_rate_err*data_rate_err+bkg_rate_err*bkg_rate_err,0.5)
            Hist_MDM_Sig_Rate[len(Hist_MDM_Sig_Rate)-1].SetBinContent(e+1,sig_rate)
            Hist_MDM_Sig_Rate[len(Hist_MDM_Sig_Rate)-1].SetBinError(e+1,sig_rate_err)
            if bkg_rate!=0:
                norm_sig_rate = sig_rate/incl_bkg_rate*Hist_MDM_InclBkg_Rate[len(Hist_MDM_InclBkg_Rate)-1].GetBinContent(e+1)
                #norm_sig_rate = sig_rate/bkg_rate*Hist_MDM_Bkg_Rate[len(Hist_MDM_InclBkg_Rate)-1].GetBinContent(e+1)
                Hist_MDM_NormSig_Rate[len(Hist_MDM_NormSig_Rate)-1].SetBinContent(e+1,norm_sig_rate)
                Hist_MDM_NormSig_Rate[len(Hist_MDM_NormSig_Rate)-1].SetBinError(e+1,sig_rate_err)
    MakeComparisonPlot(Hist_MDM_NormSig_Rate,legend_S2B,color_S2B,'E [GeV]','events/hour (normalized)','SignalFluxNorm_MDM%s_%s'%(PercentCrab,folder_tag),0.,0.,True,True)
    MakeComparisonPlot(Hist_MDM_Sig_Rate,legend_S2B,color_S2B,'E [GeV]','events/hour','SignalFlux_MDM%s_%s'%(PercentCrab,folder_tag),0.,0.,True,True)
    MakeComparisonPlot(Hist_MDM_Bkg_Rate,legend_S2B,color_S2B,'E [GeV]','events/hour','BkgdFlux_MDM%s_%s'%(PercentCrab,folder_tag),0.,0.,True,True)
    MakeComparisonPlot(Hist_MDM_InclBkg_Rate,legend_S2B,color_S2B,'E [GeV]','events/hour','InclBkgdFlux_MDM%s_%s'%(PercentCrab,folder_tag),0.,0.,True,True)

def SingleSourceSkyMap(source_name_input):

    global ErecS_lower_cut
    global ErecS_upper_cut
    global exposure_hours
    global exposure_hours_dark
    global NSB_avg
    global NSB_avg_dark

    folder_path = FileFolder[0]
    folder_tag = FileTag[0]
    folder_label = FileLabel[0]

    FilePath_Folder0 = []
    for source in range(0,len(source_list)):
        if not source_name_input==source_list[source]: continue
        FilePath_Folder0 += [[source_list[source]]]
        for elev in range(0,len(elev_range)):
            file_elev_lower = elev_range[elev][0]
            file_elev_upper = elev_range[elev][1]
            FilePath = "%s/Netflix_"%(folder_path)+source_list[source]+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(ONOFF)+".root";
            FilePath_Folder0[len(FilePath_Folder0)-1] += [FilePath]

    for source in range(0,len(FilePath_Folder0)):
        source_name = FilePath_Folder0[source][0]
        exposure_hours = 0.
        exposure_hours_dark = 0.
        NSB_avg = 0.
        NSB_avg_dark = 0.
        GetSourceInfo(FilePath_Folder0[source])
        print 'Source %s, exposure = %s, NSB = %s'%(source_name,exposure_hours,NSB_avg)
        if exposure_hours==0: continue
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
                NormalizeSkyMapHistograms(FilePath_Folder0[source][path])
                StackSkymapHistograms()

        ErecS_lower_cut = energy_list[0]
        PlotsStackedHistograms('%s%s_%s_%s'%(source_name,PercentCrab,ONOFF,folder_tag))

        Hist_Data_Skymap_SumE.Rebin2D(n_rebin,n_rebin)
        Hist_Dark_Skymap_SumE.Rebin2D(n_rebin,n_rebin)
        Hist_Bkgd_Skymap_SumE.Rebin2D(n_rebin,n_rebin)
        Hist_Bkgd_Skymap_Raw_SumE.Rebin2D(n_rebin,n_rebin)
        Hist_Highlight_Skymap.Rebin2D(n_rebin,n_rebin)
        Hist_Highlight_Bias_Skymap.Rebin2D(n_rebin,n_rebin)

        Hist_Data_Skymap_smooth = Smooth2DMap(Hist_Data_Skymap_SumE,smooth_size,False)
        Hist_Bkgd_Skymap_smooth = Smooth2DMap(Hist_Bkgd_Skymap_SumE,smooth_size,False)

        Hist_Dark_Skymap_smooth = Smooth2DMap(Hist_Dark_Skymap_SumE,smooth_size,False)
        Make2DSignificancePlot(Syst_Dark,Hist_Data_Skymap_SumE,Hist_Dark_Skymap_SumE,Hist_Data_Skymap_smooth,Hist_Dark_Skymap_smooth,'RA','Dec','Skymap_Dark_%s%s_%s_%s'%(source_name,PercentCrab,ONOFF,folder_tag))

        Hist_Bkgd_Skymap_Raw_smooth = Smooth2DMap(Hist_Bkgd_Skymap_Raw_SumE,smooth_size,False)
        Make2DSignificancePlot(Syst_MDM,Hist_Data_Skymap_SumE,Hist_Bkgd_Skymap_Raw_SumE,Hist_Data_Skymap_smooth,Hist_Bkgd_Skymap_Raw_smooth,'RA','Dec','Skymap_MDM_Raw_%s%s_%s_%s'%(source_name,PercentCrab,ONOFF,folder_tag))

        Make2DSignificancePlot(Syst_MDM,Hist_Data_Skymap_SumE,Hist_Bkgd_Skymap_SumE,Hist_Data_Skymap_smooth,Hist_Bkgd_Skymap_smooth,'RA','Dec','Skymap_MDM_%s%s_%s_%s'%(source_name,PercentCrab,ONOFF,folder_tag))

def HMS2deg(ra='', dec=''):
    print 'ra = %s'%(ra)
    print 'dec = %s'%(dec)
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)           
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC


SystDarkVsMDM() # run this to get method systematics
#SystAsFunctionOfEnergy()
#SystAsFunctionOfElevation()
#SystAsFunctionOfNSB()
#SystAsFunctionOfGalLat()
#SystAsFunctionOfSignal()
#SensitivityAsFunctionOfExposure()
#CompareDiffFolders()

#source_of_interest = 'Proton_NSB200'
#source_of_interest = 'Proton_NSB750'
#source_of_interest = 'Crab'
source_of_interest = 'Mrk421'
#source_of_interest = 'H1426'
#source_of_interest = 'PKS1424'
#source_of_interest = '3C264'
#source_of_interest = 'OJ287V6'
#source_of_interest = '1ES0229'
#source_of_interest = 'S3_1227_V6'
#source_of_interest = 'MS1221V6'
#source_of_interest = 'PKS1441V6'
#source_of_interest = 'RBS0413V6'
#source_of_interest = 'PG1553V6'
#source_of_interest = 'Segue1V6'
#source_of_interest = '1ES1011V6'
#source_of_interest = 'NGC1275V6'
#source_of_interest = '1ES0647V6'
#source_of_interest = '1ES1440V6'
#source_of_interest = '1ES1741V6'
#source_of_interest = 'IC443HotSpot'
#source_of_interest = 'RGBJ0710'
#source_of_interest = 'CasA'
#source_of_interest = 'M82'
#source_of_interest = 'G079'
#source_of_interest = 'WComaeV6'
#source_of_interest = '1ES1218V6'
#source_of_interest = 'MGRO_J1908_V6'
#source_of_interest = 'MGRO_J1908_V5'
#source_of_interest = 'Segue1V5'
#source_of_interest = 'IC443HotSpotV5'
#source_of_interest = 'GemingaV6'
#source_of_interest = 'GemingaV5'

source_idx = FindSourceIndex(source_of_interest)
source_ra = round(float(HMS2deg(sky_coord[source_idx].split('+')[0],sky_coord[source_idx].split('+')[1])[0]),3)
source_dec = round(float(HMS2deg(sky_coord[source_idx].split('+')[0],sky_coord[source_idx].split('+')[1])[1]),3)
Hist_Data_Skymap = ROOT.TH2D("Hist_Data_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Dark_Skymap = ROOT.TH2D("Hist_Dark_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Bkgd_Skymap = ROOT.TH2D("Hist_Bkgd_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Bkgd_Skymap_Raw = ROOT.TH2D("Hist_Bkgd_Skymap_Raw","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Data_Skymap_SumE = ROOT.TH2D("Hist_Data_Skymap_SumE","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Dark_Skymap_SumE = ROOT.TH2D("Hist_Dark_Skymap_SumE","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Bkgd_Skymap_SumE = ROOT.TH2D("Hist_Bkgd_Skymap_SumE","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Bkgd_Skymap_Raw_SumE = ROOT.TH2D("Hist_Bkgd_Skymap_Raw_SumE","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Highlight_Skymap = ROOT.TH2D("Hist_Highlight_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Highlight_Bias_Skymap = ROOT.TH2D("Hist_Highlight_Bias_Skymap","",150,source_ra-3,source_ra+3,150,source_dec-3,source_dec+3)
Hist_Data_CameraFoV = ROOT.TH2D("Hist_Data_CameraFoV","",150,-3,3,150,-3,3)
Hist_Bkgd_CameraFoV = ROOT.TH2D("Hist_Bkgd_CameraFoV","",150,-3,3,150,-3,3)
Hist_Bkgd_CameraFoV_Raw = ROOT.TH2D("Hist_Bkgd_CameraFoV_Raw","",150,-3,3,150,-3,3)
Hist_Data_CameraFoV_SumE = ROOT.TH2D("Hist_Data_CameraFoV_SumE","",150,-3,3,150,-3,3)
Hist_Bkgd_CameraFoV_SumE = ROOT.TH2D("Hist_Bkgd_CameraFoV_SumE","",150,-3,3,150,-3,3)
Hist_Bkgd_CameraFoV_Raw_SumE = ROOT.TH2D("Hist_Bkgd_CameraFoV_Raw_SumE","",150,-3,3,150,-3,3)

#n_rebin = 8
#smooth_size = 0.1
#n_rebin = 4
#smooth_size = 0.1
n_rebin = 2
smooth_size = 0.1

theta2_range = 1.0
highlight_threshold = 3.0

#ONOFF = "OFF"
#PercentCrab = "_Crab0"
#RadialAcceptance()

ONOFF = "ON"
#ONOFF = "OFF"

PercentCrab = "_Crab0"
#PercentCrab = "_Crab10"
#PercentCrab = "_Crab20"
#PercentCrab = "_Crab50"
#PercentCrab = "_Crab100"
#PercentCrab = "_Crab200"

SingleSourceSkyMap(source_of_interest)
SingleSourceSpectrum(source_of_interest)
