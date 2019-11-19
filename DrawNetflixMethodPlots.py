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

PlotAllEnergies = False
PlotAllSRs = False
DoSkymap = True
DoSmooth = False
UseRDBM = False
UseDark = False
UseRing = False

theta2_lower = 0.0
theta2_upper = 10.0

#file_elev_lower = 55
#file_elev_upper = 65
#file_elev_lower = 65
#file_elev_upper = 75
file_elev_lower = 75
file_elev_upper = 85

gal_lat_cut = 0.

elev_range = []
elev_range += [[75,85]]
elev_range += [[65,75]]
elev_range += [[55,65]]
elev_range += [[45,55]]
#elev_range += [[35,45]]
#elev_range += [[25,35]]

PercentCrab = ""
#PercentCrab = "_Crab0"
#PercentCrab = "_Crab50"
#file_tag = "ON"
file_tag = "OFF"

file_theta2_lower = 0.0
file_theta2_upper = 0.2

if "ON" in file_tag:
    file_theta2_lower = 0.0
    file_theta2_upper = 10.0
    UseRDBM = True
    UseDark = True
if "OFF" in file_tag:
    file_theta2_lower = 0.2
    file_theta2_upper = 10.0
    UseRDBM = True
    UseDark = True

tag = file_tag

FileFolder = 'output_unblind_4x4_3tel_tight'
tag += '_unblind_4x4_3tel_tight'
#FileFolder = 'output_unblind_4x4_3tel_medium'
#tag += '_unblind_4x4_3tel_medium'
#FileFolder = 'output_unblind_4x4_3tel_loose'
#tag += '_unblind_4x4_3tel_loose'
#FileFolder = 'output_unblind_4x4_3tel'
#tag += '_unblind_4x4_3tel'
#FileFolder = 'output_unblind_5x5_3tel_tight'
#tag += '_unblind_5x5_3tel_tight'
#FileFolder = 'output_unblind_16x16_3tel'
#tag += '_unblind_16x16_3tel'
#FileFolder = 'output_unique'
#tag += '_unique'
#FileFolder = 'output_eigen'
#tag += '_eigen'
#FileFolder = 'output_test'
#tag += '_test'

#tag += "_Data"
tag += "_MC"

target = ""
source = []
sky_coord = []
source += ['Proton_Crab0']
sky_coord += ['10 07 04 +16 04 55']
#source += ['Proton_Crab50']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['Proton_Crab100']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['Proton_Crab200']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['Proton_Crab400']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['Crab']
#sky_coord += ['05 34 31.97 +22 00 52.1']
#source += ['Mrk421']
#sky_coord += ['11 04 19 +38 11 41']
#source += ['H1426']
#sky_coord += ['14 28 32.609 +42 40 21.05']
#source += ['1ES0229']
#sky_coord += ['02 32 53.2 +20 16 21']
#source += ['PKS1424']
#sky_coord += ['14 27 00 +23 47 00']
#source += ['3C264']
#sky_coord += ['11 45 5.009 +19 36 22.74']
#source += ['OJ287V6']
#sky_coord += ['08 54 49.1 +20 05 58.89']
#source += ['S3_1227_V6']
#sky_coord += ['12 30 14.1 +25 18 07']
#source += ['MS1221V6']
#sky_coord += ['12 24 24.2 +24 36 24']
#source += ['PKS1441V6']
#sky_coord += ['14 43 56.9 +25 01 44']
#source += ['RBS0413V6']
#sky_coord += ['03 19 47 +18 45 42']
#source += ['PG1553V6']
#sky_coord += ['15 55 44.7 +11 11 41']
#source += ['Segue1V6']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['1ES1011V6']
#sky_coord += ['10 15 4.139 +49 26 0.71']
#source += ['NGC1275V6']
#sky_coord += ['03 19 48.1 +41 30 42']
#source += ['1ES0647V6']
#sky_coord += ['06 50 46.490 +25 02 59.62']
#source += ['1ES1440V6']
#sky_coord += ['14 42 48.277 +12 00 40.37']
#source += ['1ES1741V6']
#sky_coord += ['17 44 01.2 +19 32 47']
#source += ['IC443HotSpot']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#source += ['RGBJ0710']
#sky_coord += ['07 10 26.4 +59 09 00']
#source += ['CasA']
#sky_coord += ['23 23 13.8 +58 48 26']
#source += ['WComaeV6']
#sky_coord += ['12 21 31.7 +28 13 59']
#source += ['M82']
#sky_coord += ['09 55 52.7 +69 40 46']
#source += ['G079']
#sky_coord += ['20 32 28.56 +40 19 41.52']
#source += ['1ES1218V6']
#sky_coord += ['12 21 26.3 +30 11 29']
#source += ['MGRO_J1908_V6']
#sky_coord += ['19 07 54 +06 16 07']
#source += ['MGRO_J1908_V5']
#sky_coord += ['19 07 54 +06 16 07']
#source += ['Segue1V5']
#sky_coord += ['10 07 04 +16 04 55']
#source += ['IC443HotSpotV5']
#sky_coord += ['06 18 2.700 +22 39 36.00']
#source += ['Test']
#sky_coord += ['14 28 32.609 +42 40 21.05']

N_bins_for_deconv = 40
MSCW_plot_lower = -1.
MSCW_plot_upper = 3.
MSCL_plot_lower = -1.
MSCL_plot_upper = 3.

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
energy_list += [8913]

MSCW_lower_cut = -1.0
MSCW_upper_cut = 1.0
MSCW_blind_cut = 1.0
MSCL_lower_cut = -1.0
MSCL_upper_cut = 1.0
MSCL_blind_cut = 1.0
exposure_hours = 0.
total_exposure_hours = 0.

NRGBs = 5
#NCont = 512
NCont = 104
stops = [0.00,0.45,0.50,0.55,1.00]
red =   [0.00,1.00,1.00,1.00,1.00]
green = [0.00,1.00,1.00,1.00,0.00]
blue =  [1.00,1.00,1.00,1.00,0.00]
ROOT.TColor.CreateGradientColorTable(NRGBs,array('d',stops),array('d',red),array('d',green),array('d',blue),NCont)
ROOT.gStyle.SetNumberContours(NCont)

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

def AddStat2D(Hist_Data,Hist_Bkgd):

    Hist_Stat = Hist_Bkgd.Clone()
    for bx in range(1,Hist_Bkgd.GetNbinsX()+1):
        for by in range(1,Hist_Bkgd.GetNbinsY()+1):
            data_err = Hist_Data.GetBinError(bx,by)
            Hist_Stat.SetBinError(bx,by,data_err)

    return Hist_Stat

def AddSyst2D(Hist_Data,Hist_Dark,Hist_Bkgd):

    bin_lower_x = Hist_Bkgd.GetXaxis().FindBin(MSCL_lower_cut)
    bin_upper_x = Hist_Bkgd.GetXaxis().FindBin(MSCL_blind_cut)-1
    bin_lower_y = Hist_Bkgd.GetYaxis().FindBin(MSCW_lower_cut)
    bin_upper_y = Hist_Bkgd.GetYaxis().FindBin(MSCW_blind_cut)-1
    bkgd_cr_total = 0
    syst_cr = 0
    nbins_cr = 0

    Hist_sys = Hist_Bkgd.Clone()
    for bx in range(1,Hist_Bkgd.GetNbinsX()+1):
        for by in range(1,Hist_Bkgd.GetNbinsY()+1):
            data_cr_x = Hist_Data.Integral(bx,bx,bin_upper_y,Hist_Bkgd.GetNbinsY())
            bkgd_cr_x = Hist_Dark.Integral(bx,bx,bin_upper_y,Hist_Bkgd.GetNbinsY())
            data_cr_y = Hist_Data.Integral(bin_upper_x,Hist_Bkgd.GetNbinsX(),by,by)
            bkgd_cr_y = Hist_Dark.Integral(bin_upper_x,Hist_Bkgd.GetNbinsX(),by,by)
            syst_cr = 0.
            if not (abs(data_cr_x)+abs(data_cr_y)+abs(bkgd_cr_x)+abs(bkgd_cr_y))==0:
                syst_cr = 2.*(abs(data_cr_x-bkgd_cr_x)+abs(data_cr_y-bkgd_cr_y))/(abs(data_cr_x)+abs(data_cr_y)+abs(bkgd_cr_x)+abs(bkgd_cr_y))
            new_err = abs(Hist_Bkgd.GetBinContent(bx,by))*syst_cr
            old_err = Hist_Bkgd.GetBinError(bx,by)
            Hist_sys.SetBinError(bx,by,pow(new_err*new_err+old_err*old_err,0.5))

    return Hist_sys

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

def GetSkyScaleFactor(Hist_Bkg_MSCW,MSCW_cut_lower,MSCW_cut_upper,Hist_CR_theta2):

    norm_bin_low_target = Hist_Bkg_MSCW.FindBin(MSCW_cut_lower)
    norm_bin_up_target = Hist_Bkg_MSCW.FindBin(MSCW_cut_upper)-1
    bkg_total, bkg_err = IntegralAndError(Hist_Bkg_MSCW,norm_bin_low_target,norm_bin_up_target)
    #bkg_total, bkg_err = IntegralAndSystError(Hist_Bkg_MSCW,norm_bin_low_target,norm_bin_up_target,-1)
    cr_total, cr_err = IntegralAndError(Hist_CR_theta2,1,Hist_CR_theta2.GetNbinsX())
    scale = 0
    scale_err = 0
    if not cr_total==0 and not bkg_total==0:
        scale = bkg_total/cr_total
        scale_err = scale*(bkg_err/bkg_total)
    return scale, scale_err

def RaDecHistScale(Hist,Hist_Syst,scale,scale_err):

    #Hist.Scale(scale)

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
            Hist_Syst.SetBinContent(bx,by,new_error)
            Hist_Syst.SetBinError(bx,by,0)

def Theta2HistScale(Hist,Hist_Syst,scale,scale_err):

    for b in range(1,Hist.GetNbinsX()+1):
        old_content = Hist.GetBinContent(b)
        old_error = Hist.GetBinError(b)
        new_content = old_content*scale
        new_error = 0
        if old_content>0 and scale>0:
            new_error = new_content*(old_error/old_content)
        Hist.SetBinContent(b,new_content)
        Hist.SetBinError(b,new_error)
        new_error = 0
        if old_content>0 and scale>0:
            new_error = new_content*(scale_err/scale)
        Hist_Syst.SetBinContent(b,new_error)
        Hist_Syst.SetBinError(b,0)

def Event_rate(Hist_SR,time):

    err_SR = 0
    data_SR = Hist_SR.Integral()
    err_SR = pow(data_SR,0.5)
    return data_SR/time, err_SR/time

def Variation_ratio(Hist_SR, Hist_Bkg,range_lower,range_upper,syst):

    range_lower = 0.4

    norm_bin_low_target = Hist_SR.FindBin(range_lower)
    norm_bin_up_target = Hist_Bkg.FindBin(range_upper)-1
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hist_SR,norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    #predict_bkg, err_bkg = IntegralAndError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target)
    predict_bkg, err_bkg = IntegralAndSystError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target,syst)
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0 and not data_SR==0:
        sbratio = (data_SR)/(predict_bkg)
        sbratio_err = (data_SR)/(predict_bkg)*pow(pow(err_SR/data_SR,2)+pow(err_bkg/predict_bkg,2),0.5)
    return sbratio, sbratio_err

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeComparisonPlot(Hists,legends,colors,title_x,title_y,name,logx,logy):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
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

    #Hists[max_hist].SetMaximum(max_heigh)
    #Hists[max_hist].SetMinimum(min_heigh)
    Hists[max_hist].Draw("E")
    #Hists[0].Draw("E")

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
    legend.SetTextSize(0.15)
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

    c_both.SaveAs('output_plots/%s_%s.png'%(name,tag))

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
    for h in range(0,len(Hists)):
        mean += [Hists[h].GetMean()]
        rms += [Hists[h].GetRMS()]
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            if legends[h]=='Bkg (OFF)':
                Hists[h].SetLineWidth(3)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #if not 'Energy' in name:
    #    Hists[max_hist].SetMinimum(0)
    if not 'ZoomIn' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        low_end = -1.
        high_end = 3.
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    else:
        low_end = range_lower
        high_end = range_upper
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    Hists[max_hist].Draw("E")
    #Hists[0].SetMinimum(0)
    #Hists[0].Draw("E")

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
    predict_bkg, err_bkg = IntegralAndError(Hists[1],norm_bin_low_target,norm_bin_up_target)
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
    if 'Theta2' in name:
        print '+++++++++++++++++++++++++++++'
        print '%s, %0.1f hrs'%(target,exposure_hours)
        print 'Total = %0.1f'%(data_SR)
        print 'Bkg = %0.1f#pm%0.1f'%(predict_bkg,err_bkg)
        print 'Excess = %0.1f#pm%0.1f (%0.1f#sigma)'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5),Sig)
        print '$%0.0f$& $%0.0f$& $\pm%0.0f$& $%0.1f$'%(data_SR,predict_bkg,err_bkg,Sig)
        print '+++++++++++++++++++++++++++++'

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
    Hist_Band.GetYaxis().SetTitle("Data/Predict")
    Hist_Band.GetYaxis().SetTitleSize(0.10)
    Hist_Band.GetYaxis().SetNdivisions(505)
    #Hist_Band.SetMaximum(5)
    #Hist_Band.SetMinimum(-5)
    Hist_Band.SetMaximum(1.2)
    Hist_Band.SetMinimum(0.8)
    Hist_Band.GetXaxis().SetRangeUser(low_end,high_end)
    Hist_Band.Draw("e2")
    line2 = ROOT.TLine(Hist_Band.GetBinLowEdge(1),1,Hist_Band.GetBinLowEdge(Hist_Band.GetNbinsX()+1),1)
    line2.SetLineStyle(1)
    line2.SetLineColor(1)
    line2.SetLineWidth(2)
    line2.Draw("same")
    Hist_Ratio = []
    for h in range(0,len(Hists)):
        if colors[h]==1:
            Hist_Ratio += [Hists[h].Clone()]
            Hist_Ratio[h].Divide(Hist_Sum)
            Hist_Ratio[h].SetLineWidth(2)
            Hist_Ratio[h].Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_%s.pdf'%(name,target,tag))


def Make2DSkyCoordPlot(Hist_Data,xtitle,ytitle,name,doProj):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    #pad1.SetGrid()
    pad1.Draw()
    pad1.cd()
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.Draw("COL4Z")
    #Hist_Data.Draw("CONT3 same")
    #pad1.SetLogz()

    canvas.SaveAs('output_plots/%s.pdf'%(name))

def Make2DProjectionPlot(Hist_Data,xtitle,ytitle,name,doProj):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    #pad1.SetGrid()
    pad1.Draw()
    pad1.cd()
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.Draw("COL4Z")
    #Hist_Data.Draw("CONT3 same")
    bins = []
    for b in range(0,Hist_Data.GetNbinsX()+1):
        bins += [Hist_Data.GetXaxis().GetBinLowEdge(b+1)]
    Hist_1D = ROOT.TH1D("Hist_1D","",len(bins)-1,array('d',bins))
    for b in range(0,Hist_Data.GetNbinsX()):
        hist_temp = Hist_Data.ProjectionY("hist_temp",b+1,b+1)
        Hist_1D.SetBinContent(b+1,hist_temp.GetMean())
        Hist_1D.SetBinError(b+1,hist_temp.GetRMS())
    Hist_1D.SetLineColor(2)
    if doProj: Hist_1D.Draw("E same")
    #pad1.SetLogz()
    line1 = ROOT.TLine(MSCW_lower_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCW_upper_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCW_lower_cut,MSCW_lower_cut,MSCW_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()

    canvas.SaveAs('output_plots/%s.pdf'%(name))

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

    c_both.SaveAs('output_plots/%s_%s_%s.png'%(name,target,tag))

def Make2DSignificancePlot(Hist_SR,Hist_Bkg,xtitle,ytitle,name):

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

    syst = 0.02
    Hist_Data = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = pow(pow(Hist_Bkg.GetBinError(bx+1,by+1),2)+pow(syst*Hist_Bkg.GetBinContent(bx+1,by+1),2),0.5)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Data.SetBinContent(bx+1,by+1,Sig)

    Hist_Highlight = Hist_SR.Clone()
    Hist_Ratio = Hist_SR.Clone()
    Hist_Ratio.Add(Hist_Bkg,-1.)
    Hist_Ratio.Divide(Hist_Bkg)
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            Hist_Highlight.SetBinContent(bx+1,by+1,0.)
            if Hist_Data.GetBinContent(bx+1,by+1)>=5.0: 
                Hist_Highlight.SetBinContent(bx+1,by+1,5.0)
            if Hist_Data.GetBinContent(bx+1,by+1)<2.5: 
                Hist_Ratio.SetBinContent(bx+1,by+1,0.)

    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.SetMaximum(5)
    Hist_Data.SetMinimum(-5)
    Hist_Data.Draw("COL4Z")
    Hist_Highlight.Draw("CONT3 same")
    canvas.SaveAs('output_plots/SkymapSig_%s_%s_%s.png'%(name,target,tag))

    Hist_Ratio.GetYaxis().SetTitle(ytitle)
    Hist_Ratio.GetXaxis().SetTitle(xtitle)
    Hist_Ratio.SetMaximum(1.0)
    Hist_Ratio.SetMinimum(0)
    Hist_Ratio.Draw("COL4Z")
    Hist_Highlight.Draw("CONT3 same")
    canvas.SaveAs('output_plots/SkymapRatio_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    lumilab1 = ROOT.TLatex(0.15,0.90,'Integral = %0.1f'%(Hist_Excess.Integral()) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.05)
    lumilab1.Draw()
    canvas.SaveAs('output_plots/SkymapTotal_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_Bkg.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    lumilab2 = ROOT.TLatex(0.15,0.90,'Integral = %0.1f'%(Hist_Excess.Integral()) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.05)
    lumilab2.Draw()
    canvas.SaveAs('output_plots/SkymapBkg_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.Add(Hist_Bkg,-1.)
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    Hist_Highlight.Draw("CONT3 same")
    canvas.SaveAs('output_plots/SkymapExcess_%s_%s_%s.png'%(name,target,tag))

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    Hist_Sig_cut = ROOT.TH1D("Hist_Sig_cut","",65,-5,8)
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0 and not Hist_Bkg.GetBinContent(bx+1,by+1)==0:
                content = Hist_Data.GetBinContent(bx+1,by+1)
                Hist_Sig.Fill(content)
                locationx1 = Hist_Data.GetXaxis().GetBinCenter(bx)
                locationy1 = Hist_Data.GetYaxis().GetBinCenter(by)
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
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SkymapSigDist_%s'%(name))
    return Hist_Sig

def Make2DSignificancePlotShowerShape(Hist_SR,Hist_CR,Hist_gamma,xtitle,ytitle,theta2_low,theta2_up,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    Hist_Bkg = Hist_CR.Clone()
    Hist_Bkg.Add(Hist_gamma)
    if not theta2_low==0 and not theta2_up==0:
        for bx in range(0,Hist_SR.GetNbinsX()):
            for by in range(0,Hist_SR.GetNbinsY()):
                theta2 = pow(Hist_SR.GetXaxis().GetBinCenter(bx+1),2)+pow(Hist_SR.GetYaxis().GetBinCenter(by+1),2)
                if theta2<theta2_low or theta2>theta2_up: 
                    Hist_SR.SetBinContent(bx+1,by+1,0)
                    Hist_SR.SetBinError(bx+1,by+1,0)
                    Hist_Bkg.SetBinContent(bx+1,by+1,0)
                    Hist_Bkg.SetBinError(bx+1,by+1,0)
        low_end = -1.5*pow(theta2_up,0.5)
        high_end = 1.5*pow(theta2_up,0.5)
        Hist_SR.GetXaxis().SetRangeUser(low_end,high_end)
        Hist_SR.GetYaxis().SetRangeUser(low_end,high_end)
        Hist_Bkg.GetXaxis().SetRangeUser(low_end,high_end)
        Hist_Bkg.GetYaxis().SetRangeUser(low_end,high_end)

    Hist_Data = Hist_SR.Clone()
    Hist_Data_abs = Hist_SR.Clone()
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = Hist_Bkg.GetBinError(bx+1,by+1)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Sig = min(5,Sig)
            Sig = max(-5,Sig)
            Hist_Data.SetBinContent(bx+1,by+1,Sig)
            Hist_Data_abs.SetBinContent(bx+1,by+1,abs(Sig))
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.SetMaximum(5)
    Hist_Data.SetMinimum(-5)
    Hist_Data.Draw("COL4Z")
    line1 = ROOT.TLine(MSCW_lower_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCW_upper_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCW_lower_cut,MSCW_lower_cut,MSCW_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()
    canvas.SaveAs('output_plots/ShowerShapeSig_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    line1 = ROOT.TLine(MSCW_lower_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCW_upper_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCW_lower_cut,MSCW_lower_cut,MSCW_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()
    canvas.SaveAs('output_plots/ShowerShapeTotal_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_Bkg.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    line1 = ROOT.TLine(MSCW_lower_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCW_upper_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCW_lower_cut,MSCW_lower_cut,MSCW_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()
    canvas.SaveAs('output_plots/ShowerShapeBkg_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.Add(Hist_Bkg,-1.)
    Hist_Excess.Divide(Hist_SR)
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.SetMaximum(1.0)
    Hist_Excess.SetMinimum(-1.0)
    Hist_Excess.Draw("COL4Z")
    line1 = ROOT.TLine(MSCW_lower_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCW_upper_cut,MSCW_blind_cut,MSCW_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCW_lower_cut,MSCW_lower_cut,MSCW_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()
    canvas.SaveAs('output_plots/ShowerShapeExcess_%s_%s_%s.png'%(name,target,tag))

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0 and not Hist_Bkg.GetBinContent(bx+1,by+1)==0:
                content = Hist_Data.GetBinContent(bx+1,by+1)
                locationx1 = Hist_Data.GetXaxis().GetBinCenter(bx)
                locationy1 = Hist_Data.GetYaxis().GetBinCenter(by)
                if (locationx1<1 and locationy1<1):
                    Hist_Sig.Fill(content)
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
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','ShowerShapeSigDist_%s'%(name))
    return Hist_Sig

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

    c_both.SaveAs('output_plots/%s_%s.png'%(name,tag))

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

    c_both.SaveAs('output_plots/%s_%s.png'%(name,tag))


Hist_SumE_S2B = []
Hist_RDBM_S2B = []
Hist_Dark_S2B = []
Hist_Ring_S2B = []
Hist_MSCW = []
Hist_Eigenvalue_Rank = []
Hist_RDBM_S2B_GalLat = []
legend_S2B = []
color_S2B = []

Hist_Bkg_Rate = []
Hist_Sig_Rate = []
Hist_NormSig_Rate = []
Hist_S2B_Rate = []
legend2_S2B = []
color2_S2B = []

c = SkyCoord(sky_coord[:], unit=(u.hourangle, u.deg))
for s in range(0,len(source)):

    y = c.galactic.b.degree[s]
    if abs(y)<gal_lat_cut: continue

    target = source[s]
    target_label = target
    if not PercentCrab=="_Crab0":
        target_label += " + MC #gamma"

    Hist_Dark_ShowerDirection_Sum = ROOT.TH2D("Hist_Dark_ShowerDirection_Sum","",180,0,360,90,0,90)
    Hist_Data_ShowerDirection_Sum = ROOT.TH2D("Hist_Data_ShowerDirection_Sum","",180,0,360,90,0,90)
    Hist2D_Data_SumE = ROOT.TH2D("Hist2D_Data_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_Dark_SumE = ROOT.TH2D("Hist2D_Dark_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_Bkgd_SumE = ROOT.TH2D("Hist2D_Bkgd_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_GammaDark_SumE = ROOT.TH2D("Hist2D_GammaDark_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_GammaRDBM_SumE = ROOT.TH2D("Hist2D_GammaRDBM_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_TrueBkgd_SumE = ROOT.TH2D("Hist2D_TrueBkgd_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_GammaDark_MSCW_SumE = ROOT.TH1D("Hist_GammaDark_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_GammaRDBM_MSCW_SumE = ROOT.TH1D("Hist_GammaRDBM_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Data_MSCW_SumE = ROOT.TH1D("Hist_Data_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Ring_MSCW_SumE = ROOT.TH1D("Hist_Ring_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Bkgd_MSCW_SumE = ROOT.TH1D("Hist_Bkgd_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_TrueBkgd_MSCW_SumE = ROOT.TH1D("Hist_TrueBkgd_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Dark_MSCW_SumE = ROOT.TH1D("Hist_Dark_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_GammaDark_MSCL_SumE = ROOT.TH1D("Hist_GammaDark_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_GammaRDBM_MSCL_SumE = ROOT.TH1D("Hist_GammaRDBM_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Data_MSCL_SumE = ROOT.TH1D("Hist_Data_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Ring_MSCL_SumE = ROOT.TH1D("Hist_Ring_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Bkgd_MSCL_SumE = ROOT.TH1D("Hist_Bkgd_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_TrueBkgd_MSCL_SumE = ROOT.TH1D("Hist_TrueBkgd_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Dark_MSCL_SumE = ROOT.TH1D("Hist_Dark_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Data_Theta2 = ROOT.TH1D("Hist_Data_Theta2","",1024,0,10)
    Hist_TrueBkgd_Theta2 = ROOT.TH1D("Hist_TrueBkgd_Theta2","",1024,0,10)

    souce_ra = round(float(HMS2deg(sky_coord[s].split('+')[0],sky_coord[s].split('+')[1])[0]),3)
    souce_dec = round(float(HMS2deg(sky_coord[s].split('+')[0],sky_coord[s].split('+')[1])[1]),3)
    if source[s]=='Proton':
        source_ra = 0.
        source_dec = 0.

    Hist_Data_Skymap = ROOT.TH2D("Hist_Data_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Bkgd_Theta2 = ROOT.TH1D("Hist_Bkgd_Theta2","",1024,0,10)
    Hist_Bkgd_Theta2_Raw = ROOT.TH1D("Hist_Bkgd_Theta2_Raw","",1024,0,10)
    Hist_Bkgd_Skymap = ROOT.TH2D("Hist_Bkgd_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Bkgd_Skymap_Raw = ROOT.TH2D("Hist_Bkgd_Skymap_Raw","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Bkgd_Syst_Theta2 = ROOT.TH1D("Hist_Bkgd_Syst_Theta2","",1024,0,10)
    Hist_Bkgd_Syst_Skymap = ROOT.TH2D("Hist_Bkgd_Syst_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Dark_Theta2 = ROOT.TH1D("Hist_Dark_Theta2","",1024,0,10)
    Hist_Dark_Skymap = ROOT.TH2D("Hist_Dark_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Dark_Syst_Theta2 = ROOT.TH1D("Hist_Dark_Syst_Theta2","",1024,0,10)
    Hist_Dark_Syst_Skymap = ROOT.TH2D("Hist_Dark_Syst_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Ring_Theta2 = ROOT.TH1D("Hist_Ring_Theta2","",1024,0,10)
    Hist_Ring_Skymap = ROOT.TH2D("Hist_Ring_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    Hist_Ring_Syst_Theta2 = ROOT.TH1D("Hist_Ring_Syst_Theta2","",1024,0,10)
    Hist_Ring_Syst_Skymap = ROOT.TH2D("Hist_Ring_Syst_Skymap","",150,souce_ra-3,souce_ra+3,150,souce_dec-3,souce_dec+3)
    exposure_hours = 0.

    Hist_SumE_S2B += [ROOT.TH1D("Hist_SumE_S2B_%s"%(target),"",2,0,2)]
    Hist_RDBM_S2B_GalLat += [ROOT.TH1D("Hist_RDBM_S2B_GalLat_%s"%(target),"",6,0,90)]
    Hist_RDBM_S2B += [ROOT.TH1D("Hist_RDBM_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    Hist_Dark_S2B += [ROOT.TH1D("Hist_Dark_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    Hist_Ring_S2B += [ROOT.TH1D("Hist_Ring_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    Hist_Eigenvalue_Rank += [ROOT.TH1D("Hist_Eigenvalue_Rank_%s"%(target),"",N_bins_for_deconv,0,N_bins_for_deconv)]
    Hist_MSCW += [ROOT.TH1D("Hist_MSCW_%s"%(target),"",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)]

    Hist_Data_Theta2_SumElev = []
    Hist_TrueBkgd_Theta2_SumElev = []
    Hist_Bkgd_Theta2_SumElev = []
    Hist_Dark_Theta2_SumElev = []
    Hist_Bkgd_Syst_Theta2_SumElev = []
    Hist_Dark_Syst_Theta2_SumElev = []
    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]
        Hist_Data_Theta2_SumElev += [ROOT.TH1D("Hist_Data_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]
        Hist_TrueBkgd_Theta2_SumElev += [ROOT.TH1D("Hist_TrueBkgd_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]
        Hist_Bkgd_Theta2_SumElev += [ROOT.TH1D("Hist_Bkgd_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]
        Hist_Dark_Theta2_SumElev += [ROOT.TH1D("Hist_Dark_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]
        Hist_Bkgd_Syst_Theta2_SumElev += [ROOT.TH1D("Hist_Bkgd_Syst_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]
        Hist_Dark_Syst_Theta2_SumElev += [ROOT.TH1D("Hist_Dark_Syst_Theta2_SumElev_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut),"",1024,0,10)]

    for elev in range(0,len(elev_range)):

        file_elev_lower = elev_range[elev][0]
        file_elev_upper = elev_range[elev][1]

        Hist_Bkg_Rate += [ROOT.TH1D("Hist_Bkg_Rate_%s_%s"%(target,file_elev_lower),"",len(energy_list)-1,array('d',energy_list))]
        Hist_Sig_Rate += [ROOT.TH1D("Hist_Sig_Rate_%s_%s"%(target,file_elev_lower),"",len(energy_list)-1,array('d',energy_list))]
        Hist_NormSig_Rate += [ROOT.TH1D("Hist_NormSig_Rate_%s_%s"%(target,file_elev_lower),"",len(energy_list)-1,array('d',energy_list))]
        Hist_S2B_Rate += [ROOT.TH1D("Hist_S2B_Rate_%s_%s"%(target,file_elev_lower),"",len(energy_list)-1,array('d',energy_list))]
        legend2_S2B += ['%s (%0.1f hrs, elev %s-%s)'%(target,exposure_hours,file_elev_lower,file_elev_upper)]
        color2_S2B += [elev+1]

        FilePath = "%s/Netflix_"%(FileFolder)+target+PercentCrab+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_%s"%(file_tag)+".root";
        if not os.path.isfile(FilePath):continue
        InputFile = ROOT.TFile(FilePath)
        InfoTree = InputFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        exposure_hours += InfoTree.exposure_hours
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_blind_cut = InfoTree.MSCL_cut_blind

        HistName = "Hist_Dark_ShowerDirection"
        Hist_Dark_ShowerDirection = InputFile.Get(HistName)
        Hist_Dark_ShowerDirection_Sum.Add(Hist_Dark_ShowerDirection)
        HistName = "Hist_Data_ShowerDirection"
        Hist_Data_ShowerDirection = InputFile.Get(HistName)
        Hist_Data_ShowerDirection_Sum.Add(Hist_Data_ShowerDirection)

        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            HistName = "Hist_Data_Eigenvalues_real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_Eigenvalues_real = InputFile.Get(HistName)
            HistName = "Hist_Dark_Eigenvalues_real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_Eigenvalues_real = InputFile.Get(HistName)
            HistName = "Hist_Data_Eigenvalues_imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_Eigenvalues_imag = InputFile.Get(HistName)
            HistName = "Hist_Dark_Eigenvalues_imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_Eigenvalues_imag = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Data_EigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_EigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Data_InvEigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_InvEigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Dark_EigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_EigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Dark_InvEigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_InvEigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorReal_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorReal_0 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorReal_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorReal_1 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorReal_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorReal_2 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Fit_EigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_EigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorImag_0_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorImag_0 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorImag_1_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorImag_1 = InputFile.Get(HistName)
            HistName = "Hist_Fit_InvEigenvectorImag_2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Fit_InvEigenvectorImag_2 = InputFile.Get(HistName)
            HistName = "Hist_GammaDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_GammaDark = InputFile.Get(HistName)
            HistName = "Hist_GammaRDBM_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_GammaRDBM = InputFile.Get(HistName)
            HistName = "Hist_Data_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_Data = InputFile.Get(HistName)
            HistName = "Hist_Ring_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_Ring = InputFile.Get(HistName)
            HistName = "Hist_Dark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_Dark = InputFile.Get(HistName)
            HistName = "Hist_TrueBkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_TrueBkgd = InputFile.Get(HistName)
            HistName = "Hist_GammaDark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_GammaDark = InputFile.Get(HistName)
            HistName = "Hist_GammaRDBM_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_GammaRDBM = InputFile.Get(HistName)
            HistName = "Hist_Bkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist2D_Bkgd = InputFile.Get(HistName)
            Hist2D_Bkgd = AddStat2D(Hist2D_Dark,Hist2D_Bkgd)
            #Hist2D_Bkgd = AddSyst2D(Hist2D_Data,Hist2D_Dark,Hist2D_Bkgd)
            #Make2DProjectionPlot(Hist2D_Data,"MSCL","MSCW",'Data_2D_E%s'%(ErecS_lower_cut),False)
            #Make2DProjectionPlot(Hist2D_Bkgd,"MSCL","MSCW",'Bkgd_2D_E%s'%(ErecS_lower_cut),False)
            #Make2DProjectionPlot(Hist2D_Ring,"MSCL","MSCW",'Ring_2D_E%s'%(ErecS_lower_cut),False)
            #Make2DProjectionPlot(Hist2D_Dark,"MSCL","MSCW",'Dark_2D_E%s'%(ErecS_lower_cut),False)
            print HistName
            print 'Hist2D_Data.Integral() = %s'%(Hist2D_Data.Integral())
            print 'Hist2D_Bkgd.Integral() = %s'%(Hist2D_Bkgd.Integral())
            bin_lower = Hist2D_Data.GetXaxis().FindBin(MSCL_lower_cut)
            bin_upper = Hist2D_Data.GetXaxis().FindBin(MSCL_blind_cut)-1
            Hist1D_GammaDark_MSCW_SumSRs = Hist2D_GammaDark.ProjectionY("Hist1D_GammaDark_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_GammaRDBM_MSCW_SumSRs = Hist2D_GammaRDBM.ProjectionY("Hist1D_GammaRDBM_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_Data_MSCW_SumSRs = Hist2D_Data.ProjectionY("Hist1D_Data_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_Ring_MSCW_SumSRs = Hist2D_Ring.ProjectionY("Hist1D_Ring_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_Bkgd_MSCW_SumSRs = Hist2D_Bkgd.ProjectionY("Hist1D_Bkgd_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_TrueBkgd_MSCW_SumSRs = Hist2D_TrueBkgd.ProjectionY("Hist1D_TrueBkgd_MSCW_SumSRs",bin_lower,bin_upper)
            Hist1D_Dark_MSCW_SumSRs = Hist2D_Dark.ProjectionY("Hist1D_Dark_MSCW_SumSRs",bin_lower,bin_upper)
            bin_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
            bin_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1
            Hist1D_GammaDark_MSCL_SumSRs = Hist2D_GammaDark.ProjectionX("Hist1D_GammaDark_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_GammaRDBM_MSCL_SumSRs = Hist2D_GammaRDBM.ProjectionX("Hist1D_GammaRDBM_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_Data_MSCL_SumSRs = Hist2D_Data.ProjectionX("Hist1D_Data_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_Ring_MSCL_SumSRs = Hist2D_Ring.ProjectionX("Hist1D_Ring_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_Bkgd_MSCL_SumSRs = Hist2D_Bkgd.ProjectionX("Hist1D_Bkgd_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_TrueBkgd_MSCL_SumSRs = Hist2D_TrueBkgd.ProjectionX("Hist1D_TrueBkgd_MSCL_SumSRs",bin_lower,bin_upper)
            Hist1D_Dark_MSCL_SumSRs = Hist2D_Dark.ProjectionX("Hist1D_Dark_MSCL_SumSRs",bin_lower,bin_upper)


            Hist2D_Data_SumE.Add(Hist2D_Data)
            Hist2D_Dark_SumE.Add(Hist2D_Dark)
            Hist2D_Bkgd_SumE.Add(Hist2D_Bkgd)
            Hist2D_GammaDark_SumE.Add(Hist2D_GammaDark)
            Hist2D_GammaRDBM_SumE.Add(Hist2D_GammaRDBM)
            Hist2D_TrueBkgd_SumE.Add(Hist2D_TrueBkgd)
            Hist_GammaDark_MSCW_SumE.Add(Hist1D_GammaDark_MSCW_SumSRs)
            Hist_GammaRDBM_MSCW_SumE.Add(Hist1D_GammaRDBM_MSCW_SumSRs)
            Hist_Data_MSCW_SumE.Add(Hist1D_Data_MSCW_SumSRs)
            Hist_Ring_MSCW_SumE.Add(Hist1D_Ring_MSCW_SumSRs)
            Hist_Bkgd_MSCW_SumE.Add(Hist1D_Bkgd_MSCW_SumSRs)
            Hist_TrueBkgd_MSCW_SumE.Add(Hist1D_TrueBkgd_MSCW_SumSRs)
            Hist_Dark_MSCW_SumE.Add(Hist1D_Dark_MSCW_SumSRs)
            Hist_GammaDark_MSCL_SumE.Add(Hist1D_GammaDark_MSCL_SumSRs)
            Hist_GammaRDBM_MSCL_SumE.Add(Hist1D_GammaRDBM_MSCL_SumSRs)
            Hist_Data_MSCL_SumE.Add(Hist1D_Data_MSCL_SumSRs)
            Hist_Ring_MSCL_SumE.Add(Hist1D_Ring_MSCL_SumSRs)
            Hist_Bkgd_MSCL_SumE.Add(Hist1D_Bkgd_MSCL_SumSRs)
            Hist_TrueBkgd_MSCL_SumE.Add(Hist1D_TrueBkgd_MSCL_SumSRs)
            Hist_Dark_MSCL_SumE.Add(Hist1D_Dark_MSCL_SumSRs)

            bin_lower = Hist1D_Bkgd_MSCW_SumSRs.FindBin(MSCW_lower_cut)
            bin_upper = Hist1D_Bkgd_MSCW_SumSRs.FindBin(MSCW_blind_cut)-1

            HistName = "Hist_TrueBkgd_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_TrueBkgd_Theta2_Tmp = InputFile.Get(HistName)
            HistName = "Hist_Data_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_Theta2_Tmp = InputFile.Get(HistName)
            Hist_Data_Theta2.Add(Hist_Data_Theta2_Tmp)
            Hist_Data_Theta2_SumElev[e].Add(Hist_Data_Theta2_Tmp)
            HistName = "Hist_Data_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Data_Skymap.Add(InputFile.Get(HistName))

            bkg_total, bkg_err = IntegralAndError(Hist1D_Bkgd_MSCW_SumSRs,bin_lower,bin_upper)
            truebkg_total, truebkg_err = IntegralAndError(Hist1D_TrueBkgd_MSCW_SumSRs,bin_lower,bin_upper)
            HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Bkgd_Theta2_Tmp = InputFile.Get(HistName)
            HistName = "Hist_Data_CR_SelectFoV_Theta2_Raw_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Bkgd_Theta2_Raw_Tmp = InputFile.Get(HistName)
            Hist_Bkgd_Syst_Theta2_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Bkgd_Theta2_Tmp.Integral()
            old_integral_raw = Hist_Bkgd_Theta2_Raw_Tmp.Integral()
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
            print 'bkgd_scale = %s, bkgd_scale_err = %s'%(bkgd_scale,bkgd_scale_err)
            if not bkg_total==0:
                Theta2HistScale(Hist_Bkgd_Theta2_Tmp,Hist_Bkgd_Syst_Theta2_Tmp,bkgd_scale,bkgd_scale_err)
                Theta2HistScale(Hist_Bkgd_Theta2_Raw_Tmp,Hist_Bkgd_Syst_Theta2_Tmp,bkgd_scale_raw,bkgd_scale_raw_err)
            else:
                Hist_Bkgd_Theta2_Tmp.Scale(0)
                Hist_Bkgd_Theta2_Raw_Tmp.Scale(0)
            Hist_Bkgd_Theta2.Add(Hist_Bkgd_Theta2_Tmp)
            Hist_Bkgd_Theta2_Raw.Add(Hist_Bkgd_Theta2_Raw_Tmp)
            Hist_Bkgd_Theta2_SumElev[e].Add(Hist_Bkgd_Theta2_Tmp)
            Hist_TrueBkgd_Theta2.Add(Hist_TrueBkgd_Theta2_Tmp)
            Hist_TrueBkgd_Theta2_SumElev[e].Add(Hist_TrueBkgd_Theta2_Tmp)
            Hist_Bkgd_Syst_Theta2.Add(Hist_Bkgd_Syst_Theta2_Tmp)
            Hist_Bkgd_Syst_Theta2_SumElev[e].Add(Hist_Bkgd_Syst_Theta2_Tmp)

            data_rate, data_rate_err = Event_rate(Hist2D_Data,exposure_hours)
            bkg_rate, bkg_rate_err = Event_rate(Hist2D_Bkgd,exposure_hours)
            sig_rate = data_rate-bkg_rate
            sig_rate_err = pow(data_rate_err*data_rate_err+bkg_rate_err*bkg_rate_err,0.5)
            Hist_Bkg_Rate[len(Hist_Bkg_Rate)-1].SetBinContent(e+1,bkg_rate)
            Hist_Bkg_Rate[len(Hist_Bkg_Rate)-1].SetBinError(e+1,bkg_rate_err)
            Hist_Sig_Rate[len(Hist_Sig_Rate)-1].SetBinContent(e+1,sig_rate)
            Hist_Sig_Rate[len(Hist_Sig_Rate)-1].SetBinError(e+1,sig_rate_err)
            Hist_S2B_Rate[len(Hist_Sig_Rate)-1].SetBinContent(e+1,sig_rate/bkg_rate)
            ratio_err = pow(sig_rate_err/sig_rate,2)+pow(bkg_rate_err/bkg_rate,2)
            ratio_err = sig_rate/bkg_rate*pow(ratio_err,0.5)
            Hist_S2B_Rate[len(Hist_Sig_Rate)-1].SetBinError(e+1,ratio_err)
            norm_sig_rate = sig_rate/bkg_rate*Hist_Bkg_Rate[0].GetBinContent(e+1)
            Hist_NormSig_Rate[len(Hist_NormSig_Rate)-1].SetBinContent(e+1,norm_sig_rate)
            Hist_NormSig_Rate[len(Hist_NormSig_Rate)-1].SetBinError(e+1,sig_rate_err)

            HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Bkgd_Skymap_Tmp = InputFile.Get(HistName)
            HistName = "Hist_Data_CR_Skymap_Raw_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Bkgd_Skymap_Raw_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Bkgd_Skymap_Tmp.Integral()
            old_integral_raw = Hist_Bkgd_Skymap_Raw_Tmp.Integral()
            scale = 0
            scale_err = 0
            scale_raw = 0
            scale_raw_err = 0
            Hist_Bkgd_Syst_Skymap_Tmp = InputFile.Get(HistName)
            Hist_Bkgd_Syst_Skymap_Tmp.Reset()
            if not bkg_total==0 and not old_integral==0:
                scale = bkg_total/old_integral
                scale_err = scale*(bkg_err/bkg_total)
                scale_raw = bkg_total/old_integral_raw
                scale_raw_err = scale_raw*(bkg_err/bkg_total)
            else:
                scale = 0
                scale_err = 0
                scale_raw = 0
                scale_raw_err = 0
            if not bkg_total==0:
                RaDecHistScale(Hist_Bkgd_Skymap_Tmp,Hist_Bkgd_Syst_Skymap_Tmp,scale,scale_err)
                RaDecHistScale(Hist_Bkgd_Skymap_Raw_Tmp,Hist_Bkgd_Syst_Skymap_Tmp,scale_raw,scale_raw_err)
            else:
                Hist_Bkgd_Skymap_Tmp.Scale(0)
                Hist_Bkgd_Skymap_Raw_Tmp.Scale(0)
            Hist_Bkgd_Skymap.Add(Hist_Bkgd_Skymap_Tmp)
            Hist_Bkgd_Skymap_Raw.Add(Hist_Bkgd_Skymap_Raw_Tmp)
            Hist_Bkgd_Syst_Skymap.Add(Hist_Bkgd_Syst_Skymap_Tmp)

            bkg_total, bkg_err = IntegralAndError(Hist1D_Dark_MSCW_SumSRs,bin_lower,bin_upper)
            HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_Theta2_Tmp = InputFile.Get(HistName)
            Hist_Dark_Syst_Theta2_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Dark_Theta2_Tmp.Integral()
            dark_scale = 0
            dark_scale_err = 0
            if not bkg_total==0 and not old_integral==0:
                dark_scale = bkg_total/old_integral
                dark_scale_err = dark_scale*(bkg_err/bkg_total)
            else:
                dark_scale = 0
                dark_scale_err = 0
            print 'dark_scale = %s, dark_scale_err = %s'%(dark_scale,dark_scale_err)
            if not bkg_total==0:
                Theta2HistScale(Hist_Dark_Theta2_Tmp,Hist_Dark_Syst_Theta2_Tmp,dark_scale,dark_scale_err)
            else:
                Hist_Dark_Theta2_Tmp.Scale(0)
            Hist_Dark_Theta2.Add(Hist_Dark_Theta2_Tmp)
            Hist_Dark_Theta2_SumElev[e].Add(Hist_Dark_Theta2_Tmp)
            Hist_Dark_Syst_Theta2.Add(Hist_Dark_Syst_Theta2_Tmp)
            Hist_Dark_Syst_Theta2_SumElev[e].Add(Hist_Dark_Syst_Theta2_Tmp)

            HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Dark_Skymap_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Dark_Skymap_Tmp.Integral()
            scale = 0
            scale_err = 0
            Hist_Dark_Syst_Skymap_Tmp = InputFile.Get(HistName)
            Hist_Dark_Syst_Skymap_Tmp.Reset()
            if not bkg_total==0 and not old_integral==0:
                scale = bkg_total/old_integral
                scale_err = scale*(bkg_err/bkg_total)
            else:
                scale = 0
                scale_err = 0
            if not bkg_total==0:
                RaDecHistScale(Hist_Dark_Skymap_Tmp,Hist_Dark_Syst_Skymap_Tmp,scale,scale_err)
            else:
                Hist_Dark_Skymap_Tmp.Scale(0)
            Hist_Dark_Skymap.Add(Hist_Dark_Skymap_Tmp)
            Hist_Dark_Syst_Skymap.Add(Hist_Dark_Syst_Skymap_Tmp)

            bkg_total, bkg_err = IntegralAndError(Hist1D_Ring_MSCW_SumSRs,bin_lower,bin_upper)
            HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Ring_Theta2_Tmp = InputFile.Get(HistName)
            Hist_Ring_Syst_Theta2_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Ring_Theta2_Tmp.Integral()
            Ring_scale = 0
            Ring_scale_err = 0
            if not bkg_total==0 and not old_integral==0:
                Ring_scale = bkg_total/old_integral
                Ring_scale_err = Ring_scale*(bkg_err/bkg_total)
            else:
                Ring_scale = 0
                Ring_scale_err = 0
            print 'Ring_scale = %s, Ring_scale_err = %s'%(Ring_scale,Ring_scale_err)
            if not bkg_total==0:
                Theta2HistScale(Hist_Ring_Theta2_Tmp,Hist_Ring_Syst_Theta2_Tmp,Ring_scale,Ring_scale_err)
            else:
                Hist_Ring_Theta2_Tmp.Scale(0)
            Hist_Ring_Theta2.Add(Hist_Ring_Theta2_Tmp)
            Hist_Ring_Syst_Theta2.Add(Hist_Ring_Syst_Theta2_Tmp)

            Ring_syst = 0
            if Ring_scale>0: Ring_syst = Ring_scale_err/Ring_scale
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_Tmp,Hist_Ring_Theta2_Tmp,0.2,theta2_upper,Ring_syst)
            Hist_Ring_S2B[len(Hist_Ring_S2B)-1].SetBinContent(e+1,s2b)
            Hist_Ring_S2B[len(Hist_Ring_S2B)-1].SetBinError(e+1,s2b_err)

            HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
            Hist_Ring_Skymap_Tmp = InputFile.Get(HistName)
            old_integral = Hist_Ring_Skymap_Tmp.Integral()
            scale = 0
            scale_err = 0
            Hist_Ring_Syst_Skymap_Tmp = InputFile.Get(HistName)
            Hist_Ring_Syst_Skymap_Tmp.Reset()
            if not bkg_total==0 and not old_integral==0:
                scale = bkg_total/old_integral
                scale_err = scale*(bkg_err/bkg_total)
            else:
                scale = 0
                scale_err = 0
            if not bkg_total==0:
                RaDecHistScale(Hist_Ring_Skymap_Tmp,Hist_Ring_Syst_Skymap_Tmp,scale,scale_err)
            else:
                Hist_Ring_Skymap_Tmp.Scale(0)
            Hist_Ring_Skymap.Add(Hist_Ring_Skymap_Tmp)
            Hist_Ring_Syst_Skymap.Add(Hist_Ring_Syst_Skymap_Tmp)

            for binx in range(1,Hist_Data_Eigenvalues_real.GetNbinsX()+1):
                new_content = pow(Hist_Data_Eigenvalues_real.GetBinContent(binx),2)
                new_content += pow(Hist_Data_Eigenvalues_imag.GetBinContent(binx),2)
                new_content = pow(new_content,0.5)
                Hist_Data_Eigenvalues_real.SetBinContent(binx,new_content)
            for binx in range(1,Hist_Dark_Eigenvalues_real.GetNbinsX()+1):
                new_content = pow(Hist_Dark_Eigenvalues_real.GetBinContent(binx),2)
                new_content += pow(Hist_Dark_Eigenvalues_imag.GetBinContent(binx),2)
                new_content = pow(new_content,0.5)
                Hist_Dark_Eigenvalues_real.SetBinContent(binx,new_content)
            Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Add(Hist_Data_Eigenvalues_real)
            Hist_MSCW[len(Hist_MSCW)-1].Add(Hist1D_Data_MSCW_SumSRs)
            if PlotAllEnergies:
                if UseRDBM:
                    if Hist_Data_Eigenvalues_real.GetBinContent(1)>0:
                        scale = 1./Hist_Data_Eigenvalues_real.GetBinContent(1)
                        Hist_Data_Eigenvalues_real.Scale(scale)
                        Hist_Data_Eigenvalues_imag.Scale(scale)
                    if Hist_Dark_Eigenvalues_real.GetBinContent(1)>0:
                        scale = 1./Hist_Dark_Eigenvalues_real.GetBinContent(1)
                        Hist_Dark_Eigenvalues_real.Scale(scale)
                        Hist_Dark_Eigenvalues_imag.Scale(scale)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist_Data_Eigenvalues_real]
                    legends += ['proton+#gamma']
                    colors += [1]
                    Hists += [Hist_Dark_Eigenvalues_real]
                    legends += ['proton only']
                    colors += [2]
                    plotname = 'Target_Eigenvalues_E%s_%s'%(ErecS_lower_cut,target)
                    title_x = 'rank'
                    title_y = 'size of eigenvalue'
                    MakeComparisonPlot(Hists,legends,colors,title_x,title_y,plotname,False,False)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist1D_Data_MSCW_SumSRs]
                    legends += ['%s'%(target_label)]
                    colors += [1]
                    Hists += [Hist1D_TrueBkgd_MSCW_SumSRs]
                    legends += ['true bkg.']
                    colors += [4]
                    plotname = 'Target_W_True_E%s'%(ErecS_lower_cut)
                    title = 'MSCW'
                    MakeChi2Plot(Hists,legends,colors,title,plotname,True,True,MSCW_lower_cut,MSCW_blind_cut,-1)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist1D_Data_MSCW_SumSRs]
                    legends += ['%s'%(target_label)]
                    colors += [1]
                    Hists += [Hist1D_Dark_MSCW_SumSRs]
                    legends += ['predict. bkg.']
                    colors += [4]
                    plotname = 'Target_W_Dark_E%s'%(ErecS_lower_cut)
                    title = 'MSCW'
                    MakeChi2Plot(Hists,legends,colors,title,plotname,True,True,MSCW_lower_cut,MSCW_blind_cut,-1)

    ErecS_lower_cut = energy_list[0]

    if UseRDBM:
        bkgd_syst = Hist_Bkgd_Syst_Theta2.Integral()/Hist_Bkgd_Theta2.Integral()
        print 'bkgd_syst = %s'%(bkgd_syst)
        s2b = 0
        s2b_err = 0
        if not PercentCrab=="_Crab0":
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2,Hist_Bkgd_Theta2,theta2_lower,theta2_upper,bkgd_syst)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2,Hist_Bkgd_Theta2,0.2,theta2_upper,bkgd_syst)
        #s2b, s2b_err = Variation_ratio(Hist_Data_Theta2,Hist_Bkgd_Theta2,theta2_lower,theta2_upper,bkgd_syst)
        Hist_SumE_S2B[len(Hist_SumE_S2B)-1].SetBinContent(2,s2b)
        Hist_SumE_S2B[len(Hist_SumE_S2B)-1].SetBinError(2,s2b_err)
        Hist_RDBM_S2B_GalLat[len(Hist_RDBM_S2B_GalLat)-1].SetBinContent(Hist_RDBM_S2B_GalLat[0].FindBin(abs(y)),s2b)
        Hist_RDBM_S2B_GalLat[len(Hist_RDBM_S2B_GalLat)-1].SetBinError(Hist_RDBM_S2B_GalLat[0].FindBin(abs(y)),s2b_err)
    if UseDark:
        dark_syst = Hist_Dark_Syst_Theta2.Integral()/Hist_Dark_Theta2.Integral()
        print 'dark_syst = %s'%(bkgd_syst)
        s2b = 0
        s2b_err = 0
        if not PercentCrab=="_Crab0":
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2,Hist_Dark_Theta2,theta2_lower,theta2_upper,dark_syst)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2,Hist_Dark_Theta2,0.2,theta2_upper,dark_syst)
        #s2b, s2b_err = Variation_ratio(Hist_Data_Theta2,Hist_Dark_Theta2,theta2_lower,theta2_upper,dark_syst)
        Hist_SumE_S2B[len(Hist_SumE_S2B)-1].SetBinContent(1,s2b)
        Hist_SumE_S2B[len(Hist_SumE_S2B)-1].SetBinError(1,s2b_err)

    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCW_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Bkgd_MSCW_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        Hists += [Hist_GammaRDBM_MSCW_SumE]
        legends += ['pred. #gamma']
        colors += [3]
        plotname = 'Target_W_SumSRs_SumE_RDBM'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCL_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Bkgd_MSCL_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        Hists += [Hist_GammaRDBM_MSCL_SumE]
        legends += ['pred. #gamma']
        colors += [3]
        plotname = 'Target_L_SumSRs_SumE_RDBM'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
        Hist_Sig = Make2DSignificancePlotShowerShape(Hist2D_Data_SumE,Hist2D_Bkgd_SumE,Hist2D_GammaRDBM_SumE,'MSCL','MSCW',0,0,'RDBM_Sig2D')
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCW_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Ring_MSCW_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        plotname = 'Target_W_SumSRs_SumE_Ring'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCL_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Ring_MSCL_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        plotname = 'Target_L_SumSRs_SumE_Ring'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCW_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Dark_MSCW_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        Hists += [Hist_GammaDark_MSCW_SumE]
        legends += ['pred. #gamma']
        colors += [3]
        plotname = 'Target_W_SumSRs_SumE_Dark'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_MSCL_SumE]
            legends += ['true bkg.']
            colors += [2]
        Hists += [Hist_Dark_MSCL_SumE]
        legends += ['pred. bkg.']
        colors += [4]
        Hists += [Hist_GammaDark_MSCL_SumE]
        legends += ['pred. #gamma']
        colors += [3]
        plotname = 'Target_L_SumSRs_SumE_Dark'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
        Hist_Sig = Make2DSignificancePlotShowerShape(Hist2D_Data_SumE,Hist2D_Dark_SumE,Hist2D_GammaDark_SumE,'MSCL','MSCW',0,0,'Dark_Sig2D')

    n_rebin = 1
    if theta2_upper>0.25: n_rebin = 4
    if theta2_upper>1.00: n_rebin = 16
    Hist_TrueBkgd_Theta2.Rebin(n_rebin)
    Hist_Bkgd_Theta2.Rebin(n_rebin)
    Hist_Bkgd_Theta2_Raw.Rebin(n_rebin)
    Hist_Bkgd_Syst_Theta2.Rebin(n_rebin)
    Hist_Dark_Theta2.Rebin(n_rebin)
    Hist_Dark_Syst_Theta2.Rebin(n_rebin)
    Hist_Ring_Theta2.Rebin(n_rebin)
    Hist_Ring_Syst_Theta2.Rebin(n_rebin)
    Hist_Data_Theta2.Rebin(n_rebin)
    
    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_Theta2]
            legends += ['true bkg']
            colors += [2]
        Hists += [Hist_Bkgd_Theta2]
        legends += ['bkg']
        colors += [4]
        Hists += [Hist_Bkgd_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_RDBM_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_Theta2]
            legends += ['true bkg']
            colors += [2]
        Hists += [Hist_Bkgd_Theta2_Raw]
        legends += ['bkg']
        colors += [4]
        Hists += [Hist_Bkgd_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_RDBM_ZoomInTheta2_Raw_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,-1)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target_label)]
        colors += [1]
        if not PercentCrab=="_Crab0":
            Hists += [Hist_TrueBkgd_Theta2]
            legends += ['true bkg']
            colors += [2]
        Hists += [Hist_Dark_Theta2]
        legends += ['dark']
        colors += [4]
        Hists += [Hist_Dark_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_Dark_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,-1)
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target_label)]
        colors += [1]
        Hists += [Hist_Ring_Theta2]
        legends += ['ring']
        colors += [4]
        Hists += [Hist_Ring_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_Ring_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,-1)
    
    if DoSmooth:
        smooth_size = 0.07
        if UseRDBM:
            Hist_Data_Skymap_smooth = Smooth2DMap(Hist_Data_Skymap,smooth_size,False)
            Hist_Bkgd_Skymap_smooth = Smooth2DMap(Hist_Bkgd_Skymap,smooth_size,False)
            Hist_Bkgd_Skymap_Raw_smooth = Smooth2DMap(Hist_Bkgd_Skymap_Raw,smooth_size,False)
            plotname = 'Target_SRall_RDBM_Skymap_Smooth_E%s'%(ErecS_lower_cut)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap_smooth,Hist_Bkgd_Skymap_smooth,'RA','Dec',plotname)
            plotname = 'Target_SRall_RDBM_Skymap_Raw_Smooth_E%s'%(ErecS_lower_cut)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap_smooth,Hist_Bkgd_Skymap_Raw_smooth,'RA','Dec',plotname)
        if UseDark:
            Hist_Data_Skymap_smooth = Smooth2DMap(Hist_Data_Skymap,smooth_size,False)
            Hist_Dark_Skymap_smooth = Smooth2DMap(Hist_Dark_Skymap,smooth_size,False)
            plotname = 'Target_SRall_Dark_Skymap_Smooth_E%s'%(ErecS_lower_cut)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap_smooth,Hist_Dark_Skymap_smooth,'RA','Dec',plotname)
    if DoSkymap:
        smooth_size = 0.07
        n_rebin = 4
        Hist_Data_Skymap.Rebin2D(n_rebin,n_rebin)
        if UseRDBM:
            Hist_Bkgd_Skymap.Rebin2D(n_rebin,n_rebin)
            Hist_Bkgd_Skymap_Raw.Rebin2D(n_rebin,n_rebin)
            Hist_Bkgd_Syst_Skymap.Rebin2D(n_rebin,n_rebin)
            plotname = 'Target_SRall_RDBM_Skymap_E%s'%(ErecS_lower_cut)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Bkgd_Skymap,'RA','Dec',plotname)
            plotname = 'Target_SRall_RDBM_Skymap_Raw_E%s'%(ErecS_lower_cut)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Bkgd_Skymap_Raw,'RA','Dec',plotname)
        if UseDark:
            plotname = 'Target_SRall_Dark_Skymap_E%s'%(ErecS_lower_cut)
            Hist_Dark_Skymap.Rebin2D(n_rebin,n_rebin)
            Hist_Dark_Syst_Skymap.Rebin2D(n_rebin,n_rebin)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Dark_Skymap,'RA','Dec',plotname)
        if UseRing:
            plotname = 'Target_SRall_Ring_Skymap_E%s'%(ErecS_lower_cut)
            Hist_Ring_Skymap.Rebin2D(n_rebin,n_rebin)
            Hist_Ring_Syst_Skymap.Rebin2D(n_rebin,n_rebin)
            Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Ring_Skymap,'RA','Dec',theta2_lower,theta2_upper,plotname)

    for e in range(0,len(energy_list)-1):
        bkgd_syst = 0.
        if Hist_Bkgd_Theta2_SumElev[e].Integral()!=0.: bkgd_syst = Hist_Bkgd_Syst_Theta2_SumElev[e].Integral()/Hist_Bkgd_Theta2_SumElev[e].Integral()
        print 'bkgd_syst = %s'%(bkgd_syst)
        s2b = 0.
        s2b_err = 0.
        if not PercentCrab=="_Crab0":
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumElev[e],Hist_Bkgd_Theta2_SumElev[e],theta2_lower,theta2_upper,bkgd_syst)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumElev[e],Hist_Bkgd_Theta2_SumElev[e],0.2,theta2_upper,bkgd_syst)
        #s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumElev[e],Hist_Bkgd_Theta2_SumElev[e],theta2_lower,theta2_upper,bkgd_syst)
        Hist_RDBM_S2B[len(Hist_RDBM_S2B)-1].SetBinContent(e+1,s2b)
        Hist_RDBM_S2B[len(Hist_RDBM_S2B)-1].SetBinError(e+1,s2b_err)
        dark_syst = 0.
        if Hist_Dark_Theta2_SumElev[e].Integral()!=0.: dark_syst = Hist_Dark_Syst_Theta2_SumElev[e].Integral()/Hist_Dark_Theta2_SumElev[e].Integral()
        if not PercentCrab=="_Crab0":
            s2b, s2b_err = Variation_ratio(Hist_TrueBkgd_Theta2_SumElev[e],Hist_Dark_Theta2_SumElev[e],theta2_lower,theta2_upper,dark_syst)
        else:
            s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumElev[e],Hist_Dark_Theta2_SumElev[e],0.2,theta2_upper,dark_syst)
        #s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_SumElev[e],Hist_Dark_Theta2_SumElev[e],theta2_lower,theta2_upper,dark_syst)
        Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinContent(e+1,s2b)
        Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinError(e+1,s2b_err)

    legend_S2B += ['%s (%0.1f hrs)'%(target,exposure_hours)]
    color_S2B += [49-s]

    InputFile.Close()
    Make2DProjectionPlot(Hist_Dark_ShowerDirection_Sum,'Azimuth','Zenith','Dark_ShowerDirection',False)
    Make2DProjectionPlot(Hist_Data_ShowerDirection_Sum,'Azimuth','Zenith','Data_ShowerDirection',False)
    if Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)>0:
        scale = 1./Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].GetBinContent(1)
        Hist_Eigenvalue_Rank[len(Hist_Eigenvalue_Rank)-1].Scale(scale)
    if Hist_MSCW[len(Hist_MSCW)-1].Integral()>0:
        bin_low = Hist_MSCW[len(Hist_MSCW)-1].FindBin(1.)
        bin_up = Hist_MSCW[len(Hist_MSCW)-1].FindBin(3.)-1
        scale = 1./Hist_MSCW[len(Hist_MSCW)-1].Integral(bin_low,bin_up)
        Hist_MSCW[len(Hist_MSCW)-1].Scale(scale)

count = 0
for s in range(0,len(source)):
    y = c.galactic.b.degree[s]
    if abs(y)<gal_lat_cut: continue
    Hist_SumE_S2B[count].GetXaxis().SetBinLabel(1,'Initial')
    Hist_SumE_S2B[count].GetXaxis().SetBinLabel(2,'MDM')
    count += 1
if UseRDBM:
    MakeComparisonPlotSigDist(Hist_RDBM_S2B,legend_S2B,color_S2B,'E [GeV]','Mismodeling_RDBM_systematics',0.6,1.4,True,False)
    MakeComparisonPlotSigDist(Hist_RDBM_S2B_GalLat,legend_S2B,color_S2B,'gal. lat.','Mismodeling_RDBM_systematics_GalLat',0.6,1.4,False,False)
    MakeComparisonPlotTwoColumn(Hist_SumE_S2B,legend_S2B,color_S2B,'methods','Mismodeling_systematics',0.6,1.4,False,False)
if UseDark:
    MakeComparisonPlotSigDist(Hist_Dark_S2B,legend_S2B,color_S2B,'E [GeV]','Mismodeling_Dark_systematics',0.6,1.4,True,False)
MakeComparisonPlot(Hist_Bkg_Rate,legend2_S2B,color2_S2B,'E [GeV]','event rate','BackgroundRate',True,True)
MakeComparisonPlot(Hist_Sig_Rate,legend2_S2B,color2_S2B,'E [GeV]','event rate','SignalRate',True,True)
MakeComparisonPlot(Hist_NormSig_Rate,legend2_S2B,color2_S2B,'E [GeV]','event rate','NormSignalRate',True,True)
MakeComparisonPlot(Hist_S2B_Rate,legend2_S2B,color2_S2B,'E [GeV]','event rate','Sig2BkgRate',True,True)

for hist in range(0,len(Hist_Eigenvalue_Rank)):
    Hist_Eigenvalue_Rank[hist].GetXaxis().SetRangeUser(0,16)
MakeComparisonPlot(Hist_Eigenvalue_Rank,legend_S2B,color_S2B,'rank','size of eigenvalues','Eigenvalue_Rank',False,False)
MakeComparisonPlot(Hist_MSCW,legend_S2B,color_S2B,'MSCW','normalized counts','MultiSourceMSCW',False,False)

#c = SkyCoord(frame=ICRS, ra=308.119, dec=40.3282, unit=(u.degree, u.degree))
#print c.to_string('hmsdms')
Hists = []
legends = []
colors = []
Hist_SkyCoord = ROOT.TH2D("Hist_SkyCoord","",180,-180,180,90,-90,90)
Prof_GalExcess = ROOT.TProfile("Prof_GalExcess","",6,0,90,-100,100)
count = 0
for s in range(0,len(source)):
    y = c.galactic.b.degree[s]
    if abs(y)<gal_lat_cut: continue
    print c.galactic.l.degree[s]
    print c.galactic.b.degree[s]
    x = c.galactic.l.degree[s]
    if x>180.: x = x-360.
    y = c.galactic.b.degree[s]
    Hist_SkyCoord.Fill(x,y,100.*(Hist_SumE_S2B[count].GetBinContent(2)-1.))
    Prof_GalExcess.Fill(pow(y*y,0.5),100.*(Hist_SumE_S2B[count].GetBinContent(2)-1.))
    count += 1
Make2DSkyCoordPlot(Hist_SkyCoord,'x','y','SkyCoord',False)
Hists += [Prof_GalExcess]
legends += ["VERITAS extragalactic sources"]
colors += [2]
MakeComparisonPlot(Hists,legends,colors,"Gal. lat.","variance","GalExcess",False,False)
