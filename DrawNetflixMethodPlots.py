
import sys,ROOT
import array
import math
from array import *
from ROOT import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

PlotAllEnergies = False
PlotAllSRs = False
DoSkymap = True
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

#tag = "Small"
#file_theta2_lower = 0.0
#file_theta2_upper = 0.2
#UseRing = True
#UseRDBM = True
tag = "Large"
file_theta2_lower = 0.2
file_theta2_upper = 10.0
UseRDBM = True
UseDark = True

target = ""
source = []
#source += ['Everything']
#source += ['Segue1V6']
#source += ['IC443HotSpot']
#source += ['Crab']
#source += ['Mrk421']
#source += ['H1426']
#source += ['1ES0229']
source += ['PKS1424']
#source += ['3C264']
#source += ['G079']
#source += ['RGBJ0710']
#source += ['CasA']
#source += ['M82']
#source += ['WComaeV6']
#source += ['1ES1218V6']

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
energy_list += [1585]
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

NRGBs = 3
NCont = 512
stops = [0.,0.5,1.]
red = [0.2,0.5,1.]
green = [0.5,1.,0.5]
blue = [1.,0.5,0.2]
ROOT.TColor.CreateGradientColorTable(NRGBs,array('d',stops),array('d',red),array('d',green),array('d',blue),NCont)
ROOT.gStyle.SetNumberContours(NCont)

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
    for bx in range(1,Hist_Bkgd.GetNbinsX()+1):
        for by in range(1,Hist_Bkgd.GetNbinsY()+1):
            if not (bx<=bin_upper_x and by<=bin_upper_y):
                if not (bx>bin_upper_x and by>bin_upper_y):
                    nbins_cr += 1
                    data_cr = Hist_Data.GetBinContent(bx,by)
                    bkgd_cr = Hist_Dark.GetBinContent(bx,by)
                    if bkgd_cr>0:
                        syst_cr += abs(data_cr-bkgd_cr)
                        bkgd_cr_total += bkgd_cr
    syst_cr = syst_cr/bkgd_cr_total
    print 'syst_cr = %s'%(syst_cr)

    Hist_sys = Hist_Bkgd.Clone()
    for bx in range(1,Hist_Bkgd.GetNbinsX()+1):
        for by in range(1,Hist_Bkgd.GetNbinsY()+1):
            new_err = abs(Hist_Bkgd.GetBinContent(bx,by))*syst_cr
            old_err = Hist_sys.GetBinError(bx,by)
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

def Variation_ratio(Hist_SR, Hist_Bkg,range_lower,range_upper,syst):

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

def MakeComparisonPlot(Hists,legends,colors,title,name,logx,logy):
    
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
            Hists[h].GetXaxis().SetTitle(title)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
            if min_heigh > Hists[h].GetMinimum(): 
                min_heigh = Hists[h].GetMinimum()

    if 'MSCW' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        for h in range(0,len(Hists)):
            Hists[h].GetXaxis().SetRangeUser(low_end,high_end)

    Hists[max_hist].SetMaximum(max_heigh)
    Hists[max_hist].SetMinimum(min_heigh)
    Hists[max_hist].Draw("E")
    #Hists[0].Draw("E")

    for h in range(0,len(Hists)):
        #if colors[h]==1 or colors[h]==2: Hists[0].SetLineWidth(3)
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].SetLineWidth(2)
            Hists[h].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
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
            if Hists[h]!=0:
                scale = Hists[0].Integral()/Hists[h].Integral()
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

    fill_color = [38,30,46]
    if doSum:
        stack = ROOT.THStack("stack", "")
        Hist_Sum = Hists[1].Clone()
        Hist_Sum.Reset()
        Hist_Sys = Hists[1].Clone()
        Hist_Sys.Reset()
        for h in range(1,len(Hists)):
            if not legends[h]=='syst.':
                set_histStyle( Hists[h] , fill_color[h-1])
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


    Hists[0].Draw("E same")
    #for h in range(0,len(Hists)):
    #    if Hists[h]!=0:
    #        Hists[h].Draw("E same")
    Hists[0].SetLineWidth(3)
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
    predict_bkg, err_bkg = IntegralAndError(Hist_Sum,norm_bin_low_target,norm_bin_up_target)
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
    Hist_Band.GetYaxis().SetTitle("Data/BKG")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
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
    Hist_Ratio = Hists[0].Clone()
    Hist_Ratio.Divide(Hist_Sum)
    Hist_Ratio.SetLineWidth(2)
    Hist_Ratio.Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_%s.pdf'%(name,target,tag))


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
    Hist_Data.Draw("CONT3 same")
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

def Make2DSignificancePlot(Hist_SR,Hist_Bkg,xtitle,ytitle,theta2_low,theta2_up,name):

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
    for bx in range(0,Hist_SR.GetNbinsX()):
        for by in range(0,Hist_SR.GetNbinsY()):
            if Hist_Bkg.GetBinContent(bx+1,by+1)==0: continue
            NSR = Hist_SR.GetBinContent(bx+1,by+1)
            NSR_Err = Hist_SR.GetBinError(bx+1,by+1)
            NBkg = Hist_Bkg.GetBinContent(bx+1,by+1)
            NBkg_Err = Hist_Bkg.GetBinError(bx+1,by+1)
            Sig = 1.*CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Data.SetBinContent(bx+1,by+1,Sig)
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.SetMaximum(5)
    Hist_Data.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapSig_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapTotal_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_Bkg.Clone()
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    canvas.SaveAs('output_plots/SkymapBkg_%s_%s_%s.png'%(name,target,tag))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.Add(Hist_Bkg,-1.)
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
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

def Make2DSignificancePlotShowerShape(Hist_SR,Hist_Bkg,xtitle,ytitle,theta2_low,theta2_up,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
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
    Hist_Excess.SetMaximum(0.5)
    Hist_Excess.SetMinimum(-0.5)
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


for s in range(0,len(source)):

    target = source[s]

    Hist2D_Data_SumE = ROOT.TH2D("Hist2D_Data_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_Dark_SumE = ROOT.TH2D("Hist2D_Dark_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist2D_Bkgd_SumE = ROOT.TH2D("Hist2D_Bkgd_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Data_MSCW_SumE = ROOT.TH1D("Hist_Data_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Ring_MSCW_SumE = ROOT.TH1D("Hist_Ring_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Bkgd_MSCW_SumE = ROOT.TH1D("Hist_Bkgd_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Dark_MSCW_SumE = ROOT.TH1D("Hist_Dark_MSCW_SumE","",N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper)
    Hist_Data_MSCL_SumE = ROOT.TH1D("Hist_Data_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Ring_MSCL_SumE = ROOT.TH1D("Hist_Ring_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Bkgd_MSCL_SumE = ROOT.TH1D("Hist_Bkgd_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Dark_MSCL_SumE = ROOT.TH1D("Hist_Dark_MSCL_SumE","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper)
    Hist_Data_Theta2 = ROOT.TH1D("Hist_Data_Theta2","",1024,0,10)
    Hist_Data_Skymap = ROOT.TH2D("Hist_Data_Skymap","",150,-3,3,150,-3,3)
    Hist_Bkgd_Theta2 = ROOT.TH1D("Hist_Bkgd_Theta2","",1024,0,10)
    Hist_Bkgd_Skymap = ROOT.TH2D("Hist_Bkgd_Skymap","",150,-3,3,150,-3,3)
    Hist_Bkgd_Syst_Theta2 = ROOT.TH1D("Hist_Bkgd_Syst_Theta2","",1024,0,10)
    Hist_Bkgd_Syst_Skymap = ROOT.TH2D("Hist_Bkgd_Syst_Skymap","",150,-3,3,150,-3,3)
    ring_integral = 0
    dark_integral = 0
    ring_integral_err = 0
    dark_integral_err = 0
    exposure_hours = 0.

    FilePath = "output_Jul16/Netflix_"+target+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_Theta2%sto%s"%(int(file_theta2_lower*10.),int(file_theta2_upper*10.))+".root";
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours

    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]
        HistName = "Hist_FourierCoeff_Eigenvector_0_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_0_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_Eigenvector_0_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_0_Imag = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_0_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_0_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_0_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_0_Imag = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_Eigenvector_1_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_1_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_Eigenvector_1_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_1_Imag = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_1_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_1_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_1_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_1_Imag = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_Eigenvector_2_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_2_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_Eigenvector_2_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_Eigenvector_2_Imag = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_2_Real_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_2_Real = InputFile.Get(HistName)
        HistName = "Hist_FourierCoeff_InvEigenvector_2_Imag_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_FourierCoeff_InvEigenvector_2_Imag = InputFile.Get(HistName)
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
        HistName = "Hist_Data_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_Data = InputFile.Get(HistName)
        HistName = "Hist_Ring_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_Ring = InputFile.Get(HistName)
        HistName = "Hist_Dark_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_Dark = InputFile.Get(HistName)
        HistName = "Hist_Bkgd_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_Bkgd = InputFile.Get(HistName)
        Hist2D_Bkgd = AddStat2D(Hist2D_Dark,Hist2D_Bkgd)
        Hist2D_Bkgd = AddSyst2D(Hist2D_Data,Hist2D_Dark,Hist2D_Bkgd)
        #Make2DProjectionPlot(Hist2D_Data,"MSCL","MSCW",'Data_2D_E%s'%(ErecS_lower_cut),False)
        #Make2DProjectionPlot(Hist2D_Bkgd,"MSCL","MSCW",'Bkgd_2D_E%s'%(ErecS_lower_cut),False)
        #Make2DProjectionPlot(Hist2D_Ring,"MSCL","MSCW",'Ring_2D_E%s'%(ErecS_lower_cut),False)
        #Make2DProjectionPlot(Hist2D_Dark,"MSCL","MSCW",'Dark_2D_E%s'%(ErecS_lower_cut),False)
        print HistName
        print 'Hist2D_Data.Integral() = %s'%(Hist2D_Data.Integral())
        print 'Hist2D_Bkgd.Integral() = %s'%(Hist2D_Bkgd.Integral())
        bin_lower = Hist2D_Data.GetXaxis().FindBin(MSCL_lower_cut)
        bin_upper = Hist2D_Data.GetXaxis().FindBin(MSCL_blind_cut)-1
        Hist1D_Data_MSCW_SumSRs = Hist2D_Data.ProjectionY("Hist1D_Data_MSCW_SumSRs",bin_lower,bin_upper)
        Hist1D_Ring_MSCW_SumSRs = Hist2D_Ring.ProjectionY("Hist1D_Ring_MSCW_SumSRs",bin_lower,bin_upper)
        Hist1D_Bkgd_MSCW_SumSRs = Hist2D_Bkgd.ProjectionY("Hist1D_Bkgd_MSCW_SumSRs",bin_lower,bin_upper)
        Hist1D_Dark_MSCW_SumSRs = Hist2D_Dark.ProjectionY("Hist1D_Dark_MSCW_SumSRs",bin_lower,bin_upper)
        bin_lower = Hist2D_Data.GetYaxis().FindBin(MSCW_lower_cut)
        bin_upper = Hist2D_Data.GetYaxis().FindBin(MSCW_blind_cut)-1
        Hist1D_Data_MSCL_SumSRs = Hist2D_Data.ProjectionX("Hist1D_Data_MSCL_SumSRs",bin_lower,bin_upper)
        Hist1D_Ring_MSCL_SumSRs = Hist2D_Ring.ProjectionX("Hist1D_Ring_MSCL_SumSRs",bin_lower,bin_upper)
        Hist1D_Bkgd_MSCL_SumSRs = Hist2D_Bkgd.ProjectionX("Hist1D_Bkgd_MSCL_SumSRs",bin_lower,bin_upper)
        Hist1D_Dark_MSCL_SumSRs = Hist2D_Dark.ProjectionX("Hist1D_Dark_MSCL_SumSRs",bin_lower,bin_upper)

        Hist2D_Data_SumE.Add(Hist2D_Data)
        Hist2D_Dark_SumE.Add(Hist2D_Dark)
        Hist2D_Bkgd_SumE.Add(Hist2D_Bkgd)
        Hist_Data_MSCW_SumE.Add(Hist1D_Data_MSCW_SumSRs)
        Hist_Ring_MSCW_SumE.Add(Hist1D_Ring_MSCW_SumSRs)
        Hist_Bkgd_MSCW_SumE.Add(Hist1D_Bkgd_MSCW_SumSRs)
        Hist_Dark_MSCW_SumE.Add(Hist1D_Dark_MSCW_SumSRs)
        Hist_Data_MSCL_SumE.Add(Hist1D_Data_MSCL_SumSRs)
        Hist_Ring_MSCL_SumE.Add(Hist1D_Ring_MSCL_SumSRs)
        Hist_Bkgd_MSCL_SumE.Add(Hist1D_Bkgd_MSCL_SumSRs)
        Hist_Dark_MSCL_SumE.Add(Hist1D_Dark_MSCL_SumSRs)

        bin_lower = Hist1D_Bkgd_MSCW_SumSRs.FindBin(MSCW_lower_cut)
        bin_upper = Hist1D_Bkgd_MSCW_SumSRs.FindBin(MSCW_blind_cut)-1

        HistName = "Hist_Data_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Data_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Data_Theta2.Add(Hist_Data_Theta2_Tmp)
        HistName = "Hist_Data_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Data_Skymap.Add(InputFile.Get(HistName))

        bkg_total, bkg_err = IntegralAndError(Hist1D_Bkgd_MSCW_SumSRs,bin_lower,bin_upper)
        HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Bkgd_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Bkgd_Syst_Theta2_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Bkgd_Theta2_Tmp.Integral()
        bkgd_scale = 0
        bkgd_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            bkgd_scale = bkg_total/old_integral
            bkgd_scale_err = bkgd_scale*(bkg_err/bkg_total)
        else:
            bkgd_scale = 0
            bkgd_scale_err = 0
        print 'bkgd_scale = %s, bkgd_scale_err = %s'%(bkgd_scale,bkgd_scale_err)
        if not bkg_total==0:
            Theta2HistScale(Hist_Bkgd_Theta2_Tmp,Hist_Bkgd_Syst_Theta2_Tmp,bkgd_scale,bkgd_scale_err)
        else:
            Hist_Bkgd_Theta2_Tmp.Scale(0)
        Hist_Bkgd_Theta2.Add(Hist_Bkgd_Theta2_Tmp)
        Hist_Bkgd_Syst_Theta2.Add(Hist_Bkgd_Syst_Theta2_Tmp)

        HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Bkgd_Skymap_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Bkgd_Skymap_Tmp.Integral()
        scale = 0
        scale_err = 0
        Hist_Bkgd_Syst_Skymap_Tmp = InputFile.Get(HistName)
        Hist_Bkgd_Syst_Skymap_Tmp.Reset()
        if not bkg_total==0 and not old_integral==0:
            scale = bkg_total/old_integral
            scale_err = scale*(bkg_err/bkg_total)
        else:
            scale = 0
            scale_err = 0
        if not bkg_total==0:
            RaDecHistScale(Hist_Bkgd_Skymap_Tmp,Hist_Bkgd_Syst_Skymap_Tmp,scale,scale_err)
        else:
            Hist_Bkgd_Skymap_Tmp.Scale(0)
        Hist_Bkgd_Skymap.Add(Hist_Bkgd_Skymap_Tmp)
        Hist_Bkgd_Syst_Skymap.Add(Hist_Bkgd_Syst_Skymap_Tmp)

        if PlotAllEnergies:
            if UseRDBM:
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_FourierCoeff_Eigenvector_0_Real]
                legends += ['1st eigenvec sin']
                colors += [1]
                Hists += [Hist_FourierCoeff_Eigenvector_0_Imag]
                legends += ['1st eigenvec cos']
                colors += [2]
                Hists += [Hist_FourierCoeff_InvEigenvector_0_Real]
                legends += ['1st inverse eigenvec sin']
                colors += [3]
                Hists += [Hist_FourierCoeff_InvEigenvector_0_Imag]
                legends += ['1st inverse eigenvec cos']
                colors += [4]
                plotname = 'FourierCoeff_Eignvec_0_E%s'%(ErecS_lower_cut)
                title = 'Fourier mode'
                MakeComparisonPlot(Hists,legends,colors,title,plotname,False,False)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_FourierCoeff_Eigenvector_1_Real]
                legends += ['2nd eigenvec sin']
                colors += [1]
                Hists += [Hist_FourierCoeff_Eigenvector_1_Imag]
                legends += ['2nd eigenvec cos']
                colors += [2]
                Hists += [Hist_FourierCoeff_InvEigenvector_1_Real]
                legends += ['2nd inverse eigenvec sin']
                colors += [3]
                Hists += [Hist_FourierCoeff_InvEigenvector_1_Imag]
                legends += ['2nd inverse eigenvec cos']
                colors += [4]
                plotname = 'FourierCoeff_Eignvec_1_E%s'%(ErecS_lower_cut)
                title = 'Fourier mode'
                MakeComparisonPlot(Hists,legends,colors,title,plotname,False,False)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_FourierCoeff_Eigenvector_2_Real]
                legends += ['3rd eigenvec sin']
                colors += [1]
                Hists += [Hist_FourierCoeff_Eigenvector_2_Imag]
                legends += ['3rd eigenvec cos']
                colors += [2]
                Hists += [Hist_FourierCoeff_InvEigenvector_2_Real]
                legends += ['3rd inverse eigenvec sin']
                colors += [3]
                Hists += [Hist_FourierCoeff_InvEigenvector_2_Imag]
                legends += ['3rd inverse eigenvec cos']
                colors += [4]
                plotname = 'FourierCoeff_Eignvec_2_E%s'%(ErecS_lower_cut)
                title = 'Fourier mode'
                MakeComparisonPlot(Hists,legends,colors,title,plotname,False,False)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_Data_EigenvectorReal_0]
                legends += ['1st']
                colors += [1]
                Hists += [Hist_Data_EigenvectorReal_1]
                legends += ['2nd']
                colors += [2]
                Hists += [Hist_Data_EigenvectorReal_2]
                legends += ['3rd']
                colors += [3]
                plotname = 'Target_Eigenvector_E%s'%(ErecS_lower_cut)
                title = 'eigenvector'
                MakeComparisonPlot(Hists,legends,colors,title,plotname,False,False)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCW_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Bkgd_MSCW_SumSRs]
                legends += ['bkg']
                colors += [4]
                plotname = 'Target_W_SumSRs_RDBM_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCL_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Bkgd_MSCL_SumSRs]
                legends += ['bkg']
                colors += [4]
                plotname = 'Target_L_SumSRs_RDBM_E%s'%(ErecS_lower_cut)
                title = 'MSCL'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
            if UseRing:
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCW_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Ring_MSCW_SumSRs]
                legends += ['ring']
                colors += [4]
                plotname = 'Target_W_SumSRs_Ring_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCL_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Ring_MSCL_SumSRs]
                legends += ['ring']
                colors += [4]
                plotname = 'Target_L_SumSRs_Ring_E%s'%(ErecS_lower_cut)
                title = 'MSCL'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
            if UseDark:
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCW_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Dark_MSCW_SumSRs]
                legends += ['dark']
                colors += [4]
                plotname = 'Target_W_SumSRs_Dark_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_Data_MSCL_SumSRs]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_Dark_MSCL_SumSRs]
                legends += ['dark']
                colors += [4]
                plotname = 'Target_L_SumSRs_Dark_E%s'%(ErecS_lower_cut)
                title = 'MSCL'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)

    ErecS_lower_cut = energy_list[0]

    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Bkgd_MSCW_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_W_SumSRs_SumE_RDBM'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Bkgd_MSCL_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_L_SumSRs_SumE_RDBM'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
        Hist_Sig = Make2DSignificancePlotShowerShape(Hist2D_Data_SumE,Hist2D_Bkgd_SumE,'MSCL','MSCW',0,0,'RDBM_Sig2D')
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Ring_MSCW_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_W_SumSRs_SumE_Ring'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Ring_MSCL_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_L_SumSRs_SumE_Ring'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCW_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_W_SumSRs_SumE_Dark'
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCL_SumE]
        legends += ['bkg']
        colors += [4]
        plotname = 'Target_L_SumSRs_SumE_Dark'
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCL_lower_cut,MSCL_blind_cut,-1)
        Hist_Sig = Make2DSignificancePlotShowerShape(Hist2D_Data_SumE,Hist2D_Dark_SumE,'MSCL','MSCW',0,0,'Dark_Sig2D')

    n_rebin = 1
    if theta2_upper>0.25: n_rebin = 4
    if theta2_upper>1.00: n_rebin = 16
    Hist_Bkgd_Theta2.Rebin(n_rebin)
    Hist_Bkgd_Syst_Theta2.Rebin(n_rebin)
    Hist_Data_Theta2.Rebin(n_rebin)
    
    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Bkgd_Theta2]
        legends += ['bkg']
        colors += [4]
        Hists += [Hist_Bkgd_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_RDBM_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,-1)
    
    if DoSkymap and UseRDBM:
        
        plotname = 'Target_SRall_RDBM_Skymap_E%s'%(ErecS_lower_cut)
        n_rebin = 2
        Hist_Data_Skymap.Rebin2D(n_rebin,n_rebin)
        Hist_Bkgd_Skymap.Rebin2D(n_rebin,n_rebin)
        Hist_Bkgd_Syst_Skymap.Rebin2D(n_rebin,n_rebin)
        Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Bkgd_Skymap,'RA','Dec',theta2_lower,theta2_upper,plotname)
        
        smooth_size = 0.05
        Hist_Data_Skymap_smooth = Smooth2DMap(Hist_Data_Skymap,smooth_size,False)
        Hist_Bkgd_Skymap_smooth = Smooth2DMap(Hist_Bkgd_Skymap,smooth_size,False)
        Hist_Bkgd_Syst_Skymap_smooth = Smooth2DMap(Hist_Bkgd_Syst_Skymap,smooth_size,True)
        Hist_Bkgd_Skymap_smooth.Add(Hist_Bkgd_Syst_Skymap_smooth)
        plotname = 'Target_SRall_RDBM_Skymap_Smooth_E%s'%(ErecS_lower_cut)
        Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap_smooth,Hist_Bkgd_Skymap_smooth,'RA','Dec',theta2_lower,theta2_upper,plotname)

