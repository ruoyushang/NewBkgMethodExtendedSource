
import sys,ROOT
import array
import math
from array import *
from ROOT import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

PlotAllSRs = False
DoSkymap = False
UseRDBM = False
UseDark = False
UseRing = False

theta2_lower = 0.2
theta2_upper = 10.0

file_elev_lower = 70
file_elev_upper = 85

#tag = "Small"
#file_theta2_lower = 0.0
#file_theta2_upper = 0.2
#UseRing = True
tag = "Large"
file_theta2_lower = 0.0
file_theta2_upper = 10.0
UseRDBM = True
UseDark = True

target = ""
source = []
source += ['Everything']
#source += ['Segue1V6']
#source += ['IC443HotSpot']
#source += ['Crab']
#source += ['Mrk421']
#source += ['H1426']
#source += ['1ES0229']
#source += ['PKS1424']
#source += ['3C264']
#source += ['G079']
#source += ['RGBJ0710']
#source += ['CasA']
#source += ['M82']
#source += ['WComaeV6']
#source += ['1ES1218V6']

energy_list = []
energy_list += [200]
energy_list += [237]
energy_list += [282]
energy_list += [335]
energy_list += [398] # works from here
energy_list += [473]
energy_list += [562]
energy_list += [667]
energy_list += [794]
energy_list += [943]
energy_list += [1122]
energy_list += [1332]
energy_list += [1585]
energy_list += [1882]
energy_list += [2239]
energy_list += [3162]
energy_list += [4467]
energy_list += [6310]
energy_list += [8913]

MSCW_lower_cut = -1.0
MSCW_upper_cut = 1.0
MSCW_blind_cut = 1.0
exposure_hours = 0.
total_exposure_hours = 0.

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

def AddSystCR(Hist,bin1,bin2,Hist_CR_Data,Hist_CR_Bkgd):

    Hist_sys = Hist.Clone()
    for b in range(1,Hist.GetNbinsX()+1):
        #if b<bin1: continue
        #if b>bin2: continue
        syst = 0
        norm = 0.5*(Hist_CR_Data.GetBinContent(b)+Hist_CR_Bkgd.GetBinContent(b))
        if norm>0:
            syst = 2.*abs(Hist_CR_Data.GetBinContent(b)-Hist_CR_Bkgd.GetBinContent(b))/norm
        if norm>0:
            old_content = Hist.GetBinContent(b)
            new_content = Hist.GetBinContent(b)*(Hist_CR_Data.GetBinContent(b)-Hist_CR_Bkgd.GetBinContent(b))/norm
            Hist.SetBinContent(b,old_content+new_content)
        old_err = Hist.GetBinError(b)
        new_err = Hist.GetBinContent(b)*syst
        Hist_sys.SetBinContent(b,0)
        Hist_sys.SetBinError(b,new_err)
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
    #line1 = ROOT.TLine(MSCL_lower_cut,MSCW_blind_cut,MSCL_upper_cut,MSCW_blind_cut)
    #line1.SetLineStyle(1)
    #line1.SetLineColor(2)
    #line1.SetLineWidth(2)
    #line1.Draw()
    #line2 = ROOT.TLine(MSCL_upper_cut,MSCW_blind_cut,MSCL_upper_cut,MSCW_lower_cut)
    #line2.SetLineStyle(1)
    #line2.SetLineColor(2)
    #line2.SetLineWidth(2)
    #line2.Draw()
    #line3 = ROOT.TLine(MSCL_lower_cut,MSCW_lower_cut,MSCL_upper_cut,MSCW_lower_cut)
    #line3.SetLineStyle(1)
    #line3.SetLineColor(2)
    #line3.SetLineWidth(2)
    #line3.Draw()
    zenith_mean = 0
    azimuth_mean_1 = 0
    azimuth_mean_2 = 0
    azimuth_rms_1 = 0
    azimuth_rms_2 = 0
    total_weight = 0
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            x_pos = Hist_Data.GetXaxis().GetBinCenter(bx+1)
            y_pos = Hist_Data.GetYaxis().GetBinCenter(by+1)
            zenith_mean += Hist_Data.GetBinContent(bx+1,by+1)*y_pos
            azimuth_mean_1 += Hist_Data.GetBinContent(bx+1,by+1)*x_pos
            if x_pos<180.:
                azimuth_mean_2 += Hist_Data.GetBinContent(bx+1,by+1)*x_pos
            else:
                azimuth_mean_2 += Hist_Data.GetBinContent(bx+1,by+1)*(x_pos-360.)
            total_weight += Hist_Data.GetBinContent(bx+1,by+1)
    if total_weight>0:
        zenith_mean = zenith_mean/total_weight
        azimuth_mean_1 = azimuth_mean_1/total_weight
        azimuth_mean_2 = azimuth_mean_2/total_weight
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            x_pos = Hist_Data.GetXaxis().GetBinCenter(bx+1)
            y_pos = Hist_Data.GetYaxis().GetBinCenter(by+1)
            azimuth_rms_1 += Hist_Data.GetBinContent(bx+1,by+1)*pow(x_pos-azimuth_mean_1,2)
            if x_pos<180.:
                azimuth_rms_2 += Hist_Data.GetBinContent(bx+1,by+1)*pow(x_pos-azimuth_mean_2,2)
            else:
                azimuth_rms_2 += Hist_Data.GetBinContent(bx+1,by+1)*pow(x_pos-360.-azimuth_mean_2,2)
    azimuth_mean = 0
    if azimuth_rms_1<azimuth_rms_2:
        azimuth_mean = azimuth_mean_1
    else:
        azimuth_mean = azimuth_mean_2

    lumilab1 = ROOT.TLatex(0.5,0.8,'avg. zenith = %0.1f'%(zenith_mean) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.03)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.5,0.75,'avg. azimuth = %0.1f'%(azimuth_mean) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.03)
    lumilab2.Draw()

    canvas.SaveAs('output_plots/%s_%s_%s.png'%(name,target,tag))

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
        low_end = -2.
        high_end = 4.
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

    pad2.cd()
    Hist_Band = Hist_Sum.Clone()
    #Hist_Band.Divide(Hists[0])
    for b in range(0,Hist_Band.GetNbinsX()):
        Hist_Band.SetBinContent(b+1,0)
        Hist_Band.SetBinError(b+1,3)
        #Hist_Band.SetBinContent(b+1,1)
        #Hist_Band.SetBinError(b+1,0.1)
    Hist_Band.SetLineStyle(1)
    Hist_Band.SetLineColor(1)
    Hist_Band.SetLineWidth(2)
    Hist_Band.SetFillColor(1)
    Hist_Band.SetFillStyle(3004)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("#sqrt{#chi^{2}} (per bin)")
    #Hist_Band.GetYaxis().SetTitle("Data/BKG")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(5)
    Hist_Band.SetMinimum(-5)
    #Hist_Band.SetMaximum(1.2)
    #Hist_Band.SetMinimum(0.8)
    Hist_Band.GetXaxis().SetRangeUser(low_end,high_end)
    Hist_Band.Draw("e2")
    Hist_Ratio = Hists[0].Clone()
    #Hist_Ratio.Divide(Hist_Sum)
    for b in range(0,Hist_Ratio.GetNbinsX()):
        ndata = Hists[0].GetBinContent(b+1)
        edata = Hists[0].GetBinError(b+1)
        nbkg = Hist_Sum.GetBinContent(b+1)
        ebkg = Hist_Err.GetBinError(b+1)
        if not (ebkg*ebkg+edata*edata)==0 and not (ndata-nbkg)==0:
            Hist_Ratio.SetBinContent(b+1,(ndata-nbkg)/abs(ndata-nbkg)*pow((ndata-nbkg)*(ndata-nbkg)/(ebkg*ebkg+edata*edata),0.5))
        else:
            Hist_Ratio.SetBinContent(b+1,0)
        Hist_Ratio.SetBinError(b+1,0)
    Hist_Ratio.SetFillColor(1)
    Hist_Ratio.SetLineWidth(2)
    Hist_Ratio.Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_%s.png'%(name,target,tag))

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

def MakeComparisonPlotNoRatio(Hists,legends,colors,title,name,minheight,maxheight,logx,logy):
    
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
    Hists[max_hist].GetYaxis().SetTitle("R = data / bkg")
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
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")

    if logx: 
        pad1.SetLogx()

    pad2.cd()
    Hist_Sig = ROOT.TH1D("Hist_Sig","",32,-8,8)
    for h in range(0,len(Hists)):
        for b in range(0,Hists[h].GetNbinsX()):
            if not Hists[h].GetBinError(b+1)==0:
                Hist_Sig.Fill((Hists[h].GetBinContent(b+1)-1.)/Hists[h].GetBinError(b+1))
    Hist_Sig.GetXaxis().SetLabelSize(0.06)
    Hist_Sig.GetYaxis().SetLabelSize(0.06)
    Hist_Sig.GetXaxis().SetTitleSize(0.06)
    Hist_Sig.GetXaxis().SetTitle('s = (R-1)/#DeltaR')
    Hist_Sig.GetXaxis().SetTitleOffset(1.0)
    Hist_Sig.Draw("E")

    lumilab1 = ROOT.TLatex(0.15,0.80,'mean = %0.1f, RMS = %0.1f'%(Hist_Sig.GetMean(),Hist_Sig.GetRMS()) )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.10)
    lumilab1.Draw()

    c_both.SaveAs('output_plots/%s_%s.png'%(name,tag))

def MakeComparisonPlot(Hists,legends,colors,title,name,maxhight,logx,logy):
    
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

    if not maxhight==0: Hists[max_hist].SetMaximum(maxhight)
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
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))
    return Hist_Sig



Hist_RDBM_S2B = []
Hist_Dark_S2B = []
Hist_Ring_S2B = []
legend_S2B = []
color_S2B = []

for s in range(0,len(source)):

    target = source[s]

    Hist_Hadr_Skymap = ROOT.TH2D("Hist_Hadr_Skymap","",150,-3,3,150,-3,3)
    Hist_Elec_Skymap = ROOT.TH2D("Hist_Elec_Skymap","",150,-3,3,150,-3,3)
    Hist_Hadr_Syst_Skymap = ROOT.TH2D("Hist_Hadr_Syst_Skymap","",150,-3,3,150,-3,3)
    Hist_Elec_Syst_Skymap = ROOT.TH2D("Hist_Elec_Syst_Skymap","",150,-3,3,150,-3,3)
    Hist_Data_Skymap = ROOT.TH2D("Hist_Data_Skymap","",150,-3,3,150,-3,3)
    Hist_Data_Theta2 = ROOT.TH1D("Hist_Data_Theta2","",1024,0,10)
    Hist_Hadr_Theta2 = ROOT.TH1D("Hist_Hadr_Theta2","",1024,0,10)
    Hist_Hadr_Syst_Theta2 = ROOT.TH1D("Hist_Hadr_Syst_Theta2","",1024,0,10)
    Hist_Elec_Theta2 = ROOT.TH1D("Hist_Elec_Theta2","",1024,0,10)
    Hist_Bkgd_Theta2 = ROOT.TH1D("Hist_Bkgd_Theta2","",1024,0,10)
    Hist_Ring_Theta2 = ROOT.TH1D("Hist_Ring_Theta2","",1024,0,10)
    Hist_Ring_Syst_Theta2 = ROOT.TH1D("Hist_Ring_Syst_Theta2","",1024,0,10)
    Hist_Dark_Theta2 = ROOT.TH1D("Hist_Dark_Theta2","",1024,0,10)
    Hist_Dark_Syst_Theta2 = ROOT.TH1D("Hist_Dark_Syst_Theta2","",1024,0,10)
    Hist_Data_MSCW_SumE = ROOT.TH1D("Hist_Data_MSCW_SumE","",960,-30,30)
    Hist_Hadr_MSCW_SumE = ROOT.TH1D("Hist_Hadr_MSCW_SumE","",960,-30,30)
    Hist_Elec_MSCW_SumE = ROOT.TH1D("Hist_Elec_MSCW_SumE","",960,-30,30)
    Hist_CR_Data_MSCW_SumE = ROOT.TH1D("Hist_CR_Data_MSCW_SumE","",960,-30,30)
    Hist_CR_Hadr_MSCW_SumE = ROOT.TH1D("Hist_CR_Hadr_MSCW_SumE","",960,-30,30)
    Hist_CR_Elec_MSCW_SumE = ROOT.TH1D("Hist_CR_Elec_MSCW_SumE","",960,-30,30)
    Hist_Data_MSCL_SumE = ROOT.TH1D("Hist_Data_MSCL_SumE","",960,-30,30)
    Hist_Hadr_MSCL_SumE = ROOT.TH1D("Hist_Hadr_MSCL_SumE","",960,-30,30)
    Hist_Elec_MSCL_SumE = ROOT.TH1D("Hist_Elec_MSCL_SumE","",960,-30,30)
    Hist_CR_Data_MSCL_SumE = ROOT.TH1D("Hist_CR_Data_MSCL_SumE","",960,-30,30)
    Hist_CR_Hadr_MSCL_SumE = ROOT.TH1D("Hist_CR_Hadr_MSCL_SumE","",960,-30,30)
    Hist_CR_Elec_MSCL_SumE = ROOT.TH1D("Hist_CR_Elec_MSCL_SumE","",960,-30,30)
    Hist_Ring_MSCW_SumE = ROOT.TH1D("Hist_Ring_MSCW_SumE","",960,-30,30)
    Hist_Ring_Syst_MSCW_SumE = ROOT.TH1D("Hist_Ring_Syst_MSCW_SumE","",960,-30,30)
    Hist_Dark_MSCW_SumE = ROOT.TH1D("Hist_Dark_MSCW_SumE","",960,-30,30)
    Hist_Dark_MSCL_SumE = ROOT.TH1D("Hist_Dark_MSCL_SumE","",960,-30,30)
    ring_integral = 0
    dark_integral = 0
    ring_integral_err = 0
    dark_integral_err = 0
    exposure_hours = 0.

    #FilePath = "output_backup/Deconvolution_"+target+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_Theta2%sto%s"%(int(file_theta2_lower*10.),int(file_theta2_upper*10.))+".root";
    FilePath = "output_Jul16/Deconvolution_"+target+"_TelElev%sto%s"%(int(file_elev_lower),int(file_elev_upper))+"_Theta2%sto%s"%(int(file_theta2_lower*10.),int(file_theta2_upper*10.))+".root";
    InputFile = ROOT.TFile(FilePath)
    InfoTree = InputFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    exposure_hours = InfoTree.exposure_hours

    Hist_RDBM_S2B += [ROOT.TH1D("Hist_RDBM_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    Hist_Dark_S2B += [ROOT.TH1D("Hist_Dark_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    Hist_Ring_S2B += [ROOT.TH1D("Hist_Ring_S2B_%s"%(target),"",len(energy_list)-1,array('d',energy_list))]
    legend_S2B += ['%s (%0.1f hrs)'%(target,exposure_hours)]
    color_S2B += [s+1]

    HistName = "Hist_Dark_ShowerDirection"
    Hist_Dark_ShowerDirection = InputFile.Get(HistName)
    HistName = "Hist_Data_ShowerDirection"
    Hist_Data_ShowerDirection = InputFile.Get(HistName)
    Make2DProjectionPlot(Hist_Dark_ShowerDirection,'Azimuth','Zenith','Dark_ShowerDirection',False)
    Make2DProjectionPlot(Hist_Data_ShowerDirection,'Azimuth','Zenith','Data_ShowerDirection',False)

    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]
        HistName = "Hist_Target_Bkg_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Hadr = InputFile.Get(HistName)
        HistName = "Hist_Target_Ele_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Elec = InputFile.Get(HistName)
        HistName = "Hist_Data_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Data = InputFile.Get(HistName)
        HistName = "Hist_Target_Bkg_MSCWL_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_L_Hadr = InputFile.Get(HistName)
        HistName = "Hist_Target_Ele_MSCWL_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_L_Elec = InputFile.Get(HistName)
        HistName = "Hist_Data_MSCWL_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_L_Data = InputFile.Get(HistName)
        HistName = "Hist_Ring_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Ring = InputFile.Get(HistName)
        HistName = "Hist_Ring_Syst_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Ring_Syst = InputFile.Get(HistName)
        HistName = "Hist_DarkScaled_MSCLW_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_W_Dark = InputFile.Get(HistName)
        HistName = "Hist_DarkScaled_MSCWL_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist2D_L_Dark = InputFile.Get(HistName)
        print HistName
        print 'Hist2D_W_Data.Integral() = %s'%(Hist2D_W_Data.Integral())
        print 'Hist2D_L_Data.Integral() = %s'%(Hist2D_L_Data.Integral())
        Hist1D_W_Hadr_SumSRs = Hist2D_W_Hadr.ProjectionY("Hist1D_W_Hadr",1,1)
        Hist1D_W_Elec_SumSRs = Hist2D_W_Elec.ProjectionY("Hist1D_W_Elec",1,1)
        Hist1D_W_Data_SumSRs = Hist2D_W_Data.ProjectionY("Hist1D_W_Data",1,1)
        Hist1D_W_CR_Hadr_SumSRs = Hist2D_W_Hadr.ProjectionY("Hist1D_W_CR_Hadr",1,1)
        Hist1D_W_CR_Elec_SumSRs = Hist2D_W_Elec.ProjectionY("Hist1D_W_CR_Elec",1,1)
        Hist1D_W_CR_Data_SumSRs = Hist2D_W_Data.ProjectionY("Hist1D_W_CR_Data",1,1)
        Hist1D_L_Hadr_SumSRs = Hist2D_L_Hadr.ProjectionY("Hist1D_L_Hadr",1,1)
        Hist1D_L_Elec_SumSRs = Hist2D_L_Elec.ProjectionY("Hist1D_L_Elec",1,1)
        Hist1D_L_Data_SumSRs = Hist2D_L_Data.ProjectionY("Hist1D_L_Data",1,1)
        Hist1D_L_CR_Hadr_SumSRs = Hist2D_L_Hadr.ProjectionY("Hist1D_L_CR_Hadr",1,1)
        Hist1D_L_CR_Elec_SumSRs = Hist2D_L_Elec.ProjectionY("Hist1D_L_CR_Elec",1,1)
        Hist1D_L_CR_Data_SumSRs = Hist2D_L_Data.ProjectionY("Hist1D_L_CR_Data",1,1)
        Hist1D_W_Ring_SumSRs = Hist2D_W_Ring.ProjectionY("Hist1D_W_Ring",1,1)
        Hist1D_W_Ring_Syst_SumSRs = Hist2D_W_Ring_Syst.ProjectionY("Hist1D_W_Ring",1,1)
        Hist1D_W_Dark_SumSRs = Hist2D_W_Dark.ProjectionY("Hist1D_W_Dark",1,1)
        Hist1D_L_Dark_SumSRs = Hist2D_L_Dark.ProjectionY("Hist1D_L_Dark",1,1)
        Hist1D_W_Hadr_SumSRs.Reset()
        Hist1D_W_Elec_SumSRs.Reset()
        Hist1D_W_Data_SumSRs.Reset()
        Hist1D_W_CR_Hadr_SumSRs.Reset()
        Hist1D_W_CR_Elec_SumSRs.Reset()
        Hist1D_W_CR_Data_SumSRs.Reset()
        Hist1D_L_Hadr_SumSRs.Reset()
        Hist1D_L_Elec_SumSRs.Reset()
        Hist1D_L_Data_SumSRs.Reset()
        Hist1D_L_CR_Hadr_SumSRs.Reset()
        Hist1D_L_CR_Elec_SumSRs.Reset()
        Hist1D_L_CR_Data_SumSRs.Reset()
        Hist1D_W_Ring_SumSRs.Reset()
        Hist1D_W_Ring_Syst_SumSRs.Reset()
        Hist1D_W_Dark_SumSRs.Reset()
        Hist1D_L_Dark_SumSRs.Reset()
        for region in range(0,Hist2D_W_Data.GetNbinsX()-1):
            Hist1D_W_Hadr = Hist2D_W_Hadr.ProjectionY("Hist1D_W_Hadr",region+1,region+1)
            Hist1D_W_Elec = Hist2D_W_Elec.ProjectionY("Hist1D_W_Elec",region+1,region+1)
            Hist1D_W_Data = Hist2D_W_Data.ProjectionY("Hist1D_W_Data",region+1,region+1)
            Hist1D_L_Hadr = Hist2D_L_Hadr.ProjectionY("Hist1D_L_Hadr",region+1,region+1)
            Hist1D_L_Elec = Hist2D_L_Elec.ProjectionY("Hist1D_L_Elec",region+1,region+1)
            Hist1D_L_Data = Hist2D_L_Data.ProjectionY("Hist1D_L_Data",region+1,region+1)
            Hist1D_W_Ring = Hist2D_W_Ring.ProjectionY("Hist1D_W_Ring",region+1,region+1)
            Hist1D_W_Ring_Syst = Hist2D_W_Ring_Syst.ProjectionY("Hist1D_W_Ring_Syst",region+1,region+1)
            Hist1D_W_Dark = Hist2D_W_Dark.ProjectionY("Hist1D_W_Dark",region+1,region+1)
            Hist1D_L_Dark = Hist2D_L_Dark.ProjectionY("Hist1D_L_Dark",region+1,region+1)
            if region+1<=6:
                Hist1D_W_Hadr_SumSRs.Add(Hist1D_W_Hadr)
                Hist1D_W_Elec_SumSRs.Add(Hist1D_W_Elec)
                Hist1D_W_Data_SumSRs.Add(Hist1D_W_Data)
                Hist1D_L_Hadr_SumSRs.Add(Hist1D_L_Hadr)
                Hist1D_L_Elec_SumSRs.Add(Hist1D_L_Elec)
                Hist1D_L_Data_SumSRs.Add(Hist1D_L_Data)
                Hist1D_W_Ring_SumSRs.Add(Hist1D_W_Ring)
                Hist1D_W_Ring_Syst_SumSRs.Add(Hist1D_W_Ring_Syst)
                Hist1D_W_Dark_SumSRs.Add(Hist1D_W_Dark)
                Hist1D_L_Dark_SumSRs.Add(Hist1D_L_Dark)
            else:
                Hist1D_W_CR_Hadr_SumSRs.Add(Hist1D_W_Hadr)
                Hist1D_W_CR_Elec_SumSRs.Add(Hist1D_W_Elec)
                Hist1D_W_CR_Data_SumSRs.Add(Hist1D_W_Data)
                Hist1D_L_CR_Hadr_SumSRs.Add(Hist1D_L_Hadr)
                Hist1D_L_CR_Elec_SumSRs.Add(Hist1D_L_Elec)
                Hist1D_L_CR_Data_SumSRs.Add(Hist1D_L_Data)
            if PlotAllSRs:
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_W_Data]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_W_Hadr]
                legends += ['hadron']
                colors += [4]
                Hists += [Hist1D_W_Elec]
                legends += ['electron (MC)']
                colors += [3]
                plotname = 'Target_SR%s_RDBM_E%s'%(region,ErecS_lower_cut)
                title = 'MSCW'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist1D_L_Data]
                legends += ['%s'%(target)]
                colors += [1]
                Hists += [Hist1D_L_Hadr]
                legends += ['hadron']
                colors += [4]
                Hists += [Hist1D_L_Elec]
                legends += ['electron (MC)']
                colors += [3]
                plotname = 'Target_SR%s_RDBM_E%s'%(region,ErecS_lower_cut)
                title = 'MSCL'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)

        bin_lower = Hist1D_W_Hadr_SumSRs.FindBin(MSCW_blind_cut)
        bin_upper = Hist1D_W_Hadr_SumSRs.FindBin(MSCW_blind_cut*3)-1

        w_cr_data_total, w_cr_data_err = IntegralAndError(Hist1D_W_Data_SumSRs,bin_lower,bin_upper)
        w_cr_hadr_total, w_cr_data_err = IntegralAndError(Hist1D_W_Hadr_SumSRs,bin_lower,bin_upper)
        l_cr_data_total, l_cr_data_err = IntegralAndError(Hist1D_L_Data_SumSRs,bin_lower,bin_upper)
        l_cr_hadr_total, l_cr_data_err = IntegralAndError(Hist1D_L_Hadr_SumSRs,bin_lower,bin_upper)

        #w_cr2_syst = 0
        #w_cr2_syst = abs(w_cr_data_total-w_cr_hadr_total)/w_cr_data_total
        #print 'w_cr2_syst = %s'%(w_cr2_syst)
        #l_cr2_syst = 0
        #l_cr2_syst = abs(l_cr_data_total-l_cr_hadr_total)/l_cr_data_total
        #print 'l_cr2_syst = %s'%(l_cr2_syst)

        bin_lower = Hist1D_W_Hadr_SumSRs.FindBin(MSCW_lower_cut)
        bin_upper = Hist1D_W_Hadr_SumSRs.FindBin(MSCW_blind_cut)-1

        w_cr_data_total, w_cr_data_err = IntegralAndError(Hist1D_W_CR_Data_SumSRs,bin_lower,bin_upper)
        w_cr_hadr_total, w_cr_data_err = IntegralAndError(Hist1D_W_CR_Hadr_SumSRs,bin_lower,bin_upper)
        l_cr_data_total, l_cr_data_err = IntegralAndError(Hist1D_L_CR_Data_SumSRs,bin_lower,bin_upper)
        l_cr_hadr_total, l_cr_data_err = IntegralAndError(Hist1D_L_CR_Hadr_SumSRs,bin_lower,bin_upper)

        w_cr1_syst = 0
        if not w_cr_data_total==0:
            w_cr1_syst = abs(w_cr_data_total-w_cr_hadr_total)/w_cr_data_total
        print 'w_cr1_syst = %s'%(w_cr1_syst)
        Hist1D_W_Hadr_Syst_SumSRs = AddSystCR(Hist1D_W_Hadr_SumSRs,bin_lower,bin_upper,Hist1D_W_CR_Data_SumSRs,Hist1D_W_CR_Hadr_SumSRs)
        l_cr1_syst = 0
        if not l_cr_data_total==0:
            l_cr1_syst = abs(l_cr_data_total-l_cr_hadr_total)/l_cr_data_total
        print 'l_cr1_syst = %s'%(l_cr1_syst)
        Hist1D_L_Hadr_Syst_SumSRs = AddSystCR(Hist1D_L_Hadr_SumSRs,bin_lower,bin_upper,Hist1D_L_CR_Data_SumSRs,Hist1D_L_CR_Hadr_SumSRs)

        Hist_Data_MSCW_SumE.Add(Hist1D_W_Data_SumSRs)
        Hist_Hadr_MSCW_SumE.Add(Hist1D_W_Hadr_SumSRs)
        Hist_Hadr_MSCW_SumE.Add(Hist1D_W_Hadr_Syst_SumSRs)
        Hist_Elec_MSCW_SumE.Add(Hist1D_W_Elec_SumSRs)
        Hist_CR_Data_MSCW_SumE.Add(Hist1D_W_CR_Data_SumSRs)
        Hist_CR_Hadr_MSCW_SumE.Add(Hist1D_W_CR_Hadr_SumSRs)
        Hist_CR_Elec_MSCW_SumE.Add(Hist1D_W_CR_Elec_SumSRs)
        Hist_Data_MSCL_SumE.Add(Hist1D_L_Data_SumSRs)
        Hist_Hadr_MSCL_SumE.Add(Hist1D_L_Hadr_SumSRs)
        Hist_Hadr_MSCL_SumE.Add(Hist1D_L_Hadr_Syst_SumSRs)
        Hist_Elec_MSCL_SumE.Add(Hist1D_L_Elec_SumSRs)
        Hist_CR_Data_MSCL_SumE.Add(Hist1D_L_CR_Data_SumSRs)
        Hist_CR_Hadr_MSCL_SumE.Add(Hist1D_L_CR_Hadr_SumSRs)
        Hist_CR_Elec_MSCL_SumE.Add(Hist1D_L_CR_Elec_SumSRs)
        Hist_Ring_MSCW_SumE.Add(Hist1D_W_Ring_SumSRs)
        Hist_Ring_Syst_MSCW_SumE.Add(Hist1D_W_Ring_Syst_SumSRs)
        Hist_Dark_MSCW_SumE.Add(Hist1D_W_Dark_SumSRs)
        Hist_Dark_MSCL_SumE.Add(Hist1D_L_Dark_SumSRs)


        w_bkg_total, w_bkg_err = IntegralAndError(Hist1D_W_Hadr_SumSRs,bin_lower,bin_upper)
        l_bkg_total, l_bkg_err = IntegralAndError(Hist1D_L_Hadr_SumSRs,bin_lower,bin_upper)
        w_sys_total, w_sys_err = IntegralAndError(Hist1D_W_Hadr_Syst_SumSRs,bin_lower,bin_upper)
        l_sys_total, l_sys_err = IntegralAndError(Hist1D_L_Hadr_Syst_SumSRs,bin_lower,bin_upper)
        w_weight = 1
        if (w_cr1_syst>0): w_weight = pow(1./w_cr1_syst,1)
        l_weight = 1
        if (l_cr1_syst>0): l_weight = pow(1./l_cr1_syst,1)
        bkg_total = w_bkg_total*w_weight + l_bkg_total*l_weight
        if (w_weight+l_weight)>0:
            bkg_total = bkg_total/(w_weight+l_weight)
        else:
            bkg_total = 0
        bkg_err = 0
        if (w_cr1_syst>0 and w_bkg_total>0):
            bkg_err += w_weight*(pow(w_bkg_err,2)+pow(w_sys_err,2))
        if (l_cr1_syst>0 and l_bkg_total>0):
            bkg_err += l_weight*(pow(l_bkg_err,2)+pow(l_sys_err,2))
        bkg_err = pow(bkg_err/(w_weight+l_weight),0.5)

        HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Hadr_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Hadr_Syst_Theta2_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Hadr_Theta2_Tmp.Integral()
        hadr_scale = 0
        hadr_scale_err = 0
        print 'w_bkg_total = %s, w_bkg_err = %s'%(w_bkg_total,w_bkg_err)
        print 'l_bkg_total = %s, l_bkg_err = %s'%(l_bkg_total,l_bkg_err)
        print 'w_sys_total = %s, w_sys_err = %s'%(w_sys_total,w_sys_err)
        print 'l_sys_total = %s, l_sys_err = %s'%(l_sys_total,l_sys_err)
        print 'bkg_total = %s, bkg_err = %s'%(bkg_total,bkg_err)
        if not bkg_total==0 and not old_integral==0:
            hadr_scale = bkg_total/old_integral
            hadr_scale_err = hadr_scale*(bkg_err/bkg_total)
        else:
            hadr_scale = 0
            hadr_scale_err = 0
        print 'hadr_scale = %s, hadr_scale_err = %s'%(hadr_scale,hadr_scale_err)
        if not bkg_total==0:
            Theta2HistScale(Hist_Hadr_Theta2_Tmp,Hist_Hadr_Syst_Theta2_Tmp,hadr_scale,hadr_scale_err)
        else:
            Hist_Hadr_Theta2_Tmp.Scale(0)
        Hist_Hadr_Theta2.Add(Hist_Hadr_Theta2_Tmp)
        Hist_Hadr_Syst_Theta2.Add(Hist_Hadr_Syst_Theta2_Tmp)
        HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Hadr_Skymap_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Hadr_Skymap_Tmp.Integral()
        scale = 0
        scale_err = 0
        Hist_Hadr_Syst_Skymap_Tmp = InputFile.Get(HistName)
        Hist_Hadr_Syst_Skymap_Tmp.Reset()
        if not bkg_total==0 and not old_integral==0:
            scale = bkg_total/old_integral
            scale_err = scale*(bkg_err/bkg_total)
        else:
            scale = 0
            scale_err = 0
        if not bkg_total==0:
            RaDecHistScale(Hist_Hadr_Skymap_Tmp,Hist_Hadr_Syst_Skymap_Tmp,scale,scale_err)
        else:
            Hist_Hadr_Skymap_Tmp.Scale(0)
        Hist_Hadr_Skymap.Add(Hist_Hadr_Skymap_Tmp)
        Hist_Hadr_Syst_Skymap.Add(Hist_Hadr_Syst_Skymap_Tmp)

        w_bkg_total, w_bkg_err = IntegralAndError(Hist1D_W_Elec_SumSRs,bin_lower,bin_upper)
        l_bkg_total, l_bkg_err = IntegralAndError(Hist1D_L_Elec_SumSRs,bin_lower,bin_upper)
        w_weight = 1
        if (w_cr1_syst>0): w_weight = pow(1./w_cr1_syst,1)
        l_weight = 1
        if (l_cr1_syst>0): l_weight = pow(1./l_cr1_syst,1)
        bkg_total = w_bkg_total*w_weight + l_bkg_total*l_weight
        if (w_weight+l_weight)>0:
            bkg_total = bkg_total/(w_weight+l_weight)
        else:
            bkg_total = 0
        bkg_err = pow(0.5*(w_bkg_err*w_bkg_err+l_bkg_err*l_bkg_err),0.5)

        HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Elec_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Elec_Syst_Theta2_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Elec_Theta2_Tmp.Integral()
        elec_scale = 0
        elec_scale_err = 0
        if not bkg_total==0 and not old_integral==0:
            elec_scale = bkg_total/old_integral
            elec_scale_err = elec_scale*(bkg_err/bkg_total)
        else:
            elec_scale = 0
            elec_scale_err = 0
        if not bkg_total==0:
            Theta2HistScale(Hist_Elec_Theta2_Tmp,Hist_Elec_Syst_Theta2_Tmp,elec_scale,elec_scale_err)
        else:
            Hist_Elec_Theta2_Tmp.Scale(0)
        Hist_Elec_Theta2.Add(Hist_Elec_Theta2_Tmp)
        HistName = "Hist_Data_CR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Elec_Skymap_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Elec_Skymap_Tmp.Integral()
        scale = 0
        scale_err = 0
        Hist_Elec_Syst_Skymap_Tmp = InputFile.Get(HistName)
        Hist_Elec_Syst_Skymap_Tmp.Reset()
        if not bkg_total==0 and not old_integral==0:
            scale = bkg_total/old_integral
            scale_err = scale*(bkg_err/bkg_total)
        else:
            scale = 0
            scale_err = 0
        if not bkg_total==0:
            RaDecHistScale(Hist_Elec_Skymap_Tmp,Hist_Elec_Syst_Skymap_Tmp,scale,scale_err)
        else:
            Hist_Elec_Skymap_Tmp.Scale(0)
        Hist_Elec_Skymap.Add(Hist_Elec_Skymap_Tmp)
        Hist_Elec_Syst_Skymap.Add(Hist_Elec_Syst_Skymap_Tmp)

        HistName = "Hist_Data_SR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Data_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Data_Theta2.Add(Hist_Data_Theta2_Tmp)
        HistName = "Hist_Data_SR_Skymap_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Data_Skymap.Add(InputFile.Get(HistName))

        bkg_total, bkg_err = IntegralAndError(Hist1D_W_Ring_SumSRs,bin_lower,bin_upper)
        sys_total, sys_err = IntegralAndError(Hist1D_W_Ring_Syst_SumSRs,bin_lower,bin_upper)
        bkg_err = pow(bkg_err*bkg_err+sys_total*sys_total,0.5)
        ring_integral += bkg_total
        ring_integral_err += bkg_err*bkg_err
        HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Ring_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Ring_Syst_Theta2_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Ring_Theta2_Tmp.Integral()
        ring_scale = 0
        ring_scale_err = 0
        if old_integral>0 and bkg_total>0:
          ring_scale = bkg_total/old_integral
          ring_scale_err = ring_scale*(bkg_err/bkg_total)
        Theta2HistScale(Hist_Ring_Theta2_Tmp,Hist_Ring_Syst_Theta2_Tmp,ring_scale,ring_scale_err)
        Hist_Ring_Theta2.Add(Hist_Ring_Theta2_Tmp)
        Hist_Ring_Syst_Theta2.Add(Hist_Ring_Syst_Theta2_Tmp)

        bkg_total, bkg_err = IntegralAndError(Hist1D_W_Dark_SumSRs,bin_lower,bin_upper)
        dark_integral += bkg_total
        dark_integral_err += bkg_err*bkg_err
        HistName = "Hist_Data_CR_SelectFoV_Theta2_ErecS%sto%s"%(ErecS_lower_cut,ErecS_upper_cut)
        Hist_Dark_Theta2_Tmp = InputFile.Get(HistName)
        Hist_Dark_Syst_Theta2_Tmp = InputFile.Get(HistName)
        old_integral = Hist_Dark_Theta2_Tmp.Integral()
        dark_scale = 0
        dark_scale_err = 0
        if old_integral>0:
          dark_scale = bkg_total/old_integral
          dark_scale_err = dark_scale*(bkg_err/bkg_total)
        Theta2HistScale(Hist_Dark_Theta2_Tmp,Hist_Dark_Syst_Theta2_Tmp,dark_scale,dark_scale_err)
        Hist_Dark_Theta2.Add(Hist_Dark_Theta2_Tmp)
        Hist_Dark_Syst_Theta2.Add(Hist_Dark_Syst_Theta2_Tmp)

        Hist_Bkgd_Theta2.Reset()
        Hist_Bkgd_Theta2.Add(Hist_Hadr_Theta2_Tmp)
        Hist_Bkgd_Theta2.Add(Hist_Elec_Theta2_Tmp)
        hadr_syst = 0
        if hadr_scale>0: hadr_syst = hadr_scale_err/hadr_scale
        print 'hadr_syst = %s'%(hadr_syst)
        s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_Tmp,Hist_Bkgd_Theta2,theta2_lower,theta2_upper,hadr_syst)
        Hist_RDBM_S2B[len(Hist_RDBM_S2B)-1].SetBinContent(e+1,s2b)
        Hist_RDBM_S2B[len(Hist_RDBM_S2B)-1].SetBinError(e+1,s2b_err)

        Hist_Bkgd_Theta2.Reset()
        Hist_Bkgd_Theta2.Add(Hist_Dark_Theta2_Tmp)
        dark_syst = 0
        if dark_scale>0: dark_syst = dark_scale_err/dark_scale
        s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_Tmp,Hist_Bkgd_Theta2,theta2_lower,theta2_upper,dark_syst)
        Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinContent(e+1,s2b)
        Hist_Dark_S2B[len(Hist_Dark_S2B)-1].SetBinError(e+1,s2b_err)

        Hist_Bkgd_Theta2.Reset()
        Hist_Bkgd_Theta2.Add(Hist_Ring_Theta2_Tmp)
        ring_syst = 0
        if ring_scale>0: ring_syst = ring_scale_err/ring_scale
        s2b, s2b_err = Variation_ratio(Hist_Data_Theta2_Tmp,Hist_Bkgd_Theta2,theta2_lower,theta2_upper,ring_syst)
        Hist_Ring_S2B[len(Hist_Ring_S2B)-1].SetBinContent(e+1,s2b)
        Hist_Ring_S2B[len(Hist_Ring_S2B)-1].SetBinError(e+1,s2b_err)

    ErecS_lower_cut = energy_list[0]
    
    ring_integral_err = pow(ring_integral_err,0.5)
    dark_integral_err = pow(dark_integral_err,0.5)
    
    n_rebin = 1
    if theta2_upper>0.25: n_rebin = 4
    if theta2_upper>1.00: n_rebin = 16
    Hist_Hadr_Theta2.Rebin(n_rebin)
    Hist_Hadr_Syst_Theta2.Rebin(n_rebin)
    Hist_Elec_Theta2.Rebin(n_rebin)
    Hist_Ring_Theta2.Rebin(n_rebin)
    Hist_Ring_Syst_Theta2.Rebin(n_rebin)
    Hist_Dark_Theta2.Rebin(n_rebin)
    Hist_Dark_Syst_Theta2.Rebin(n_rebin)
    Hist_Data_Theta2.Rebin(n_rebin)
    
    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Hadr_Theta2]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_Elec_Theta2]
        legends += ['electron']
        colors += [3]
        Hists += [Hist_Hadr_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_RDBM_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,hadr_syst)
    
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Ring_Theta2]
        legends += ['ring']
        colors += [4]
        Hists += [Hist_Ring_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_Ring_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,ring_syst)
    
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_Theta2]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_Theta2]
        legends += ['Segue1']
        colors += [4]
        Hists += [Hist_Dark_Syst_Theta2]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_Dark_ZoomInTheta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,theta2_lower,theta2_upper,dark_syst)
    
    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_CR_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_CR_Hadr_MSCW_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_CR_Elec_MSCW_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_CRall_RDBM_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_CR_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_CR_Hadr_MSCL_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_CR_Elec_MSCL_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_CRall_RDBM_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Hadr_MSCW_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_Elec_MSCW_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_SRall_RDBM_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Hadr_MSCL_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_Elec_MSCL_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_SRall_RDBM_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Ring_MSCW_SumE]
        legends += ['Ring']
        colors += [4]
        Hists += [Hist_Ring_Syst_MSCW_SumE]
        legends += ['syst.']
        colors += [0]
        plotname = 'Target_SRall_Ring_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCW_SumE]
        legends += ['Segue1']
        colors += [4]
        plotname = 'Target_SRall_Dark_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCL_SumE]
        legends += ['Segue1']
        colors += [4]
        plotname = 'Target_SRall_Dark_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseRDBM:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_CR_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_CR_Hadr_MSCW_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_CR_Elec_MSCW_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_CRall_ZoomInRDBM_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_CR_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_CR_Hadr_MSCL_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_CR_Elec_MSCL_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_CRall_ZoomInRDBM_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Hadr_MSCW_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_Elec_MSCW_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_SRall_ZoomInRDBM_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Hadr_MSCL_SumE]
        legends += ['hadron']
        colors += [4]
        Hists += [Hist_Elec_MSCL_SumE]
        legends += ['electron (MC)']
        colors += [3]
        plotname = 'Target_SRall_ZoomInRDBM_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseRing:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Ring_MSCW_SumE]
        legends += ['Ring']
        colors += [4]
        Hists += [Hist_Ring_Syst_MSCW_SumE]
        legends += ['syst.']
        colors += [4]
        plotname = 'Target_SRall_ZoomInRing_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    if UseDark:
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCW_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCW_SumE]
        legends += ['Segue1']
        colors += [4]
        plotname = 'Target_SRall_ZoomInDark_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Data_MSCL_SumE]
        legends += ['%s'%(target)]
        colors += [1]
        Hists += [Hist_Dark_MSCL_SumE]
        legends += ['Segue1']
        colors += [4]
        plotname = 'Target_SRall_ZoomInDark_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,0)
    
    if DoSkymap and UseRDBM:
        
        
        plotname = 'Target_SRall_RDBM_Skymap_E%s'%(ErecS_lower_cut)
        Hist_Bkg_Skymap = Hist_Hadr_Skymap.Clone()
        Hist_Bkg_Skymap.Add(Hist_Elec_Skymap)
        Hist_Bkg_Syst_Skymap = Hist_Hadr_Syst_Skymap.Clone()
        Hist_Bkg_Syst_Skymap.Add(Hist_Elec_Syst_Skymap)
        Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap,Hist_Bkg_Skymap,'RA','Dec',theta2_lower,theta2_upper,plotname)
        
        smooth_size = 0.05
        Hist_Data_Skymap_smooth = Smooth2DMap(Hist_Data_Skymap,smooth_size,False)
        Hist_Bkg_Skymap_smooth = Smooth2DMap(Hist_Bkg_Skymap,smooth_size,False)
        Hist_Bkg_Syst_Skymap_smooth = Smooth2DMap(Hist_Bkg_Syst_Skymap,smooth_size,True)
        Hist_Bkg_Skymap_smooth.Add(Hist_Bkg_Syst_Skymap_smooth)
        plotname = 'Target_SRall_RDBM_Skymap_Smooth_E%s'%(ErecS_lower_cut)
        Hist_Sig = Make2DSignificancePlot(Hist_Data_Skymap_smooth,Hist_Bkg_Skymap_smooth,'RA','Dec',theta2_lower,theta2_upper,plotname)

if UseRDBM:
    MakeComparisonPlotNoRatio(Hist_RDBM_S2B,legend_S2B,color_S2B,'E [GeV]','Mismodeling_RDBM_systematics',0.6,1.4,True,False)
if UseDark:
    MakeComparisonPlotNoRatio(Hist_Dark_S2B,legend_S2B,color_S2B,'E [GeV]','Mismodeling_Dark_systematics',0.6,1.4,True,False)
if UseRing:
    MakeComparisonPlotNoRatio(Hist_Ring_S2B,legend_S2B,color_S2B,'E [GeV]','Mismodeling_Ring_systematics',0.6,1.4,True,False)
