
import sys,ROOT
import array
import math
from array import *
from ROOT import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output_Apr23'
blindness = 'Deconvolution'
#blindness = 'FitMethod'
converge = ''
#converge = '_NoConverge'
NTel = "_Ntel3to4"

e2p_source = 'DarkField'
e2p_folder = 'output_Apr23'
doTheta2 = True
doRaDec = False
doMSCW = True
doONOFF = False

source = ''
source_list = []
#source_list  += ['ComaV4']
#source_list  += ['CrabV4']
#source_list  += ['3C58']
#source_list  += ['BrandonValidation']
#source_list  += ['Segue1V6']
#source_list  += ['PKS1424']
source_list  += ['Crab']
#source_list  += ['H1426']
#source_list  += ['3C264']
#source_list  += ['Ton599']
#source_list  += ['IC443HotSpot']
#source_list  += ['DarkField']
#source_list  += ['MGRO_J1908_V5']
#source_list  += ['Segue1V5']
#source_list  += ['VA_Segue1']
#source_list  += ['VA_Geminga']

Region = 'SR'
#Region = 'VR'

Elev_lower_list = []
Elev_upper_list = []
Azim_lower_list = []
Azim_upper_list = []
Elev_lower_list += [50]
Elev_upper_list += [90]
Azim_lower_list += [0]
Azim_upper_list += [360]
#Elev_lower_list += [50]
#Elev_upper_list += [70]
#Azim_lower_list += [0]
#Azim_upper_list += [360]
#Elev_lower_list += [70]
#Elev_upper_list += [90]
#Azim_lower_list += [0]
#Azim_upper_list += [360]
#Elev_lower_list += [50]
#Elev_upper_list += [90]
#Azim_lower_list += [0]
#Azim_upper_list += [180]
#Elev_lower_list += [50]
#Elev_upper_list += [90]
#Azim_lower_list += [180]
#Azim_upper_list += [360]

SRs_included = [1,2,3,4,5,6]

#Elev_lower_cut = 20
#Elev_upper_cut = 50
Elev_lower_cut = 50
Elev_upper_cut = 70
#Elev_lower_cut = 70
#Elev_upper_cut = 90
Azim_lower_cut = 0
Azim_upper_cut = 360

Theta2_lower_cut = 0
Theta2_upper_cut = 100
mscw_cut = 'MSCWCut10'
mscw_blind = 'MSCWBlind20'

ErecS_lower_cut = 100
ErecS_upper_cut = 1e10 
MSCL_lower_cut = 0
MSCL_upper_cut = 0
MSCW_lower_cut = -1.0
MSCW_upper_cut = 1.0
MSCW_blind_cut = 1.0
exposure_hours = 0.

energy_list = []
#energy_list += [200]
#energy_list += [282]
#energy_list += [398]
#energy_list += [562]
#energy_list += [794]
energy_list += [1122]
energy_list += [1585]
#energy_list += [2239]
#energy_list += [3162]
#energy_list += [4467]
#energy_list += [6310]
#energy_list += [8913]

Variable = ''
xtitle = ''
variable_list = []
xtitle_list = []
variable_list += ['_ErecS_']
xtitle_list += ['E [GeV]']
variable_list += ['_Theta2_']
xtitle_list += ['#theta^{2}']
variable_list += ['_Theta2ZoomIn_']
xtitle_list += ['#theta^{2}']
variable_list += ['_Xoff_vs_Yoff_']
xtitle_list += ['']

def SelectDiagnosticaHistograms(folder,which_source,var,isOFF):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

    theta2_lower = Theta2_lower_cut
    theta2_upper = Theta2_upper_cut
    FilePath = '%s/%s_%s%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s_%s%s.root'%(folder,blindness,which_source,NTel,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,theta2_lower,theta2_upper,mscw_cut,mscw_blind,converge)
    #if isOFF:
    #    FilePath = '%s/%s_%s%s_Elev%sto%s_Azim%sto%s_Theta22to%s_%s_%s%s.root'%(folder,blindness,which_source,NTel,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,theta2_upper,mscw_cut,mscw_blind,converge)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    #print 'HistList.GetSize() = %s'%(HistList.GetSize())
    histogram_found = False
    for h in range(0,HistList.GetSize()):
        name = HistList.At(h).GetName()
        hist = InputFile.Get(name)
        ErecS_lower = 0
        ErecS_upper = 0
        if not var in name: continue
        for e in range(0,len(name.split("_"))):
            element = name.split("_")[e]
            if 'ErecS' in element and 'to' in element:
                element = element.strip('ErecS')
                ErecS_lower = float(element.split('to')[0])
                ErecS_upper = float(element.split('to')[1])
        #print 'ErecS_lower_cut = %s'%(ErecS_lower_cut)
        #print 'ErecS_upper_cut = %s'%(ErecS_upper_cut)
        #print 'ErecS_lower = %s'%(ErecS_lower)
        #print 'ErecS_upper = %s'%(ErecS_upper)
        if ErecS_lower<int(ErecS_lower_cut): continue
        if ErecS_upper>int(ErecS_upper_cut): continue
        #print 'found histogram'
        histogram_found = True
        if Hist_Data.Integral()==0:
            Hist_Data.Delete()
            Hist_Data = hist.Clone()
        else:
            Hist_Data.Add(hist)

    if not histogram_found: print 'Looking for %s, histogram not found!!'%(var)
    return Hist_Data

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

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
            legend.AddEntry(Hists[h],'%s, mean = %.5f, RMS = %.5f'%(legends[h],Hists[h].GetMean(),Hists[h].GetRMS()),"pl")
    legend.Draw("SAME")

    pad1.SetLogy()

    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

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

    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
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

    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeComparisonPlot(Hists,legends,colors,title,name,maxhight,logx,logy):
    
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
    else: Hists[max_hist].SetMinimum(0.1)
    if not maxhight==0: Hists[max_hist].SetMaximum(maxhight)
    Hists[max_hist].Draw("E")

    Hists[0].SetLineWidth(3)
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
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

    pad2.cd()
    Hist_Band = Hists[0].Clone()
    for b in range(0,Hist_Band.GetNbinsX()):
        Hist_Band.SetBinContent(b+1,1)
        Hist_Band.SetBinError(b+1,0.0)
    Hist_Band.SetMarkerSize(0)
    Hist_Band.GetXaxis().SetTitle(title)
    Hist_Band.GetXaxis().SetTitleOffset(1.1)
    Hist_Band.GetXaxis().SetTitleSize(0.13)
    Hist_Band.GetXaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetLabelSize(0.1)
    Hist_Band.GetYaxis().SetTitleOffset(0.3)
    Hist_Band.GetYaxis().SetTitle("ratio")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(1.5)
    Hist_Band.SetMinimum(0.5)
    if logy:
        Hist_Band.SetMaximum(2.0)
        Hist_Band.SetMinimum(0.0)
    Hist_Band.Draw("e2")
    Hist_Ratio = []
    for h in range(0,len(Hists)):
        Hist_Ratio += [Hists[h].Clone()]
    for h in range(1,len(Hists)):
        Hist_Ratio[h].Divide(Hists[0])
        if logy: 
            for b in range(0,Hist_Ratio[h].GetNbinsX()):
                Hist_Ratio[h].SetBinError(b+1,0.)
        Hist_Ratio[h].Draw("E same")

    if logx: 
        pad1.SetLogx()
        pad2.SetLogx()

    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeGaussianPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
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

    if doNorm:
        for h in range(1,len(Hists)):
            if Hists[h]!=0:
                scale = Hists[0].Integral()/Hists[h].Integral()
                Hists[h].Scale(scale)

    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,1):
        mean += [0]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            Hists[h].GetXaxis().SetRangeUser(-2.0,Hists[0].GetMean()+3.*Hists[0].GetRMS())
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")

    power2 = 2
    func = ROOT.TF1("func","[2]*exp(-0.5*pow((x-[1])/[0],%s))"%(power2), -1.0, 20.5)
    func.SetParameter(0,Hists[0].GetRMS())
    func.SetParameter(1,Hists[0].GetMean())
    func.SetParameter(2,30.)
    func.SetParName(0,"RMS")
    func.SetParName(1,"shift")
    func.SetParName(2,"Amp")

    for h in range(0,1):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    Hists[0].Fit("func","","",-1.0,20.0)
    rms[0] = func.GetParameter(0)
    mean[0] = func.GetParameter(1)
    amp[0] = func.GetParameter(2)

    Hists[0].SetLineWidth(3)
    Hists[0].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,1):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f,%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeBlindReyleighPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
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
    pad1.SetLogy()

    pad1.cd()

    if doNorm:
        for h in range(1,len(Hists)):
            if Hists[h]!=0:
                scale = Hists[0].Integral()/Hists[h].Integral()
                Hists[h].Scale(scale)

    max_heigh = 0
    max_hist = 0
    mean = []
    rms = []
    amp = []
    for h in range(0,len(Hists)):
        mean += [-1]
        rms += [0]
        amp += [0]
        if Hists[h]!=0:
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            Hists[h].GetXaxis().SetRangeUser(-2.0,Hists[0].GetMean()+3.*Hists[0].GetRMS())
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")

    power1 = 2
    power2 = 2
    func = ROOT.TF1("func","[1]*pow((x+1)/[0],%s)*exp(-0.5*pow((x+1)/[0],%s))"%(power1,power2), -1.0, 20.5)
    func.SetParameter(0,3.)
    func.SetParameter(1,100.)
    func.SetParName(0,"RMS")
    func.SetParName(1,"Amp")

    for h in range(0,2):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    for h in range(0,len(Hists)):
        Hists[h].Fit("func","","",2.0,20.0)
        rms[h] = func.GetParameter(0)
        amp[h] = func.GetParameter(1)

    Hists[0].SetLineWidth(3)
    Hists[0].Draw("E")
    func_result = []
    for h in range(0,len(Hists)):
        Hists[h].SetLineWidth(3)
        Hists[h].Draw("E same")
        func_result += [ROOT.TF1("func_result_%s"%(h),"%0.3f*pow((x+1)/%0.3f,%s)*exp(-0.5*pow((x+1)/%0.3f,%s))"%(amp[h],rms[h],power1,rms[h],power2), -1.0, 20.5)]
        func_result[h].Draw("same")

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
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f,%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/BlindFit_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeReyleighPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
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
    pad1.SetLogy()

    pad1.cd()

    if doNorm:
        for h in range(1,len(Hists)):
            if Hists[h]!=0:
                scale = Hists[0].Integral()/Hists[h].Integral()
                Hists[h].Scale(scale)

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
            Hists[h].SetLineColor(colors[h])
            Hists[h].GetXaxis().SetTitle(title)
            Hists[h].GetXaxis().SetRangeUser(-2.0,Hists[0].GetMean()+3.*Hists[0].GetRMS())
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")

    power1 = 2
    power2 = 2
    func = ROOT.TF1("func","[2]*pow((x-[1])/[0],%s)*exp(-0.5*pow((x-[1])/[0],%s))"%(power1,power2), -1.0, 20.5)
    func.SetParameter(0,3.)
    func.SetParameter(1,-0.9)
    func.SetParameter(2,100.)
    func.SetParName(0,"RMS")
    func.SetParName(1,"shift")
    func.SetParName(2,"Amp")

    for h in range(0,2):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    for h in range(0,len(Hists)):
        Hists[h].Fit("func","","",-1.0,20.0)
        rms[h] = func.GetParameter(0)
        mean[h] = func.GetParameter(1)
        amp[h] = func.GetParameter(2)

    Hists[0].SetLineWidth(3)
    Hists[0].Draw("E")
    for h in range(0,len(Hists)):
        Hists[h].SetLineWidth(3)
        Hists[h].Draw("E same")
        func_result = ROOT.TF1("func_result","[0]*pow((x-%0.3f)/%0.3f,%s)*exp(-0.5*pow((x-%0.3f)/%0.3f,%s))"%(mean[h],rms[h],power1,mean[h],rms[h],power2), -1.0, 20.5)
        #Hists[3].Fit("func_result","","",2.0,20.0)
        func_result.Draw("same")

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
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f,%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/Fit_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))


def IntegralAndError(Hist,bin1,bin2):
    
    integral = 0
    error = 0
    for b in range(1,Hist.GetNbinsX()+1):
        if b<bin1: continue
        if b>bin2: continue
        integral += Hist.GetBinContent(b)
        error += pow(Hist.GetBinError(b),2)
    error = pow(error,0.5)
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
    sum_stat_err = pow(sum_stat_err,0.5)
    error = pow(sum_syst_err*sum_syst_err+sum_stat_err*sum_stat_err,0.5)
    return integral, error

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

def RaDecHistScale(Hist,scale,scale_err):

    for bx in range(1,Hist.GetNbinsX()+1):
        for by in range(1,Hist.GetNbinsY()+1):
            old_content = Hist.GetBinContent(bx,by)
            old_error = Hist.GetBinError(bx,by)
            new_content = old_content*scale
            if not old_content==0:
                new_error = new_content*pow(pow(old_error/old_content,2)+pow(scale_err/scale,2),0.5)
                Hist.SetBinContent(bx,by,new_content)
                Hist.SetBinError(bx,by,pow(new_error*new_error+old_error*old_error,0.5))

def Theta2HistScale(Hist,scale,scale_err):

    for b in range(1,Hist.GetNbinsX()+1):
        old_content = Hist.GetBinContent(b)
        old_error = Hist.GetBinError(b)
        new_content = old_content*scale
        if not old_content==0:
            new_error = new_content*pow(pow(old_error/old_content,2)+pow(scale_err/scale,2),0.5)
            Hist.SetBinContent(b,new_content)
            Hist.SetBinError(b,pow(new_error*new_error+old_error*old_error,0.5))

def S2B_ratio(Hist_SR, Hist_Bkg,range_lower,range_upper,syst):

    norm_bin_low_target = Hist_SR.FindBin(range_lower)
    norm_bin_up_target = Hist_Bkg.FindBin(range_upper)-1
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hist_SR,norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target)
    #predict_bkg, err_bkg = IntegralAndSystError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target,syst)
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0 and not data_SR==0:
        sbratio = (data_SR-predict_bkg)/(predict_bkg)
        sbratio_err = data_SR/predict_bkg*pow(pow(err_SR/data_SR,2)+pow(err_bkg/predict_bkg,2),0.5)
        #sbratio = (data_SR-predict_bkg)/(data_SR)
        #sbratio_err = predict_bkg/data_SR*pow(pow(err_SR/data_SR,2)+pow(err_bkg/predict_bkg,2),0.5)
    return sbratio, sbratio_err

def Variation_ratio(Hist_SR, Hist_Bkg,range_lower,range_upper,syst,e2p_ratio,e2p_error):

    norm_bin_low_target = Hist_SR.FindBin(range_lower)
    norm_bin_up_target = Hist_Bkg.FindBin(range_upper)-1
    err_SR = 0
    data_SR = 0
    data_SR, err_SR = IntegralAndError(Hist_SR,norm_bin_low_target,norm_bin_up_target)
    err_bkg = 0
    predict_bkg = 0
    predict_bkg, err_bkg = IntegralAndError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target)
    #predict_bkg, err_bkg = IntegralAndSystError(Hist_Bkg,norm_bin_low_target,norm_bin_up_target,syst)
    predict_bkg = predict_bkg*(1.+e2p_ratio)
    err_scale = predict_bkg*e2p_error
    err_bkg = pow(err_bkg*err_bkg+err_scale*err_scale,0.5)
    sbratio = 0
    sbratio_err = 0
    if not predict_bkg==0 and not data_SR==0:
        sbratio = (predict_bkg)/(data_SR)
        sbratio_err = predict_bkg/data_SR*pow(pow(err_SR/data_SR,2)+pow(err_bkg/predict_bkg,2),0.5)
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
    if 'MSCW' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    elif 'MSCL' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    else:
        low_end = range_lower
        high_end = range_upper
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    Hists[max_hist].Draw("E")
    #Hists[0].SetMinimum(0)
    #Hists[0].Draw("E")

    if doSum:
        Hist_Sum = Hists[1].Clone()
        #Hist_Sum.Add(Hists[2])
        if not 'Energy' in name:
            Hist_Sum.SetMinimum(0)
        set_histStyle( Hist_Sum , 38)
        stack = ROOT.THStack("stack", "")
        stack.Add( Hist_Sum )
        stack.Draw("hist same")
        Hist_Err = Hist_Sum.Clone()
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
        if Hists[h]!=0:
            if 'NoElec' in name and legends[h]=='electron': continue
            Hists[h].Draw("E same")
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
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            if 'NoElec' in name and legends[h]=='electron': continue
            #legend.AddEntry(Hists[h],legends[h]+"(%0.2f#pm%0.2f)"%(mean[h],rms[h]),"pl")
            legend.AddEntry(Hists[h],legends[h],"pl")
    #if doSum:
    #    legend.AddEntry(Hist_Sum,'total bkg',"f")
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
    #predict_bkg, err_bkg = IntegralAndSystError(Hists[1],norm_bin_low_target,norm_bin_up_target,syst)
    lumilab2 = ROOT.TLatex(0.15,0.60,'Excess (RDBM) = %0.1f#pm%0.1f'%(data_SR-predict_bkg,pow(err_SR*err_SR+err_bkg*err_bkg,0.5)) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab3 = ROOT.TLatex(0.15,0.40,'Bkg (RDBM) = %0.1f#pm%0.1f'%(predict_bkg,err_bkg) )
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
    #Hist_Band.Divide(Hist_All)
    for b in range(0,Hist_Band.GetNbinsX()):
        Hist_Band.SetBinContent(b+1,0)
        Hist_Band.SetBinError(b+1,3)
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
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(5)
    Hist_Band.SetMinimum(-5)
    Hist_Band.GetXaxis().SetRangeUser(low_end,high_end)
    Hist_Band.Draw("e2")
    Hist_Ratio = Hists[0].Clone()
    for b in range(0,Hist_Ratio.GetNbinsX()):
        ndata = Hists[0].GetBinContent(b+1)
        edata = Hists[0].GetBinError(b+1)
        nbkg = Hist_Sum.GetBinContent(b+1)
        ebkg = Hist_Sum.GetBinError(b+1)
        if not (ebkg*ebkg+edata*edata)==0 and not (ndata-nbkg)==0:
            Hist_Ratio.SetBinContent(b+1,(ndata-nbkg)/abs(ndata-nbkg)*pow((ndata-nbkg)*(ndata-nbkg)/(ebkg*ebkg+edata*edata),0.5))
        else:
            Hist_Ratio.SetBinContent(b+1,0)
        Hist_Ratio.SetBinError(b+1,0)
    Hist_Ratio.SetFillColor(1)
    Hist_Ratio.Draw("B same")

    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeDiagnosticPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
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
            if legends[h]=='Ring':
                Hists[h].SetLineWidth(3)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    #if not 'Energy' in name:
    #    Hists[max_hist].SetMinimum(0)
    if 'MSCW' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    if 'MSCL' in name:
        low_end = Hists[0].GetMean()-3.*Hists[0].GetRMS()
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(low_end,high_end)
    Hists[max_hist].Draw("E")

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
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            if 'NoElec' in name and legends[h]=='electron': continue
            #legend.AddEntry(Hists[h],legends[h]+"(%0.2f#pm%0.2f)"%(mean[h],rms[h]),"pl")
            legend.AddEntry(Hists[h],legends[h],"pl")
    if doSum:
        legend.AddEntry(Hist_Sum,'total bkg',"f")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab5 = ROOT.TLatex(0.15,0.20,'E >%0.1f GeV'%(ErecS_lower_cut) )
    lumilab5.SetNDC()
    lumilab5.SetTextSize(0.15)
    lumilab5.Draw()
    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def BuildTheta2(Hist_2D):

    Hist_1D = ROOT.TH1D("Hist_1D","",100,0,10)
    for bx in range(1,Hist_2D.GetNbinsX()+1):
        for by in range(1,Hist_2D.GetNbinsY()+1):
            bin_content = Hist_2D.GetBinContent(bx,by)
            locationx = Hist_2D.GetXaxis().GetBinCenter(bx)
            locationy = Hist_2D.GetYaxis().GetBinCenter(by)
            theta2 = locationx*locationx+locationy*locationy
            Hist_1D.Fill(theta2,bin_content)
    return Hist_1D

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

def RooFitSmooth2DMap(Hist_Old,smooth_size):

    Hist_Smooth = Hist_Old.Clone()
    #x = ROOT.RooFit.RooRealVar("x","x",-5,5)
    #y = ROOT.RooFit.RooRealVar("y","y",-5,5)
    return Hist_Smooth

def CalculateSignificance(s,b,err):
    if (b*b+(s+b)*err*err)==0.: return 0.
    if (s+b)*(b+err*err)==0.: return 0.
    first_term = (s+b)*math.log((s+b)*(b+err*err)/(b*b+(s+b)*err*err))
    second_term = b*b/(err*err)*math.log(1.+err*err*s/(b*(b+err*err)))
    result = 0.
    if first_term>second_term: result = pow(2*(first_term-second_term),0.5)
    else: result = pow(2*(-first_term+second_term),0.5)
    if s>0: return result
    else: return -1.*result

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
    canvas.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

    Hist_Excess = Hist_SR.Clone()
    Hist_Excess.Add(Hist_Bkg,-1.)
    Hist_Excess.GetYaxis().SetTitle(ytitle)
    Hist_Excess.GetXaxis().SetTitle(xtitle)
    Hist_Excess.Draw("COL4Z")
    canvas.SaveAs('output_plots/Excess_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

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
    Hist_Model.FillRandom("func",int(Hist_Sig.GetEntries()))
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Normal']
    color_list += [2]
    Hist_list += [Hist_Sig]
    legend_list += ['Data']
    color_list += [4]
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(name))
    return Hist_Sig

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
    line1 = ROOT.TLine(MSCL_lower_cut,MSCW_blind_cut,MSCL_upper_cut,MSCW_blind_cut)
    line1.SetLineStyle(1)
    line1.SetLineColor(2)
    line1.SetLineWidth(2)
    line1.Draw()
    line2 = ROOT.TLine(MSCL_upper_cut,MSCW_blind_cut,MSCL_upper_cut,MSCW_lower_cut)
    line2.SetLineStyle(1)
    line2.SetLineColor(2)
    line2.SetLineWidth(2)
    line2.Draw()
    line3 = ROOT.TLine(MSCL_lower_cut,MSCW_lower_cut,MSCL_upper_cut,MSCW_lower_cut)
    line3.SetLineStyle(1)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw()
    if 'ErecS' in name:
        f_mean = ROOT.TF1("f_mean","5.+log2(pow(x/0.08,0.4))",100,pow(10,4))
        f_mean.SetLineColor(2)
        f_mean.SetLineStyle(2)
        f_mean.Draw("same")
        lumilab = ROOT.TLatex(0.5,0.8,'d(E) = 5+log2[(E/80MeV)^{0.4}]' )
        lumilab.SetNDC()
        lumilab.SetTextSize(0.03)
        lumilab.Draw()
        pad1.SetLogx()
    #pad1.SetLogy()
    canvas.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def Make2DTrajectoryPlot(Hist_1,Hist_2,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
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
    Hist_1.GetYaxis().SetTitle(ytitle)
    Hist_1.GetXaxis().SetTitle(xtitle)
    Hist_2.GetYaxis().SetTitle(ytitle)
    Hist_2.GetXaxis().SetTitle(xtitle)
    Hist_1.SetLineColor(4)
    Hist_1.Draw("CONT3")
    Hist_2.SetLineColor(2)
    Hist_2.Draw("CONT3 same")
    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_1,'Crab 2017-18',"pl")
    legend.AddEntry(Hist_2,source,"pl")
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/%s_%s.pdf'%(name,source))

    #canvas2 = ROOT.TCanvas("canvas2","canvas2", 200, 10, 600, 600)
    #Hist_1.RebinX(2)
    #Hist_1.RebinY(2)
    #Hist_1.Draw("COL4Z")
    #canvas2.SaveAs('output_plots/%s_%s_Lego.pdf'%(name,'Crab'))
    #Hist_2.Draw("COL4Z")
    #Hist_2.RebinX(2)
    #Hist_2.RebinY(2)
    #canvas2.SaveAs('output_plots/%s_%s_Lego.pdf'%(name,source))

Elev_lower_cut = Elev_lower_list[0]
Elev_upper_cut = Elev_upper_list[0]
Azim_lower_cut = Azim_lower_list[0]
Azim_upper_cut = Azim_upper_list[0]
FilePath = '%s/%s_%s%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s_%s%s.root'%(e2p_folder,blindness,e2p_source,NTel,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,mscw_cut,mscw_blind,converge)
TargetFile=ROOT.TFile(FilePath)
InfoTree = TargetFile.Get("InfoTree")
InfoTree.GetEntry(0)
dark_exposure_hours = InfoTree.exposure_hours

Hist_e2p = []
legend_e2p = []
color_e2p = []

Hist_ruo_predict = []
legend_ruo_predict = []
color_ruo_predict = []

Hist_dark_predict = []
legend_dark_predict = []
color_dark_predict = []

for s in range(0,len(source_list)):

    for elev in range(0,len(Elev_lower_list)):

        Elev_lower_cut = Elev_lower_list[elev]
        Elev_upper_cut = Elev_upper_list[elev]
        Azim_lower_cut = Azim_lower_list[elev]
        Azim_upper_cut = Azim_upper_list[elev]

        tele_pointing = ''
        if Elev_lower_cut>=70: tele_pointing = '(up)'
        if Elev_upper_cut<=70: tele_pointing = '(down)'
        if Azim_lower_cut>=180: tele_pointing = '(west)'
        if Azim_upper_cut<=180: tele_pointing = '(east)'
        source = source_list[s]

        ErecS_lower_cut = 0
        ErecS_upper_cut = 1e10
        FilePath = '%s/%s_%s%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s_%s%s.root'%(folder,blindness,source,NTel,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,mscw_cut,mscw_blind,converge)
        TargetFile=ROOT.TFile(FilePath)
        InfoTree = TargetFile.Get("InfoTree")
        InfoTree.GetEntry(0)
        MSCW_lower_cut = InfoTree.MSCW_cut_lower
        MSCW_upper_cut = InfoTree.MSCW_cut_upper
        MSCW_blind_cut = InfoTree.MSCW_cut_blind
        MSCL_upper_cut = InfoTree.MSCL_cut_upper
        MSCL_lower_cut = InfoTree.MSCL_cut_lower
        Theta2_upper_limit = InfoTree.Theta2_upper_limit
        Number_of_CR = InfoTree.Number_of_CR
        Number_of_SR = InfoTree.Number_of_SR
        used_runs = InfoTree.used_runs
        energy_vec = InfoTree.energy_vec
        exposure_hours = InfoTree.exposure_hours

        Hist_e2p += [ROOT.TH1D("Hist_e2p_%s"%(source),"",len(energy_list)-1,array('d',energy_list))]
        legend_e2p += ['%s %s (%0.1f hrs)'%(source,tele_pointing,exposure_hours)]
        Hist_ruo_predict += [ROOT.TH1D("Hist_ruo_predict_%s"%(source),"",len(energy_list)-1,array('d',energy_list))]
        legend_ruo_predict += ['%s %s (%0.1f hrs)'%(source,tele_pointing,exposure_hours)]
        Hist_dark_predict += [ROOT.TH1D("Hist_dark_predict_%s"%(source),"",len(energy_list)-1,array('d',energy_list))]
        legend_dark_predict += ['%s %s (%0.1f hrs)'%(source,tele_pointing,exposure_hours)]
        if s+elev<4:
            color_e2p += [s+elev+1]
            color_ruo_predict += [s+elev+1]
            color_dark_predict += [s+elev+1]
        else:
            color_e2p += [s+elev+2]
            color_ruo_predict += [s+elev+2]
            color_dark_predict += [s+elev+2]

        Hist_Target_TelElevAzim = SelectDiagnosticaHistograms(folder,source,'Target_TelElevAzim',False)
        Make2DProjectionPlot(Hist_Target_TelElevAzim,'Elev','Azim','Target_TelElevAzim',False)

        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]

            evec_match = 0
            for evec in range(0,len(energy_vec)-1):
                if energy_list[e]==energy_vec[evec]: evec_match = evec

            which_method = 'MSCW'

            if doMSCW:
                Hist_Target_SR_MSCL = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCL',False)
                Hist_Target_Deconv_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_Deconv_MSCW',False)
                Hist_Target_MSCLW = SelectDiagnosticaHistograms(folder,source,'Target_MSCLW',False)
                Make2DProjectionPlot(Hist_Target_MSCLW,'MSCL','MSCW','Target_MSCLW_E%s'%(ErecS_lower_cut),False)
                #Make2DProjectionPlot(Hist_Target_MSCLW,'MSL','MSW','Target_MSCLW_E%s'%(ErecS_lower_cut),False)

            if doMSCW:
                Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_SR1_MSCW_SumRuns',False)
                Nbins = Hist_Target_SR_MSCW.GetNbinsX()
                lower_end = Hist_Target_SR_MSCW.GetBinLowEdge(1)
                upper_end = Hist_Target_SR_MSCW.GetBinLowEdge(Nbins+1)
                Hist_Target_SRall_MSCW = ROOT.TH1D("Hist_Target_SRall_MSCW","",Nbins,lower_end,upper_end)
                Hist_Target_BkgSRall_MSCW = ROOT.TH1D("Hist_Target_BkgSRall_MSCW","",Nbins,lower_end,upper_end)
                for sr in SRs_included:
                    Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_SR%s_MSCW_SumRuns'%(sr),False)
                    Hist_Target_BkgSR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),False)
                    Hist_Target_SRall_MSCW.Add(Hist_Target_SR_MSCW)
                    Hist_Target_BkgSRall_MSCW.Add(Hist_Target_BkgSR_MSCW)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist_Target_SR_MSCW]
                    legends += ['%s %s'%(source,tele_pointing)]
                    colors += [1]
                    Hists += [Hist_Target_BkgSR_MSCW]
                    legends += ['Bkg (RDBM)']
                    colors += [4]
                    plotname = 'Target_SR%s_MSCW_E%s'%(sr,ErecS_lower_cut)
                    title = 'MSCW'
                    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_Target_SRall_MSCW]
                legends += ['%s %s'%(source,tele_pointing)]
                colors += [1]
                Hists += [Hist_Target_BkgSRall_MSCW]
                legends += ['Bkg (RDBM)']
                colors += [4]
                plotname = 'Target_SRall_MSCW_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,MSCW_blind_cut,-1)
                for sr in range(1,Number_of_CR):
                    Hist_Target_CR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_CR%s_MSCW_SumRuns'%(sr),False)
                    Hist_Target_BkgCR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgCR%s_MSCW_SumRuns'%(sr),False)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist_Target_CR_MSCW]
                    legends += ['%s %s'%(source,tele_pointing)]
                    colors += [1]
                    Hists += [Hist_Target_BkgCR_MSCW]
                    legends += ['Bkg']
                    colors += [4]
                    plotname = 'Target_CR%s_MSCW_E%s'%(sr,ErecS_lower_cut)
                    title = 'MSCW'
                    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,MSCW_lower_cut,1000,-1)
                Hists = []
                legends = []
                colors = []
                count = 0
                plotname = 'Target_MultiSR_MSCW_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                for sr in SRs_included:
                    Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_SR%s_MSCW_SumRuns'%(sr),False)
                    Hist_Target_SR_MSCW.Rebin(2)
                    Hists += [Hist_Target_SR_MSCW]
                    count += 1
                    legends += ['L-%s'%(count-1)]
                    colors += [count]
                MakeComparisonPlot(Hists,legends,colors,title,plotname,0,False,True)
                MakeReyleighPlot(Hists,legends,colors,title,plotname,False,False)
                MakeBlindReyleighPlot(Hists,legends,colors,title,plotname,False,False)
                Hists = []
                legends = []
                colors = []
                count = 0
                plotname = 'Target_MultiBkgSR_MSCW_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                for sr in SRs_included:
                    Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),False)
                    Hists += [Hist_Target_SR_MSCW]
                    count += 1
                    legends += ['L-%s'%(count-1)]
                    colors += [count]
                MakeComparisonPlot(Hists,legends,colors,title,plotname,0,False,True)

            if doONOFF:
                temp_theta2_lower = Theta2_lower_cut
                temp_theta2_upper = Theta2_upper_cut
                Theta2_lower_cut = 0
                Theta2_upper_cut = 100
                Hist_Target_SR_MSCW_ON = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCW_SumRuns_SumSRs',False)
                Hist_Target_BkgSR_MSCW_ON = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR_MSCW_SumRuns_SumSRs',False)
                Theta2_lower_cut = 2
                Theta2_upper_cut = 100
                Hist_Target_SR_MSCW_OFF = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCW_SumRuns_SumSRs',False)
                Hist_Target_BkgSR_MSCW_OFF = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR_MSCW_SumRuns_SumSRs',False)
                Theta2_lower_cut = temp_theta2_lower 
                Theta2_upper_cut = temp_theta2_upper 
                Hists = []
                legends = []
                colors = []
                plotname = 'Target_ONOFF_MSCW_E%s'%(ErecS_lower_cut)
                title = 'MSCW'
                Hists += [Hist_Target_SR_MSCW_ON]
                legends += ['SR ON']
                colors += [1]
                Hists += [Hist_Target_BkgSR_MSCW_ON]
                legends += ['Bkg ON']
                colors += [2]
                Hists += [Hist_Target_SR_MSCW_OFF]
                legends += ['SR OFF']
                colors += [3]
                Hists += [Hist_Target_BkgSR_MSCW_OFF]
                legends += ['Bkg OFF']
                colors += [4]
                MakeComparisonPlot(Hists,legends,colors,title,plotname,0,False,False)

            #for sr in SRs_included:
            #    Hist_Target_SR_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_SR%s_theta2'%(sr),False)
            #    Hist_Target_Bkg_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_CR%s_theta2'%(sr),False)
            #    Hist_Target_BkgSR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),False)
            #    scalesky, scalesky_err = GetSkyScaleFactor(Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,Hist_Target_Bkg_theta2)
            #    e2p_file = open("e2p_ratio_%s.txt"%(mscw_cut),"read")
            #    energy_string = '%s-%s'%(ErecS_lower_cut,ErecS_upper_cut)
            #    print energy_string
            #    e2p_ratio = 0
            #    e2p_error = 0
            #    for line in e2p_file:
            #        if energy_string in line:
            #            e2p_ratio = float(line.split()[4])
            #            e2p_error = float(line.split()[6])
            #    print e2p_ratio
            #    print e2p_error
            #    ideal_nbins = float(Hist_Target_SR_theta2.Integral())/10.
            #    n_merge = 1
            #    while Hist_Target_SR_theta2.GetNbinsX()>ideal_nbins and n_merge<8:
            #        Hist_Target_Bkg_theta2.Rebin(2)
            #        Hist_Target_SR_theta2.Rebin(2)
            #        n_merge = n_merge*2
            #    Theta2HistScale(Hist_Target_Bkg_theta2,scalesky, scalesky_err)
            #    Theta2HistScale(Hist_Target_Bkg_theta2,e2p_ratio+1.,e2p_error)
            #    Hists = []
            #    legends = []
            #    colors = []
            #    Hists += [Hist_Target_SR_theta2]
            #    legends += ['%s %s'%(source,tele_pointing)]
            #    colors += [1]
            #    Hists += [Hist_Target_Bkg_theta2]
            #    legends += ['Bkg']
            #    colors += [4]
            #    plotname = 'Target_SR%s_fine_theta2_E%s'%(sr,ErecS_lower_cut)
            #    title = 'theta2'
            #    syst = (scalesky*e2p_error+(e2p_ratio+1.)*scalesky_err)/(scalesky*(e2p_ratio+1.))
            #    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,0,0.5,syst)

            if doTheta2:
                avg_events_per_bin = 1000.
                Hist_Target_SRall_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCW_SumRuns_SumSRs',False)
                Hist_Target_SR_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_SR1_theta2',False)
                Nbins = Hist_Target_SR_theta2.GetNbinsX()
                lower_end = Hist_Target_SR_theta2.GetBinLowEdge(1)
                upper_end = Hist_Target_SR_theta2.GetBinLowEdge(Nbins+1)
                Hist_Target_SRall_theta2 = ROOT.TH1D("Hist_Target_SRall_theta2","",Nbins,lower_end,upper_end)
                Hist_Target_BkgSRall_theta2 = ROOT.TH1D("Hist_Target_BkgSRall_theta2","",Nbins,lower_end,upper_end)
                Hist_Target_LZABkgSRall_theta2 = ROOT.TH1D("Hist_Target_LZABkgSRall_theta2","",Nbins,lower_end,upper_end)
                Hist_Target_SZABkgSRall_theta2 = ROOT.TH1D("Hist_Target_SZABkgSRall_theta2","",Nbins,lower_end,upper_end)
                ideal_nbins = float(Hist_Target_SRall_MSCW.Integral())/avg_events_per_bin
                n_merge = 1
                while (Hist_Target_SRall_theta2.GetNbinsX()>ideal_nbins and n_merge<64) or (n_merge<8):
                    Hist_Target_BkgSRall_theta2.Rebin(2)
                    Hist_Target_LZABkgSRall_theta2.Rebin(2)
                    Hist_Target_SZABkgSRall_theta2.Rebin(2)
                    Hist_Target_SRall_theta2.Rebin(2)
                    n_merge = n_merge*2
                for sr in SRs_included:
                    Hist_Target_SR_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_SR%s_theta2'%(sr),False)
                    Hist_Target_Bkg_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_CR%s_theta2'%(sr),False)
                    Hist_Target_LZABkg_theta2 = SelectDiagnosticaHistograms(folder,source,'TargetLZA_CR%s_theta2'%(sr),False)
                    Hist_Target_SZABkg_theta2 = SelectDiagnosticaHistograms(folder,source,'TargetSZA_CR%s_theta2'%(sr),False)
                    Hist_Target_BkgSR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),False)
                    scalesky, scalesky_err = GetSkyScaleFactor(Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,Hist_Target_Bkg_theta2)
                    scalesky_LZA, scalesky_err_LZA = GetSkyScaleFactor(Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,Hist_Target_LZABkg_theta2)
                    scalesky_SZA, scalesky_err_SZA = GetSkyScaleFactor(Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,Hist_Target_SZABkg_theta2)
                    Hist_e2p_SR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_SR%s_MSCW_SumRuns'%(sr),True)
                    Hist_e2p_BkgSR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),True)
                    e2p_ratio = 0
                    e2p_error = 0
                    e2p_ratio, e2p_error = S2B_ratio(Hist_e2p_SR_MSCW,Hist_e2p_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,-1)
                    if e2p_ratio<0.:
                        e2p_ratio = 0
                        e2p_error = 0
                    #if 1.0*e2p_error>=e2p_ratio:
                    #    e2p_ratio = 0
                    #    e2p_error = 0

                    ideal_nbins = float(Hist_Target_SRall_MSCW.Integral())/avg_events_per_bin
                    n_merge = 1
                    while (Hist_Target_SR_theta2.GetNbinsX()>ideal_nbins and n_merge<64) or (n_merge<8):
                        Hist_Target_Bkg_theta2.Rebin(2)
                        Hist_Target_LZABkg_theta2.Rebin(2)
                        Hist_Target_SZABkg_theta2.Rebin(2)
                        Hist_Target_SR_theta2.Rebin(2)
                        n_merge = n_merge*2
                    Theta2HistScale(Hist_Target_Bkg_theta2,scalesky,scalesky_err)
                    Theta2HistScale(Hist_Target_Bkg_theta2,e2p_ratio+1.,e2p_error)
                    Theta2HistScale(Hist_Target_LZABkg_theta2,scalesky_LZA,scalesky_err_LZA)
                    Theta2HistScale(Hist_Target_LZABkg_theta2,e2p_ratio+1.,e2p_error)
                    Theta2HistScale(Hist_Target_SZABkg_theta2,scalesky_SZA,scalesky_err_SZA)
                    Theta2HistScale(Hist_Target_SZABkg_theta2,e2p_ratio+1.,e2p_error)
                    Hist_Target_SRall_theta2.Add(Hist_Target_SR_theta2)
                    Hist_Target_BkgSRall_theta2.Add(Hist_Target_Bkg_theta2)
                    Hist_Target_LZABkgSRall_theta2.Add(Hist_Target_LZABkg_theta2)
                    Hist_Target_SZABkgSRall_theta2.Add(Hist_Target_SZABkg_theta2)
                    Hists = []
                    legends = []
                    colors = []
                    Hists += [Hist_Target_SR_theta2]
                    legends += ['%s %s'%(source,tele_pointing)]
                    colors += [1]
                    Hists += [Hist_Target_Bkg_theta2]
                    legends += ['Bkg']
                    colors += [4]
                    plotname = 'Target_SR%s_theta2_E%s'%(sr,ErecS_lower_cut)
                    title = 'theta2'
                    syst = (scalesky*e2p_error+(e2p_ratio+1.)*scalesky_err)/(scalesky*(e2p_ratio+1.))
                    MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,0,Theta2_upper_limit,syst)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_Target_SRall_theta2]
                legends += ['%s %s'%(source,tele_pointing)]
                colors += [1]
                Hists += [Hist_Target_BkgSRall_theta2]
                legends += ['Bkg']
                colors += [4]
                plotname = 'Target_SRall_theta2_E%s'%(ErecS_lower_cut)
                title = 'theta2'
                syst = 0
                MakeChi2Plot(Hists,legends,colors,title,plotname,True,False,0,Theta2_upper_limit,syst)
                Hists = []
                legends = []
                colors = []
                Hists += [Hist_Target_SRall_theta2]
                legends += ['data']
                colors += [1]
                Hists += [Hist_Target_BkgSRall_theta2]
                legends += ['elev. 50-90']
                colors += [2]
                Hists += [Hist_Target_LZABkgSRall_theta2]
                legends += ['elev. 50-70']
                colors += [3]
                Hists += [Hist_Target_SZABkgSRall_theta2]
                legends += ['elev. 70-90']
                colors += [4]
                plotname = 'Target_diff_theta2_E%s'%(ErecS_lower_cut)
                title = 'theta2'
                MakeComparisonPlot(Hists,legends,colors,title,plotname,0,False,False)


        for run in range(0,len(used_runs)):
            print 'using run %s'%(used_runs[run])
        print '%s runs.'%(len(used_runs))
        print 'exposure_hours = %s'%(exposure_hours)

        for e in range(0,len(energy_list)-1):
            ErecS_lower_cut = energy_list[e]
            ErecS_upper_cut = energy_list[e+1]
            Hist_Dark_SR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_SR_MSCW_SumRuns_SumSRs',True)
            Hist_Dark_SR_MSCW.Scale(exposure_hours/dark_exposure_hours)
            Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCW_SumRuns_SumSRs',False)
            Hist_Target_BkgSR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR_MSCW_SumRuns_SumSRs',False)
            s2b = 0.
            s2b_err = 0.
            s2b, s2b_err = S2B_ratio(Hist_Target_SR_MSCW, Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,-1)
            Hist_e2p_SR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_SR_MSCW_SumRuns_SumSRs',True)
            Hist_e2p_BkgSR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_BkgSR_MSCW_SumRuns_SumSRs',True)
            e2p_ratio = 0
            e2p_error = 0
            e2p_ratio, e2p_error = S2B_ratio(Hist_e2p_SR_MSCW,Hist_e2p_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,-1)
            if e2p_ratio<0.:
                e2p_ratio = 0
                e2p_error = 0
            #if 1.0*e2p_error>=e2p_ratio:
            #    e2p_ratio = 0
            #    e2p_error = 0
            r2b, r2b_err = Variation_ratio(Hist_Target_SR_MSCW, Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,0,e2p_ratio,e2p_error)
            d2b, d2b_err = Variation_ratio(Hist_Target_SR_MSCW, Hist_Dark_SR_MSCW,MSCW_lower_cut,MSCW_upper_cut,0,0,0)
            print 'E %s-%s, S/B = %0.3f +/- %0.3f'%(ErecS_lower_cut,ErecS_upper_cut,s2b,s2b_err)
            Hist_e2p[len(Hist_e2p)-1].SetBinContent(e+1,s2b)
            Hist_e2p[len(Hist_e2p)-1].SetBinError(e+1,s2b_err)
            Hist_ruo_predict[len(Hist_ruo_predict)-1].SetBinContent(e+1,r2b)
            Hist_ruo_predict[len(Hist_ruo_predict)-1].SetBinError(e+1,r2b_err)
            Hist_dark_predict[len(Hist_dark_predict)-1].SetBinContent(e+1,d2b)
            Hist_dark_predict[len(Hist_dark_predict)-1].SetBinError(e+1,d2b_err)

MakeComparisonPlot(Hist_e2p,legend_e2p,color_e2p,'E [GeV]','Target_e2p_ratio',0.5,True,False)
MakeComparisonPlotNoRatio(Hist_ruo_predict,legend_ruo_predict,color_ruo_predict,'E [GeV]','Ruo_predict_ratio',0.6,1.4,True,False)
MakeComparisonPlotNoRatio(Hist_dark_predict,legend_dark_predict,color_dark_predict,'E [GeV]','Dark_predict_ratio',0.6,1.4,True,False)

for e in range(0,len(energy_list)-1):
    Hist_mscw = []
    legend_mscw = []
    color_mscw = []
    ErecS_lower_cut = energy_list[e]
    ErecS_upper_cut = energy_list[e+1]
    for s in range(0,len(source_list)):
        for elev in range(0,len(Elev_lower_list)):
            Elev_lower_cut = Elev_lower_list[elev]
            Elev_upper_cut = Elev_upper_list[elev]
            Azim_lower_cut = Azim_lower_list[elev]
            Azim_upper_cut = Azim_upper_list[elev]
            tele_pointing = ''
            if Elev_lower_cut>=70: tele_pointing = '(up)'
            if Elev_upper_cut<=70: tele_pointing = '(down)'
            if Azim_lower_cut>=180: tele_pointing = '(west)'
            if Azim_upper_cut<=180: tele_pointing = '(east)'
            source = source_list[s]
            legend_mscw += ['%s %s'%(source,tele_pointing)]
            if s+elev<4:
                color_mscw += [s+elev+1]
            else:
                color_mscw += [s+elev+2]
            FilePath = '%s/%s_%s%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s_%s%s.root'%(folder,blindness,source,NTel,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,mscw_cut,mscw_blind,converge)
            TargetFile=ROOT.TFile(FilePath)
            InfoTree = TargetFile.Get("InfoTree")
            InfoTree.GetEntry(0)
            MSCW_lower_cut = InfoTree.MSCW_cut_lower
            MSCW_upper_cut = InfoTree.MSCW_cut_upper
            MSCW_blind_cut = InfoTree.MSCW_cut_blind
            MSCL_upper_cut = InfoTree.MSCL_cut_upper
            MSCL_lower_cut = InfoTree.MSCL_cut_lower
            Theta2_upper_limit = InfoTree.Theta2_upper_limit
            Number_of_CR = InfoTree.Number_of_CR
            Number_of_SR = InfoTree.Number_of_SR
            used_runs = InfoTree.used_runs
            energy_vec = InfoTree.energy_vec
            exposure_hours = InfoTree.exposure_hours

            Hist_Target_SR_MSCW_Combined = SelectDiagnosticaHistograms(folder,source,'Target_SR_MSCW_SumRuns_SumSRs',False)
            Hist_Target_SR_MSCW_Combined.Scale(1./exposure_hours)
            Hist_mscw += [Hist_Target_SR_MSCW_Combined]
    plotname = 'Target_MultiSource_MSCW_E%s'%(ErecS_lower_cut)
    MakeComparisonPlot(Hist_mscw,legend_mscw,color_mscw,'MSCW',plotname,0,False,True)

if doRaDec:
    ErecS_lower_cut = energy_list[0]
    ErecS_upper_cut = energy_list[1]
    Hist_Target_SRall_RaDec = SelectDiagnosticaHistograms(folder,source,'Target_SR1_RaDec',False)
    Hist_Target_BkgSRall_RaDec = SelectDiagnosticaHistograms(folder,source,'Target_CR1_RaDec',False)
    Hist_Sig_all = ROOT.TH1D("Hist_Sig_all","",65,-5,8)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    for s in range(0,len(source_list)):
        source = source_list[s]
        for elev in range(0,len(Elev_lower_list)):
            Hist_Target_SRall_RaDec.Reset()
            Hist_Target_BkgSRall_RaDec.Reset()
            Elev_lower_cut = Elev_lower_list[elev]
            Elev_upper_cut = Elev_upper_list[elev]
            Azim_lower_cut = Azim_lower_list[elev]
            Azim_upper_cut = Azim_upper_list[elev]
            tele_pointing = ''
            if Elev_lower_cut>=70: tele_pointing = '(up)'
            if Elev_upper_cut<=70: tele_pointing = '(down)'
            if Azim_lower_cut>=180: tele_pointing = '(west)'
            if Azim_upper_cut<=180: tele_pointing = '(east)'
            smooth_size = 0.1
            ideal_nbins = 32.*(2./smooth_size)
            plotname = 'Target_SR_RaDec'
            while Hist_Target_SRall_RaDec.GetNbinsX()>ideal_nbins:
                Hist_Target_SRall_RaDec.Rebin2D(2,2)
                Hist_Target_BkgSRall_RaDec.Rebin2D(2,2)
            for e in range(0,len(energy_list)-1):
                ErecS_lower_cut = energy_list[e]
                ErecS_upper_cut = energy_list[e+1]
                for sr in SRs_included:
                    Hist_Target_Bkg_theta2 = SelectDiagnosticaHistograms(folder,source,'Target_CR%s_theta2'%(sr),False)
                    Hist_Target_BkgSR_MSCW = SelectDiagnosticaHistograms(folder,source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),False)
                    scalesky, scalesky_err = GetSkyScaleFactor(Hist_Target_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,Hist_Target_Bkg_theta2)
                    Hist_Target_SR_RaDec = SelectDiagnosticaHistograms(folder,source,'Target_SR%s_RaDec'%(sr),False)
                    Hist_Target_Bkg_RaDec = SelectDiagnosticaHistograms(folder,source,'Target_CR%s_RaDec'%(sr),False)
                    while Hist_Target_SR_RaDec.GetNbinsX()>ideal_nbins:
                        Hist_Target_SR_RaDec.Rebin2D(2,2)
                        Hist_Target_Bkg_RaDec.Rebin2D(2,2)
                    Hist_e2p_SR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_SR%s_MSCW_SumRuns'%(sr),True)
                    Hist_e2p_BkgSR_MSCW = SelectDiagnosticaHistograms(e2p_folder,e2p_source,'Target_BkgSR%s_MSCW_SumRuns'%(sr),True)
                    e2p_ratio = 0
                    e2p_error = 0
                    e2p_ratio, e2p_error = S2B_ratio(Hist_e2p_SR_MSCW,Hist_e2p_BkgSR_MSCW,MSCW_lower_cut,MSCW_upper_cut,-1)
                    if e2p_ratio<0.:
                        e2p_ratio = 0
                        e2p_error = 0
                    #if 1.0*e2p_error>=e2p_ratio:
                    #    e2p_ratio = 0
                    #    e2p_error = 0
                    RaDecHistScale(Hist_Target_Bkg_RaDec,scalesky,scalesky_err)
                    RaDecHistScale(Hist_Target_Bkg_RaDec,e2p_ratio+1.,e2p_error)
                    Hist_Target_SRall_RaDec.Add(Hist_Target_SR_RaDec)
                    Hist_Target_BkgSRall_RaDec.Add(Hist_Target_Bkg_RaDec)
            Hist_Target_SRall_RaDec_smooth = Smooth2DMap(Hist_Target_SRall_RaDec,smooth_size,False)
            Hist_Target_BkgSRall_RaDec_smooth = Smooth2DMap(Hist_Target_BkgSRall_RaDec,smooth_size,False)
            Hist_Sig = Make2DSignificancePlot(Hist_Target_SRall_RaDec_smooth,Hist_Target_BkgSRall_RaDec_smooth,'RA','Dec',plotname)
            Hist_Sig_all.Add(Hist_Sig)
    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Model = ROOT.TH1D("Hist_Model","",65,-5,8)
    Hist_Model.FillRandom("func",int(Hist_Sig_all.GetEntries()))
    Hist_list = []
    legend_list = []
    color_list = []
    Hist_list += [Hist_Model]
    legend_list += ['Normal']
    color_list += [2]
    Hist_list += [Hist_Sig_all]
    legend_list += ['Data']
    color_list += [4]
    plotname = 'Mmultisource_SR_RaDec'
    MakeGaussComparisonPlot(Hist_list,legend_list,color_list,'significance','SigDist_%s'%(plotname))

