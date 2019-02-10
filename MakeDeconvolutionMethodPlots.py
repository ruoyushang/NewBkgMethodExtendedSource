
import sys,ROOT
import array
import math
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'

source = ''
source_list = []
#source_list  += ['ComaV4']
#source_list  += ['3C58']
#source_list  += ['BrandonValidation']
#source_list  += ['PKS1424']
#source_list  += ['Crab']
#source_list  += ['CrabV4']
#source_list  += ['H1426']
source_list  += ['3C264']
#source_list  += ['Ton599']
#source_list  += ['IC443']

Region = 'SR'

global MSCW_lower_cut
global MSCW_upper_cut
global MSCL_lower_cut
global MSCL_upper_cut
global Depth_cut_width

Elev_lower_cut = 60
Elev_upper_cut = 70
Azim_lower_cut = 100
Azim_upper_cut = 150

Theta2_lower_cut = 0
Theta2_upper_cut = 10

UseMethod1 = False
UseMethod2 = True
UseMethod3 = True

Depth_cut_width = 10
ErecS_lower_cut = 100
ErecS_upper_cut = 1e10 
MSCL_lower_cut = -0.5
MSCL_upper_cut = 0.6
MSCW_lower_cut = -0.5
MSCW_upper_cut = 0.6

energy_list = []
#energy_list += [200]
energy_list += [400]
energy_list += [500]
energy_list += [600]
energy_list += [700]
energy_list += [800]
energy_list += [1000]
energy_list += [1200]
energy_list += [1500]
energy_list += [2000]
energy_list += [3000]
energy_list += [5000]
energy_list += [1e10]

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

def SelectDiagnosticaHistograms(folder,method,isSR,var):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

    FilePath = '%s/Deconvolution_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region)
    #print 'Read %s'%(FilePath)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    #print 'HistList.GetSize() = %s'%(HistList.GetSize())
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
        if Hist_Data.Integral()==0:
            Hist_Data.Delete()
            Hist_Data = hist.Clone()
        else:
            Hist_Data.Add(hist)

    return Hist_Data

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def MakeGaussianPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
    global MSCW_lower_cut
    global MSCW_upper_cut
    global MSCL_lower_cut
    global MSCL_upper_cut
    global Depth_cut_width

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
    func.SetParameter(0,10.)
    func.SetParameter(1,0.)
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

def MakeReyleighPlot(Hists,legends,colors,title,name,doSum,doNorm):
    
    global MSCW_lower_cut
    global MSCW_upper_cut
    global MSCL_lower_cut
    global MSCL_upper_cut
    global Depth_cut_width

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
    for h in range(0,2):
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

    power1 = 2
    power2 = 2
    func = ROOT.TF1("func","[2]*pow(x-[1],%s)*exp(-0.5*pow((x-[1])/[0],%s))"%(power1,power2), -1.0, 20.5)
    func.SetParameter(0,3.)
    func.SetParameter(1,0.)
    func.SetParameter(2,100.)
    func.SetParName(0,"RMS")
    func.SetParName(1,"shift")
    func.SetParName(2,"Amp")

    for h in range(0,2):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    Hists[0].Fit("func","","",-1.0,20.0)
    rms[0] = func.GetParameter(0)
    mean[0] = func.GetParameter(1)
    amp[0] = func.GetParameter(2)
    Hists[1].Fit("func","","",-1.0,20.0)
    rms[1] = func.GetParameter(0)
    mean[1] = func.GetParameter(1)
    amp[1] = func.GetParameter(2)

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
    for h in range(0,2):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f,%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/Aux_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

    pad1.cd()
    max_heigh = 0
    max_hist = 0
    for h in range(2,4):
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

    func = ROOT.TF1("func","[2]*pow(x-[1],%s)*exp(-0.5*pow((x-[1])/[0],%s))"%(power1,power2), -1.0, 20.5)
    func.SetParameter(0,3.)
    func.SetParameter(1,0.)
    func.SetParameter(2,100.)
    func.SetParName(0,"RMS")
    func.SetParName(1,"shift")
    func.SetParName(2,"Amp")

    for h in range(2,4):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    Hists[2].Fit("func","","",-1.0,20.0)
    rms[2] = func.GetParameter(0)
    mean[2] = func.GetParameter(1)
    amp[2] = func.GetParameter(2)

    rms[3] = rms[1]-rms[0]+rms[2]
    mean[3] = mean[1]-mean[0]+mean[2]
    amp[3] = amp[1]/amp[0]*amp[2]
    func_result = ROOT.TF1("func_result","[0]*pow(x-%0.3f,%s)*exp(-0.5*pow((x-%0.3f)/%0.3f,%s))"%(mean[3],power1,mean[3],rms[3],power2), -1.0, 20.5)
    Hists[3].Fit("func_result","","",2.0,20.0)
    func_result.Draw("same")

    Hists[3].SetLineWidth(3)
    Hists[3].Draw("E same")

    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.15)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(2,4):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f,pm%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/Sig_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def MakeChi2Plot(Hists,legends,colors,title,name,doSum,doNorm):
    
    global MSCW_lower_cut
    global MSCW_upper_cut
    global MSCL_lower_cut
    global MSCL_upper_cut
    global Depth_cut_width

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
    #pad1.SetLogy()

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
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(-2.0,high_end)
    if 'MSCL' in name:
        high_end = Hists[0].GetMean()+3.*Hists[0].GetRMS()
        Hists[max_hist].GetXaxis().SetRangeUser(-2.0,high_end)
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

    Hist_Excess = Hists[0].Clone()
    Hist_Excess.Add(Hists[1],-1.)
    Hist_Excess.SetFillColor(0)
    Hist_Excess.SetFillStyle(0)
    Hist_Excess.SetLineColor(2)
    Hist_Excess.Draw("e2 same")

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
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    #if not title=='MSCW':
    #    lumilab2 = ROOT.TLatex(0.15,0.65,'%0.2f < MSCW < %0.2f'%(MSCW_lower_cut,MSCW_upper_cut) )
    #    lumilab2.SetNDC()
    #    lumilab2.SetTextSize(0.15)
    #    lumilab2.Draw()
    #if not title=='MSCL':
    #    lumilab5 = ROOT.TLatex(0.15,0.50,'%0.2f < MSCL < %0.2f'%(MSCL_lower_cut,MSCL_upper_cut) )
    #    lumilab5.SetNDC()
    #    lumilab5.SetTextSize(0.15)
    #    lumilab5.Draw()
    #if not title=='Depth':
    #    lumilab4 = ROOT.TLatex(0.15,0.35,'Slant depth (central #pm %0.2f)'%(Depth_cut_width) )
    #    lumilab4.SetNDC()
    #    lumilab4.SetTextSize(0.15)
    #    lumilab4.Draw()
    lumilab5 = ROOT.TLatex(0.15,0.20,'E >%0.1f GeV'%(ErecS_lower_cut) )
    lumilab5.SetNDC()
    lumilab5.SetTextSize(0.15)
    lumilab5.Draw()

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
    Hist_Band.SetMinimum(0)
    if 'MSCW' in name:
        Hist_Band.GetXaxis().SetRangeUser(-2.0,high_end)
    if 'MSCL' in name:
        Hist_Band.GetXaxis().SetRangeUser(-2.0,high_end)
    Hist_Band.Draw("e2")
    Hist_Ratio = Hists[0].Clone()
    for b in range(0,Hist_Ratio.GetNbinsX()):
        ndata = Hists[0].GetBinContent(b+1)
        edata = Hists[0].GetBinError(b+1)
        nbkg = Hist_Sum.GetBinContent(b+1)
        ebkg = Hist_Sum.GetBinError(b+1)
        if not (ebkg*ebkg+edata*edata)==0:
            Hist_Ratio.SetBinContent(b+1,pow((ndata-nbkg)*(ndata-nbkg)/(ebkg*ebkg+edata*edata),0.5))
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
    
    global MSCW_lower_cut
    global MSCW_upper_cut
    global MSCL_lower_cut
    global MSCL_upper_cut
    global Depth_cut_width

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

    if not 'Energy' in name:
        Hists[max_hist].SetMinimum(0)
    if 'MSCW' in name:
        Hists[max_hist].GetXaxis().SetRangeUser(-2.0,Hists[0].GetMean()+3.*Hists[0].GetRMS())
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
    if doSum:
        legend.AddEntry(Hist_Sum,'total bkg',"f")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    #if not title=='MSCW':
    #    lumilab2 = ROOT.TLatex(0.15,0.65,'%0.2f < MSCW < %0.2f'%(MSCW_lower_cut,MSCW_upper_cut) )
    #    lumilab2.SetNDC()
    #    lumilab2.SetTextSize(0.15)
    #    lumilab2.Draw()
    #if not title=='MSCL':
    #    lumilab5 = ROOT.TLatex(0.15,0.50,'%0.2f < MSCL < %0.2f'%(MSCL_lower_cut,MSCL_upper_cut) )
    #    lumilab5.SetNDC()
    #    lumilab5.SetTextSize(0.15)
    #    lumilab5.Draw()
    #if not title=='Depth':
    #    lumilab4 = ROOT.TLatex(0.15,0.35,'Slant depth (central #pm %0.2f)'%(Depth_cut_width) )
    #    lumilab4.SetNDC()
    #    lumilab4.SetTextSize(0.15)
    #    lumilab4.Draw()
    lumilab5 = ROOT.TLatex(0.15,0.20,'E >%0.1f GeV'%(ErecS_lower_cut) )
    lumilab5.SetNDC()
    lumilab5.SetTextSize(0.15)
    lumilab5.Draw()
    if 'Energy' in name:
        pad1.SetLogy()
        pad1.SetLogx()
    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def Smooth2DMap(Hist_Old,smooth_size):

    Hist_Smooth = Hist_Old.Clone()
    bin_size = Hist_Old.GetXaxis().GetBinCenter(2)-Hist_Old.GetXaxis().GetBinCenter(1)
    nbin_smooth = int(smooth_size/bin_size) + 1
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
                            if distance<smooth_size:
                                bin_content += Hist_Old.GetBinContent(bx2,by2)
                                bin_error += pow(Hist_Old.GetBinError(bx2,by2),2)
            Hist_Smooth.SetBinContent(bx1,by1,bin_content)
            Hist_Smooth.SetBinError(bx1,by1,pow(bin_error,0.5))
    return Hist_Smooth

def CalculateSignificance(s,b,err):
    if (s+b)*(b+err*err)/(b*b+(s+b)*err*err)==0.: return 0.
    first_term = (s+b)*math.log((s+b)*(b+err*err)/(b*b+(s+b)*err*err))
    second_term = b*b/(err*err)*math.log(1.+err*err*s/(b*(b+err*err)))
    result = pow(2*(first_term-second_term),0.5)
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
            #Sig = (NSR-NBkg)/(pow(NSR_Err*NSR_Err+NBkg_Err*NBkg_Err,0.5))
            Sig = CalculateSignificance(NSR-NBkg,NBkg,NBkg_Err)
            Hist_Data.SetBinContent(bx+1,by+1,Sig)
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.SetMaximum(5)
    Hist_Data.Draw("COL4Z")
    canvas.SaveAs('output_plots/%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

    func = ROOT.TF1("func","gaus", -5, 8)
    func.SetParameters(10.,0.,1.0)
    Hist_Sig = ROOT.TH1D("Hist_Sig","",65,-5,8)
    Hist_Sig_cut = ROOT.TH1D("Hist_Sig_cut","",65,-5,8)
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0:
                content = Hist_Data.GetBinContent(bx+1,by+1)
                Hist_Sig.Fill(content)
                locationx1 = Hist_Data.GetXaxis().GetBinCenter(bx)
                locationy1 = Hist_Data.GetYaxis().GetBinCenter(by)
                if (pow(locationx1*locationx1+locationy1*locationy1,0.5)>1.0):
                    Hist_Sig_cut.Fill(content)
    Hist_Model = ROOT.TH1D("Hist_Model","",65,-5,8)
    Hist_Model.FillRandom("func",int(Hist_Sig.GetEntries()))
    pad1.SetLogy()
    Hist_Sig.Draw()
    Hist_Sig_cut.SetLineColor(3)
    Hist_Sig_cut.Draw("same")
    Hist_Model.SetLineColor(2)
    Hist_Model.Draw("same")
    legend = ROOT.TLegend(0.6,0.65,0.9,0.85)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    legend.AddEntry(Hist_Model,'true gaus',"pl")
    legend.AddEntry(Hist_Sig,'#theta^{2} #in [0,10]',"pl")
    legend.AddEntry(Hist_Sig_cut,'#theta^{2} #in [1,10]',"pl")
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/Sig_%s_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region))

def Make2DProjectionPlot(Hist_Data,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.SetGrid()
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
    #Hist_1D.Draw("E same")
    #pad1.SetLogz()
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
    if source=='2ndCrab':
        legend.AddEntry(Hist_2,'Crab 2016',"pl")
    else:
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


for s in source_list:

    source = s
    ErecS_lower_cut = 0
    ErecS_upper_cut = 1e10
    FilePath = '%s/Deconvolution_%s_Elev%sto%s_Azim%sto%s_Theta2%sto%s_%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut,Azim_lower_cut,Azim_upper_cut,Theta2_lower_cut,Theta2_upper_cut,Region)
    #print 'Read %s'%(FilePath)
    TargetFile=ROOT.TFile(FilePath)
    InfoTree = TargetFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    MSCW_lower_cut = InfoTree.MSCW_cut_lower
    MSCW_upper_cut = InfoTree.MSCW_cut_upper
    MSCL_lower_cut = InfoTree.MSCL_signal_cut_lower
    MSCL_upper_cut = InfoTree.MSCL_signal_cut_upper

    #Hist_Target_SR_Energy = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_Energy')
    #Hist_Target_Bkg_Energy = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_Energy')
    #Hists = []
    #legends = []
    #colors = []
    #Hists += [Hist_Target_SR_Energy]
    #legends += ['SR']
    #colors += [4]
    #Hists += [Hist_Target_Bkg_Energy]
    #legends += ['Bkg']
    #colors += [2]
    #plotname = 'Target_Energy'
    #title = 'E [GeV]'
    #MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)

    #Hist_Dark_TelRaDec = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TelRaDec')
    #Make2DProjectionPlot(Hist_Dark_TelRaDec,'RA','Dec','Dark_TelRaDec')
    Hist_Target_TelRaDec = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TelRaDec_AfterCut')
    Make2DProjectionPlot(Hist_Target_TelRaDec,'RA','Dec','Target_TelRaDec')

    #Hist_Dark_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TelElevAzim')
    #Make2DProjectionPlot(Hist_Dark_TelElevAzim,'Elev','Azim','Dark_TelElevAzim')
    Hist_Target_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TelElevAzim')
    Make2DProjectionPlot(Hist_Target_TelElevAzim,'Elev','Azim','Target_TelElevAzim')

    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]

        which_method = 'MSCW'

        Hist_Target_SR_theta2 = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_theta2')
        Hist_Target_Bkg_theta2 = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_theta2_Sum')
        Hist_Target_SR_RaDec = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_RaDec')
        Hist_Target_Bkg_RaDec = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_RaDec_Sum')
        Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCW_Sum')
        Hist_Target_SR_MSCL = SelectDiagnosticaHistograms(folder,'MSCL','SR','Target_SR_MSCL')
        Hist_Target_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCW_Sum')
        Hist_Target_BkgCR1_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_BkgCR1_MSCW_Sum')
        Hist_Target_Deconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Deconv_MSCW')
        Hist_Target_TrueDeconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TrueDeconv_MSCW')
        Hist_Target_CR1_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR1_MSCW_Sum')
        Hist_Target_Ring_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Ring_MSCW')
        Hist_Target_Ring_MSCL = SelectDiagnosticaHistograms(folder,'MSCL','SR','Target_Ring_MSCL')
        Hist_Dark_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCW')
        Hist_Dark_CR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','CR','Dark_CR_MSCW')
        Hist_Dark_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Bkg_MSCW')
        Hist_Dark_Deconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Deconv_MSCW')
        Hist_Dark_TrueDeconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TrueDeconv_MSCW')

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_SR_MSCW]
        Hists[0].Rebin(2)
        if source=='2ndCrab':
            legends += ['Crab 2016 ON']
        else:
            legends += ['%s'%(source)]
        colors += [1]
        Hists += [Hist_Target_Bkg_MSCW]
        Hists[1].Rebin(2)
        legends += ['Bkg']
        colors += [4]
        plotname = 'Target_SR_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False)

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_CR1_MSCW]
        Hists[0].Rebin(2)
        if source=='2ndCrab':
            legends += ['Crab 2016 ON']
        else:
            legends += ['%s'%(source)]
        colors += [1]
        Hists += [Hist_Target_BkgCR1_MSCW]
        Hists[1].Rebin(2)
        legends += ['Bkg']
        colors += [4]
        plotname = 'Target_CR1_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False)

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_SR_theta2]
        if source=='2ndCrab':
            legends += ['Crab 2016 ON']
        else:
            legends += ['%s'%(source)]
        colors += [1]
        Hists += [Hist_Target_Bkg_theta2]
        legends += ['Bkg']
        colors += [4]
        plotname = 'Target_SR_theta2_E%s'%(ErecS_lower_cut)
        title = 'theta2'
        Hists[0].Rebin(2)
        Hists[1].Rebin(2)
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False)

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_SR_MSCL]
        if source=='2ndCrab':
            legends += ['Crab 2016 ON']
        else:
            legends += ['%s'%(source)]
        colors += [1]
        Hists += [Hist_Target_Ring_MSCL]
        legends += ['Ring']
        colors += [2]
        plotname = 'Target_MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeChi2Plot(Hists,legends,colors,title,plotname,True,False)


        #Hists = []
        #legends = []
        #colors = []
        #Hists += [Hist_Dark_SR_MSCW]
        #legends += ['SR']
        #colors += [1]
        #Hists += [Hist_Dark_Bkg_MSCW]
        #legends += ['Bkg']
        #colors += [4]
        ##Hists += [Hist_Dark_Elec_MSCW]
        ##legends += ['residual']
        ##colors += [3]
        #plotname = 'Dark_Deconv_MSCW_E%s'%(ErecS_lower_cut)
        #title = 'MSCW'
        ##MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)
        #MakeChi2Plot(Hists,legends,colors,title,plotname,True,False)

        #Hists = []
        #legends = []
        #colors = []
        #Hists += [Hist_Dark_SR_MSCW]
        #legends += ['MSCL < 0.5']
        #colors += [1]
        #Hists += [Hist_Dark_CR_MSCW]
        #legends += ['MSCL > 1.0']
        #colors += [2]
        #plotname = 'Dark_SRCR_MSCW_E%s'%(ErecS_lower_cut)
        #title = 'MSCW'
        #MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)


        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_TrueDeconv_MSCW]
        legends += ['true response']
        colors += [4]
        plotname = 'Target_Response_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeGaussianPlot(Hists,legends,colors,title,plotname,False,False)

        #Hists = []
        #legends = []
        #colors = []
        #Hists += [Hist_Dark_TrueDeconv_MSCW]
        #legends += ['true response']
        #colors += [4]
        #plotname = 'Dark_Response_MSCW_E%s'%(ErecS_lower_cut)
        #title = 'MSCW'
        #MakeGaussianPlot(Hists,legends,colors,title,plotname,False,False)

        #plotname = 'Target_SR_RaDec_E%s'%(ErecS_lower_cut)
        #smooth_size = 0.1
        #if energy_list[e]>=3000: smooth_size = 0.5
        #Hist_Target_SR_RaDec.Rebin2D(5,5)
        #Hist_Target_Bkg_RaDec.Rebin2D(5,5)
        #Hist_Target_SR_RaDec_Smooth = Smooth2DMap(Hist_Target_SR_RaDec,smooth_size)
        #Hist_Target_Bkg_RaDec_Smooth = Smooth2DMap(Hist_Target_Bkg_RaDec,smooth_size)
        #Make2DSignificancePlot(Hist_Target_SR_RaDec_Smooth,Hist_Target_Bkg_RaDec_Smooth,'RA','Dec',plotname)
