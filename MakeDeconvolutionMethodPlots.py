
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'

UseDarkRun = True

source = ''
source_list = []
#source_list  += ['Crab']
#source_list  += ['2ndCrab']
#source_list  += ['PKS1424']
#source_list  += ['3C264']
source_list  += ['H1426']
#source_list  += ['Ton599']
#source_list  += ['IC443']

Elev_lower_cut = 75
Elev_upper_cut = 85
#Elev_lower_cut = 25
#Elev_upper_cut = 55

Depth_cut_width = 10
Azim_lower_cut = 0
Azim_upper_cut = 360
ErecS_lower_cut = 100
ErecS_upper_cut = 1e10 
MSCL_lower_cut = -0.5
MSCL_upper_cut = 0.6
MSCW_lower_cut = -0.5
MSCW_upper_cut = 0.6

energy_list = []
energy_list += [100]
energy_list += [200]
energy_list += [500]
energy_list += [700]
energy_list += [1000]
energy_list += [2000]
energy_list += [4000]
energy_list += [10000]
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

    FilePath = '%s/Deconvolution_%s_Elev%sto%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut)
    print 'Read %s'%(FilePath)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    print 'HistList.GetSize() = %s'%(HistList.GetSize())
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
        if ErecS_lower<ErecS_lower_cut: continue
        if ErecS_upper>ErecS_upper_cut: continue
        if Hist_Data.Integral()==0:
            Hist_Data.Delete()
            Hist_Data = hist.Clone()
        else:
            Hist_Data.Add(hist)

    return Hist_Data

def MakeDiagnosticPlot(Hists,legends,title,name):
    
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
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].SetLineColor(h+1)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h
            Hists[h].GetXaxis().SetTitle(title)
            Hists[h].GetXaxis().SetRangeUser(-1,10)
            if title=='Depth': Hists[h].GetXaxis().SetRangeUser(0,20)

    Hists[max_hist].Draw("E")
    Hists[max_hist].SetMinimum(0)
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            Hists[h].Draw("E same")
    Hists[0].SetLineWidth(3)
    Hists[0].Draw("E same")
    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    for h in range(0,len(Hists)):
        if Hists[h]!=0:
            legend.AddEntry(Hists[h],legends[h],"pl")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    lumilab2 = ROOT.TLatex(0.15,0.65,'%0.2f < Tel elev. < %0.2f'%(Elev_lower_cut,Elev_upper_cut) )
    lumilab2.SetNDC()
    lumilab2.SetTextSize(0.15)
    lumilab2.Draw()
    lumilab5 = ROOT.TLatex(0.15,0.50,'%0.2f < Tel azim. < %0.2f'%(Azim_lower_cut,Azim_upper_cut) )
    lumilab5.SetNDC()
    lumilab5.SetTextSize(0.15)
    lumilab5.Draw()
    #pad1.SetLogy()
    c_both.SaveAs('output_plots/%s_%s.pdf'%(name,source))

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

    canvas2 = ROOT.TCanvas("canvas2","canvas2", 200, 10, 600, 600)
    Hist_1.RebinX(2)
    Hist_1.RebinY(2)
    Hist_1.Draw("COL4Z")
    canvas2.SaveAs('output_plots/%s_%s_Lego.pdf'%(name,'Crab'))
    Hist_2.Draw("COL4Z")
    Hist_2.RebinX(2)
    Hist_2.RebinY(2)
    canvas2.SaveAs('output_plots/%s_%s_Lego.pdf'%(name,source))


for s in source_list:
    source = s
    ErecS_lower_cut = 0
    ErecS_upper_cut = 1e10
    FilePath = '%s/Deconvolution_%s_Elev%sto%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut)
    print 'Read %s'%(FilePath)
    TargetFile=ROOT.TFile(FilePath)
    InfoTree = TargetFile.Get("InfoTree")
    Depth_cut_width = InfoTree.Depth_cut_width
    MSCW_lower_cut = InfoTree.MSCW_cut_lower
    MSCW_upper_cut = InfoTree.MSCW_cut_upper
    MSCL_lower_cut = InfoTree.MSCL_cut_lower
    MSCL_upper_cut = InfoTree.MSCL_cut_upper
    Hist_Target_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TelElevAzim')
    Hist_Dark_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TelElevAzim')
    Make2DTrajectoryPlot(Hist_Dark_TelElevAzim,Hist_Target_TelElevAzim,'Tel. elev.','Tel. azim.','TelElevAzim')
    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]

        which_method = 'MSCW'

        Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCW')
        Hist_Dark_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Bkg_MSCW')
        Hist_Dark_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCW')
        Hist_Target_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCW')
        Hist_Target_CR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_MSCW')

        Hist_Target_SR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCL')
        Hist_Dark_Bkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Bkg_MSCL')
        Hist_Dark_SR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCL')
        Hist_Target_Bkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCL')
        Hist_Target_CR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_MSCL')

        Hist_Target_SR_Depth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_Depth')
        Hist_Dark_Bkg_Depth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Bkg_Depth')
        Hist_Dark_SR_Depth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_Depth')
        Hist_Target_Bkg_Depth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_Depth')
        Hist_Target_CR_Depth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_Depth')

        Hists = []
        legends = []
        Hists += [Hist_Target_SR_MSCW]
        legends += ['Target SR']
        Hists += [Hist_Dark_SR_MSCW]
        legends += ['method 1']
        Hists += [Hist_Dark_Bkg_MSCW]
        legends += ['method 2']
        Hists += [Hist_Target_Bkg_MSCW]
        legends += ['method 3']
        plotname = 'MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeDiagnosticPlot(Hists,legends,title,plotname)

        Hists = []
        legends = []
        Hists += [Hist_Target_SR_MSCL]
        legends += ['Target SR']
        Hists += [Hist_Dark_SR_MSCL]
        legends += ['method 1']
        Hists += [Hist_Dark_Bkg_MSCL]
        legends += ['method 2']
        Hists += [Hist_Target_Bkg_MSCL]
        legends += ['method 3']
        plotname = 'MSCL_E%s'%(ErecS_lower_cut)
        title = 'MSCL'
        MakeDiagnosticPlot(Hists,legends,title,plotname)

        Hists = []
        legends = []
        Hists += [Hist_Target_SR_Depth]
        legends += ['Target SR']
        Hists += [Hist_Dark_SR_Depth]
        legends += ['method 1']
        Hists += [Hist_Dark_Bkg_Depth]
        legends += ['method 2']
        Hists += [Hist_Target_Bkg_Depth]
        legends += ['method 3']
        plotname = 'Depth_E%s'%(ErecS_lower_cut)
        title = 'Depth'
        MakeDiagnosticPlot(Hists,legends,title,plotname)
