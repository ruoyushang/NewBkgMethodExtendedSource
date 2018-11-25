
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'

UseDarkRun = False

source = ''
source_list = []
#source_list  += ['Crab']
#source_list  += ['2ndCrab']
#source_list  += ['PKS1424']
#source_list  += ['3C264']
#source_list  += ['H1426']
#source_list  += ['Ton599']
source_list  += ['IC443']

Elev_lower_cut = 65
Elev_upper_cut = 70
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

    #FilePath = '%s/Histograms_%s_%s_%s_Elev%sto%s.root'%(folder,source,method,isSR,Elev_lower_cut,Elev_upper_cut)
    FilePath = '%s/Histograms_%s_%s.root'%(folder,source,method)
    print 'Read %s'%(FilePath)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    print 'HistList.GetSize() = %s'%(HistList.GetSize())
    for h in range(0,HistList.GetSize()):
        name = HistList.At(h).GetName()
        hist = InputFile.Get(name)
        if not var in name: continue
        for e in range(0,len(name.split("_"))):
            element = name.split("_")[e]
            if 'ErecS' in element and 'to' in element:
                element = element.strip('ErecS')
                ErecS_lower = float(element.split('to')[0])
                ErecS_upper = float(element.split('to')[1])
        if ErecS_lower<ErecS_lower_cut: continue
        if ErecS_upper>ErecS_upper_cut: continue
        print 'found histogram.'
        if Hist_Data.Integral()==0:
            Hist_Data.Delete()
            Hist_Data = hist.Clone()
        else:
            Hist_Data.Add(hist)

    return Hist_Data

def MakeDiagnosticPlot(Hist_Crab_On,Hist_Crab_Off,Hist_Target_On,title,name):
    
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
    Hist_Crab_On.SetLineColor(2)
    Hist_Crab_On.SetLineWidth(3)
    Hist_Crab_Off.SetLineColor(4)
    Hist_Crab_Off.SetLineWidth(3)
    Hist_Target_On.SetLineColor(3)
    Hist_Target_On.SetLineWidth(3)
    if Hist_Crab_On.GetMaximum()>Hist_Crab_Off.GetMaximum():
        Hist_Crab_On.GetXaxis().SetTitle(title)
        Hist_Crab_On.GetXaxis().SetRangeUser(-1,10)
        Hist_Crab_On.Draw("E")
    else:
        Hist_Crab_Off.GetXaxis().SetTitle(title)
        Hist_Crab_Off.GetXaxis().SetRangeUser(-1,10)
        Hist_Crab_Off.Draw("E")
    Hist_Crab_On.Draw("E same")
    Hist_Crab_Off.Draw("E same")
    Hist_Target_On.Draw("E same")
    pad3.cd()
    legend = ROOT.TLegend(0.55,0.1,0.94,0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.2)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.Clear()
    mean = Hist_Crab_On.GetMean()
    rms = Hist_Crab_On.GetRMS()
    legend.AddEntry(Hist_Crab_On,'Crab ON (%.2f#pm%.2f)'%(mean,rms),"pl")
    mean = Hist_Crab_Off.GetMean()
    rms = Hist_Crab_Off.GetRMS()
    legend.AddEntry(Hist_Crab_Off,'Crab OFF (%.2f#pm%.2f)'%(mean,rms),"pl")
    mean = Hist_Target_On.GetMean()
    rms = Hist_Target_On.GetRMS()
    legend.AddEntry(Hist_Target_On,'%s (%.2f#pm%.2f)'%(source,mean,rms),"pl")
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

for s in source_list:
    source = s
    for e in range(0,len(energy_list)):
        ErecS_lower_cut = energy_list[e]

        which_method = 'MSCW'

        if not UseDarkRun:
            Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCW')
            Hist_Target_CR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_MSCW')
            Hist_Target_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCW')
            MakeDiagnosticPlot(Hist_Target_SR_MSCW,Hist_Target_CR_MSCW,Hist_Target_Bkg_MSCW,'MSCW','Target_Bkg_MSCW_E%s'%(ErecS_lower_cut))
            Hist_Target_ASR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ASR_MSCW')
            Hist_Target_ACR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ACR_MSCW')
            Hist_Target_ABkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ABkg_MSCW')
            MakeDiagnosticPlot(Hist_Target_ASR_MSCW,Hist_Target_ACR_MSCW,Hist_Target_ABkg_MSCW,'MSCW','Target_ABkg_MSCW_E%s'%(ErecS_lower_cut))

            #Hist_Target_SR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCL')
            #Hist_Target_CR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_MSCL')
            #Hist_Target_Bkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCL')
            #MakeDiagnosticPlot(Hist_Target_SR_MSCL,Hist_Target_CR_MSCL,Hist_Target_Bkg_MSCL,'MSCL','Target_Bkg_MSCL_E%s'%(ErecS_lower_cut))
            #Hist_Target_ASR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ASR_MSCL')
            #Hist_Target_ACR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ACR_MSCL')
            #Hist_Target_ABkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ABkg_MSCL')
            #MakeDiagnosticPlot(Hist_Target_ASR_MSCL,Hist_Target_ACR_MSCL,Hist_Target_ABkg_MSCL,'MSCL','Target_ABkg_MSCL_E%s'%(ErecS_lower_cut))

        if UseDarkRun:
            Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCW')
            Hist_Target_CR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCW')
            Hist_Target_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCW')
            MakeDiagnosticPlot(Hist_Target_SR_MSCW,Hist_Target_CR_MSCW,Hist_Target_Bkg_MSCW,'MSCW','Target_Bkg_MSCW_E%s'%(ErecS_lower_cut))
            Hist_Target_ASR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ASR_MSCW')
            Hist_Target_ACR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_ASR_MSCW')
            Hist_Target_ABkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ABkg_MSCW')
            MakeDiagnosticPlot(Hist_Target_ASR_MSCW,Hist_Target_ACR_MSCW,Hist_Target_ABkg_MSCW,'MSCW','Target_ABkg_MSCW_E%s'%(ErecS_lower_cut))

            #Hist_Target_SR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCL')
            #Hist_Target_CR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCL')
            #Hist_Target_Bkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCL')
            #MakeDiagnosticPlot(Hist_Target_SR_MSCL,Hist_Target_CR_MSCL,Hist_Target_Bkg_MSCL,'MSCL','Target_Bkg_MSCL_E%s'%(ErecS_lower_cut))
            #Hist_Target_ASR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ASR_MSCL')
            #Hist_Target_ACR_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_ASR_MSCL')
            #Hist_Target_ABkg_MSCL = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ABkg_MSCL')
            #MakeDiagnosticPlot(Hist_Target_ASR_MSCL,Hist_Target_ACR_MSCL,Hist_Target_ABkg_MSCL,'MSCL','Target_ABkg_MSCL_E%s'%(ErecS_lower_cut))

