
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'
MSCW_lower_cut = -1.0
MSCW_upper_cut = 1.0
#MSCW_lower_cut = 2
#MSCW_upper_cut = 3

source = ''
source_list = []
source_list  += ['2ndCrab']
source_list  += ['PKS1424']
source_list  += ['3C264']
source_list  += ['H1426']
source_list  += ['Ton599']
source_list  += ['IC443']

Elev_lower_cut = 55
Elev_upper_cut = 85
Azim_lower_cut = 0
Azim_upper_cut = 360
ErecS_lower_cut = 100
ErecS_upper_cut = 1e10 
MSCL_lower_cut = -1.0
MSCL_upper_cut = 1.0

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

AddSyst = True
#AddSyst = False

#UseLooseControlRegions = True
UseLooseControlRegions = False

#UseAuxRegion = True
UseAuxRegion = False

IncludeHeightMethod = False
IncludeDepthUpperMethod = False
IncludeDepthLowerMethod = False
IncludeMSCWMethod = False
IncludeMSCLMethod = False

#IncludeHeightMethod = True
IncludeDepthUpperMethod = True
#IncludeDepthLowerMethod = True
IncludeMSCWMethod = True
IncludeMSCLMethod = True

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def Make2DSignificancePlot(Hist_Data,Hist_CR,title,name):

    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    Hist_Sig = Hist_Data.Clone()
    for bx in range(0,Hist_Sig.GetNbinsX()):
        for by in range(0,Hist_Sig.GetNbinsY()):
            this_data = Hist_Data.GetBinContent(bx+1,by+1)
            this_data_err = Hist_Data.GetBinError(bx+1,by+1)
            this_bkg = Hist_CR.GetBinContent(bx+1,by+1)
            this_bkg_err = Hist_CR.GetBinError(bx+1,by+1)
            value = 0
            if not this_bkg_err==0:
                value = (this_data-this_bkg)/pow(this_bkg_err*this_bkg_err+this_data_err*this_data_err,0.5)
            Hist_Sig.SetBinContent(bx+1,by+1,value)
    Hist_Sig.GetYaxis().SetTitle('Yoff')
    Hist_Sig.GetXaxis().SetTitle('Xoff')
    Hist_Sig.GetZaxis().SetTitle('Significance')
    Hist_Sig.SetMaximum(5)
    Hist_Sig.SetMinimum(-3)
    Hist_Sig.Draw("COL4Z")
    if UseLooseControlRegions:
        c_both.SaveAs('output_plots/%s.pdf'%(name))
    else:
        c_both.SaveAs('output_plots/%s_Tight.pdf'%(name))

def Make1DProjectionPlot(Hist_Data,xtitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    hist_temp = []
    merge_bins = 7
    for b in range(0,Hist_Data.GetNbinsX()/merge_bins):
        hist_temp += [Hist_Data.ProjectionY("hist_temp_%s"%(b),b*merge_bins+1,b*merge_bins+merge_bins)]
    hist_temp[0].GetXaxis().SetTitle(xtitle)
    hist_temp[0].Rebin(2)
    #hist_temp[0].SetLineColor(1)
    hist_temp[0].Draw("HIST")
    for b in range(0,Hist_Data.GetNbinsX()/merge_bins):
        hist_temp[b].Rebin(2)
        #hist_temp[b].SetLineColor(b+2)
        hist_temp[b].Draw("HIST same")
    pad1.SetLogy()
    if not 'Crab' in name:
        canvas.SaveAs('output_plots/%s_%s.pdf'%(name,source))
    else:
        canvas.SaveAs('output_plots/%s.pdf'%(name))

def Make2DProjectionPlot(Hist_Data,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
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
    Hist_1D.Draw("E same")
    pad1.SetLogz()
    if 'Energy' in name:
        pad1.SetLogx()
    if not 'Crab' in name:
        canvas.SaveAs('output_plots/%s_%s.pdf'%(name,source))
    else:
        canvas.SaveAs('output_plots/%s.pdf'%(name))

def Make2DTrajectoryPlot(Hist_1,Hist_2,xtitle,ytitle,name):

    canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
    pad1.SetBottomMargin(0.15)
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.15)
    pad1.SetTopMargin(0.15)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad1.cd()
    Hist_1.GetYaxis().SetTitle(ytitle)
    Hist_1.GetXaxis().SetTitle(xtitle)
    Hist_1.SetLineColor(4)
    Hist_1.Draw("CONT3")
    Hist_2.SetLineColor(2)
    Hist_2.Draw("CONT3 same")
    canvas.SaveAs('output_plots/%s_%s.pdf'%(name,source))

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
        Hist_Crab_On.Draw("E")
    else:
        Hist_Crab_Off.GetXaxis().SetTitle(title)
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

def MakeStackPlot(Hist_Data,Hist_CR,title,name):
    
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
    stack = ROOT.THStack("stack", "")
    stack.Add( Hist_CR )
    #if Hist_Data.GetMaximum()>stack.GetMaximum():
    #    if Hist_Data.GetMinimum()>stack.GetMinimum():
    #        Hist_Data.SetMinimum(max(1.,stack.GetMinimum()/2.))
    #    Hist_Data.Draw("E")
    #else:
    #    if Hist_Data.GetMinimum()<stack.GetMinimum():
    #        stack.SetMinimum(max(1.,Hist_Data.GetMinimum()/2.))
    #    stack.Draw("hist")
    Hist_Data.SetMaximum(Hist_Data.GetMaximum()*1.2)
    Hist_Data.Draw("E")
    stack.Draw("hist same")
    Hist_Data.Draw("E same")
    Hist_All = Hist_CR.Clone()
    Hist_All.SetFillColor(1)
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
    legend.AddEntry(Hist_Data,'Data',"pl")
    legend.AddEntry(Hist_CR,'Background',"f")
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
    lumilab3 = ROOT.TLatex(0.15,0.35,'%0.2f < MSCL < %0.2f'%(MSCL_lower_cut,MSCL_upper_cut) )
    lumilab3.SetNDC()
    lumilab3.SetTextSize(0.15)
    lumilab3.Draw()
    lumilab4 = ROOT.TLatex(0.15,0.20,'%0.2f < MSCW < %0.2f'%(MSCW_lower_cut,MSCW_upper_cut) )
    lumilab4.SetNDC()
    lumilab4.SetTextSize(0.15)
    lumilab4.Draw()
    pad2.cd()
    Hist_Band = Hist_All.Clone()
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
    Hist_Band.GetYaxis().SetTitle("#sigma(data vs. bg)")
    #    #Hist_Band.GetYaxis().SetTitle("Data/Bkg")
    #    Hist_Band.GetYaxis().SetTitle("#sigma(data vs. bg)")
    Hist_Band.GetYaxis().SetTitleSize(0.13)
    Hist_Band.GetYaxis().SetNdivisions(505)
    Hist_Band.SetMaximum(5)
    Hist_Band.SetMinimum(-5)
    Hist_Band.Draw("e2")
    Hist_Ratio = Hist_Data.Clone()
    #    Hist_Ratio = Hist_Data.Clone()
    #Hist_Ratio.Divide(Hist_All)
    Hist_Ratio.Add(Hist_All,-1.)
    for b in range(0,Hist_Ratio.GetNbinsX()):
        unc = Hist_All.GetBinError(b+1)*Hist_All.GetBinError(b+1) + Hist_Data.GetBinError(b+1)*Hist_Data.GetBinError(b+1)
        #    unc = Hist_All.GetBinError(b+1)*Hist_All.GetBinError(b+1) + Hist_Data.GetBinError(b+1)*Hist_Data.GetBinError(b+1)
        unc = pow(unc,0.5)
        if not unc==0:
            Hist_Ratio.SetBinContent(b+1,Hist_Ratio.GetBinContent(b+1)/unc)
        else:
            Hist_Ratio.SetBinContent(b+1,0)
        Hist_Ratio.SetBinError(b+1,0)
    Hist_Ratio.SetFillColor(1)
    Hist_Ratio.Draw("B same")
    if title=='E [GeV]':
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()
    if UseLooseControlRegions:
        c_both.SaveAs('output_plots/%s.pdf'%(name))
    else:
        c_both.SaveAs('output_plots/%s_Tight.pdf'%(name))
    pad1.SetLogy(0)
    pad1.SetLogx(0)
    pad2.SetLogx(0)


def SelectHistograms(folder,method,isData,isSR):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

    if UseLooseControlRegions:
        FilePath = '%s/Histograms_%s_%s_%s.root'%(folder,source,method,isSR)
    else:
        FilePath = '%s/Histograms_%s_%s_%s_Tight.root'%(folder,source,method,isSR)
    print 'Read %s'%(FilePath)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    print 'HistList.GetSize() = %s'%(HistList.GetSize())
    for h in range(0,HistList.GetSize()):
        name = HistList.At(h).GetName()
        hist = InputFile.Get(name)
        Elev_lower = 0
        Elev_upper = 0
        Azim_lower = 0
        Azim_upper = 0
        ErecS_lower = 0
        ErecS_upper = 0
        if not Variable in name: continue
        for e in range(0,len(name.split("_"))):
            element = name.split("_")[e]
            if 'ErecS' in element and 'to' in element:
                element = element.strip('ErecS')
                ErecS_lower = float(element.split('to')[0])
                ErecS_upper = float(element.split('to')[1])
        if ErecS_lower<ErecS_lower_cut: continue
        if ErecS_upper>ErecS_upper_cut: continue
        if isData in name: 
            if Hist_Data.Integral()==0:
                Hist_Data.Delete()
                Hist_Data = hist.Clone()
            else:
                Hist_Data.Add(hist)

    return Hist_Data

def SelectDiagnosticaHistograms(folder,method,isSR,var):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

    if UseLooseControlRegions:
        FilePath = '%s/Histograms_%s_%s_%s.root'%(folder,source,method,isSR)
    else:
        FilePath = '%s/Histograms_%s_%s_%s_Tight.root'%(folder,source,method,isSR)
    print 'Read %s'%(FilePath)
    InputFile=ROOT.TFile(FilePath)
    
    HistList = InputFile.GetListOfKeys()
    
    if '2D' in var:
        var = 'Hist'+var
    else:
        var = 'Hist_'+var
    print 'HistList.GetSize() = %s'%(HistList.GetSize())
    for h in range(0,HistList.GetSize()):
        name = HistList.At(h).GetName()
        hist = InputFile.Get(name)
        if not var == name: continue
        Hist_Data = hist.Clone()

    return Hist_Data


for s in source_list:
    source = s
    for v in range(0,len(variable_list)):
        Variable = variable_list[v]
        xtitle = xtitle_list[v]
        if IncludeDepthUpperMethod:
            Hist2D_TelElev_vs_SlantDepth_Target = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','2D_TelElev_vs_SlantDepth_Target')
            Make2DProjectionPlot(Hist2D_TelElev_vs_SlantDepth_Target,'Tel. elev','Slant depth [1000kg/m^{2}]','TelElev_vs_SlantDepth')
            Hist2D_Energy_vs_SlantDepth_Target = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','2D_Energy_vs_SlantDepth_Target')
            Make2DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Target,'E [GeV]','Slant depth [1000kg/m^{2}]','Energy_vs_SlantDepth')
            Make1DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Target,'Slant depth [1000kg/m^{2}]','SlantDepth_as_FofE')
            Hist2D_Energy_vs_SlantDepth_Crab_On = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','2D_Energy_vs_SlantDepth_Crab_On')
            Make2DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Crab_On,'E [GeV]','Slant depth [1000kg/m^{2}]','Energy_vs_SlantDepth_Crab_On')
            Make1DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Crab_On,'Slant depth [1000kg/m^{2}]','SlantDepth_as_FofE_Crab_On')
            Hist2D_Energy_vs_SlantDepth_Crab_Off = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','2D_Energy_vs_SlantDepth_Crab_Off')
            Make2DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Crab_Off,'E [GeV]','Slant depth [1000kg/m^{2}]','Energy_vs_SlantDepth_Crab_Off')
            Make1DProjectionPlot(Hist2D_Energy_vs_SlantDepth_Crab_Off,'Slant depth [1000kg/m^{2}]','SlantDepth_as_FofE_Crab_Off')
            Hist2D_TelElev_vs_EmissionHeight_Target = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','2D_TelElev_vs_EmissionHeight_Target')
            Make2DProjectionPlot(Hist2D_TelElev_vs_EmissionHeight_Target,'Tel. elev','Emission height [km]','TelElev_vs_EmissionHeight')
            Hist_TelElevAzim_Counts_Target = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','TelElevAzim_Counts_Target')
            Hist_TelElevAzim_Counts_Crab = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','TelElevAzim_Counts_Crab')
            Make2DTrajectoryPlot(Hist_TelElevAzim_Counts_Crab,Hist_TelElevAzim_Counts_Target,'Tel. elev.','Tel. azim.','TelElevAzim')
            Hist_Loose_SlantDepth_Crab_On = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Crab_On')
            Hist_Loose_SlantDepth_Crab_Off = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Crab_Off')
            Hist_Loose_SlantDepth_Target_On = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Target_On')
            MakeDiagnosticPlot(Hist_Loose_SlantDepth_Crab_On,Hist_Loose_SlantDepth_Crab_Off,Hist_Loose_SlantDepth_Target_On,'Slant Depth [1000 kg/m^{2}]','Loose_SlantDepth_Upper')
            Hist_Loose_SlantDepth_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Crab_On_TeV')
            Hist_Loose_SlantDepth_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Crab_Off_TeV')
            Hist_Loose_SlantDepth_Target_On_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Loose_SlantDepth_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Loose_SlantDepth_Crab_On_TeV,Hist_Loose_SlantDepth_Crab_Off_TeV,Hist_Loose_SlantDepth_Target_On_TeV,'Slant Depth [1000 kg/m^{2}]','Loose_SlantDepth_TeV_Upper')
            Hist_Tight_SlantDepth_Crab_On = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Crab_On')
            Hist_Tight_SlantDepth_Crab_Off = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Crab_Off')
            Hist_Tight_SlantDepth_Target_On = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Target_On')
            MakeDiagnosticPlot(Hist_Tight_SlantDepth_Crab_On,Hist_Tight_SlantDepth_Crab_Off,Hist_Tight_SlantDepth_Target_On,'Slant Depth [1000 kg/m^{2}]','Tight_SlantDepth_Upper')
            Hist_Tight_SlantDepth_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Crab_On_TeV')
            Hist_Tight_SlantDepth_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Crab_Off_TeV')
            Hist_Tight_SlantDepth_Target_On_TeV = SelectDiagnosticaHistograms(folder,'DepthUpper','SR','Tight_SlantDepth_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Tight_SlantDepth_Crab_On_TeV,Hist_Tight_SlantDepth_Crab_Off_TeV,Hist_Tight_SlantDepth_Target_On_TeV,'Slant Depth [1000 kg/m^{2}]','Tight_SlantDepth_TeV_Upper')
        if IncludeDepthLowerMethod:
            Hist_Loose_SlantDepth_Crab_On = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Crab_On')
            Hist_Loose_SlantDepth_Crab_Off = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Crab_Off')
            Hist_Loose_SlantDepth_Target_On = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Target_On')
            MakeDiagnosticPlot(Hist_Loose_SlantDepth_Crab_On,Hist_Loose_SlantDepth_Crab_Off,Hist_Loose_SlantDepth_Target_On,'Slant Depth [1000 kg/m^{2}]','Loose_SlantDepth_Lower')
            Hist_Loose_SlantDepth_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Crab_On_TeV')
            Hist_Loose_SlantDepth_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Crab_Off_TeV')
            Hist_Loose_SlantDepth_Target_On_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Loose_SlantDepth_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Loose_SlantDepth_Crab_On_TeV,Hist_Loose_SlantDepth_Crab_Off_TeV,Hist_Loose_SlantDepth_Target_On_TeV,'Slant Depth [1000 kg/m^{2}]','Loose_SlantDepth_TeV_Lower')
            Hist_Tight_SlantDepth_Crab_On = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Crab_On')
            Hist_Tight_SlantDepth_Crab_Off = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Crab_Off')
            Hist_Tight_SlantDepth_Target_On = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Target_On')
            MakeDiagnosticPlot(Hist_Tight_SlantDepth_Crab_On,Hist_Tight_SlantDepth_Crab_Off,Hist_Tight_SlantDepth_Target_On,'Slant Depth [1000 kg/m^{2}]','Tight_SlantDepth_Lower')
            Hist_Tight_SlantDepth_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Crab_On_TeV')
            Hist_Tight_SlantDepth_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Crab_Off_TeV')
            Hist_Tight_SlantDepth_Target_On_TeV = SelectDiagnosticaHistograms(folder,'DepthLower','SR','Tight_SlantDepth_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Tight_SlantDepth_Crab_On_TeV,Hist_Tight_SlantDepth_Crab_Off_TeV,Hist_Tight_SlantDepth_Target_On_TeV,'Slant Depth [1000 kg/m^{2}]','Tight_SlantDepth_TeV_Lower')
        if IncludeMSCWMethod:
            Hist_Loose_MSCW_Crab_On = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Crab_On')
            Hist_Loose_MSCW_Crab_Off = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Crab_Off')
            Hist_Loose_MSCW_Target_On = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Target_On')
            MakeDiagnosticPlot(Hist_Loose_MSCW_Crab_On,Hist_Loose_MSCW_Crab_Off,Hist_Loose_MSCW_Target_On,'MSCW','Loose_MSCW')
            Hist_Loose_MSCW_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Crab_On_TeV')
            Hist_Loose_MSCW_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Crab_Off_TeV')
            Hist_Loose_MSCW_Target_On_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Loose_MSCW_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Loose_MSCW_Crab_On_TeV,Hist_Loose_MSCW_Crab_Off_TeV,Hist_Loose_MSCW_Target_On_TeV,'MSCW','Loose_MSCW_TeV')
            Hist_Tight_MSCW_Crab_On = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Crab_On')
            Hist_Tight_MSCW_Crab_Off = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Crab_Off')
            Hist_Tight_MSCW_Target_On = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Target_On')
            MakeDiagnosticPlot(Hist_Tight_MSCW_Crab_On,Hist_Tight_MSCW_Crab_Off,Hist_Tight_MSCW_Target_On,'MSCW','Tight_MSCW')
            Hist_Tight_MSCW_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Crab_On_TeV')
            Hist_Tight_MSCW_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Crab_Off_TeV')
            Hist_Tight_MSCW_Target_On_TeV = SelectDiagnosticaHistograms(folder,'MSCW','SR','Tight_MSCW_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Tight_MSCW_Crab_On_TeV,Hist_Tight_MSCW_Crab_Off_TeV,Hist_Tight_MSCW_Target_On_TeV,'MSCW','Tight_MSCW_TeV')
        
        if IncludeMSCLMethod:
            Hist_Loose_MSCL_Crab_On = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Crab_On')
            Hist_Loose_MSCL_Crab_Off = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Crab_Off')
            Hist_Loose_MSCL_Target_On = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Target_On')
            MakeDiagnosticPlot(Hist_Loose_MSCL_Crab_On,Hist_Loose_MSCL_Crab_Off,Hist_Loose_MSCL_Target_On,'MSCL','Loose_MSCL')
            Hist_Loose_MSCL_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Crab_On_TeV')
            Hist_Loose_MSCL_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Crab_Off_TeV')
            Hist_Loose_MSCL_Target_On_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Loose_MSCL_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Loose_MSCL_Crab_On_TeV,Hist_Loose_MSCL_Crab_Off_TeV,Hist_Loose_MSCL_Target_On_TeV,'MSCL','Loose_MSCL_TeV')
            Hist_Tight_MSCL_Crab_On = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Crab_On')
            Hist_Tight_MSCL_Crab_Off = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Crab_Off')
            Hist_Tight_MSCL_Target_On = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Target_On')
            MakeDiagnosticPlot(Hist_Tight_MSCL_Crab_On,Hist_Tight_MSCL_Crab_Off,Hist_Tight_MSCL_Target_On,'MSCL','Tight_MSCL')
            Hist_Tight_MSCL_Crab_On_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Crab_On_TeV')
            Hist_Tight_MSCL_Crab_Off_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Crab_Off_TeV')
            Hist_Tight_MSCL_Target_On_TeV = SelectDiagnosticaHistograms(folder,'MSCL','SR','Tight_MSCL_Target_On_TeV')
            MakeDiagnosticPlot(Hist_Tight_MSCL_Crab_On_TeV,Hist_Tight_MSCL_Crab_Off_TeV,Hist_Tight_MSCL_Target_On_TeV,'MSCL','Tight_MSCL_TeV')
        
        if not UseAuxRegion:
            if IncludeHeightMethod:
                Hist_Data_Height = SelectHistograms(folder,'Height','SR','SR')
                Hist_CR_Height = SelectHistograms(folder,'Height','Bkg','SR')
                Hist_Temp = Hist_Data_Height.Clone()
            if IncludeMSCWMethod:
                Hist_Data_MSCW = SelectHistograms(folder,'MSCW','SR','SR')
                Hist_CR_MSCW = SelectHistograms(folder,'MSCW','Bkg','SR')
                Hist_Temp = Hist_Data_MSCW.Clone()
            if IncludeMSCLMethod:
                Hist_Data_MSCL = SelectHistograms(folder,'MSCL','SR','SR')
                Hist_CR_MSCL = SelectHistograms(folder,'MSCL','Bkg','SR')
                Hist_Temp = Hist_Data_MSCL.Clone()
            if IncludeDepthUpperMethod:
                Hist_Data_DepthUpper = SelectHistograms(folder,'DepthUpper','SR','SR')
                Hist_CR_DepthUpper = SelectHistograms(folder,'DepthUpper','Bkg','SR')
                Hist_Temp = Hist_Data_DepthUpper.Clone()
            if IncludeDepthLowerMethod:
                Hist_Data_DepthLower = SelectHistograms(folder,'DepthLower','SR','SR')
                Hist_CR_DepthLower = SelectHistograms(folder,'DepthLower','Bkg','SR')
                Hist_Temp = Hist_Data_DepthLower.Clone()
        else:
            if IncludeHeightMethod:
                Hist_Data_Height = SelectHistograms(folder,'Height','SR','VR')
                Hist_CR_Height = SelectHistograms(folder,'Height','Bkg','VR')
                Hist_Temp = Hist_Data_Height.Clone()
            if IncludeMSCWMethod:
                Hist_Data_MSCW = SelectHistograms(folder,'MSCW','SR','VR')
                Hist_CR_MSCW = SelectHistograms(folder,'MSCW','Bkg','VR')
                Hist_Temp = Hist_Data_MSCW.Clone()
            if IncludeMSCLMethod:
                Hist_Data_MSCL = SelectHistograms(folder,'MSCL','SR','VR')
                Hist_CR_MSCL = SelectHistograms(folder,'MSCL','Bkg','VR')
                Hist_Temp = Hist_Data_MSCL.Clone()
            if IncludeDepthUpperMethod:
                Hist_Data_DepthUpper = SelectHistograms(folder,'DepthUpper','SR','VR')
                Hist_CR_DepthUpper = SelectHistograms(folder,'DepthUpper','Bkg','VR')
                Hist_Temp = Hist_Data_DepthUpper.Clone()
            if IncludeDepthLowerMethod:
                Hist_Data_DepthLower = SelectHistograms(folder,'DepthLower','SR','VR')
                Hist_CR_DepthLower = SelectHistograms(folder,'DepthLower','Bkg','VR')
                Hist_Temp = Hist_Data_DepthLower.Clone()
        
        if AddSyst:
            if IncludeHeightMethod:
                Hist_VR_Data_Height = SelectHistograms(folder,'Height','SR','VR')
                Hist_VR_CR_Height = SelectHistograms(folder,'Height','Bkg','VR')
                Hist_Syst_Height = Hist_VR_Data_Height.Clone()
                Hist_Syst_Height.Add(Hist_VR_CR_Height,-1)
                Hist_Syst_Height.Divide(Hist_VR_CR_Height)
            if IncludeDepthUpperMethod:
                Hist_VR_Data_DepthUpper = SelectHistograms(folder,'DepthUpper','SR','VR')
                Hist_VR_CR_DepthUpper = SelectHistograms(folder,'DepthUpper','Bkg','VR')
                Hist_Syst_DepthUpper = Hist_VR_Data_DepthUpper.Clone()
                Hist_Syst_DepthUpper.Add(Hist_VR_CR_DepthUpper,-1)
                Hist_Syst_DepthUpper.Divide(Hist_VR_CR_DepthUpper)
            if IncludeDepthLowerMethod:
                Hist_VR_Data_DepthLower = SelectHistograms(folder,'DepthLower','SR','VR')
                Hist_VR_CR_DepthLower = SelectHistograms(folder,'DepthLower','Bkg','VR')
                Hist_Syst_DepthLower = Hist_VR_Data_DepthLower.Clone()
                Hist_Syst_DepthLower.Add(Hist_VR_CR_DepthLower,-1)
                Hist_Syst_DepthLower.Divide(Hist_VR_CR_DepthLower)
            if IncludeMSCWMethod:
                Hist_VR_Data_MSCW = SelectHistograms(folder,'MSCW','SR','VR')
                Hist_VR_CR_MSCW = SelectHistograms(folder,'MSCW','Bkg','VR')
                Hist_Syst1_MSCW = Hist_VR_Data_MSCW.Clone()
                Hist_Syst1_MSCW.Add(Hist_VR_CR_MSCW,-1)
                Hist_Syst1_MSCW.Divide(Hist_VR_CR_MSCW)
                Hist_VR2_Data_MSCW = SelectHistograms(folder,'MSCW','SR','VR2')
                Hist_VR2_CR_MSCW = SelectHistograms(folder,'MSCW','Bkg','VR2')
                Hist_Syst2_MSCW = Hist_VR2_Data_MSCW.Clone()
                Hist_Syst2_MSCW.Add(Hist_VR2_CR_MSCW,-1)
                Hist_Syst2_MSCW.Divide(Hist_VR2_CR_MSCW)
                Hist_Syst_MSCW = Hist_VR_Data_MSCW.Clone()
                for b in range(0,Hist_Syst1_MSCW.GetNbinsX()):
                    Hist_Syst_MSCW.SetBinContent(b+1,pow(pow(Hist_Syst1_MSCW.GetBinContent(b+1),2)+pow(Hist_Syst2_MSCW.GetBinContent(b+1),2),0.5))
            if IncludeMSCLMethod:
                Hist_VR_Data_MSCL = SelectHistograms(folder,'MSCL','SR','VR')
                Hist_VR_CR_MSCL = SelectHistograms(folder,'MSCL','Bkg','VR')
                Hist_Syst1_MSCL = Hist_VR_Data_MSCL.Clone()
                Hist_Syst1_MSCL.Add(Hist_VR_CR_MSCL,-1)
                Hist_Syst1_MSCL.Divide(Hist_VR_CR_MSCL)
                Hist_VR2_Data_MSCL = SelectHistograms(folder,'MSCL','SR','VR2')
                Hist_VR2_CR_MSCL = SelectHistograms(folder,'MSCL','Bkg','VR2')
                Hist_Syst2_MSCL = Hist_VR2_Data_MSCL.Clone()
                Hist_Syst2_MSCL.Add(Hist_VR2_CR_MSCL,-1)
                Hist_Syst2_MSCL.Divide(Hist_VR2_CR_MSCL)
                Hist_Syst_MSCL = Hist_VR_Data_MSCL.Clone()
                for b in range(0,Hist_Syst1_MSCL.GetNbinsX()):
                    Hist_Syst_MSCL.SetBinContent(b+1,pow(pow(Hist_Syst1_MSCL.GetBinContent(b+1),2)+pow(Hist_Syst2_MSCL.GetBinContent(b+1),2),0.5))
        
        if IncludeHeightMethod:
            for b in range(0,Hist_CR_Height.GetNbinsX()):
                CR_Height = Hist_CR_Height.GetBinContent(b+1)
                ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
                if AddSyst:
                    SystCR_Height = (Hist_Syst_Height.GetBinContent(b+1)*Hist_CR_Height.GetBinContent(b+1))
                    SystCR_Err_Height = (Hist_Syst_Height.GetBinError(b+1)*Hist_CR_Height.GetBinContent(b+1))
                    SystCR_Height = pow(SystCR_Height*SystCR_Height+SystCR_Err_Height*SystCR_Err_Height,0.5)
                    ErrCR_Height = pow(ErrCR_Height*ErrCR_Height+SystCR_Height*SystCR_Height,0.5)
                    Hist_CR_Height.SetBinError(b+1,ErrCR_Height)
                    #Hist_CR_Height.SetBinContent(b+1,SystCR_Height+Hist_CR_Height.GetBinContent(b+1))
        if IncludeDepthUpperMethod:
            for b in range(0,Hist_CR_DepthUpper.GetNbinsX()):
                CR_DepthUpper = Hist_CR_DepthUpper.GetBinContent(b+1)
                ErrCR_DepthUpper = Hist_CR_DepthUpper.GetBinError(b+1)
                if AddSyst:
                    SystCR_DepthUpper = (Hist_Syst_DepthUpper.GetBinContent(b+1)*Hist_CR_DepthUpper.GetBinContent(b+1))
                    if not Hist_VR_CR_DepthUpper.GetBinContent(b+1)==0:
                        SystCR_Err_DepthUpper = Hist_VR_Data_DepthUpper.GetBinError(b+1)/Hist_VR_CR_DepthUpper.GetBinContent(b+1)*Hist_CR_DepthUpper.GetBinContent(b+1)
                    else:
                        SystCR_Err_DepthUpper = Hist_CR_DepthUpper.GetBinContent(b+1)
                    SystCR_DepthUpper = pow(SystCR_DepthUpper*SystCR_DepthUpper+SystCR_Err_DepthUpper*SystCR_Err_DepthUpper,0.5)
                    ErrCR_DepthUpper = pow(ErrCR_DepthUpper*ErrCR_DepthUpper+SystCR_DepthUpper*SystCR_DepthUpper,0.5)
                    Hist_CR_DepthUpper.SetBinError(b+1,ErrCR_DepthUpper)
                    #Hist_CR_DepthUpper.SetBinContent(b+1,SystCR_DepthUpper+Hist_CR_DepthUpper.GetBinContent(b+1))
        if IncludeDepthLowerMethod:
            for b in range(0,Hist_CR_DepthLower.GetNbinsX()):
                CR_DepthLower = Hist_CR_DepthLower.GetBinContent(b+1)
                ErrCR_DepthLower = Hist_CR_DepthLower.GetBinError(b+1)
                if AddSyst:
                    SystCR_DepthLower = (Hist_Syst_DepthLower.GetBinContent(b+1)*Hist_CR_DepthLower.GetBinContent(b+1))
                    if not Hist_VR_CR_DepthLower.GetBinContent(b+1)==0:
                        SystCR_Err_DepthLower = Hist_VR_Data_DepthLower.GetBinError(b+1)/Hist_VR_CR_DepthLower.GetBinContent(b+1)*Hist_CR_DepthLower.GetBinContent(b+1)
                    else:
                        SystCR_Err_DepthLower = Hist_CR_DepthLower.GetBinContent(b+1)
                    SystCR_DepthLower = pow(SystCR_DepthLower*SystCR_DepthLower+SystCR_Err_DepthLower*SystCR_Err_DepthLower,0.5)
                    ErrCR_DepthLower = pow(ErrCR_DepthLower*ErrCR_DepthLower+SystCR_DepthLower*SystCR_DepthLower,0.5)
                    Hist_CR_DepthLower.SetBinError(b+1,ErrCR_DepthLower)
                    #Hist_CR_DepthLower.SetBinContent(b+1,SystCR_DepthLower+Hist_CR_DepthLower.GetBinContent(b+1))
        if IncludeMSCWMethod:
            for b in range(0,Hist_CR_MSCW.GetNbinsX()):
                CR_MSCW = Hist_CR_MSCW.GetBinContent(b+1)
                ErrCR_MSCW = Hist_CR_MSCW.GetBinError(b+1)
                if AddSyst:
                    SystCR_MSCW = (Hist_Syst_MSCW.GetBinContent(b+1)*Hist_CR_MSCW.GetBinContent(b+1))
                    #SystCR_Err_MSCW = 0
                    #if not Hist_VR_CR_MSCW.GetBinContent(b+1)==0:
                    #    SystCR_Err_MSCW += Hist_VR_Data_MSCW.GetBinError(b+1)/Hist_VR_CR_MSCW.GetBinContent(b+1)*Hist_CR_MSCW.GetBinContent(b+1)
                    #else:
                    #    SystCR_Err_MSCW += Hist_CR_MSCW.GetBinContent(b+1)
                    #if not Hist_VR2_CR_MSCW.GetBinContent(b+1)==0:
                    #    SystCR_Err_MSCW += Hist_VR2_Data_MSCW.GetBinError(b+1)/Hist_VR2_CR_MSCW.GetBinContent(b+1)*Hist_CR_MSCW.GetBinContent(b+1)
                    #else:
                    #    SystCR_Err_MSCW += Hist_CR_MSCW.GetBinContent(b+1)
                    #SystCR_MSCW = pow(SystCR_MSCW*SystCR_MSCW+SystCR_Err_MSCW*SystCR_Err_MSCW,0.5)
                    ErrCR_MSCW = pow(ErrCR_MSCW*ErrCR_MSCW+SystCR_MSCW*SystCR_MSCW,0.5)
                    Hist_CR_MSCW.SetBinError(b+1,ErrCR_MSCW)
                    #Hist_CR_MSCW.SetBinContent(b+1,SystCR_MSCW+Hist_CR_MSCW.GetBinContent(b+1))
        if IncludeMSCLMethod:
            for b in range(0,Hist_CR_MSCL.GetNbinsX()):
                CR_MSCL = Hist_CR_MSCL.GetBinContent(b+1)
                ErrCR_MSCL = Hist_CR_MSCL.GetBinError(b+1)
                if AddSyst:
                    SystCR_MSCL = (Hist_Syst_MSCL.GetBinContent(b+1)*Hist_CR_MSCL.GetBinContent(b+1))
                    #SystCR_Err_MSCL = 0
                    #if not Hist_VR_CR_MSCL.GetBinContent(b+1)==0:
                    #    SystCR_Err_MSCL += Hist_VR_Data_MSCL.GetBinError(b+1)/Hist_VR_CR_MSCL.GetBinContent(b+1)*Hist_CR_MSCL.GetBinContent(b+1)
                    #else:
                    #    SystCR_Err_MSCL += Hist_CR_MSCL.GetBinContent(b+1)
                    #if not Hist_VR2_CR_MSCL.GetBinContent(b+1)==0:
                    #    SystCR_Err_MSCL += Hist_VR2_Data_MSCL.GetBinError(b+1)/Hist_VR2_CR_MSCL.GetBinContent(b+1)*Hist_CR_MSCL.GetBinContent(b+1)
                    #else:
                    #    SystCR_Err_MSCL += Hist_CR_MSCL.GetBinContent(b+1)
                    #SystCR_MSCL = pow(SystCR_MSCL*SystCR_MSCL+SystCR_Err_MSCL*SystCR_Err_MSCL,0.5)
                    ErrCR_MSCL = pow(ErrCR_MSCL*ErrCR_MSCL+SystCR_MSCL*SystCR_MSCL,0.5)
                    Hist_CR_MSCL.SetBinError(b+1,ErrCR_MSCL)
                    #Hist_CR_MSCL.SetBinContent(b+1,SystCR_MSCL+Hist_CR_MSCL.GetBinContent(b+1))
        
        if IncludeHeightMethod:
            print 'Hist_Data_Height.Integral() = %s'%(Hist_Data_Height.Integral())
        if IncludeDepthUpperMethod:
            print 'Hist_Data_DepthUpper.Integral() = %s'%(Hist_Data_DepthUpper.Integral())
        if IncludeDepthLowerMethod:
            print 'Hist_Data_DepthLower.Integral() = %s'%(Hist_Data_DepthLower.Integral())
        if IncludeMSCWMethod:
            print 'Hist_Data_MSCW.Integral() = %s'%(Hist_Data_MSCW.Integral())
        if IncludeMSCLMethod:
            print 'Hist_Data_MSCL.Integral() = %s'%(Hist_Data_MSCL.Integral())
        
        if IncludeHeightMethod:
            tag = Variable.strip('_')
            tag += '_%s'%(source)
            tag += '_%s'%('EmissionHeight')
            if UseAuxRegion:
                tag += '_Aux'
            if AddSyst:
                tag += '_Syst'
            if Variable == '_XoffVsYoff_':
                Make2DSignificancePlot(Hist_Data_Height,Hist_CR_Height,xtitle,tag)
            else:
                MakeStackPlot(Hist_Data_Height,Hist_CR_Height,xtitle,tag)
            
        if IncludeDepthUpperMethod:
            tag = Variable.strip('_')
            tag += '_%s'%(source)
            tag += '_%s'%('SlantDepthUpper')
            if UseAuxRegion:
                tag += '_Aux'
            if AddSyst:
                tag += '_Syst'
            if Variable == '_XoffVsYoff_':
                Make2DSignificancePlot(Hist_Data_DepthUpper,Hist_CR_DepthUpper,xtitle,tag)
            else:
                MakeStackPlot(Hist_Data_DepthUpper,Hist_CR_DepthUpper,xtitle,tag)
            
        if IncludeDepthLowerMethod:
            tag = Variable.strip('_')
            tag += '_%s'%(source)
            tag += '_%s'%('SlantDepthLower')
            if UseAuxRegion:
                tag += '_Aux'
            if AddSyst:
                tag += '_Syst'
            if Variable == '_XoffVsYoff_':
                Make2DSignificancePlot(Hist_Data_DepthLower,Hist_CR_DepthLower,xtitle,tag)
            else:
                MakeStackPlot(Hist_Data_DepthLower,Hist_CR_DepthLower,xtitle,tag)
            
        if IncludeMSCWMethod:
            tag = Variable.strip('_')
            tag += '_%s'%(source)
            tag += '_%s'%('MSCW')
            if UseAuxRegion:
                tag += '_Aux'
            if AddSyst:
                tag += '_Syst'
            if Variable == '_XoffVsYoff_':
                Make2DSignificancePlot(Hist_Data_MSCW,Hist_CR_MSCW,xtitle,tag)
            else:
                MakeStackPlot(Hist_Data_MSCW,Hist_CR_MSCW,xtitle,tag)
            
        if IncludeMSCLMethod:
            tag = Variable.strip('_')
            tag += '_%s'%(source)
            tag += '_%s'%('MSCL')
            if UseAuxRegion:
                tag += '_Aux'
            if AddSyst:
                tag += '_Syst'
            if Variable == '_XoffVsYoff_':
                Make2DSignificancePlot(Hist_Data_MSCL,Hist_CR_MSCL,xtitle,tag)
            else:
                MakeStackPlot(Hist_Data_MSCL,Hist_CR_MSCL,xtitle,tag)
        
        Hist_CR4 = Hist_Temp.Clone()
        for b in range(0,Hist_CR4.GetNbinsX()):
            if IncludeHeightMethod:
                CR_Height = Hist_CR_Height.GetBinContent(b+1)
                ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
            if IncludeDepthUpperMethod:
                CR_DepthUpper = Hist_CR_DepthUpper.GetBinContent(b+1)
                ErrCR_DepthUpper = Hist_CR_DepthUpper.GetBinError(b+1)
            if IncludeDepthLowerMethod:
                CR_DepthLower = Hist_CR_DepthLower.GetBinContent(b+1)
                ErrCR_DepthLower = Hist_CR_DepthLower.GetBinError(b+1)
            if IncludeMSCWMethod:
                CR_MSCW = Hist_CR_MSCW.GetBinContent(b+1)
                ErrCR_MSCW = Hist_CR_MSCW.GetBinError(b+1)
            if IncludeMSCLMethod:
                CR_MSCL = Hist_CR_MSCL.GetBinContent(b+1)
                ErrCR_MSCL = Hist_CR_MSCL.GetBinError(b+1)
            CR4 = 0
            Norm4 = 0
            ErrCR4 = 0
            if IncludeHeightMethod:
                if not ErrCR_Height==0:
                    CR4 += CR_Height/(ErrCR_Height*ErrCR_Height)
                    Norm4 += 1./(ErrCR_Height*ErrCR_Height)
                    ErrCR4 += 1./(ErrCR_Height*ErrCR_Height)
            if IncludeDepthUpperMethod:
                if not ErrCR_DepthUpper==0:
                    CR4 += CR_DepthUpper/(ErrCR_DepthUpper*ErrCR_DepthUpper)
                    Norm4 += 1./(ErrCR_DepthUpper*ErrCR_DepthUpper)
                    ErrCR4 += 1./(ErrCR_DepthUpper*ErrCR_DepthUpper)
            if IncludeDepthLowerMethod:
                if not ErrCR_DepthLower==0:
                    CR4 += CR_DepthLower/(ErrCR_DepthLower*ErrCR_DepthLower)
                    Norm4 += 1./(ErrCR_DepthLower*ErrCR_DepthLower)
                    ErrCR4 += 1./(ErrCR_DepthLower*ErrCR_DepthLower)
            if IncludeMSCWMethod:
                if not ErrCR_MSCW==0:
                    CR4 += CR_MSCW/(ErrCR_MSCW*ErrCR_MSCW)
                    Norm4 += 1./(ErrCR_MSCW*ErrCR_MSCW)
                    ErrCR4 += 1./(ErrCR_MSCW*ErrCR_MSCW)
            if IncludeMSCLMethod:
                if not ErrCR_MSCL==0:
                    CR4 += CR_MSCL/(ErrCR_MSCL*ErrCR_MSCL)
                    Norm4 += 1./(ErrCR_MSCL*ErrCR_MSCL)
                    ErrCR4 += 1./(ErrCR_MSCL*ErrCR_MSCL)
            if not Norm4==0: 
                CR4 = CR4/Norm4
                ErrCR4 = 1./pow(ErrCR4,0.5)
            Hist_CR4.SetBinContent(b+1,CR4)
            Hist_CR4.SetBinError(b+1,ErrCR4)
        tag = Variable.strip('_')
        tag += '_%s'%(source)
        tag += '_%s'%('Combined')
        if UseAuxRegion:
            tag += '_Aux'
        if AddSyst:
            tag += '_Syst'
        if Variable == '_XoffVsYoff_':
            Make2DSignificancePlot(Hist_Temp,Hist_CR4,xtitle,tag)
        else:
            MakeStackPlot(Hist_Temp,Hist_CR4,xtitle,tag)
        
        if not Variable == '_XoffVsYoff_':
            bins = []
            for b in range(0,Hist_Temp.GetNbinsX()+1):
                bins += [Hist_Temp.GetBinLowEdge(b+1)]
            Hist_Sensitivity_Height = ROOT.TH1D("Hist_Sensitivity_Height","",len(bins)-1,array('d',bins))
            Hist_Sensitivity_DepthUpper = ROOT.TH1D("Hist_Sensitivity_DepthUpper","",len(bins)-1,array('d',bins))
            Hist_Sensitivity_DepthLower = ROOT.TH1D("Hist_Sensitivity_DepthLower","",len(bins)-1,array('d',bins))
            Hist_Sensitivity_MSCW = ROOT.TH1D("Hist_Sensitivity_MSCW","",len(bins)-1,array('d',bins))
            Hist_Sensitivity_MSCL = ROOT.TH1D("Hist_Sensitivity_MSCL","",len(bins)-1,array('d',bins))
            Hist_Sensitivity4 = ROOT.TH1D("Hist_Sensitivity4","",len(bins)-1,array('d',bins))
            Hist_Sensitivity5 = ROOT.TH1D("Hist_Sensitivity5","",len(bins)-1,array('d',bins))
            for b in range(0,Hist_CR4.GetNbinsX()):
                if IncludeHeightMethod:
                    CR_Height = Hist_CR_Height.GetBinContent(b+1)
                    ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
                if IncludeDepthUpperMethod:
                    CR_DepthUpper = Hist_CR_DepthUpper.GetBinContent(b+1)
                    ErrCR_DepthUpper = Hist_CR_DepthUpper.GetBinError(b+1)
                if IncludeDepthLowerMethod:
                    CR_DepthLower = Hist_CR_DepthLower.GetBinContent(b+1)
                    ErrCR_DepthLower = Hist_CR_DepthLower.GetBinError(b+1)
                if IncludeMSCWMethod:
                    CR_MSCW = Hist_CR_MSCW.GetBinContent(b+1)
                    ErrCR_MSCW = Hist_CR_MSCW.GetBinError(b+1)
                if IncludeMSCLMethod:
                    CR_MSCL = Hist_CR_MSCL.GetBinContent(b+1)
                    ErrCR_MSCL = Hist_CR_MSCL.GetBinError(b+1)
                CR4 = 0
                Norm4 = 0
                ErrCR4 = 0
                if IncludeHeightMethod:
                    if not ErrCR_Height==0:
                        CR4 += CR_Height/(ErrCR_Height*ErrCR_Height)
                        Norm4 += 1./(ErrCR_Height*ErrCR_Height)
                        ErrCR4 += 1./(ErrCR_Height*ErrCR_Height)
                if IncludeDepthUpperMethod:
                    if not ErrCR_DepthUpper==0:
                        CR4 += CR_DepthUpper/(ErrCR_DepthUpper*ErrCR_DepthUpper)
                        Norm4 += 1./(ErrCR_DepthUpper*ErrCR_DepthUpper)
                        ErrCR4 += 1./(ErrCR_DepthUpper*ErrCR_DepthUpper)
                if IncludeDepthLowerMethod:
                    if not ErrCR_DepthLower==0:
                        CR4 += CR_DepthLower/(ErrCR_DepthLower*ErrCR_DepthLower)
                        Norm4 += 1./(ErrCR_DepthLower*ErrCR_DepthLower)
                        ErrCR4 += 1./(ErrCR_DepthLower*ErrCR_DepthLower)
                if IncludeMSCWMethod:
                    if not ErrCR_MSCW==0:
                        CR4 += CR_MSCW/(ErrCR_MSCW*ErrCR_MSCW)
                        Norm4 += 1./(ErrCR_MSCW*ErrCR_MSCW)
                        ErrCR4 += 1./(ErrCR_MSCW*ErrCR_MSCW)
                if IncludeMSCLMethod:
                    if not ErrCR_MSCL==0:
                        CR4 += CR_MSCL/(ErrCR_MSCL*ErrCR_MSCL)
                        Norm4 += 1./(ErrCR_MSCL*ErrCR_MSCL)
                        ErrCR4 += 1./(ErrCR_MSCL*ErrCR_MSCL)
                #CR4 = 0
                #ErrCR4 = 0
                #if not Norm4==0:
                #    CR4 = CR4/Norm4
                #if not ErrCR4==0:
                #    ErrCR4 = 1./pow(ErrCR4,0.5)
                CR4 = CR4/Norm4
                ErrCR4 = 1./pow(ErrCR4,0.5)
                if IncludeHeightMethod:
                    if not CR_Height==0:
                        Hist_Sensitivity_Height.SetBinContent(b+1,5.*pow(ErrCR_Height*ErrCR_Height+CR_Height,0.5)/CR_Height)
                if IncludeDepthUpperMethod:
                    if not CR_DepthUpper==0:
                        Hist_Sensitivity_DepthUpper.SetBinContent(b+1,5.*pow(ErrCR_DepthUpper*ErrCR_DepthUpper+CR_DepthUpper,0.5)/CR_DepthUpper)
                if IncludeDepthLowerMethod:
                    if not CR_DepthLower==0:
                        Hist_Sensitivity_DepthLower.SetBinContent(b+1,5.*pow(ErrCR_DepthLower*ErrCR_DepthLower+CR_DepthLower,0.5)/CR_DepthLower)
                if IncludeMSCWMethod:
                    if not CR_MSCW==0:
                        Hist_Sensitivity_MSCW.SetBinContent(b+1,5.*pow(ErrCR_MSCW*ErrCR_MSCW+CR_MSCW,0.5)/CR_MSCW)
                if IncludeMSCLMethod:
                    if not CR_MSCL==0:
                        Hist_Sensitivity_MSCL.SetBinContent(b+1,5.*pow(ErrCR_MSCL*ErrCR_MSCL+CR_MSCL,0.5)/CR_MSCL)
                if not CR4==0:
                    Hist_Sensitivity4.SetBinContent(b+1,5.*pow(ErrCR4*ErrCR4+CR4,0.5)/CR4)
                Data_Depth = Hist_Temp.GetBinContent(b+1)
                ErrData_Depth = Hist_Temp.GetBinError(b+1)
                if not Data_Depth==0:
                    #Hist_Sensitivity5.SetBinContent(b+1,5.*pow(ErrData_Depth*ErrData_Depth+Data_Depth,0.5)/Data_Depth)
                    Hist_Sensitivity5.SetBinContent(b+1,5.*pow(Data_Depth,0.5)/Data_Depth)
            
            canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
            pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
            pad1.SetBottomMargin(0.15)
            pad1.SetRightMargin(0.15)
            pad1.SetLeftMargin(0.15)
            pad1.SetTopMargin(0.0)
            pad1.SetBorderMode(0)
            pad2 = ROOT.TPad("pad2","pad2",0,0.8,1,1)
            pad2.SetBottomMargin(0.0)
            pad2.SetRightMargin(0.15)
            pad2.SetLeftMargin(0.15)
            pad2.SetTopMargin(0.15)
            pad2.SetBorderMode(0)
            pad1.Draw()
            pad2.Draw()
            pad1.cd()
            Hist_Sensitivity4.GetYaxis().SetTitle('Excluded signal flux (Bg unit.)')
            Hist_Sensitivity4.GetXaxis().SetTitle(xtitle)
            Hist_Sensitivity4.SetLineColor(1)
            Hist_Sensitivity4.SetLineWidth(2)
            Hist_Sensitivity5.SetLineColor(6)
            if Variable == '_ErecS_':
                Hist_Sensitivity5.Draw()
                Hist_Sensitivity4.Draw("same")
            else:
                Hist_Sensitivity4.SetMaximum(1)
                Hist_Sensitivity4.SetMinimum(0)
                Hist_Sensitivity4.Draw()
                Hist_Sensitivity5.Draw("same")
            if IncludeDepthUpperMethod:
                Hist_Sensitivity_DepthUpper.SetLineColor(2)
                Hist_Sensitivity_DepthUpper.Draw("same")
            if IncludeDepthLowerMethod:
                Hist_Sensitivity_DepthLower.SetLineColor(5)
                Hist_Sensitivity_DepthLower.Draw("same")
            if IncludeMSCWMethod:
                Hist_Sensitivity_MSCW.SetLineColor(4)
                Hist_Sensitivity_MSCW.Draw("same")
            if IncludeMSCLMethod:
                Hist_Sensitivity_MSCL.SetLineColor(3)
                Hist_Sensitivity_MSCL.Draw("same")
            Hist_Sensitivity5.Draw("same")
            Hist_Sensitivity4.Draw("same")
            if Variable == '_ErecS_':
                pad1.SetLogx()
                pad1.SetLogy()
            pad2.cd()
            legend = ROOT.TLegend(0.1,0.1,0.6,0.9)
            legend.SetTextFont(42)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.2)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetLineColor(0)
            legend.Clear()
            if IncludeHeightMethod:
                legend.AddEntry(Hist_Sensitivity_Height,'Height method',"pl")
            if IncludeDepthUpperMethod:
                legend.AddEntry(Hist_Sensitivity_DepthUpper,'Depth (upper) method',"pl")
            if IncludeDepthLowerMethod:
                legend.AddEntry(Hist_Sensitivity_DepthLower,'Depth (lower) method',"pl")
            if IncludeMSCWMethod:
                legend.AddEntry(Hist_Sensitivity_MSCW,'MSCW method',"pl")
            if IncludeMSCLMethod:
                legend.AddEntry(Hist_Sensitivity_MSCL,'MSCL method',"pl")
            legend.AddEntry(Hist_Sensitivity4,'Combined method',"pl")
            legend.AddEntry(Hist_Sensitivity5,'Some perfect method',"pl")
            legend.Draw("SAME")
            if UseLooseControlRegions:
                canvas.SaveAs('output_plots/Sensitivity_%s.pdf'%(tag))
            else:
                canvas.SaveAs('output_plots/Sensitivity_%s_Tight.pdf'%(tag))
