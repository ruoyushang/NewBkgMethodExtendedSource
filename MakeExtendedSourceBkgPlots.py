
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'
MSCW_lower_cut = -2
MSCW_upper_cut = 0.7
#MSCW_lower_cut = 2
#MSCW_upper_cut = 3

#source  = '2ndCrab'
#source  = 'PKS1424'
#source  = '3C264'
#source  = 'H1426'
source  = 'IC443'
#source  = 'Ton599'

Elev_lower_cut = 60
Elev_upper_cut = 80
Azim_lower_cut = 0
Azim_upper_cut = 360
ErecS_lower_cut = 100
ErecS_upper_cut = 10000 
MSCL_lower_cut = -2
MSCL_upper_cut = 0.7

#Variable = '_Energy_'
#xtitle = 'E [TeV]'
Variable = '_theta2_'
xtitle = '#theta^{2}'
#Variable = '_theta2zoomin_'
#xtitle = '#theta^{2}'
#Variable = '_XoffVsYoff_'
#xtitle = ''

#AddSyst = True
AddSyst = False

IncludeHeightMethod = False
IncludeDepthMethod = True
IncludeMSCWMethod = False
IncludeMSCLMethod = False

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

def Make2DPlot(Hist_Data,Hist_CR,title,name):

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
    c_both.SaveAs('output_plots/%s.pdf'%(name))

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
    if Hist_Data.GetMaximum()>stack.GetMaximum():
        if Hist_Data.GetMinimum()>stack.GetMinimum():
            Hist_Data.SetMinimum(max(1.,stack.GetMinimum()/2.))
        Hist_Data.Draw("E")
    else:
        if Hist_Data.GetMinimum()<stack.GetMinimum():
            stack.SetMinimum(max(1.,Hist_Data.GetMinimum()/2.))
        stack.Draw("hist")
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
    if title=='E [TeV]':
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()
    c_both.SaveAs('output_plots/%s.pdf'%(name))
    pad1.SetLogy(0)
    pad1.SetLogx(0)
    pad2.SetLogx(0)


def SelectHistograms(folder,method,isData,isSR):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

    if isSR=='VR':
        FilePath = '%s/Histograms_%s_%s_%s.root'%(folder,source,method,isSR)
    else:
        FilePath = '%s/Histograms_%s_%s.root'%(folder,source,method)
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
            if 'Elev' in element:
                element = element.strip('Elev')
                Elev_lower = float(element.split('to')[0])
                Elev_upper = float(element.split('to')[1])
            if 'Azim' in element:
                element = element.strip('Azim')
                Azim_lower = float(element.split('to')[0])
                Azim_upper = float(element.split('to')[1])
            if 'ErecS' in element and 'to' in element:
                element = element.strip('ErecS')
                ErecS_lower = float(element.split('to')[0])
                ErecS_upper = float(element.split('to')[1])
        if Elev_lower<Elev_lower_cut: continue
        if Elev_upper>Elev_upper_cut: continue
        if Azim_lower<Azim_lower_cut: continue
        if Azim_upper>Azim_upper_cut: continue
        if ErecS_lower<ErecS_lower_cut: continue
        if ErecS_upper>ErecS_upper_cut: continue
        if isData in name: 
            if Hist_Data.Integral()==0:
                Hist_Data.Delete()
                Hist_Data = hist.Clone()
            else:
                Hist_Data.Add(hist)

    return Hist_Data

if MSCW_upper_cut<1:
    if IncludeHeightMethod:
        Hist_Data_Height = SelectHistograms(folder,'Height','Data','SR')
        Hist_CR_Height = SelectHistograms(folder,'Height','CR','SR')
    if IncludeMSCWMethod:
        Hist_Data_MSCW = SelectHistograms(folder,'MSCW','Data','SR')
        Hist_CR_MSCW = SelectHistograms(folder,'MSCW','CR','SR')
    if IncludeMSCLMethod:
        Hist_Data_MSCL = SelectHistograms(folder,'MSCL','Data','SR')
        Hist_CR_MSCL = SelectHistograms(folder,'MSCL','CR','SR')
    if IncludeDepthMethod:
        Hist_Data_Depth = SelectHistograms(folder,'Depth','Data','SR')
        Hist_CR_Depth = SelectHistograms(folder,'Depth','CR','SR')
else:
    if IncludeHeightMethod:
        Hist_Data_Height = SelectHistograms(folder,'Height','Data','VR')
        Hist_CR_Height = SelectHistograms(folder,'Height','CR','VR')
    if IncludeMSCWMethod:
        Hist_Data_MSCW = SelectHistograms(folder,'MSCW','Data','VR')
        Hist_CR_MSCW = SelectHistograms(folder,'MSCW','CR','VR')
    if IncludeMSCLMethod:
        Hist_Data_MSCL = SelectHistograms(folder,'MSCL','Data','VR')
        Hist_CR_MSCL = SelectHistograms(folder,'MSCL','CR','VR')
    if IncludeDepthMethod:
        Hist_Data_Depth = SelectHistograms(folder,'Depth','Data','VR')
        Hist_CR_Depth = SelectHistograms(folder,'Depth','CR','VR')

if AddSyst:
    if IncludeHeightMethod:
        Hist_VR_Data_Height = SelectHistograms(folder,'Height','Data','VR')
        Hist_VR_CR_Height = SelectHistograms(folder,'Height','CR','VR')
        Hist_Syst_Height = Hist_VR_Data_Height.Clone()
        Hist_Syst_Height.Add(Hist_VR_CR_Height,-1)
        Hist_Syst_Height.Divide(Hist_VR_CR_Height)
    if IncludeDepthMethod:
        Hist_VR_Data_Depth = SelectHistograms(folder,'Depth','Data','VR')
        Hist_VR_CR_Depth = SelectHistograms(folder,'Depth','CR','VR')
        Hist_Syst_Depth = Hist_VR_Data_Depth.Clone()
        Hist_Syst_Depth.Add(Hist_VR_CR_Depth,-1)
        Hist_Syst_Depth.Divide(Hist_VR_CR_Depth)
    if IncludeMSCWMethod:
        Hist_VR_Data_MSCW = SelectHistograms(folder,'MSCW','Data','VR')
        Hist_VR_CR_MSCW = SelectHistograms(folder,'MSCW','CR','VR')
        Hist_Syst_MSCW = Hist_VR_Data_MSCW.Clone()
        Hist_Syst_MSCW.Add(Hist_VR_CR_MSCW,-1)
        Hist_Syst_MSCW.Divide(Hist_VR_CR_MSCW)
        #Hist_Syst_MSCW.Scale(2)
    if IncludeMSCLMethod:
        Hist_VR_Data_MSCL = SelectHistograms(folder,'MSCL','Data','VR')
        Hist_VR_CR_MSCL = SelectHistograms(folder,'MSCL','CR','VR')
        Hist_Syst_MSCL = Hist_VR_Data_MSCL.Clone()
        Hist_Syst_MSCL.Add(Hist_VR_CR_MSCL,-1)
        Hist_Syst_MSCL.Divide(Hist_VR_CR_MSCL)
        #Hist_Syst_MSCL.Scale(2)

if IncludeHeightMethod:
    for b in range(0,Hist_CR_Height.GetNbinsX()):
        CR_Height = Hist_CR_Height.GetBinContent(b+1)
        ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
        if AddSyst:
            SystCR_Height = abs(Hist_Syst_Height.GetBinContent(b+1)*Hist_CR_Height.GetBinContent(b+1))
            ErrCR_Height = pow(ErrCR_Height*ErrCR_Height+SystCR_Height*SystCR_Height,0.5)
            Hist_CR_Height.SetBinError(b+1,ErrCR_Height)
if IncludeDepthMethod:
    for b in range(0,Hist_CR_Depth.GetNbinsX()):
        CR_Depth = Hist_CR_Depth.GetBinContent(b+1)
        ErrCR_Depth = Hist_CR_Depth.GetBinError(b+1)
        if AddSyst:
            SystCR_Depth = abs(Hist_Syst_Depth.GetBinContent(b+1)*Hist_CR_Depth.GetBinContent(b+1))
            ErrCR_Depth = pow(ErrCR_Depth*ErrCR_Depth+SystCR_Depth*SystCR_Depth,0.5)
            Hist_CR_Depth.SetBinError(b+1,ErrCR_Depth)
if IncludeMSCWMethod:
    for b in range(0,Hist_CR_MSCW.GetNbinsX()):
        CR_MSCW = Hist_CR_MSCW.GetBinContent(b+1)
        ErrCR_MSCW = Hist_CR_MSCW.GetBinError(b+1)
        if AddSyst:
            SystCR_MSCW = abs(Hist_Syst_MSCW.GetBinContent(b+1)*Hist_CR_MSCW.GetBinContent(b+1))
            ErrCR_MSCW = pow(ErrCR_MSCW*ErrCR_MSCW+SystCR_MSCW*SystCR_MSCW,0.5)
            Hist_CR_MSCW.SetBinError(b+1,ErrCR_MSCW)
if IncludeMSCLMethod:
    for b in range(0,Hist_CR_MSCL.GetNbinsX()):
        CR_MSCL = Hist_CR_MSCL.GetBinContent(b+1)
        ErrCR_MSCL = Hist_CR_MSCL.GetBinError(b+1)
        if AddSyst:
            SystCR_MSCL = abs(Hist_Syst_MSCL.GetBinContent(b+1)*Hist_CR_MSCL.GetBinContent(b+1))
            ErrCR_MSCL = pow(ErrCR_MSCL*ErrCR_MSCL+SystCR_MSCL*SystCR_MSCL,0.5)
            Hist_CR_MSCL.SetBinError(b+1,ErrCR_MSCL)

if IncludeHeightMethod:
    print 'Hist_Data_Height.Integral() = %s'%(Hist_Data_Height.Integral())
if IncludeDepthMethod:
    print 'Hist_Data_Depth.Integral() = %s'%(Hist_Data_Depth.Integral())
if IncludeMSCWMethod:
    print 'Hist_Data_MSCW.Integral() = %s'%(Hist_Data_MSCW.Integral())
if IncludeMSCLMethod:
    print 'Hist_Data_MSCL.Integral() = %s'%(Hist_Data_MSCL.Integral())

if IncludeHeightMethod:
    tag = Variable.strip('_')
    tag += '_%s'%(source)
    tag += '_%s'%('EmissionHeight')
    if MSCW_lower_cut>1:
        tag += '_Validation'
    if AddSyst:
        tag += '_Syst'
    if Variable == '_XoffVsYoff_':
        Make2DPlot(Hist_Data_Height,Hist_CR_Height,xtitle,tag)
    else:
        MakeStackPlot(Hist_Data_Height,Hist_CR_Height,xtitle,tag)
    
if IncludeDepthMethod:
    tag = Variable.strip('_')
    tag += '_%s'%(source)
    tag += '_%s'%('SlantDepth')
    if MSCW_lower_cut>1:
        tag += '_Validation'
    if AddSyst:
        tag += '_Syst'
    if Variable == '_XoffVsYoff_':
        Make2DPlot(Hist_Data_Depth,Hist_CR_Depth,xtitle,tag)
    else:
        MakeStackPlot(Hist_Data_Depth,Hist_CR_Depth,xtitle,tag)
    
if IncludeMSCWMethod:
    tag = Variable.strip('_')
    tag += '_%s'%(source)
    tag += '_%s'%('MSCW')
    if MSCW_lower_cut>1:
        tag += '_Validation'
    if AddSyst:
        tag += '_Syst'
    if Variable == '_XoffVsYoff_':
        Make2DPlot(Hist_Data_MSCW,Hist_CR_MSCW,xtitle,tag)
    else:
        MakeStackPlot(Hist_Data_MSCW,Hist_CR_MSCW,xtitle,tag)
    
if IncludeMSCLMethod:
    tag = Variable.strip('_')
    tag += '_%s'%(source)
    tag += '_%s'%('MSCL')
    if MSCW_lower_cut>1:
        tag += '_Validation'
    if AddSyst:
        tag += '_Syst'
    if Variable == '_XoffVsYoff_':
        Make2DPlot(Hist_Data_MSCL,Hist_CR_MSCL,xtitle,tag)
    else:
        MakeStackPlot(Hist_Data_MSCL,Hist_CR_MSCL,xtitle,tag)

Hist_CR4 = Hist_CR_Depth.Clone()
for b in range(0,Hist_CR4.GetNbinsX()):
    if IncludeHeightMethod:
        CR_Height = Hist_CR_Height.GetBinContent(b+1)
        ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
    if IncludeDepthMethod:
        CR_Depth = Hist_CR_Depth.GetBinContent(b+1)
        ErrCR_Depth = Hist_CR_Depth.GetBinError(b+1)
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
    if IncludeDepthMethod:
        if not ErrCR_Depth==0:
            CR4 += CR_Depth/(ErrCR_Depth*ErrCR_Depth)
            Norm4 += 1./(ErrCR_Depth*ErrCR_Depth)
            ErrCR4 += 1./(ErrCR_Depth*ErrCR_Depth)
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
    CR4 = CR4/Norm4
    ErrCR4 = 1./pow(ErrCR4,0.5)
    Hist_CR4.SetBinContent(b+1,CR4)
    Hist_CR4.SetBinError(b+1,ErrCR4)
tag = Variable.strip('_')
tag += '_%s'%(source)
tag += '_%s'%('Combined')
if MSCW_lower_cut>1:
    tag += '_Validation'
if AddSyst:
    tag += '_Syst'
if Variable == '_XoffVsYoff_':
    Make2DPlot(Hist_Data_Depth,Hist_CR4,xtitle,tag)
else:
    MakeStackPlot(Hist_Data_Depth,Hist_CR4,xtitle,tag)

if not Variable == '_XoffVsYoff_':
    bins = []
    for b in range(0,Hist_Data_Depth.GetNbinsX()+1):
        bins += [Hist_Data_Depth.GetBinLowEdge(b+1)]
    Hist_Sensitivity1 = ROOT.TH1D("Hist_Sensitivity1","",len(bins)-1,array('d',bins))
    Hist_Sensitivity2 = ROOT.TH1D("Hist_Sensitivity2","",len(bins)-1,array('d',bins))
    Hist_Sensitivity3 = ROOT.TH1D("Hist_Sensitivity3","",len(bins)-1,array('d',bins))
    Hist_Sensitivity4 = ROOT.TH1D("Hist_Sensitivity4","",len(bins)-1,array('d',bins))
    Hist_Sensitivity5 = ROOT.TH1D("Hist_Sensitivity5","",len(bins)-1,array('d',bins))
    for b in range(0,Hist_CR4.GetNbinsX()):
        if IncludeHeightMethod:
            CR_Height = Hist_CR_Height.GetBinContent(b+1)
            ErrCR_Height = Hist_CR_Height.GetBinError(b+1)
        if IncludeDepthMethod:
            CR_Depth = Hist_CR_Depth.GetBinContent(b+1)
            ErrCR_Depth = Hist_CR_Depth.GetBinError(b+1)
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
        if IncludeDepthMethod:
            if not ErrCR_Depth==0:
                CR4 += CR_Depth/(ErrCR_Depth*ErrCR_Depth)
                Norm4 += 1./(ErrCR_Depth*ErrCR_Depth)
                ErrCR4 += 1./(ErrCR_Depth*ErrCR_Depth)
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
        CR4 = CR4/Norm4
        ErrCR4 = 1./pow(ErrCR4,0.5)
        if IncludeDepthMethod:
            if not CR_Depth==0:
                Hist_Sensitivity1.SetBinContent(b+1,5.*pow(ErrCR_Depth*ErrCR_Depth+CR_Depth,0.5)/CR_Depth)
        if IncludeDepthMethod:
            if not CR_Depth==0:
                Hist_Sensitivity1.SetBinContent(b+1,5.*pow(ErrCR_Depth*ErrCR_Depth+CR_Depth,0.5)/CR_Depth)
        if IncludeMSCWMethod:
            if not CR_MSCW==0:
                Hist_Sensitivity2.SetBinContent(b+1,5.*pow(ErrCR_MSCW*ErrCR_MSCW+CR_MSCW,0.5)/CR_MSCW)
        if IncludeMSCLMethod:
            if not CR_MSCL==0:
                Hist_Sensitivity3.SetBinContent(b+1,5.*pow(ErrCR_MSCL*ErrCR_MSCL+CR_MSCL,0.5)/CR_MSCL)
        if not CR4==0:
            Hist_Sensitivity4.SetBinContent(b+1,5.*pow(ErrCR4*ErrCR4+CR4,0.5)/CR4)
        Data_Depth = Hist_Data_Depth.GetBinContent(b+1)
        ErrData_Depth = Hist_Data_Depth.GetBinError(b+1)
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
    Hist_Sensitivity4.SetLineWidth(2)
    Hist_Sensitivity4.SetMinimum(0)
    Hist_Sensitivity4.Draw()
    if IncludeDepthMethod:
        Hist_Sensitivity1.SetLineColor(2)
        Hist_Sensitivity1.Draw("same")
    if IncludeMSCWMethod:
        Hist_Sensitivity2.SetLineColor(3)
        Hist_Sensitivity2.Draw("same")
    if IncludeMSCLMethod:
        Hist_Sensitivity3.SetLineColor(4)
        Hist_Sensitivity3.Draw("same")
    Hist_Sensitivity5.SetLineColor(6)
    Hist_Sensitivity5.Draw("same")
    Hist_Sensitivity4.Draw("same")
    if Variable == '_Energy_':
        pad1.SetLogx()
        #pad1.SetLogy()
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
        legend.AddEntry(Hist_Sensitivity1,'Height method',"pl")
    if IncludeDepthMethod:
        legend.AddEntry(Hist_Sensitivity1,'Depth method',"pl")
    if IncludeMSCWMethod:
        legend.AddEntry(Hist_Sensitivity2,'MSCW method',"pl")
    if IncludeMSCLMethod:
        legend.AddEntry(Hist_Sensitivity3,'MSCL method',"pl")
    legend.AddEntry(Hist_Sensitivity4,'Combined method',"pl")
    legend.AddEntry(Hist_Sensitivity5,'Some perfect method',"pl")
    legend.Draw("SAME")
    canvas.SaveAs('output_plots/Sensitivity_%s.pdf'%(tag))
