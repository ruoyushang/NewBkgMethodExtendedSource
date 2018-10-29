
import sys,ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output_Oct28'
MSCW_lower_cut = -2
MSCW_upper_cut = 0.7
#folder = 'output_Oct29'
#MSCW_lower_cut = 2
#MSCW_upper_cut = 3

source  = '2ndCrab'
#source  = 'PKS1424'
#source  = '3C264'
#source  = 'H1426'
#source  = 'IC443'
#source  = 'Ton599'

Elev_lower_cut = 60
Elev_upper_cut = 90
Azim_lower_cut = 0
Azim_upper_cut = 80
Erec_lower_cut = 100
Erec_upper_cut = 10000
MSCL_lower_cut = -2
MSCL_upper_cut = 0.7

#Variable = '_Energy_'
#xtitle = 'E [TeV]'
#Variable = '_theta2_'
#xtitle = '#theta^{2}'
Variable = '_theta2zoomin_'
xtitle = '#theta^{2}'
#Variable = '_XoffVsYoff_'

def set_histStyle( hist , color):
    hist.SetFillColor(color)
    #hist.SetLineColor(1)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetFillStyle(1001)
    pass

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


def SelectHistograms(method,isData):

    Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)

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
        Erec_lower = 0
        Erec_upper = 0
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
            if 'Erec' in element and 'to' in element:
                element = element.strip('Erec')
                Erec_lower = float(element.split('to')[0])
                Erec_upper = float(element.split('to')[1])
        if Elev_lower<Elev_lower_cut: continue
        if Elev_upper>Elev_upper_cut: continue
        if Azim_lower<Azim_lower_cut: continue
        if Azim_upper>Azim_upper_cut: continue
        if Erec_lower<Erec_lower_cut: continue
        if Erec_upper>Erec_upper_cut: continue
        if isData in name: 
            if Hist_Data.Integral()==0:
                Hist_Data.Delete()
                Hist_Data = hist.Clone()
            else:
                Hist_Data.Add(hist)

    return Hist_Data

Hist_Data1 = SelectHistograms('Height','Data')
Hist_CR1 = SelectHistograms('Height','CR')
Hist_Data2 = SelectHistograms('MSCW','Data')
Hist_CR2 = SelectHistograms('MSCW','CR')

print 'Hist_Data1.Integral() = %s'%(Hist_Data1.Integral())
print 'Hist_Data2.Integral() = %s'%(Hist_Data2.Integral())

tag = Variable.strip('_')
tag += '_%s'%(source)
tag += '_%s'%('EmissionHeight')
if MSCW_lower_cut>1:
    tag += '_Validation'
MakeStackPlot(Hist_Data1,Hist_CR1,xtitle,tag)

tag = Variable.strip('_')
tag += '_%s'%(source)
tag += '_%s'%('MSCW')
if MSCW_lower_cut>1:
    tag += '_Validation'
MakeStackPlot(Hist_Data2,Hist_CR2,xtitle,tag)

Hist_CR3 = Hist_CR1.Clone()
for b in range(0,Hist_CR3.GetNbinsX()):
    CR1 = Hist_CR1.GetBinContent(b+1)
    ErrCR1 = Hist_CR1.GetBinError(b+1)
    CR2 = Hist_CR2.GetBinContent(b+1)
    ErrCR2 = Hist_CR2.GetBinError(b+1)
    CR3 = 0
    ErrCR3 = 0
    if not ErrCR1==0 and not ErrCR2==0:
        CR3 = (CR1/(ErrCR1*ErrCR1) + CR2/(ErrCR2*ErrCR2))/(1./(ErrCR1*ErrCR1) + 1./(ErrCR2*ErrCR2))
        ErrCR3 = 1./pow((pow(1./ErrCR1,2) + pow(1./ErrCR2,2)),0.5)
        Hist_CR3.SetBinContent(b+1,CR3)
        Hist_CR3.SetBinError(b+1,ErrCR3)
tag = Variable.strip('_')
tag += '_%s'%(source)
tag += '_%s'%('Combined')
if MSCW_lower_cut>1:
    tag += '_Validation'
MakeStackPlot(Hist_Data1,Hist_CR3,xtitle,tag)


Hist_Sensitivity1 = Hist_Data1.Clone()
Hist_Sensitivity2 = Hist_Data1.Clone()
Hist_Sensitivity3 = Hist_Data1.Clone()
Hist_Sensitivity1.Reset()
Hist_Sensitivity2.Reset()
Hist_Sensitivity3.Reset()
for b in range(0,Hist_CR3.GetNbinsX()):
    CR1 = Hist_CR1.GetBinContent(b+1)
    ErrCR1 = Hist_CR1.GetBinError(b+1)
    CR2 = Hist_CR2.GetBinContent(b+1)
    ErrCR2 = Hist_CR2.GetBinError(b+1)
    CR3 = 0
    ErrCR3 = 0
    if not ErrCR1==0 and not ErrCR2==0:
        CR3 = (CR1/(ErrCR1*ErrCR1) + CR2/(ErrCR2*ErrCR2))/(1./(ErrCR1*ErrCR1) + 1./(ErrCR2*ErrCR2))
        ErrCR3 = 1./pow((pow(1./ErrCR1,2) + pow(1./ErrCR2,2)),0.5)
    if not CR1==0:
        Hist_Sensitivity1.SetBinContent(b+1,5.*pow(ErrCR1*ErrCR1+CR1,0.5)/CR1)
    if not CR2==0:
        Hist_Sensitivity2.SetBinContent(b+1,5.*pow(ErrCR2*ErrCR2+CR2,0.5)/CR2)
    if not CR3==0:
        Hist_Sensitivity3.SetBinContent(b+1,5.*pow(ErrCR3*ErrCR3+CR3,0.5)/CR3)

canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
pad1.SetBottomMargin(0.15)
pad1.SetRightMargin(0.15)
pad1.SetLeftMargin(0.15)
pad1.SetTopMargin(0.15)
pad1.SetBorderMode(0)
pad1.Draw()
pad1.cd()
Hist_Sensitivity3.GetYaxis().SetTitle('Excluded signal flux (Bg unit.)')
Hist_Sensitivity3.GetXaxis().SetTitle(xtitle)
Hist_Sensitivity3.Draw()
Hist_Sensitivity2.SetLineColor(2)
Hist_Sensitivity2.Draw("same")
Hist_Sensitivity1.SetLineColor(3)
Hist_Sensitivity1.Draw("same")
Hist_Sensitivity3.Draw("same")
if Variable == '_Energy_':
    pad1.SetLogx()
canvas.SaveAs('output_plots/Sensitivity_%s.pdf'%(tag))
