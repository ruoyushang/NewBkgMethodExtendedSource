
import sys,ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetPaintTextFormat("0.3f")

source  = 'IC443'
#source  = 'PKS1424'
#method = 'EmissionHeight'
method = 'MSCW'

FilePath = 'output/Histograms_%s_%s.root'%(source,method)
print 'Read %s'%(FilePath)
InputFile=ROOT.TFile(FilePath)

Elev_lower_cut = 70
Elev_upper_cut = 80
Azim_lower_cut = 40
Azim_upper_cut = 60
Erec_lower_cut = 500
Erec_upper_cut = 1000
MSCL_lower_cut = -2
MSCL_upper_cut = 0.7
MSCW_lower_cut = -2
MSCW_upper_cut = 0.7

Variable = '_theta2_'

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
    if title=='Erec':
        pad1.SetLogy()
        pad1.SetLogx()
        pad2.SetLogx()
    c_both.SaveAs('output_plots/%s.pdf'%(name))
    pad1.SetLogy(0)
    pad1.SetLogx(0)
    pad2.SetLogx(0)


HistList = InputFile.GetListOfKeys()

tag = Variable.strip('_')

print 'HistList.GetSize() = %s'%(HistList.GetSize())
Hist_Data = ROOT.TH1D("Hist_Data","",1,0,1)
Hist_CR = ROOT.TH1D("Hist_CR","",1,0,1)
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
    if 'Data' in name: 
        if Hist_Data.Integral()==0:
            Hist_Data.Delete()
            Hist_Data = hist.Clone()
        else:
            Hist_Data.Add(hist)
    else: 
        if Hist_CR.Integral()==0:
            Hist_CR.Delete()
            Hist_CR = hist.Clone()
        else:
            Hist_CR.Add(hist)

MakeStackPlot(Hist_Data,Hist_CR,'E [TeV]',tag)

