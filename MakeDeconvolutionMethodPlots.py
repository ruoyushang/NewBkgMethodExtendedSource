
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")

folder = 'output'

source = ''
source_list = []
#source_list  += ['Crab']
#source_list  += ['2ndCrab']
source_list  += ['PKS1424']
#source_list  += ['H1426']
#source_list  += ['3C264']
#source_list  += ['Ton599']
#source_list  += ['IC443']
#source_list  += ['BrandonValidation']

Region = 'SR'

global MSCW_lower_cut
global MSCW_upper_cut
global MSCL_lower_cut
global MSCL_upper_cut
global Depth_cut_width

Elev_lower_cut = 55
Elev_upper_cut = 85
#Elev_lower_cut = 75
#Elev_upper_cut = 85
#Elev_lower_cut = 35
#Elev_upper_cut = 55
#Elev_lower_cut = 35
#Elev_upper_cut = 85

UseMethod1 = False
UseMethod2 = True
UseMethod3 = True

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
#energy_list += [150]
energy_list += [200]
energy_list += [250]
energy_list += [300]
energy_list += [400]
energy_list += [600]
energy_list += [800]
energy_list += [1000]
energy_list += [1200]
energy_list += [1600]
energy_list += [2000]
energy_list += [3000]
energy_list += [4000]
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

    FilePath = '%s/Deconvolution_%s_Elev%sto%s_%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut,Region)
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
            Hists[h].GetXaxis().SetRangeUser(-1.0,20.5)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")

    power1 = 2
    power2 = 2
    func = ROOT.TF1("func","[2]*pow(x+[1],%s)*exp(-pow((x+[1])/[0],%s))"%(power1,power2), -1.0, 20.5)
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
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f#pm%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/Aux_%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))

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
            Hists[h].GetXaxis().SetRangeUser(-1.0,20.5)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")

    func = ROOT.TF1("func","[2]*pow(x+[1],%s)*exp(-pow((x+[1])/[0],%s))"%(power1,power2), -1.0, 20.5)
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
    func_result = ROOT.TF1("func_result","[0]*pow(x+%0.3f,%s)*exp(-pow((x+%0.3f)/%0.3f,%s))"%(mean[3],power1,mean[3],rms[3],power2), -1.0, 20.5)
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
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f#pm%0.2f)"%(mean[h],rms[h]),"pl")
    legend.Draw("SAME")
    c_both.SaveAs('output_plots/Sig_%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))

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
            #Hists[h].GetXaxis().SetRangeUser(-0.5,2.5)
            #Hists[h].GetXaxis().SetRangeUser(-0.5,3.5)
            #Hists[h].GetXaxis().SetRangeUser(-0.5,5.5)
            #Hists[h].GetXaxis().SetRangeUser(-0.5,10.5)
            Hists[h].GetXaxis().SetRangeUser(-1.0,9.0)
            if title=='Depth': 
                Hists[h].GetXaxis().SetRangeUser(0,20)
            if title=='EmissionHeight': 
                Hists[h].Rebin(4)
                Hists[h].GetXaxis().SetRangeUser(5,15)
            if legends[h]=='Ring':
                Hists[h].SetLineWidth(3)
            if max_heigh < Hists[h].GetMaximum(): 
                max_heigh = Hists[h].GetMaximum()
                max_hist = h

    Hists[max_hist].SetMinimum(0)
    Hists[max_hist].Draw("E")
    #Hists[0].SetMinimum(0)
    #Hists[0].Draw("E")

    if doSum:
        Hist_Sum = Hists[2].Clone()
        Hist_Sum.Add(Hists[1])
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
            legend.AddEntry(Hists[h],legends[h]+"(%0.2f#pm%0.2f)"%(mean[h],rms[h]),"pl")
            #legend.AddEntry(Hists[h],legends[h],"pl")
    if doSum:
        legend.AddEntry(Hist_Sum,'total bkg',"f")
    legend.Draw("SAME")
    lumilab1 = ROOT.TLatex(0.15,0.75,'   ' )
    lumilab1.SetNDC()
    lumilab1.SetTextSize(0.15)
    lumilab1.Draw()
    if not title=='MSCW':
        lumilab2 = ROOT.TLatex(0.15,0.65,'%0.2f < MSCW < %0.2f'%(MSCW_lower_cut,MSCW_upper_cut) )
        lumilab2.SetNDC()
        lumilab2.SetTextSize(0.15)
        lumilab2.Draw()
    if not title=='MSCL':
        lumilab5 = ROOT.TLatex(0.15,0.50,'%0.2f < MSCL < %0.2f'%(MSCL_lower_cut,MSCL_upper_cut) )
        lumilab5.SetNDC()
        lumilab5.SetTextSize(0.15)
        lumilab5.Draw()
    if not title=='Depth':
        lumilab4 = ROOT.TLatex(0.15,0.35,'Slant depth (central #pm %0.2f)'%(Depth_cut_width) )
        lumilab4.SetNDC()
        lumilab4.SetTextSize(0.15)
        lumilab4.Draw()
    lumilab5 = ROOT.TLatex(0.15,0.20,'E >%0.1f GeV'%(ErecS_lower_cut) )
    lumilab5.SetNDC()
    lumilab5.SetTextSize(0.15)
    lumilab5.Draw()
    #pad1.SetLogy()
    c_both.SaveAs('output_plots/%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))

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
            Sig = (NSR-NBkg)/(pow(NSR_Err*NSR_Err+NBkg_Err*NBkg_Err,0.5))
            Hist_Data.SetBinContent(bx+1,by+1,Sig)
    Hist_Data.GetYaxis().SetTitle(ytitle)
    Hist_Data.GetXaxis().SetTitle(xtitle)
    Hist_Data.Draw("COL4Z")
    canvas.SaveAs('output_plots/%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))
    Hist_Sig = ROOT.TH1D("Hist_Sig","",100,-5,5)
    for bx in range(0,Hist_Data.GetNbinsX()):
        for by in range(0,Hist_Data.GetNbinsY()):
            if not Hist_SR.GetBinContent(bx+1,by+1)==0:
                content = Hist_Data.GetBinContent(bx+1,by+1)/Hist_SR.GetBinError(bx+1,by+1)
                Hist_Sig.Fill(content)
    Hist_Sig.Draw()
    canvas.SaveAs('output_plots/Sig_%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))

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
    canvas.SaveAs('output_plots/%s_%s_Elev%sto%s_%s.pdf'%(name,source,Elev_lower_cut,Elev_upper_cut,Region))

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
    FilePath = '%s/Deconvolution_%s_Elev%sto%s_%s.root'%(folder,source,Elev_lower_cut,Elev_upper_cut,Region)
    #print 'Read %s'%(FilePath)
    TargetFile=ROOT.TFile(FilePath)
    InfoTree = TargetFile.Get("InfoTree")
    InfoTree.GetEntry(0)
    Depth_cut_width = InfoTree.Depth_cut_width
    MSCW_lower_cut = InfoTree.MSCW_cut_lower
    MSCW_upper_cut = InfoTree.MSCW_cut_upper
    MSCL_lower_cut = InfoTree.MSCL_cut_lower
    MSCL_upper_cut = InfoTree.MSCL_cut_upper
    #Hist_Target_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TelElevAzim')
    #Hist_Dark_TelElevAzim = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TelElevAzim')
    #Make2DTrajectoryPlot(Hist_Dark_TelElevAzim,Hist_Target_TelElevAzim,'Tel. elev.','Tel. azim.','TelElevAzim')
    Hist_Target_ErecSDepth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ON_ErecSDepth')
    #Make2DProjectionPlot(Hist_Target_ErecSDepth,'E [GeV]','Slant depth [37g/cm^{2}]','ErecS_vs_SlantDepth')
    Hist_Target_ElevDepth = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ON_ElevDepth')
    #Make2DProjectionPlot(Hist_Target_ElevDepth,'Tel. elev.','Slant depth - d(E) [37g/cm^{2}]','Elev_vs_SlantDepth')
    for e in range(0,len(energy_list)-1):
        ErecS_lower_cut = energy_list[e]
        ErecS_upper_cut = energy_list[e+1]

        which_method = 'MSCW'

        Hist_Target_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_MSCW')
        Hist_Target_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_MSCW')
        Hist_Target_TotalBkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TotalBkg_MSCW')
        Hist_Target_ASR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ASR_MSCW')
        Hist_Target_ACR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ACR_MSCW')
        Hist_Target_ABkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_ABkg_MSCW')
        Hist_Target_Deconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Deconv_MSCW')
        Hist_Target_TrueDeconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_TrueDeconv_MSCW')
        Hist_Target_CR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_CR_MSCW')
        Hist_Target_Elec_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Elec_MSCW')
        Hist_Target_Ring_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Ring_MSCW')
        Hist_Dark_SR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_MSCW')
        Hist_Dark_Bkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Bkg_MSCW')
        Hist_Dark_ASR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_ASR_MSCW')
        Hist_Dark_ACR_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_ACR_MSCW')
        Hist_Dark_ABkg_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_ABkg_MSCW')
        Hist_Dark_Deconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Deconv_MSCW')
        Hist_Dark_TrueDeconv_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_TrueDeconv_MSCW')
        Hist_Dark_Elec_MSCW = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_Elec_MSCW')
        Hist_Target_SR_EmissionHeight = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_EmissionHeight')
        Hist_Dark_SR_EmissionHeight = SelectDiagnosticaHistograms(folder,'MSCW','SR','Dark_SR_EmissionHeight')

        Hist_Target_SR_SkyMap = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_SkyMap')
        Hist_Target_Bkg_SkyMap = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_SkyMap')
        Make2DSignificancePlot(Hist_Target_SR_SkyMap,Hist_Target_Bkg_SkyMap,"RA","Dec","SkyMap_E%s"%(ErecS_lower_cut))

        Hist_Target_SR_Theta2 = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_SR_Theta2')
        Hist_Target_Bkg_Theta2 = SelectDiagnosticaHistograms(folder,'MSCW','SR','Target_Bkg_Theta2')

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_SR_MSCW]
        if source=='2ndCrab':
            legends += ['Crab 2016 ON']
        else:
            legends += ['%s'%(source)]
        colors += [1]
        Hists += [Hist_Target_Bkg_MSCW]
        legends += ['deconvolution']
        colors += [4]
        Hists += [Hist_Target_Elec_MSCW]
        legends += ['residual']
        colors += [3]
        Hists += [Hist_Dark_SR_MSCW]
        legends += ['Crab 2017-18 OFF']
        colors += [6]
        Hists += [Hist_Target_Ring_MSCW]
        legends += ['Ring']
        colors += [2]
        #plotname = 'Target_Dark_MSCW_E%s'%(ErecS_lower_cut)
        #plotname = 'Target_NoElec_MSCW_E%s'%(ErecS_lower_cut)
        plotname = 'Target_Deconv_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        #MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)
        MakeDiagnosticPlot(Hists,legends,colors,title,plotname,True,False)

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_ASR_MSCW]
        legends += ['CR1 (MSCL<1)']
        colors += [2]
        Hists += [Hist_Target_ACR_MSCW]
        legends += ['CR2 (MSCL>1)']
        colors += [3]
        Hists += [Hist_Target_ABkg_MSCW]
        legends += ['Deconv.']
        colors += [4]
        plotname = 'Target_ASR_MSCW_E%s'%(ErecS_lower_cut)
        title = 'MSCW'
        MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)

        #Hists = []
        #legends = []
        #colors = []
        #Hists += [Hist_Target_SR_MSCW]
        #legends += ['SR (MSCL<1)']
        #colors += [2]
        #Hists += [Hist_Target_CR_MSCW]
        #legends += ['CR3 (MSCL>1)']
        #colors += [3]
        #Hists += [Hist_Target_Bkg_MSCW]
        #legends += ['Deconv.']
        #colors += [4]
        #plotname = 'Target_SR_MSCW_E%s'%(ErecS_lower_cut)
        #title = 'MSCW'
        #MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)

        Hists = []
        legends = []
        colors = []
        Hists += [Hist_Target_SR_Theta2]
        legends += ['SR']
        colors += [4]
        Hists += [Hist_Target_Bkg_Theta2]
        legends += ['Bkg']
        colors += [2]
        plotname = 'Target_Theta2_E%s'%(ErecS_lower_cut)
        title = '#theta^{2}'
        MakeDiagnosticPlot(Hists,legends,colors,title,plotname,False,False)


        #Hists = []
        #legends = []
        #colors = []
        #Hists += [Hist_Target_ACR_MSCW]
        #legends += ['CR2 (MSCL>1)']
        #colors += [4]
        #Hists += [Hist_Target_ASR_MSCW]
        #legends += ['CR1 (MSCL<1)']
        #colors += [2]
        #Hists += [Hist_Target_CR_MSCW]
        #legends += ['CR3 (MSCL>1)']
        #colors += [4]
        #Hists += [Hist_Target_SR_MSCW]
        #legends += ['SR (MSCL<1)']
        #colors += [2]
        #plotname = 'Target_SR_MSCW_E%s'%(ErecS_lower_cut)
        #title = 'MSCW'
        #MakeReyleighPlot(Hists,legends,colors,title,plotname,False,False)
