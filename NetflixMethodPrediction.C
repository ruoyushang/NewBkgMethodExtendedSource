//
// 
// This code reads output from EventDisplay and creates a root file which stores histograms of MSCW and MSCL.
// What histograms are created by this code:
//      1. 2D histograms of data events in MSCL and MSCW parameter space (as function of energy).
//      2. 2D histograms of MC electron events (normalized to electron flux) in MSCL and MSCW parameter space.
//      3. 1D theta2 histograms of data events (gamma-ray region and CR region).
//      4. 2D skymap histograms of data events (gamma-ray region and CR region).
// Author: Ruo Yu Shang
//
//
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <utility>

#include <algorithm>
#include <functional>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TBranch.h"
#include "TSpline.h"

ClassImp(TSplinePoly);
ClassImp(TSplinePoly3);
ClassImp(TSplinePoly5);
ClassImp(TSpline3);
ClassImp(TSpline5);
ClassImp(TSpline);

#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;
//using Eigen::MatrixXd;

#include "NetflixMethodGetShowerImage.C"

void NormalizeHist2D(TH2D* hist_ref, TH2D* hist_sub)
{
    int binx_lower = hist_ref->GetXaxis()->FindBin(MSCW_cut_lower);
    int binx_upper = hist_ref->GetXaxis()->FindBin(MSCW_cut_blind);
    int biny_lower = hist_ref->GetYaxis()->FindBin(MSCW_cut_blind);
    int biny_upper = hist_ref->GetYaxis()->FindBin(3.*MSCW_cut_blind);
    double integral_ref = double(hist_ref->Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    double integral_sub = double(hist_sub->Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    double scale = 0.;
    if (integral_sub>0.) scale = integral_ref/integral_sub;
    hist_sub->Scale(scale);
}
MatrixXcd NoiseRemove(MatrixXcd mtx_data)
{
    //return mtx_data;
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            //if (mtx_data(i,j).real()*mtx_data(i,j).real()+mtx_data(i,j).imag()*mtx_data(i,j).imag()>1) 
            if (i>=mtx_data.cols()-2) 
            //if (i>=mtx_data.cols()-30) 
            {
                matrix(i,j) = mtx_data(i,j);
            }
            else 
            {
                matrix(i,j) = 0;
                //matrix(i,j) = mtx_data(i,j);
            }
        }
    }
    return matrix;
}
MatrixXcd replaceEigenVectorMatrix(MatrixXcd mtx_data, MatrixXcd mtx_dark, int binx_blind)
{
    MatrixXcd matrix(mtx_dark.cols(),mtx_dark.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        if (i<mtx_data.cols()-binx_blind)
        {
            matrix.col(i) = mtx_dark.col(i);
        }
        else
        {
            matrix.col(i) = mtx_data.col(i);
        }
    }
    return matrix;
}
MatrixXcd replaceMatrix(MatrixXcd mtx_data_blind, MatrixXcd mtx_data, int binx_lower, int biny_lower, int binx_blind, int biny_blind)
{
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            //matrix(i,j) = mtx_data(i,j);
            if (i>binx_lower && j>biny_lower && i<binx_blind && j<biny_blind)
            {
                matrix(i,j) = mtx_data_blind(i,j);
            }
            else
            {
                matrix(i,j) = mtx_data(i,j);
            }
        }
    }
    return matrix;
}
double GetChi2(MatrixXcd mtx_data, MatrixXcd mtx_bkgd, int binx_blind, int biny_blind)
{
    double chi2 = 0.;
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            if (i<=binx_blind && j<=biny_blind) continue;
            //double data_err = pow(mtx_data(i,j).real(),0.5);
            //double bkgd_err = pow(abs(mtx_bkgd(i,j).real()),0.5);
            //if (data_err*data_err+bkgd_err*bkgd_err>0)
            //{
            //    chi2 += pow(mtx_data(i,j).real()-mtx_bkgd(i,j).real(),2)/(data_err*data_err+bkgd_err*bkgd_err);
            //}
            chi2 += pow(mtx_data(i,j).real()-mtx_bkgd(i,j).real(),2);
        }
    }
    chi2 = 1./chi2;
    return chi2;
}
MatrixXcd MatrixBlind(MatrixXcd mtx_data, int binx_lower, int biny_lower, int binx_blind, int biny_blind)
{
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            matrix(i,j) = mtx_data(i,j);
            if (i>=binx_lower && j>=biny_lower && i<=binx_blind && j<=biny_blind)
            {
                matrix(i,j) = 0;
            }
        }
    }
    return matrix;
}
MatrixXcd MatrixRatio(MatrixXcd mtx_data, MatrixXcd mtx_dark)
{
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            double data_content = pow(mtx_data(i,j).real()*mtx_data(i,j).real()+mtx_data(i,j).imag()*mtx_data(i,j).imag(),0.5);
            double dark_content = pow(mtx_dark(i,j).real()*mtx_dark(i,j).real()+mtx_dark(i,j).imag()*mtx_dark(i,j).imag(),0.5);
            matrix(i,j) = 1;
            //if (dark_content>0)
            //{
            //    matrix(i,j) = data_content/dark_content;
            //}
            if (data_content>0)
            {
                matrix(i,j) = dark_content/data_content;
            }
        }
    }
    return matrix;
}
MatrixXcd MatrixScale(MatrixXcd mtx_ratio, MatrixXcd mtx_bkgd)
{
    MatrixXcd matrix(mtx_ratio.cols(),mtx_ratio.rows());
    for (int i=0;i<mtx_ratio.cols();i++)
    {
        for (int j=0;j<mtx_ratio.rows();j++)
        {
            double ratio_content = mtx_ratio(i,j).real();
            double bkgd_content = mtx_bkgd(i,j).real();
            matrix(i,j) = 0;
            if (ratio_content>0) matrix(i,j) = bkgd_content*ratio_content;
        }
    }
    return matrix;
}
MatrixXcd fillMatrix(TH2D* hist)
{
    MatrixXcd matrix(hist->GetNbinsX(),hist->GetNbinsY());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            matrix(binx,biny) = hist->GetBinContent(binx+1,biny+1);
        }
    }
    return matrix;
}
void fillHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,abs(mtx(binx,biny).real()));
        }
    }
}
MatrixXcd MakeSmoothSplineFunction(MatrixXcd mtx_origin)
{
    TH2D hist_2d = TH2D("hist_2d","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    TH2D hist_2d_smooth_col = TH2D("hist_2d_smooth_col","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    TH2D hist_2d_smooth_row = TH2D("hist_2d_smooth_row","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    fillHistogram(&hist_2d,mtx_origin);
    for (int col=0;col<mtx_origin.cols();col++)
    {
        TH1D* hist = hist_2d.ProjectionY("hist",col+1,col+1);
        double bin_width = hist->GetBinLowEdge(2)-hist->GetBinLowEdge(1);
        double nbins_in_RMS = hist->GetRMS()/bin_width;
        int knot_size = 5;
        TH1D hist_rebin = TH1D("hist_rebin","",mtx_origin.rows(),0,mtx_origin.rows());
        hist_rebin.Reset();
        hist_rebin.Add(hist);
        hist_rebin.Rebin(knot_size);
        TSpline3 spline = TSpline3(&hist_rebin);
        for (int i=0;i<hist->GetNbinsX();i++) {
            double x = hist->GetBinCenter(i+1);
            hist_2d_smooth_col.SetBinContent(col+1,i+1,spline.Eval(x)/double(knot_size));
        }
    }
    for (int row=0;row<mtx_origin.rows();row++)
    {
        TH1D* hist = hist_2d.ProjectionX("hist",row+1,row+1);
        double bin_width = hist->GetBinLowEdge(2)-hist->GetBinLowEdge(1);
        double nbins_in_RMS = hist->GetRMS()/bin_width;
        int knot_size = 5;
        TH1D hist_rebin = TH1D("hist_rebin","",mtx_origin.cols(),0,mtx_origin.cols());
        hist_rebin.Reset();
        hist_rebin.Add(hist);
        hist_rebin.Rebin(knot_size);
        TSpline3 spline = TSpline3(&hist_rebin);
        for (int i=0;i<hist->GetNbinsX();i++) {
            double x = hist->GetBinCenter(i+1);
            hist_2d_smooth_row.SetBinContent(i+1,row+1,spline.Eval(x)/double(knot_size));
        }
    }
    for (int col=0;col<mtx_origin.cols();col++)
    {
        for (int row=0;row<mtx_origin.rows();row++)
        {
            hist_2d.SetBinContent(col+1,row+1,0.5*(hist_2d_smooth_col.GetBinContent(col+1,row+1)+hist_2d_smooth_row.GetBinContent(col+1,row+1)));
        }
    }
    MatrixXcd matrix(mtx_origin.cols(),mtx_origin.rows());
    matrix = fillMatrix(&hist_2d);
    //std::cout << "mtx_origin.col(0) = " << mtx_origin.col(0) << std::endl;
    //std::cout << "matrix.col(0) = " << matrix.col(0) << std::endl;
    return matrix;
}
MatrixXcd EigenvaluePrediction(MatrixXcd mtx_data,MatrixXcd mtx_data_blind, MatrixXcd mtx_dark, int binx_lower, int biny_lower, int binx_blind, int biny_blind)
{

    MatrixXcd mtx_bkgd(mtx_dark.cols(),mtx_dark.rows());

    mtx_data_blind = replaceMatrix(mtx_data_blind,mtx_data,binx_lower,biny_lower,binx_blind,biny_blind);
    //mtx_data_blind = MakeSmoothSplineFunction(mtx_data_blind);

    ComplexEigenSolver<MatrixXcd> eigensolver_data_blind(mtx_data_blind);
    ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);

    MatrixXcd mtx_eigenvalue_data_blind(mtx_dark.cols(),mtx_dark.rows());
    MatrixXcd mtx_eigenvalue_dark(mtx_dark.cols(),mtx_dark.rows());
    MatrixXcd mtx_eigenvector_data_blind(mtx_dark.cols(),mtx_dark.rows());
    MatrixXcd mtx_eigenvector_dark(mtx_dark.cols(),mtx_dark.rows());
    for (int i=0;i<mtx_dark.cols();i++)
    {
        for (int j=0;j<mtx_dark.rows();j++)
        {
            if (i==j) 
            {
                mtx_eigenvalue_data_blind(i,j) = eigensolver_data_blind.eigenvalues()(i);
                mtx_eigenvalue_dark(i,j) = eigensolver_dark.eigenvalues()(i);
            }
            else
            {
                mtx_eigenvalue_data_blind(i,j) = 0;
                mtx_eigenvalue_dark(i,j) = 0;
            }
        }
    }
    mtx_eigenvector_data_blind = eigensolver_data_blind.eigenvectors();
    mtx_eigenvector_dark = eigensolver_dark.eigenvectors();

    mtx_eigenvalue_data_blind = NoiseRemove(mtx_eigenvalue_data_blind);
    mtx_eigenvalue_dark = NoiseRemove(mtx_eigenvalue_dark);
    
    MatrixXcd mtx_eigenvalue_dark_modify(mtx_dark.cols(),mtx_dark.rows());
    mtx_eigenvalue_dark_modify = replaceEigenVectorMatrix(mtx_eigenvalue_data_blind, mtx_eigenvalue_dark, 2);
    MatrixXcd mtx_eigenvector_dark_modify(mtx_dark.cols(),mtx_dark.rows());
    mtx_eigenvector_dark_modify = replaceEigenVectorMatrix(mtx_eigenvector_data_blind, mtx_eigenvector_dark, 2);

    //mtx_bkgd = mtx_eigenvector_dark_modify*mtx_eigenvalue_dark*mtx_eigenvector_dark_modify.inverse();
    mtx_bkgd = mtx_eigenvector_data_blind*mtx_eigenvalue_dark*mtx_eigenvector_data_blind.inverse();

    TH2D Hist_Data_Temp("Hist_Data_Temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D Hist_Bkgd_Temp("Hist_Bkgd_Temp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fillHistogram(&Hist_Data_Temp,mtx_data);
    fillHistogram(&Hist_Bkgd_Temp,mtx_bkgd);
    NormalizeHist2D(&Hist_Data_Temp,&Hist_Bkgd_Temp);
    mtx_bkgd = fillMatrix(&Hist_Bkgd_Temp);

    //mtx_bkgd = MakeSmoothSplineFunction(mtx_bkgd);

    return mtx_bkgd;
}
void NetflixMethodPrediction(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, double theta2_cut_lower_input, double theta2_cut_upper_input)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    Theta2_cut_lower = theta2_cut_lower_input;
    Theta2_cut_upper = theta2_cut_upper_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;

    vector<TH2D> Hist_Bkgd_MSCLW;
    TFile InputDataFile("output_Jul16/Netflix_"+TString(target_data)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root");
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        TString filename_data  = "Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_dark  = "Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
        Hist_Bkgd_MSCLW.push_back(TH2D("Hist_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));

        //if (energy_bins[e]<300.) continue;
        //if (energy_bins[e]>350.) continue;

        MatrixXcd mtx_data(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_data_blind(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark_blind(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_data_bkgd(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark_bkgd(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        mtx_data = fillMatrix(Hist_Data);
        mtx_dark = fillMatrix(Hist_Dark);
        int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCL_cut_lower);
        int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_lower);
        int binx_blind = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int biny_blind = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
        mtx_data_blind = MatrixBlind(mtx_data, binx_lower, biny_lower, binx_blind, biny_blind);
        mtx_dark_blind = MatrixBlind(mtx_dark, binx_lower, biny_lower, binx_blind, biny_blind);

        //mtx_dark = MakeSmoothSplineFunction(mtx_dark);

        //initial condition
        //mtx_data_bkgd = mtx_data_blind;
        mtx_data_bkgd = mtx_dark;

        std::cout << "Hist_Dark->Integral() = " << Hist_Dark->Integral() << std::endl;
        for (int iter=0;iter<5;iter++)
        {
            mtx_data_bkgd = EigenvaluePrediction(mtx_data,mtx_data_bkgd,mtx_dark,binx_lower,biny_lower,binx_blind,biny_blind);
            fillHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
            std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;
        }

        fillHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);

    }
    InputDataFile.Close();


    TFile OutputFile("output_Jul16/Netflix_"+TString(target_data)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Bkgd_MSCLW.at(e).Write();
    }
    OutputFile.Close();

}
