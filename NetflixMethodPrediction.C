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


#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;
//using Eigen::MatrixXd;

#include "NetflixMethodGetShowerImage.C"

void NormalizeHist2D(TH2D* hist_ref, TH2D* hist_sub)
{
    int binx_lower = hist_ref->GetXaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_ref->GetXaxis()->FindBin(3.*MSCW_cut_blind);
    int biny_lower = hist_ref->GetYaxis()->FindBin(MSCW_cut_blind);
    int biny_upper = hist_ref->GetYaxis()->FindBin(3.*MSCW_cut_blind);
    double integral_ref = double(hist_ref->Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    double integral_sub = double(hist_sub->Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    double scale = 0.;
    if (integral_sub>0.) scale = integral_ref/integral_sub;
    hist_sub->Scale(scale);
}
MatrixXcd replaceMatrix(MatrixXcd mtx_data, MatrixXcd mtx_dark, int binx_blind, int biny_blind)
{
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            //matrix(i,j) = mtx_data(i,j);
            if (i<binx_blind && j<biny_blind)
            {
                matrix(i,j) = mtx_dark(i,j);
            }
            else
            {
                matrix(i,j) = mtx_data(i,j);
            }
        }
    }
    return matrix;
}
MatrixXcd MatrixRatio(MatrixXcd mtx_data, MatrixXcd mtx_bkgd, int binx_blind, int biny_blind)
{
    MatrixXcd matrix(mtx_data.cols(),mtx_data.rows());
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            double data_content = mtx_data(i,j).real();
            double bkgd_content = mtx_bkgd(i,j).real();
            matrix(i,j) = 0;
            if (i<binx_blind && j<biny_blind)
            {
                if (data_content>0) matrix(i,j) = bkgd_content/data_content;
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
MatrixXcd EigenvaluePrediction(MatrixXcd mtx_data, MatrixXcd mtx_bkgd, MatrixXcd mtx_eigenvalue, int binx_lower, int biny_lower)
{
    MatrixXcd mtx_ratio(mtx_data.cols(),mtx_data.rows());
    mtx_bkgd = replaceMatrix(mtx_data,mtx_bkgd,binx_lower,biny_lower);
    ComplexEigenSolver<MatrixXcd> eigensolver_bkgd(mtx_bkgd);
    mtx_bkgd = eigensolver_bkgd.eigenvectors()*mtx_eigenvalue*eigensolver_bkgd.eigenvectors().inverse();
    //for (int i=0;i<mtx_data.cols();i++)
    //{
    //    for (int j=0;j<mtx_data.rows();j++)
    //    {
    //        if (mtx_bkgd(i,j).real()<0) mtx_bkgd(i,j) = 0;
    //    }
    //}
    mtx_ratio = MatrixRatio(mtx_data,mtx_bkgd,binx_lower,biny_lower);
    return mtx_bkgd;
}
MatrixXcd EigenvectorPrediction(MatrixXcd mtx_data, MatrixXcd mtx_bkgd, MatrixXcd mtx_eigenvalue, MatrixXcd mtx_eigenvector, int binx_lower, int biny_lower)
{
    mtx_bkgd = replaceMatrix(mtx_data,mtx_bkgd,binx_lower,biny_lower);
    std::cout << "mtx_data :\n" << mtx_data.real() << std::endl;
    std::cout << "mtx_bkgd :\n" << mtx_bkgd.real() << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_bkgd(mtx_bkgd);
    for (int i=0;i<mtx_data.cols();i++)
    {
        for (int j=0;j<mtx_data.rows();j++)
        {
            if (i==j) mtx_eigenvalue(i,j) = eigensolver_bkgd.eigenvalues()(i);
            else mtx_eigenvalue(i,j) = 0;
        }
    }
    mtx_bkgd = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector.inverse();
    std::cout << "mtx_bkgd :\n" << mtx_bkgd.real() << std::endl;
    return mtx_bkgd;
}
void NetflixMethodPrediction(string target_data, double theta2_cut_lower_input, double theta2_cut_upper_input)
{

    Theta2_cut_lower = theta2_cut_lower_input;
    Theta2_cut_upper = theta2_cut_upper_input;

    vector<TH2D> Hist_Data_MSCLW;
    vector<TH2D> Hist_Dark_MSCLW;
    TFile InputDataFile("output_Jul05/Netflix_"+TString(target_data)+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root");
    for (int e=0;e<N_energy_bins;e++) 
    {
        if (energy_bins[e]<200.) continue;
        if (energy_bins[e]>250.) continue;
        std::cout << "energy = " << energy_bins[e] << std::endl;
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        TString filename_data  = "Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_dark  = "Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);

        NormalizeHist2D(Hist_Data,Hist_Dark);
        MatrixXcd mtx_data(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_data_bkgd(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark_bkgd(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark_corr(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_data_corr(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_eigenvalue(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_eigenvector(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        mtx_data = fillMatrix(Hist_Data);
        mtx_dark = fillMatrix(Hist_Dark);
        int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;

        ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);
        std::cout << "The eigenvalues of mtx_dark are:\n" << eigensolver_dark.eigenvalues() << std::endl;
        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j) mtx_eigenvalue(i,j) = eigensolver_dark.eigenvalues()(i);
                else mtx_eigenvalue(i,j) = 0;
                mtx_data_bkgd(i,j) = 0;
                mtx_dark_bkgd(i,j) = 0;
                mtx_data_corr(i,j) = 0;
                mtx_dark_corr(i,j) = 0;
            }
        }
        mtx_eigenvector = eigensolver_dark.eigenvectors();

        //mtx_data_bkgd = mtx_dark;

        for (int iter=0;iter<5;iter++)
        {
            std::cout << "+++++++++++++++++++ iteration " << iter << std::endl;
            mtx_dark_bkgd = EigenvaluePrediction(mtx_dark,mtx_dark_bkgd,mtx_eigenvalue,binx_lower,biny_lower);
        }
        mtx_dark_corr = MatrixRatio(mtx_dark_bkgd,mtx_dark,binx_lower,biny_lower);
        for (int iter=0;iter<5;iter++)
        {
            std::cout << "+++++++++++++++++++ iteration " << iter << std::endl;
            mtx_data_bkgd = EigenvaluePrediction(mtx_data,mtx_data_bkgd,mtx_eigenvalue,binx_lower,biny_lower);
        }
        mtx_data_bkgd = MatrixScale(mtx_dark_corr,mtx_data_bkgd);
        std::cout << "mtx_data :\n" << mtx_data.block(0,0,binx_lower,biny_lower).real() << std::endl;
        std::cout << "mtx_data_bkgd :\n" << mtx_data_bkgd.block(0,0,binx_lower,biny_lower).real() << std::endl;
        mtx_data_corr = MatrixRatio(mtx_data,mtx_data_bkgd,binx_lower,biny_lower);
        std::cout << "mtx_dark_corr :\n" << mtx_dark_corr.block(0,0,binx_lower,biny_lower).real() << std::endl;
        std::cout << "mtx_data_corr :\n" << mtx_data_corr.block(0,0,binx_lower,biny_lower).real() << std::endl;


    }
    InputDataFile.Close();


    //TFile OutputFile("output_Jul05/Netflix_"+TString(target)+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","update");
    //for (int e=0;e<N_energy_bins;e++)
    //{
    //    Hist_Data_MSCLW.at(e).Write();
    //    Hist_Dark_MSCLW.at(e).Write();
    //}
    //OutputFile.Close();

}
