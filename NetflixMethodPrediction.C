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
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

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
//using Eigen::MatrixXcd;

#include "NetflixMethodGetShowerImage.C"

MatrixXcd mtx_data(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark(N_bins_for_deconv,N_bins_for_deconv);

void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,abs(mtx(binx,biny)));
        }
    }
}
void fill1DHistogram(TH1D* hist,VectorXcd vtr)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        hist->SetBinContent(binx+1,vtr(binx).real());
        hist->SetBinError(binx+1,0.1*abs(vtr(binx).real()));
    }
}
double BlindedChi2(TH2D* hist_data, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>=binx_blind || by>=biny_blind)
            {
                double data = hist_data->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                if (data*data+model*model>0)
                {
                    chi2 += pow(data-model,2)/(data*data+model*model);
                }
            }
        }
    }
    return chi2;
}

double Chi2Function(const double *par)
{

    // build eigenvector matrix
    MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
    MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
    MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }

    int n_fourier_modes = 6;
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    // 1st eigenvector
    col_fix = N_bins_for_deconv-1;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 0*(1+2*n_fourier_modes);
        mtx_eigenvector(row,col_fix) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    row_fix = N_bins_for_deconv-1;
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 1*(1+2*n_fourier_modes);
        mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    // 2nd eigenvector
    col_fix = N_bins_for_deconv-2;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 2*(1+2*n_fourier_modes);
        mtx_eigenvector(row,col_fix) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    row_fix = N_bins_for_deconv-2;
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 3*(1+2*n_fourier_modes);
        mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    // 3rd eigenvector
    col_fix = N_bins_for_deconv-3;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 4*(1+2*n_fourier_modes);
        mtx_eigenvector(row,col_fix) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    row_fix = N_bins_for_deconv-3;
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
        first_index = 5*(1+2*n_fourier_modes);
        mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
            mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
        }
    }
    // build eigenvalue matrix
    first_index = 6*(1+2*n_fourier_modes);
    mtx_eigenvalue(N_bins_for_deconv-1,N_bins_for_deconv-1) = par[first_index];
    mtx_eigenvalue(N_bins_for_deconv-2,N_bins_for_deconv-2) = par[first_index+1];
    mtx_eigenvalue(N_bins_for_deconv-3,N_bins_for_deconv-3) = par[first_index+2];
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_model = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&hist_data,mtx_data);
    fill2DHistogram(&hist_model,mtx_model);
    double chi2 = BlindedChi2(&hist_data,&hist_model);
    return chi2;

}

double EigenvectorFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    double func = 0.;
    func += par[0]/2.;
    func += par[1]*cos(1.*xx);
    func += par[2]*sin(1.*xx);
    func += par[3]*cos(2.*xx);
    func += par[4]*sin(2.*xx);
    func += par[5]*cos(3.*xx);
    func += par[6]*sin(3.*xx);
    func += par[7]*cos(4.*xx);
    func += par[8]*sin(4.*xx);
    func += par[9]*cos(5.*xx);
    func += par[10]*sin(5.*xx);
    func += par[11]*cos(6.*xx);
    func += par[12]*sin(6.*xx);
    return func;
}
void FitEigenvector(TH1* Hist_vtr, TH1* Hist_fit)
{
    TF1 *func = new TF1("func",EigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,13);
    func->SetParameter(0,0.);
    func->SetParameter(1,1.);
    func->SetParameter(2,1.);
    func->SetParameter(3,1.);
    func->SetParameter(4,1.);
    func->SetParameter(5,1.);
    func->SetParameter(6,1.);
    func->SetParameter(7,1.);
    func->SetParameter(8,1.);
    func->SetParameter(9,1.);
    func->SetParameter(10,1.);
    func->SetParameter(11,1.);
    func->SetParameter(12,1.);
    Hist_vtr->Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
    for (int bx=0;bx<Hist_fit->GetNbinsX();bx++)
    {
        Hist_fit->SetBinContent(bx+1,func->Eval(Hist_fit->GetBinCenter(bx+1)));
    }
}
void SetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_input);
    TF1 *func = new TF1("func",EigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,13);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    int n_fourier_modes = 6;
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;

    // 1st eigenvector
    
    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_input.cols()-1));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 0*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }

    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_input.cols()-1));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 1*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }

    // 2nd eigenvector
    
    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_input.cols()-2));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 2*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }

    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_input.cols()-2));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 3*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }

    // 3rd eigenvector
    
    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_input.cols()-3));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 4*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }

    fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_input.cols()-3));
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);

    first_index = 5*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.001);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.001);
        Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.001);
    }


    // eigenvalues
    first_index = 6*(1+2*n_fourier_modes);
    Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", eigensolver_data.eigenvalues()(N_bins_for_deconv-1).real(), 0.001);
    Chi2Minimizer->SetVariable(first_index+1, "par["+std::to_string(int(first_index+1))+"]", eigensolver_data.eigenvalues()(N_bins_for_deconv-2).real(), 0.001);
    Chi2Minimizer->SetVariable(first_index+2, "par["+std::to_string(int(first_index+2))+"]", eigensolver_data.eigenvalues()(N_bins_for_deconv-3).real(), 0.001);

}
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
MatrixXcd NoiseRemove(MatrixXcd mtx_input)
{
    //return mtx_input;
    MatrixXcd matrix(mtx_input.cols(),mtx_input.rows());
    for (int i=0;i<mtx_input.cols();i++)
    {
        for (int j=0;j<mtx_input.rows();j++)
        {
            //if (mtx_input(i,j).real()*mtx_input(i,j).real()+mtx_input(i,j).imag()*mtx_input(i,j).imag()>1) 
            //if (i>=mtx_input.cols()-4) 
            if (i>=mtx_input.cols()-3) 
            //if (i>=mtx_input.cols()-2) 
            {
                matrix(i,j) = mtx_input(i,j);
            }
            else 
            {
                matrix(i,j) = 0;
                //matrix(i,j) = mtx_input(i,j);
            }
        }
    }
    return matrix;
}
MatrixXcd MatrixBlind(MatrixXcd mtx_input, int binx_lower, int biny_lower, int binx_blind, int biny_blind)
{
    MatrixXcd matrix(mtx_input.cols(),mtx_input.rows());
    for (int i=0;i<mtx_input.cols();i++)
    {
        for (int j=0;j<mtx_input.rows();j++)
        {
            matrix(i,j) = mtx_input(i,j);
            if (i>=binx_lower && j>=biny_lower && i<=binx_blind && j<=biny_blind)
            {
                matrix(i,j) = 0;
            }
        }
    }
    return matrix;
}
VectorXcd fillVector(TH1D* hist)
{
    VectorXcd vtr(hist->GetNbinsX());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        vtr(binx) = hist->GetBinContent(binx+1);
    }
    return vtr;
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
MatrixXcd MakeSmoothSplineFunction(MatrixXcd mtx_origin)
{
    TH2D hist_2d = TH2D("hist_2d","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    TH2D hist_2d_smooth_col = TH2D("hist_2d_smooth_col","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    TH2D hist_2d_smooth_row = TH2D("hist_2d_smooth_row","",mtx_origin.cols(),0,mtx_origin.cols(),mtx_origin.rows(),0,mtx_origin.rows());
    fill2DHistogram(&hist_2d,mtx_origin);
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
void NetflixMethodPrediction(string target_data, double tel_elev_lower_input, double tel_elev_upper_input, double theta2_cut_lower_input, double theta2_cut_upper_input)
{

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    Theta2_cut_lower = theta2_cut_lower_input;
    Theta2_cut_upper = theta2_cut_upper_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;

    vector<TH2D> Hist_Bkgd_MSCLW;
    vector<TH1D> Hist_Data_Eigenvector_0;
    vector<TH1D> Hist_Fit_Eigenvector_0;
    vector<TH1D> Hist_Dark_Eigenvector_0;
    vector<TH1D> Hist_Data_Eigenvector_1;
    vector<TH1D> Hist_Fit_Eigenvector_1;
    vector<TH1D> Hist_Dark_Eigenvector_1;
    vector<TH1D> Hist_Data_Eigenvector_2;
    vector<TH1D> Hist_Fit_Eigenvector_2;
    vector<TH1D> Hist_Dark_Eigenvector_2;
    vector<TH1D> Hist_Data_InvEigenvector_0;
    vector<TH1D> Hist_Fit_InvEigenvector_0;
    vector<TH1D> Hist_Dark_InvEigenvector_0;
    vector<TH1D> Hist_Data_InvEigenvector_1;
    vector<TH1D> Hist_Fit_InvEigenvector_1;
    vector<TH1D> Hist_Dark_InvEigenvector_1;
    vector<TH1D> Hist_Data_InvEigenvector_2;
    vector<TH1D> Hist_Fit_InvEigenvector_2;
    vector<TH1D> Hist_Dark_InvEigenvector_2;
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
        Hist_Data_Eigenvector_0.push_back(TH1D("Hist_Data_Eigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_Eigenvector_0.push_back(TH1D("Hist_Fit_Eigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_Eigenvector_0.push_back(TH1D("Hist_Dark_Eigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_Eigenvector_1.push_back(TH1D("Hist_Data_Eigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_Eigenvector_1.push_back(TH1D("Hist_Fit_Eigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_Eigenvector_1.push_back(TH1D("Hist_Dark_Eigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_Eigenvector_2.push_back(TH1D("Hist_Data_Eigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_Eigenvector_2.push_back(TH1D("Hist_Fit_Eigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_Eigenvector_2.push_back(TH1D("Hist_Dark_Eigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvector_0.push_back(TH1D("Hist_Data_InvEigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvector_0.push_back(TH1D("Hist_Fit_InvEigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvector_0.push_back(TH1D("Hist_Dark_InvEigenvector_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvector_1.push_back(TH1D("Hist_Data_InvEigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvector_1.push_back(TH1D("Hist_Fit_InvEigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvector_1.push_back(TH1D("Hist_Dark_InvEigenvector_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvector_2.push_back(TH1D("Hist_Data_InvEigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvector_2.push_back(TH1D("Hist_Fit_InvEigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvector_2.push_back(TH1D("Hist_Dark_InvEigenvector_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));

        if (energy_bins[e]<237) continue;
        if (energy_bins[e]>=282) continue;

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

        ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_data);
        fill1DHistogram(&Hist_Data_Eigenvector_0.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-1));
        FitEigenvector(&Hist_Data_Eigenvector_0.at(e),&Hist_Fit_Eigenvector_0.at(e));
        fill1DHistogram(&Hist_Data_Eigenvector_1.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-2));
        FitEigenvector(&Hist_Data_Eigenvector_1.at(e),&Hist_Fit_Eigenvector_1.at(e));
        fill1DHistogram(&Hist_Data_Eigenvector_2.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-3));
        FitEigenvector(&Hist_Data_Eigenvector_2.at(e),&Hist_Fit_Eigenvector_2.at(e));
        fill1DHistogram(&Hist_Data_InvEigenvector_0.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1));
        FitEigenvector(&Hist_Data_InvEigenvector_0.at(e),&Hist_Fit_InvEigenvector_0.at(e));
        fill1DHistogram(&Hist_Data_InvEigenvector_1.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2));
        FitEigenvector(&Hist_Data_InvEigenvector_1.at(e),&Hist_Fit_InvEigenvector_1.at(e));
        fill1DHistogram(&Hist_Data_InvEigenvector_2.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3));
        FitEigenvector(&Hist_Data_InvEigenvector_2.at(e),&Hist_Fit_InvEigenvector_2.at(e));
        ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);
        fill1DHistogram(&Hist_Dark_Eigenvector_0.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Dark_Eigenvector_1.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Dark_Eigenvector_2.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Dark_InvEigenvector_0.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Dark_InvEigenvector_1.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Dark_InvEigenvector_2.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-3));

        //mtx_dark = MakeSmoothSplineFunction(mtx_dark);

        //initial condition
        //mtx_data_bkgd = mtx_data_blind;
        mtx_data_bkgd = mtx_dark;

        //std::cout << "Hist_Dark->Integral() = " << Hist_Dark->Integral() << std::endl;
        //for (int iter=0;iter<5;iter++)
        //{
        //    mtx_data_bkgd = EigenvaluePrediction(mtx_data,mtx_data_bkgd,mtx_dark,binx_lower,biny_lower,binx_blind,biny_blind);
        //    fill2DHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        //    std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;
        //}
        //mtx_data_bkgd = EigenvaluePrediction(mtx_data,mtx_data_bkgd,mtx_dark,binx_lower,biny_lower,binx_blind,biny_blind);
        
        
        // Choose method upon creation between:
        // kConjugateFR, kConjugatePR, kVectorBFGS,
        // kVectorBFGS2, kSteepestDescent
        ROOT::Math::GSLMinimizer Chi2Minimizer( ROOT::Math::kVectorBFGS );
        Chi2Minimizer.SetMaxFunctionCalls(1000000);
        Chi2Minimizer.SetMaxIterations(100000);
        Chi2Minimizer.SetTolerance(0.001);
        ROOT::Math::Functor Chi2Func(&Chi2Function,42); 
        Chi2Minimizer.SetFunction(Chi2Func);
        // Set the free variables to be minimized!
        SetInitialVariables(&Chi2Minimizer, mtx_data);
        Chi2Minimizer.Minimize();

        const double *par = Chi2Minimizer.X();
        // build eigenvector matrix
        MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
        MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
        MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenvector(row,col) = 0.;
                mtx_eigenvalue(row,col) = 0.;
            }
        }

        int n_fourier_modes = 6;
        int col_fix = 0;
        int row_fix = 0;
        int first_index = 0;
        // 1st eigenvector
        col_fix = N_bins_for_deconv-1;
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 0*(1+2*n_fourier_modes);
            mtx_eigenvector(row,col_fix) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        row_fix = N_bins_for_deconv-1;
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 1*(1+2*n_fourier_modes);
            mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        // 2nd eigenvector
        col_fix = N_bins_for_deconv-2;
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 2*(1+2*n_fourier_modes);
            mtx_eigenvector(row,col_fix) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        row_fix = N_bins_for_deconv-2;
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 3*(1+2*n_fourier_modes);
            mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        // 3rd eigenvector
        col_fix = N_bins_for_deconv-3;
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 4*(1+2*n_fourier_modes);
            mtx_eigenvector(row,col_fix) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        row_fix = N_bins_for_deconv-3;
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            first_index = 5*(1+2*n_fourier_modes);
            mtx_eigenvector_inv(row_fix,col) = par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x);
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x);
            }
        }
        // build eigenvalue matrix
        first_index = 6*(1+2*n_fourier_modes);
        mtx_eigenvalue(N_bins_for_deconv-1,N_bins_for_deconv-1) = par[first_index];
        mtx_eigenvalue(N_bins_for_deconv-2,N_bins_for_deconv-2) = par[first_index+1];
        mtx_eigenvalue(N_bins_for_deconv-3,N_bins_for_deconv-3) = par[first_index+2];
        std::cout << "Minimum: (" << par[first_index] << "," << par[first_index+1] << "," << par[first_index+2] << ") " << std::endl;


        std::cout << "eigensolver_data.eigenvalues()(N_bins_for_deconv-1) = " << eigensolver_data.eigenvalues()(N_bins_for_deconv-1) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(N_bins_for_deconv-2) = " << eigensolver_data.eigenvalues()(N_bins_for_deconv-2) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(N_bins_for_deconv-3) = " << eigensolver_data.eigenvalues()(N_bins_for_deconv-3) << std::endl;
        std::cout << "mtx_eigenvalue = " << std::endl;
        std::cout << mtx_eigenvalue.block(N_bins_for_deconv-3,N_bins_for_deconv-3,3,3) << std::endl;
        //std::cout << "eigensolver_data.eigenvectors() = " << std::endl;
        //std::cout << eigensolver_data.eigenvectors().block(N_bins_for_deconv-3,N_bins_for_deconv-3,3,3) << std::endl;
        //std::cout << "mtx_eigenvector = " << std::endl;
        //std::cout << mtx_eigenvector.block(N_bins_for_deconv-3,N_bins_for_deconv-3,3,3) << std::endl;
        std::cout << "eigensolver_data.eigenvectors() = " << std::endl;
        std::cout << eigensolver_data.eigenvectors().block(1,N_bins_for_deconv-3,N_bins_for_deconv,3) << std::endl;
        std::cout << "mtx_eigenvector = " << std::endl;
        std::cout << mtx_eigenvector.block(1,N_bins_for_deconv-3,N_bins_for_deconv,3) << std::endl;
        //std::cout << "eigensolver_data.eigenvectors().inverse() = " << std::endl;
        //std::cout << eigensolver_data.eigenvectors().inverse().block(N_bins_for_deconv-3,N_bins_for_deconv-3,3,3) << std::endl;
        //std::cout << "mtx_eigenvector_inv = " << std::endl;
        //std::cout << mtx_eigenvector_inv.block(N_bins_for_deconv-3,N_bins_for_deconv-3,3,3) << std::endl;
        std::cout << "eigensolver_data.eigenvectors().inverse() = " << std::endl;
        std::cout << eigensolver_data.eigenvectors().inverse().block(N_bins_for_deconv-3,1,3,N_bins_for_deconv).transpose() << std::endl;
        std::cout << "mtx_eigenvector_inv = " << std::endl;
        std::cout << mtx_eigenvector_inv.block(N_bins_for_deconv-3,1,3,N_bins_for_deconv).transpose() << std::endl;

        fill1DHistogram(&Hist_Fit_Eigenvector_0.at(e),mtx_eigenvector.col(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Fit_Eigenvector_1.at(e),mtx_eigenvector.col(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Fit_Eigenvector_2.at(e),mtx_eigenvector.col(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Fit_InvEigenvector_0.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Fit_InvEigenvector_1.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Fit_InvEigenvector_2.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-3));
        //mtx_eigenvector.col(mtx_data.cols()-1) = fillVector(&Hist_Data_Eigenvector_0.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-2) = fillVector(&Hist_Data_Eigenvector_1.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-3) = fillVector(&Hist_Data_Eigenvector_2.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-1) = eigensolver_data.eigenvectors().col(mtx_data.cols()-1);
        //mtx_eigenvector.col(mtx_data.cols()-2) = eigensolver_data.eigenvectors().col(mtx_data.cols()-2);
        //mtx_eigenvector.col(mtx_data.cols()-3) = eigensolver_data.eigenvectors().col(mtx_data.cols()-3);
        //for (int i=0;i<mtx_dark.cols();i++)
        //{
        //    for (int j=0;j<mtx_dark.rows();j++)
        //    {
        //        if (i==j) 
        //        {
        //            mtx_eigenvalue(i,j) = eigensolver_data.eigenvalues()(i);
        //        }
        //        else
        //        {
        //            mtx_eigenvalue(i,j) = 0;
        //        }
        //    }
        //}
        mtx_data_bkgd = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;
        fill2DHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;
        //mtx_data_bkgd = eigensolver_data.eigenvectors()*mtx_eigenvalue*eigensolver_data.eigenvectors().inverse();
        //fill2DHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        //std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;

    }
    InputDataFile.Close();


    TFile OutputFile("output_Jul16/Netflix_"+TString(target_data)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Bkgd_MSCLW.at(e).Write();
        Hist_Data_Eigenvector_0.at(e).Write();
        Hist_Fit_Eigenvector_0.at(e).Write();
        Hist_Dark_Eigenvector_0.at(e).Write();
        Hist_Data_Eigenvector_1.at(e).Write();
        Hist_Fit_Eigenvector_1.at(e).Write();
        Hist_Dark_Eigenvector_1.at(e).Write();
        Hist_Data_Eigenvector_2.at(e).Write();
        Hist_Fit_Eigenvector_2.at(e).Write();
        Hist_Dark_Eigenvector_2.at(e).Write();
        Hist_Data_InvEigenvector_0.at(e).Write();
        Hist_Fit_InvEigenvector_0.at(e).Write();
        Hist_Dark_InvEigenvector_0.at(e).Write();
        Hist_Data_InvEigenvector_1.at(e).Write();
        Hist_Fit_InvEigenvector_1.at(e).Write();
        Hist_Dark_InvEigenvector_1.at(e).Write();
        Hist_Data_InvEigenvector_2.at(e).Write();
        Hist_Fit_InvEigenvector_2.at(e).Write();
        Hist_Dark_InvEigenvector_2.at(e).Write();
    }
    OutputFile.Close();

}
