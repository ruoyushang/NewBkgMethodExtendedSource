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
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

ClassImp(TSplinePoly);
ClassImp(TSplinePoly3);
ClassImp(TSplinePoly5);
ClassImp(TSpline3);
ClassImp(TSpline5);
ClassImp(TSpline);

#include <math.h>
#include <complex>
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
#include <complex>
#include <cmath>
using namespace Eigen;
//using Eigen::MatrixXcd;

#include "NetflixMethodGetShowerImage.C"

int NumberOfEigenvectors = 0;
int n_fourier_modes = 6;
int n_taylor_modes = 6;
MatrixXcd mtx_data(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);

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
        hist->SetBinError(binx+1,0.01*abs(vtr(binx).real()));
    }
}
double BlindedChi2(TH2D* hist_data, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind*2);
    int biny_upper = hist_data->GetYaxis()->FindBin(MSCW_cut_blind*2);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            double weight = 1.;
            double dx = hist_data->GetXaxis()->GetBinCenter(bx);
            double dy = hist_data->GetYaxis()->GetBinCenter(by);
            double radius = pow(dx*dx+dy*dy,0.5);
            weight = exp(-radius*radius);
            //if (bx<=binx_upper && by<=biny_upper) weight = 100.;
            if (bx>=binx_blind || by>=biny_blind)
            {
                double data = hist_data->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                chi2 += weight*pow(data-model,2);
                //if (data*data>0)
                //{
                //    chi2 += weight*pow(data-model,2)/(data*data);
                //}
                //else
                //{
                //    chi2 += weight*pow(data-model,2);
                //}
            }
        }
    }
    return chi2;
}

void TaylorParametrizeEigenvectors(const double *par)
{
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }
    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-4)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
            double x_norm = (x-MSCL_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
            mtx_eigenvector(row,col_fix) += par[first_index];
            for (int mode=1;mode<=n_taylor_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*pow(x_norm-par[first_index+2*mode],mode);
            }
        }
        first_index = (4*NthEigenvector-3)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
            double x_norm = (x-MSCL_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
            mtx_eigenvector(row,col_fix) += If*par[first_index];
            for (int mode=1;mode<=n_taylor_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += If*par[first_index+2*mode-1]*pow(x_norm-par[first_index+2*mode],mode);
            }
        }
        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-2)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            double x_norm = (x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
            mtx_eigenvector_inv(row_fix,col) += par[first_index];
            for (int mode=1;mode<=n_taylor_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*pow(x_norm-par[first_index+2*mode],mode);
            }
        }
        first_index = (4*NthEigenvector-1)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            double x_norm = (x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
            mtx_eigenvector_inv(row_fix,col) += If*par[first_index];
            for (int mode=1;mode<=n_taylor_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += If*par[first_index+2*mode-1]*pow(x_norm-par[first_index+2*mode],mode);
            }
        }
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        mtx_eigenvalue(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvector) = par[first_index];
    }
}

void ParametrizeEigenvectors(const double *par)
{
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }
    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-4)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col_fix) += par[first_index+row];
        }
        first_index = (4*NthEigenvector-3)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col_fix) += If*par[first_index+row];
        }
        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-2)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv(row_fix,col) += par[first_index+col];
        }
        first_index = (4*NthEigenvector-1)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv(row_fix,col) += If*par[first_index+col];
        }
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*(N_bins_for_deconv)+(NthEigenvector-1);
        mtx_eigenvalue(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvector) = par[first_index];
    }
}
void FourierParametrizeEigenvectors(const double *par)
{
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }
    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-4)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
            double x_norm = 1*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
            mtx_eigenvector(row,col_fix) += par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode-1]*cos(double(mode)*x_norm);
                mtx_eigenvector(row,col_fix) += par[first_index+2*mode]*sin(double(mode)*x_norm);
            }
        }
        first_index = (4*NthEigenvector-3)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
            double x_norm = 1*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
            mtx_eigenvector(row,col_fix) += If*par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector(row,col_fix) += If*par[first_index+2*mode-1]*cos(double(mode)*x_norm);
                mtx_eigenvector(row,col_fix) += If*par[first_index+2*mode]*sin(double(mode)*x_norm);
            }
        }
        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-2)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            double x_norm = 1*M_PI*(x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
            mtx_eigenvector_inv(row_fix,col) += par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode-1]*cos(double(mode)*x_norm);
                mtx_eigenvector_inv(row_fix,col) += par[first_index+2*mode]*sin(double(mode)*x_norm);
            }
        }
        first_index = (4*NthEigenvector-1)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
            double x_norm = 1*M_PI*(x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
            mtx_eigenvector_inv(row_fix,col) += If*par[first_index]/2.;
            for (int mode=1;mode<=n_fourier_modes;mode++)
            {
                mtx_eigenvector_inv(row_fix,col) += If*par[first_index+2*mode-1]*cos(double(mode)*x_norm);
                mtx_eigenvector_inv(row_fix,col) += If*par[first_index+2*mode]*sin(double(mode)*x_norm);
            }
        }
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        mtx_eigenvalue(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvector) = par[first_index];
    }
}

double TaylorChi2Function(const double *par)
{

    TaylorParametrizeEigenvectors(par);
    
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

double Chi2Function(const double *par)
{

    ParametrizeEigenvectors(par);
    
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

double FourierChi2Function(const double *par)
{

    FourierParametrizeEigenvectors(par);
    
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

double TaylorEigenvectorFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    double func = 0.;
    func += par[0];
    double x_norm = (xx-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
    for (int mode=1;mode<=n_taylor_modes;mode++)
    {
        func += par[2*mode-1]*pow(x_norm-par[2*mode],mode);
    }
    return func;
}
double FourierEigenvectorFunction(Double_t *x, Double_t *par) {
    double xx =x[0];
    double func = 0.;
    func += par[0]/2.;
    double x_norm = 1*M_PI*(xx-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func += par[2*mode-1]*cos(double(mode)*x_norm);
        func += par[2*mode]*sin(double(mode)*x_norm);
    }
    return func;
}
void TaylorFitEigenvector(TH1* Hist_vtr, TH1* Hist_fit)
{
    TF1 *func = new TF1("func",TaylorEigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,1+2*n_taylor_modes);
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_taylor_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_vtr->Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
    for (int bx=0;bx<Hist_fit->GetNbinsX();bx++)
    {
        Hist_fit->SetBinContent(bx+1,func->Eval(Hist_fit->GetBinCenter(bx+1)));
    }
}
void FourierFitEigenvector(TH1* Hist_vtr, TH1* Hist_fit)
{
    TF1 *func = new TF1("func",FourierEigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,1+2*n_fourier_modes);
    func->SetParameter(0,0.);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        func->SetParameter(2*mode-1,1.);
        func->SetParameter(2*mode,1.);
    }
    Hist_vtr->Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
    for (int bx=0;bx<Hist_fit->GetNbinsX();bx++)
    {
        Hist_fit->SetBinContent(bx+1,func->Eval(Hist_fit->GetBinCenter(bx+1)));
    }
}
void TaylorSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_input);
    TF1 *func = new TF1("func",TaylorEigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,1+2*n_taylor_modes);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;

    // 1st eigenvector
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_input.cols()-NthEigenvector).real());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            func->SetParameter(mode,1.);
        }
        Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-4)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.01);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.01);
        }

        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_input.cols()-NthEigenvector).imag());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            func->SetParameter(mode,1.);
        }
        Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-3)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.01);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.01);
        }

        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_input.cols()-NthEigenvector).real());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            func->SetParameter(mode,1.);
        }
        Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-2)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.01);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.01);
        }

        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_input.cols()-NthEigenvector).imag());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            func->SetParameter(mode,1.);
        }
        Hist_Data_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-1)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",func->GetParameter(0),0.01);
        for (int mode=1;mode<=n_taylor_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", func->GetParameter(2*mode), 0.01);
        }

        // eigenvalues
        first_index = (4*NthEigenvector-0)*(1+2*n_taylor_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", eigensolver_data.eigenvalues()(N_bins_for_deconv-NthEigenvector).real(), 0.01);
        Chi2Minimizer->FixVariable(first_index);
    }

}
void SetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_dark, MatrixXcd mtx_data)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);
    ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_data);
    TH1D Hist_Dark_Eigenvector = TH1D("Hist_Dark_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    const std::complex<double> If(0.0, 1.0);
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale_real = 0.;
    double scale_imag = 0.;

    // 1st eigenvector
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        scale_real = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        scale_imag = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());

        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-4)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double initial_par = scale_real*eigensolver_dark.eigenvectors().col(col_fix)(row).real();
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",initial_par,0.01);
            Chi2Minimizer->SetVariableLimits(first_index+row,initial_par*0.2,initial_par*5.0);
        }
        first_index = (4*NthEigenvector-3)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double initial_par = scale_imag*eigensolver_dark.eigenvectors().col(col_fix)(row).imag();
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",initial_par,0.01);
            Chi2Minimizer->SetVariableLimits(first_index+row,initial_par*0.2,initial_par*5.0);
        }
        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-2)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double initial_par = scale_real*eigensolver_dark.eigenvectors().inverse().row(row_fix)(col).real();
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",initial_par,0.01);
            Chi2Minimizer->SetVariableLimits(first_index+col,initial_par*0.2,initial_par*5.0);
        }
        first_index = (4*NthEigenvector-1)*(N_bins_for_deconv)+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double initial_par = scale_imag*eigensolver_dark.eigenvectors().inverse().row(row_fix)(col).imag();
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",initial_par,0.01);
            Chi2Minimizer->SetVariableLimits(first_index+col,initial_par*0.2,initial_par*5.0);
        }
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*(N_bins_for_deconv)+(NthEigenvector-1);
        double initial_eigenvalue = eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real();
        Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", initial_eigenvalue, 0.01);
        Chi2Minimizer->SetVariableLimits(first_index,initial_eigenvalue*0.8,initial_eigenvalue*1.2);
        //Chi2Minimizer->FixVariable(first_index);
    }

}
void FourierSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_dark, MatrixXcd mtx_data)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);
    ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_data);
    TF1 *func = new TF1("func",FourierEigenvectorFunction,MSCL_plot_lower,MSCL_plot_upper,1+2*n_fourier_modes);
    TH1D Hist_Dark_Eigenvector = TH1D("Hist_Dark_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale = 0.;

    // 1st eigenvector
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        scale = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            func->SetParameter(2*mode-1,1.);
            func->SetParameter(2*mode,1.);
        }
        Hist_Dark_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-4)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",scale*func->GetParameter(0),0.01);
        Chi2Minimizer->SetVariableLimits(first_index,scale*func->GetParameter(0)*0.8,scale*func->GetParameter(0)*1.2);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", scale*func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,scale*func->GetParameter(2*mode-1)*0.8,scale*func->GetParameter(2*mode-1)*1.2);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", scale*func->GetParameter(2*mode), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,scale*func->GetParameter(2*mode)*0.8,scale*func->GetParameter(2*mode)*1.2);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        scale = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            func->SetParameter(2*mode-1,1.);
            func->SetParameter(2*mode,1.);
        }
        Hist_Dark_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-3)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",scale*func->GetParameter(0),0.01);
        Chi2Minimizer->SetVariableLimits(first_index,scale*func->GetParameter(0)*0.8,scale*func->GetParameter(0)*1.2);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", scale*func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,scale*func->GetParameter(2*mode-1)*0.8,scale*func->GetParameter(2*mode-1)*1.2);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", scale*func->GetParameter(2*mode), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,scale*func->GetParameter(2*mode)*0.8,scale*func->GetParameter(2*mode)*1.2);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).real());
        scale = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            func->SetParameter(2*mode-1,1.);
            func->SetParameter(2*mode,1.);
        }
        Hist_Dark_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-2)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",scale*func->GetParameter(0),0.01);
        Chi2Minimizer->SetVariableLimits(first_index,scale*func->GetParameter(0)*0.8,scale*func->GetParameter(0)*1.2);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", scale*func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,scale*func->GetParameter(2*mode-1)*0.8,scale*func->GetParameter(2*mode-1)*1.2);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", scale*func->GetParameter(2*mode), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,scale*func->GetParameter(2*mode)*0.8,scale*func->GetParameter(2*mode)*1.2);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).imag());
        scale = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        func->SetParameter(0,0.);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            func->SetParameter(2*mode-1,1.);
            func->SetParameter(2*mode,1.);
        }
        Hist_Dark_Eigenvector.Fit("func","","",MSCL_plot_lower,MSCL_plot_upper);
        first_index = (4*NthEigenvector-1)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",scale*func->GetParameter(0),0.01);
        Chi2Minimizer->SetVariableLimits(first_index,scale*func->GetParameter(0)*0.8,scale*func->GetParameter(0)*1.2);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]", scale*func->GetParameter(2*mode-1), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,scale*func->GetParameter(2*mode-1)*0.8,scale*func->GetParameter(2*mode-1)*1.2);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]", scale*func->GetParameter(2*mode), 0.01);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,scale*func->GetParameter(2*mode)*0.8,scale*func->GetParameter(2*mode)*1.2);
        }

        // eigenvalues
        first_index = (4*NthEigenvector-0)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real(), 0.01);
        Chi2Minimizer->SetVariableLimits(first_index,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real()*0.8,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real()*1.2);
        //Chi2Minimizer->FixVariable(first_index);
    }

}
void NormalizeHist2D(TH2D* hist_ref, TH2D* hist_sub)
{
    int binx_lower = hist_ref->GetXaxis()->FindBin(MSCW_cut_lower);
    int binx_upper = hist_ref->GetXaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_lower = hist_ref->GetYaxis()->FindBin(MSCW_cut_blind);
    int biny_upper = hist_ref->GetYaxis()->FindBin(3.*MSCW_cut_blind)-1;
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
VectorXd fillVector(TH1D* hist)
{
    VectorXd vtr(hist->GetNbinsX());
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
    vector<TH1D> Hist_Data_EigenvectorReal_0;
    vector<TH1D> Hist_Fit_EigenvectorReal_0;
    vector<TH1D> Hist_Dark_EigenvectorReal_0;
    vector<TH1D> Hist_Data_EigenvectorReal_1;
    vector<TH1D> Hist_Fit_EigenvectorReal_1;
    vector<TH1D> Hist_Dark_EigenvectorReal_1;
    vector<TH1D> Hist_Data_EigenvectorReal_2;
    vector<TH1D> Hist_Fit_EigenvectorReal_2;
    vector<TH1D> Hist_Dark_EigenvectorReal_2;
    vector<TH1D> Hist_Data_InvEigenvectorReal_0;
    vector<TH1D> Hist_Fit_InvEigenvectorReal_0;
    vector<TH1D> Hist_Dark_InvEigenvectorReal_0;
    vector<TH1D> Hist_Data_InvEigenvectorReal_1;
    vector<TH1D> Hist_Fit_InvEigenvectorReal_1;
    vector<TH1D> Hist_Dark_InvEigenvectorReal_1;
    vector<TH1D> Hist_Data_InvEigenvectorReal_2;
    vector<TH1D> Hist_Fit_InvEigenvectorReal_2;
    vector<TH1D> Hist_Dark_InvEigenvectorReal_2;
    vector<TH1D> Hist_Data_EigenvectorImag_0;
    vector<TH1D> Hist_Fit_EigenvectorImag_0;
    vector<TH1D> Hist_Dark_EigenvectorImag_0;
    vector<TH1D> Hist_Data_EigenvectorImag_1;
    vector<TH1D> Hist_Fit_EigenvectorImag_1;
    vector<TH1D> Hist_Dark_EigenvectorImag_1;
    vector<TH1D> Hist_Data_EigenvectorImag_2;
    vector<TH1D> Hist_Fit_EigenvectorImag_2;
    vector<TH1D> Hist_Dark_EigenvectorImag_2;
    vector<TH1D> Hist_Data_InvEigenvectorImag_0;
    vector<TH1D> Hist_Fit_InvEigenvectorImag_0;
    vector<TH1D> Hist_Dark_InvEigenvectorImag_0;
    vector<TH1D> Hist_Data_InvEigenvectorImag_1;
    vector<TH1D> Hist_Fit_InvEigenvectorImag_1;
    vector<TH1D> Hist_Dark_InvEigenvectorImag_1;
    vector<TH1D> Hist_Data_InvEigenvectorImag_2;
    vector<TH1D> Hist_Fit_InvEigenvectorImag_2;
    vector<TH1D> Hist_Dark_InvEigenvectorImag_2;
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
        Hist_Data_EigenvectorReal_0.push_back(TH1D("Hist_Data_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorReal_0.push_back(TH1D("Hist_Fit_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorReal_0.push_back(TH1D("Hist_Dark_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_EigenvectorReal_1.push_back(TH1D("Hist_Data_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorReal_1.push_back(TH1D("Hist_Fit_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorReal_1.push_back(TH1D("Hist_Dark_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_EigenvectorReal_2.push_back(TH1D("Hist_Data_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorReal_2.push_back(TH1D("Hist_Fit_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorReal_2.push_back(TH1D("Hist_Dark_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorReal_0.push_back(TH1D("Hist_Data_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorReal_0.push_back(TH1D("Hist_Fit_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorReal_0.push_back(TH1D("Hist_Dark_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorReal_1.push_back(TH1D("Hist_Data_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorReal_1.push_back(TH1D("Hist_Fit_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorReal_1.push_back(TH1D("Hist_Dark_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorReal_2.push_back(TH1D("Hist_Data_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorReal_2.push_back(TH1D("Hist_Fit_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorReal_2.push_back(TH1D("Hist_Dark_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_EigenvectorImag_0.push_back(TH1D("Hist_Data_EigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorImag_0.push_back(TH1D("Hist_Fit_EigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorImag_0.push_back(TH1D("Hist_Dark_EigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_EigenvectorImag_1.push_back(TH1D("Hist_Data_EigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorImag_1.push_back(TH1D("Hist_Fit_EigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorImag_1.push_back(TH1D("Hist_Dark_EigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_EigenvectorImag_2.push_back(TH1D("Hist_Data_EigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_EigenvectorImag_2.push_back(TH1D("Hist_Fit_EigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_EigenvectorImag_2.push_back(TH1D("Hist_Dark_EigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorImag_0.push_back(TH1D("Hist_Data_InvEigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorImag_0.push_back(TH1D("Hist_Fit_InvEigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorImag_0.push_back(TH1D("Hist_Dark_InvEigenvectorImag_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorImag_1.push_back(TH1D("Hist_Data_InvEigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorImag_1.push_back(TH1D("Hist_Fit_InvEigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorImag_1.push_back(TH1D("Hist_Dark_InvEigenvectorImag_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Data_InvEigenvectorImag_2.push_back(TH1D("Hist_Data_InvEigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Fit_InvEigenvectorImag_2.push_back(TH1D("Hist_Fit_InvEigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));
        Hist_Dark_InvEigenvectorImag_2.push_back(TH1D("Hist_Dark_InvEigenvectorImag_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper));

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

        n_fourier_modes = 6;
        n_taylor_modes = 6;
        ComplexEigenSolver<MatrixXcd> eigensolver_data(mtx_data);
        ComplexEigenSolver<MatrixXcd> eigensolver_dark(mtx_dark);
        fill1DHistogram(&Hist_Data_EigenvectorReal_0.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Data_EigenvectorReal_1.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Data_EigenvectorReal_2.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_0.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_1.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_2.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3).real());
        FourierFitEigenvector(&Hist_Data_EigenvectorReal_0.at(e),&Hist_Fit_EigenvectorReal_0.at(e));
        FourierFitEigenvector(&Hist_Data_EigenvectorReal_1.at(e),&Hist_Fit_EigenvectorReal_1.at(e));
        FourierFitEigenvector(&Hist_Data_EigenvectorReal_2.at(e),&Hist_Fit_EigenvectorReal_2.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_0.at(e),&Hist_Fit_InvEigenvectorReal_0.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_1.at(e),&Hist_Fit_InvEigenvectorReal_1.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_2.at(e),&Hist_Fit_InvEigenvectorReal_2.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorReal_0.at(e),&Hist_Fit_EigenvectorReal_0.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorReal_1.at(e),&Hist_Fit_EigenvectorReal_1.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorReal_2.at(e),&Hist_Fit_EigenvectorReal_2.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorReal_0.at(e),&Hist_Fit_InvEigenvectorReal_0.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorReal_1.at(e),&Hist_Fit_InvEigenvectorReal_1.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorReal_2.at(e),&Hist_Fit_InvEigenvectorReal_2.at(e));
        fill1DHistogram(&Hist_Dark_EigenvectorReal_0.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Dark_EigenvectorReal_1.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Dark_EigenvectorReal_2.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_0.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_1.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_2.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Data_EigenvectorImag_0.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Data_EigenvectorImag_1.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Data_EigenvectorImag_2.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-3).imag());
        fill1DHistogram(&Hist_Data_InvEigenvectorImag_0.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Data_InvEigenvectorImag_1.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Data_InvEigenvectorImag_2.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3).imag());
        FourierFitEigenvector(&Hist_Data_EigenvectorImag_0.at(e),&Hist_Fit_EigenvectorImag_0.at(e));
        FourierFitEigenvector(&Hist_Data_EigenvectorImag_1.at(e),&Hist_Fit_EigenvectorImag_1.at(e));
        FourierFitEigenvector(&Hist_Data_EigenvectorImag_2.at(e),&Hist_Fit_EigenvectorImag_2.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_0.at(e),&Hist_Fit_InvEigenvectorImag_0.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_1.at(e),&Hist_Fit_InvEigenvectorImag_1.at(e));
        FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_2.at(e),&Hist_Fit_InvEigenvectorImag_2.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorImag_0.at(e),&Hist_Fit_EigenvectorImag_0.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorImag_1.at(e),&Hist_Fit_EigenvectorImag_1.at(e));
        //TaylorFitEigenvector(&Hist_Data_EigenvectorImag_2.at(e),&Hist_Fit_EigenvectorImag_2.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorImag_0.at(e),&Hist_Fit_InvEigenvectorImag_0.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorImag_1.at(e),&Hist_Fit_InvEigenvectorImag_1.at(e));
        //TaylorFitEigenvector(&Hist_Data_InvEigenvectorImag_2.at(e),&Hist_Fit_InvEigenvectorImag_2.at(e));
        fill1DHistogram(&Hist_Dark_EigenvectorImag_0.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Dark_EigenvectorImag_1.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Dark_EigenvectorImag_2.at(e),eigensolver_dark.eigenvectors().col(mtx_data.cols()-3).imag());
        fill1DHistogram(&Hist_Dark_InvEigenvectorImag_0.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Dark_InvEigenvectorImag_1.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Dark_InvEigenvectorImag_2.at(e),eigensolver_dark.eigenvectors().inverse().row(mtx_data.cols()-3).imag());

        mtx_data_bkgd = mtx_dark;

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenvector(row,col) = 0.;
                mtx_eigenvector_inv(row,col) = 0.;
                mtx_eigenvalue(row,col) = 0.;
            }
        }
        
        int col_fix = 0;
        int row_fix = 0;
        int first_index = 0;

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenvector(row,col) = 0.;
                mtx_eigenvector_inv(row,col) = 0.;
                mtx_eigenvalue(row,col) = 0.;
            }
        }

        // Choose method upon creation between:
        // kConjugateFR, kConjugatePR, kVectorBFGS,
        // kVectorBFGS2, kSteepestDescent

        std::cout << "++++++++++++++++++++ Searching for 1st eigenvector... ++++++++++++++++++++" << std::endl;
        NumberOfEigenvectors = 3;
        ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kSteepestDescent );
        //ROOT::Minuit2::Minuit2Minimizer Chi2Minimizer_1st( ROOT::Minuit2::kMigrad );
        //Chi2Minimizer_1st.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
        //Chi2Minimizer_1st.SetMaxIterations(100000); // for GSL
        //Chi2Minimizer_1st.SetTolerance(0.001);
        Chi2Minimizer_1st.SetMaxFunctionCalls(200); // for Minuit/Minuit2
        Chi2Minimizer_1st.SetMaxIterations(200); // for GSL
        Chi2Minimizer_1st.SetTolerance(0.01);

        ROOT::Math::Functor Chi2Func_1st(&Chi2Function,4*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors);        
        Chi2Minimizer_1st.SetFunction(Chi2Func_1st);
        SetInitialVariables(&Chi2Minimizer_1st, mtx_dark, mtx_data);
        Chi2Minimizer_1st.Minimize();
        const double *par = Chi2Minimizer_1st.X();
        ParametrizeEigenvectors(par);

        //n_fourier_modes = 20;
        //ROOT::Math::Functor Chi2Func_1st(&FourierChi2Function,4*NumberOfEigenvectors*(1+2*n_fourier_modes)+NumberOfEigenvectors);        Chi2Minimizer_1st.SetFunction(Chi2Func_1st);
        //FourierSetInitialVariables(&Chi2Minimizer_1st, mtx_dark, mtx_data);
        ////Chi2Minimizer_1st.Minimize();
        //const double *par = Chi2Minimizer_1st.X();
        //FourierParametrizeEigenvectors(par);

        //n_taylor_modes = 8;
        //ROOT::Math::Functor Chi2Func_1st(&TaylorChi2Function,4*NumberOfEigenvectors*(1+2*n_taylor_modes)+NumberOfEigenvectors);        
        //Chi2Minimizer_1st.SetFunction(Chi2Func_1st);
        //TaylorSetInitialVariables(&Chi2Minimizer_1st, mtx_dark);
        ////Chi2Minimizer_1st.Minimize();
        //const double *par = Chi2Minimizer_1st.X();
        //TaylorParametrizeEigenvectors(par);

        //std::cout << "eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1).transpose()" << std::endl;
        //std::cout << eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1).transpose() << std::endl;
        //std::cout << "mtx_eigenvector_inv.row(mtx_data.cols()-1).transpose()" << std::endl;
        //std::cout << mtx_eigenvector_inv.row(mtx_data.cols()-1).transpose() << std::endl;
        //std::cout << "Hist_Data_InvEigenvectorReal_0.at(e).Print():" << std::endl;
        //Hist_Data_InvEigenvectorReal_0.at(e).Print("All");

        fill1DHistogram(&Hist_Fit_EigenvectorReal_0.at(e),mtx_eigenvector.col(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Fit_EigenvectorReal_1.at(e),mtx_eigenvector.col(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Fit_EigenvectorReal_2.at(e),mtx_eigenvector.col(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_0.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_1.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_2.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Fit_EigenvectorImag_0.at(e),mtx_eigenvector.col(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Fit_EigenvectorImag_1.at(e),mtx_eigenvector.col(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Fit_EigenvectorImag_2.at(e),mtx_eigenvector.col(mtx_data.cols()-3).imag());
        fill1DHistogram(&Hist_Fit_InvEigenvectorImag_0.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-1).imag());
        fill1DHistogram(&Hist_Fit_InvEigenvectorImag_1.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-2).imag());
        fill1DHistogram(&Hist_Fit_InvEigenvectorImag_2.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-3).imag());
        //
        //mtx_eigenvector.col(mtx_data.cols()-1).real() = fillVector(&Hist_Data_EigenvectorReal_0.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-2).real() = fillVector(&Hist_Data_EigenvectorReal_1.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-3).real() = fillVector(&Hist_Data_EigenvectorReal_2.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-1).real() = fillVector(&Hist_Data_InvEigenvectorReal_0.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-2).real() = fillVector(&Hist_Data_InvEigenvectorReal_1.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-3).real() = fillVector(&Hist_Data_InvEigenvectorReal_2.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-1).imag() = fillVector(&Hist_Data_EigenvectorImag_0.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-2).imag() = fillVector(&Hist_Data_EigenvectorImag_1.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-3).imag() = fillVector(&Hist_Data_EigenvectorImag_2.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-1).imag() = fillVector(&Hist_Data_InvEigenvectorImag_0.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-2).imag() = fillVector(&Hist_Data_InvEigenvectorImag_1.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-3).imag() = fillVector(&Hist_Data_InvEigenvectorImag_2.at(e));
        //std::cout << "mtx_eigenvector_inv.row(mtx_data.cols()-1)" << std::endl;
        //std::cout << mtx_eigenvector_inv.row(mtx_data.cols()-1) << std::endl;
        //
        //mtx_eigenvector.col(mtx_data.cols()-1).real() = fillVector(&Hist_Fit_EigenvectorReal_0.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-2).real() = fillVector(&Hist_Fit_EigenvectorReal_1.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-3).real() = fillVector(&Hist_Fit_EigenvectorReal_2.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-1).real() = fillVector(&Hist_Fit_InvEigenvectorReal_0.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-2).real() = fillVector(&Hist_Fit_InvEigenvectorReal_1.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-3).real() = fillVector(&Hist_Fit_InvEigenvectorReal_2.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-1).imag() = fillVector(&Hist_Fit_EigenvectorImag_0.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-2).imag() = fillVector(&Hist_Fit_EigenvectorImag_1.at(e));
        //mtx_eigenvector.col(mtx_data.cols()-3).imag() = fillVector(&Hist_Fit_EigenvectorImag_2.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-1).imag() = fillVector(&Hist_Fit_InvEigenvectorImag_0.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-2).imag() = fillVector(&Hist_Fit_InvEigenvectorImag_1.at(e));
        //mtx_eigenvector_inv.row(mtx_data.cols()-3).imag() = fillVector(&Hist_Fit_InvEigenvectorImag_2.at(e));
        //
        //for (int col=0;col<N_bins_for_deconv;col++)
        //{
        //    for (int row=0;row<N_bins_for_deconv;row++)
        //    {
        //        mtx_eigenvector(row,col) = 0.;
        //        mtx_eigenvector_inv(row,col) = 0.;
        //    }
        //}
        //mtx_eigenvector.col(mtx_data.cols()-1) = eigensolver_data.eigenvectors().col(mtx_data.cols()-1);
        //mtx_eigenvector.col(mtx_data.cols()-2) = eigensolver_data.eigenvectors().col(mtx_data.cols()-2);
        //mtx_eigenvector.col(mtx_data.cols()-3) = eigensolver_data.eigenvectors().col(mtx_data.cols()-3);
        //mtx_eigenvector_inv.row(mtx_data.cols()-1) = eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1);
        //mtx_eigenvector_inv.row(mtx_data.cols()-2) = eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2);
        //mtx_eigenvector_inv.row(mtx_data.cols()-3) = eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3);
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
        //NormalizeHist2D(Hist_Data,&Hist_Bkgd_MSCLW.at(e));
        std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-1):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-1) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-2):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-2) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-3):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-3) << std::endl;
        std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-1):" << std::endl;
        std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-1) << std::endl;
        std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-2):" << std::endl;
        std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-2) << std::endl;
        std::cout << "eigensolver_dark.eigenvalues()(mtx_dark.cols()-3):" << std::endl;
        std::cout << eigensolver_dark.eigenvalues()(mtx_dark.cols()-3) << std::endl;

    }
    InputDataFile.Close();


    TFile OutputFile("output_Jul16/Netflix_"+TString(target_data)+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_Theta2"+std::to_string(int(10.*Theta2_cut_lower))+"to"+std::to_string(int(10.*Theta2_cut_upper))+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_Bkgd_MSCLW.at(e).Write();
        Hist_Data_EigenvectorReal_0.at(e).Write();
        Hist_Fit_EigenvectorReal_0.at(e).Write();
        Hist_Dark_EigenvectorReal_0.at(e).Write();
        Hist_Data_EigenvectorReal_1.at(e).Write();
        Hist_Fit_EigenvectorReal_1.at(e).Write();
        Hist_Dark_EigenvectorReal_1.at(e).Write();
        Hist_Data_EigenvectorReal_2.at(e).Write();
        Hist_Fit_EigenvectorReal_2.at(e).Write();
        Hist_Dark_EigenvectorReal_2.at(e).Write();
        Hist_Data_InvEigenvectorReal_0.at(e).Write();
        Hist_Fit_InvEigenvectorReal_0.at(e).Write();
        Hist_Dark_InvEigenvectorReal_0.at(e).Write();
        Hist_Data_InvEigenvectorReal_1.at(e).Write();
        Hist_Fit_InvEigenvectorReal_1.at(e).Write();
        Hist_Dark_InvEigenvectorReal_1.at(e).Write();
        Hist_Data_InvEigenvectorReal_2.at(e).Write();
        Hist_Fit_InvEigenvectorReal_2.at(e).Write();
        Hist_Dark_InvEigenvectorReal_2.at(e).Write();
        Hist_Data_EigenvectorImag_0.at(e).Write();
        Hist_Fit_EigenvectorImag_0.at(e).Write();
        Hist_Dark_EigenvectorImag_0.at(e).Write();
        Hist_Data_EigenvectorImag_1.at(e).Write();
        Hist_Fit_EigenvectorImag_1.at(e).Write();
        Hist_Dark_EigenvectorImag_1.at(e).Write();
        Hist_Data_EigenvectorImag_2.at(e).Write();
        Hist_Fit_EigenvectorImag_2.at(e).Write();
        Hist_Dark_EigenvectorImag_2.at(e).Write();
        Hist_Data_InvEigenvectorImag_0.at(e).Write();
        Hist_Fit_InvEigenvectorImag_0.at(e).Write();
        Hist_Dark_InvEigenvectorImag_0.at(e).Write();
        Hist_Data_InvEigenvectorImag_1.at(e).Write();
        Hist_Fit_InvEigenvectorImag_1.at(e).Write();
        Hist_Dark_InvEigenvectorImag_1.at(e).Write();
        Hist_Data_InvEigenvectorImag_2.at(e).Write();
        Hist_Fit_InvEigenvectorImag_2.at(e).Write();
        Hist_Dark_InvEigenvectorImag_2.at(e).Write();
    }
    OutputFile.Close();

}
