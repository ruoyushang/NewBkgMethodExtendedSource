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
double CurrentEnergy = 0.;
MatrixXcd mtx_data(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data_redu(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
vector<TSpline3> spline_eigenvec_real;
vector<TSpline3> spline_eigenvec_imag;
vector<TSpline3> spline_eigenvec_inv_real;
vector<TSpline3> spline_eigenvec_inv_imag;

void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,max(0.,mtx(binx,biny).real()));
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
VectorXd fillVector(TH1D* hist)
{
    VectorXd vtr(hist->GetNbinsX());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        vtr(binx) = hist->GetBinContent(binx+1);
    }
    return vtr;
}
void RestoreEdge(TH2D* Hist_Origin, TH2D* Hist_Modify, double edge_x, double edge_y)
{
    int edge_binx = Hist_Origin->GetXaxis()->FindBin(edge_x);
    int edge_biny = Hist_Origin->GetYaxis()->FindBin(edge_y);
    std::cout << "edge_binx = " << edge_binx <<std::endl;
    std::cout << "edge_biny = " << edge_biny <<std::endl;
    for (int binx=1;binx<=Hist_Origin->GetNbinsX();binx++)
    {
        for (int biny=1;biny<=Hist_Origin->GetNbinsY();biny++)
        {
            if (binx<edge_binx || biny<edge_biny)
            {
                Hist_Modify->SetBinContent(binx,biny,Hist_Origin->GetBinContent(binx,biny));
            }
        }
    }
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
void SmoothEigenvectors(MatrixXcd mtx, MatrixXcd mtx_inv)
{
    spline_eigenvec_real.clear();
    spline_eigenvec_imag.clear();
    spline_eigenvec_inv_real.clear();
    spline_eigenvec_inv_imag.clear();
    TH1D hist_eigenvec_real = TH1D("hist_eigenvec_real","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D hist_eigenvec_imag = TH1D("hist_eigenvec_imag","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D hist_eigenvec_inv_real = TH1D("hist_eigenvec_inv_real","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D hist_eigenvec_inv_imag = TH1D("hist_eigenvec_inv_imag","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        fill1DHistogram(&hist_eigenvec_real,mtx.col(N_bins_for_deconv-NthEigenvector).real());
        fill1DHistogram(&hist_eigenvec_imag,mtx.col(N_bins_for_deconv-NthEigenvector).imag());
        fill1DHistogram(&hist_eigenvec_inv_real,mtx_inv.row(N_bins_for_deconv-NthEigenvector).real());
        fill1DHistogram(&hist_eigenvec_inv_imag,mtx_inv.row(N_bins_for_deconv-NthEigenvector).imag());
        spline_eigenvec_real.push_back(TSpline3(&hist_eigenvec_real));
        spline_eigenvec_imag.push_back(TSpline3(&hist_eigenvec_imag));
        spline_eigenvec_inv_real.push_back(TSpline3(&hist_eigenvec_inv_real));
        spline_eigenvec_inv_imag.push_back(TSpline3(&hist_eigenvec_inv_imag));
    }
}
double CalculateSignificance(double s,double b,double err)
{
    if ((b*b+(s+b)*err*err)==0.) return 0.;
    if ((s+b)*(b+err*err)==0.) return 0.;
    double first_term = (s+b)*log((s+b)*(b+err*err)/(b*b+(s+b)*err*err));
    double second_term = 0.;
    if (err>0. && b>0)
    {
        second_term = b*b/(err*err)*log(1.+err*err*s/(b*(b+err*err)));
    }
    else 
    {
        second_term = 0.;
    }
    double result = 0.;
    if (first_term>second_term) result = pow(2*(first_term-second_term),0.5);
    else result = pow(2*(-first_term+second_term),0.5);
    if (s>0) return result;
    else return -1.*result;
}

double UnblindedChi2(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind*2);
    int biny_upper = hist_data->GetYaxis()->FindBin(MSCW_cut_blind*2);
    int binx_lower = hist_data->GetXaxis()->FindBin(-0.5);
    int biny_lower = hist_data->GetYaxis()->FindBin(-0.5);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            double data = hist_data->GetBinContent(bx,by);
            double dark = hist_dark->GetBinContent(bx,by);
            double model = hist_model->GetBinContent(bx,by);
            double weight = 1.;
            double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(-0.5);
            double dy = hist_data->GetYaxis()->GetBinCenter(by)-(-0.5);
            double width = 1.0;
            weight = exp(-0.5*dx*dx/(1.0*width*1.0*width))*exp(-0.5*dy*dy/(1.0*width*1.0*width));
            chi2 += weight*pow(data-model,2);
        }
    }
    return chi2;
}

double BlindedChi2(TH2D* hist_data, TH2D* hist_data_redu, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind*2);
    int biny_upper = hist_data->GetYaxis()->FindBin(MSCW_cut_blind*2);
    int binx_lower = hist_data->GetXaxis()->FindBin(-0.5);
    int biny_lower = hist_data->GetYaxis()->FindBin(-0.5);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            double data = hist_data->GetBinContent(bx,by);
            double data_redu = hist_data_redu->GetBinContent(bx,by);
            double dark = hist_dark->GetBinContent(bx,by);
            double model = hist_model->GetBinContent(bx,by);
            double weight = 1.;
            double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(-1.);
            double dy = hist_data->GetYaxis()->GetBinCenter(by)-(-1.);
            double width = 1.0;
            weight = exp(-0.5*dx*dx/(1.0*width*1.0*width))*exp(-0.5*dy*dy/(1.0*width*1.0*width));
            if (abs(data-data_redu)/(data+1.)>10.) weight = 0.;
            if (bx>=binx_blind || by>=biny_blind)
            {
                chi2 += weight*pow(data_redu-model,2);
            }
            else if (bx<binx_lower || by<biny_lower)
            {
                chi2 += weight*pow(dark-model,2);
            }
            //else
            //{
            //    chi2 += 0.01*weight*pow(dark-model,2);
            //}
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
    TH1D hist_eigenvec = TH1D("hist_eigenvec","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
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
        first_index = (4*NthEigenvector-4)*2+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double xx = hist_eigenvec.GetBinCenter(row+1);
            if (xx+par[first_index+1]<MSCL_plot_lower) continue;
            if (xx+par[first_index+1]>MSCL_plot_upper) continue;
            mtx_eigenvector(row,col_fix) += par[first_index+0]*spline_eigenvec_real.at(NthEigenvector-1).Eval(xx+par[first_index+1]);
        }
        first_index = (4*NthEigenvector-3)*2+(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double xx = hist_eigenvec.GetBinCenter(row+1);
            if (xx+par[first_index+1]<MSCL_plot_lower) continue;
            if (xx+par[first_index+1]>MSCL_plot_upper) continue;
            mtx_eigenvector(row,col_fix) += If*par[first_index+0]*spline_eigenvec_imag.at(NthEigenvector-1).Eval(xx+par[first_index+1]);
        }
        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-2)*2+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double xx = hist_eigenvec.GetBinCenter(col+1);
            if (xx+par[first_index+1]<MSCL_plot_lower) continue;
            if (xx+par[first_index+1]>MSCL_plot_upper) continue;
            mtx_eigenvector_inv(row_fix,col) += par[first_index+0]*spline_eigenvec_inv_real.at(NthEigenvector-1).Eval(xx+par[first_index+1]);
        }
        first_index = (4*NthEigenvector-1)*2+(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double xx = hist_eigenvec.GetBinCenter(col+1);
            if (xx+par[first_index+1]<MSCL_plot_lower) continue;
            if (xx+par[first_index+1]>MSCL_plot_upper) continue;
            mtx_eigenvector_inv(row_fix,col) += If*par[first_index+0]*spline_eigenvec_inv_imag.at(NthEigenvector-1).Eval(xx+par[first_index+1]);
        }
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*2+(NthEigenvector-1);
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

    //for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    for (int NthEigenvector=1;NthEigenvector<=N_bins_for_deconv;NthEigenvector++)
    {
        mtx_eigenvector.col(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector);
        mtx_eigenvector_inv.row(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector);
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
            double x_norm = par[first_index]*2.0*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
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
            double x_norm = par[first_index]*2.0*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
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
            double x_norm = par[first_index]*2.0*M_PI*(x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
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
            double x_norm = par[first_index]*2.0*M_PI*(x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
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

void FillFourierParameters(const double *par, int NthEigenvector, TH1D* Hist_EigenvecReal, TH1D* Hist_EigenvecImag, TH1D* Hist_InvEigenvecReal, TH1D* Hist_InvEigenvecImag)
{

    // build eigenvector matrix
    int first_index = 0;
    first_index = (4*NthEigenvector-4)*(1+2*n_fourier_modes)+(NthEigenvector-1);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Hist_EigenvecReal->SetBinContent(mode,pow(pow(par[first_index+2*mode-1],2)+pow(par[first_index+2*mode],2),0.5));
    }
    first_index = (4*NthEigenvector-3)*(1+2*n_fourier_modes)+(NthEigenvector-1);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Hist_EigenvecImag->SetBinContent(mode,pow(pow(par[first_index+2*mode-1],2)+pow(par[first_index+2*mode],2),0.5));
    }
    first_index = (4*NthEigenvector-2)*(1+2*n_fourier_modes)+(NthEigenvector-1);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Hist_InvEigenvecReal->SetBinContent(mode,pow(pow(par[first_index+2*mode-1],2)+pow(par[first_index+2*mode],2),0.5));
    }
    first_index = (4*NthEigenvector-1)*(1+2*n_fourier_modes)+(NthEigenvector-1);
    for (int mode=1;mode<=n_fourier_modes;mode++)
    {
        Hist_InvEigenvecImag->SetBinContent(mode,pow(pow(par[first_index+2*mode-1],2)+pow(par[first_index+2*mode],2),0.5));
    }

}



double FourierChi2Function(const double *par)
{

    FourierParametrizeEigenvectors(par);
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_model = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_data_redu = TH2D("hist_data_redu","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark = TH2D("hist_dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    fill2DHistogram(&hist_data,mtx_data);
    fill2DHistogram(&hist_data_redu,mtx_data_redu);
    fill2DHistogram(&hist_dark,mtx_dark);
    fill2DHistogram(&hist_model,mtx_model);

    double chi2 = BlindedChi2(&hist_data,&hist_data_redu,&hist_dark,&hist_model);
    //double chi2 = UnblindedChi2(&hist_data,&hist_dark,&hist_model);

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
    double x_norm = 2.0*M_PI*(xx-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
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
void SetInitialEigenvectors()
{
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    TH1D Hist_Dark_Eigenvector = TH1D("Hist_Dark_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    const std::complex<double> If(0.0, 1.0);
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale_real = 1.;
    double scale_imag = 1.;
    double scale_inv_real = 1.;
    double scale_inv_imag = 1.;

    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //for (int NthEigenvector=1;NthEigenvector<=N_bins_for_deconv;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        scale_real = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        scale_imag = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).real());
        scale_inv_real = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).imag());
        scale_inv_imag = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());

        std::cout << "Nth Eigenvec = " << NthEigenvector << std::endl;
        std::cout << "scale_real = " << scale_real << std::endl;
        std::cout << "scale_imag = " << scale_imag << std::endl;
        std::cout << "scale_inv_real = " << scale_inv_real << std::endl;
        std::cout << "scale_inv_imag = " << scale_inv_imag << std::endl;

        scale_real = 1.;
        scale_imag = 1.;
        scale_inv_real = 1.;
        scale_inv_imag = 1.;

        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = scale_real*eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).real();
            mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) += If*scale_imag*eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).imag();
        }
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = scale_inv_real*eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).real();
            mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) += If*scale_inv_imag*eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).imag();
        }
    }
}
void SetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_dark, MatrixXcd mtx_data)
{
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    TH1D Hist_Dark_Eigenvector = TH1D("Hist_Dark_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    const std::complex<double> If(0.0, 1.0);
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale_real = 1.;
    double scale_imag = 1.;
    double scale_inv_real = 1.;
    double scale_inv_imag = 1.;
    double step_size = 0.001;

    // 1st eigenvector
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        scale_real = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        double eigenvec_real_range = 0.1*double(Hist_Dark_Eigenvector.Integral());
        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        scale_imag = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        double eigenvec_imag_range = 0.1*double(Hist_Dark_Eigenvector.Integral());

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).real());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).real());
        scale_inv_real = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        double eigenvec_inv_real_range = 0.1*double(Hist_Dark_Eigenvector.Integral());
        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).imag());
        fill1DHistogram(&Hist_Data_Eigenvector,eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).imag());
        scale_inv_imag = double(Hist_Data_Eigenvector.Integral())/double(Hist_Dark_Eigenvector.Integral());
        double eigenvec_inv_imag_range = 0.1*double(Hist_Dark_Eigenvector.Integral());

        std::cout << "Nth Eigenvec = " << NthEigenvector << std::endl;
        std::cout << "scale_real = " << scale_real << std::endl;
        std::cout << "scale_imag = " << scale_imag << std::endl;
        std::cout << "scale_inv_real = " << scale_inv_real << std::endl;
        std::cout << "scale_inv_imag = " << scale_inv_imag << std::endl;
        scale_real = 1.;
        scale_imag = 1.;
        scale_inv_real = 1.;
        scale_inv_imag = 1.;

        double eigenval_upper_limit = 1.2;
        double eigenval_lower_limit = 0.8;

        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = (4*NthEigenvector-4)*2+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index+0,"par["+std::to_string(int(first_index+0))+"]",1,step_size);
        Chi2Minimizer->SetVariableLimits(first_index+0,1,1);
        Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",0,step_size);
        first_index = (4*NthEigenvector-3)*2+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index+0,"par["+std::to_string(int(first_index+0))+"]",1,step_size);
        Chi2Minimizer->SetVariableLimits(first_index+0,1,1);
        Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",0,step_size);
        first_index = (4*NthEigenvector-2)*2+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index+0,"par["+std::to_string(int(first_index+0))+"]",1,step_size);
        Chi2Minimizer->SetVariableLimits(first_index+0,1,1);
        Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",0,step_size);
        first_index = (4*NthEigenvector-1)*2+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index+0,"par["+std::to_string(int(first_index+0))+"]",1,step_size);
        Chi2Minimizer->SetVariableLimits(first_index+0,1,1);
        Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",0,step_size);
        // build eigenvalue matrix
        first_index = (4*NthEigenvector-0)*2+(NthEigenvector-1);
        double initial_eigenvalue = eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real();
        Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", initial_eigenvalue,step_size);
        //Chi2Minimizer->SetVariableLimits(first_index,initial_eigenvalue*eigenval_lower_limit,initial_eigenvalue*eigenval_upper_limit);
    }

}
void FourierSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, MatrixXcd mtx_dark, MatrixXcd mtx_data)
{
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    TH1D Hist_Dark_Eigenvector = TH1D("Hist_Dark_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    TH1D Hist_Data_Eigenvector = TH1D("Hist_Data_Eigenvector","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale = 0.;
    double limit = 0.5; // 0.0-1.0

    // 1st eigenvector
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).real());
        scale = double(Hist_Dark_Eigenvector.Integral());
        first_index = (4*NthEigenvector-4)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",1.,0.001);
        Chi2Minimizer->SetVariableLimits(first_index,0.5,2.0);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            if (mode<NthEigenvector) limit = 0.;
            else limit = 0.5*scale;
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]",0.,0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,-limit,limit);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]",0.,0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,-limit,limit);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).imag());
        scale = double(Hist_Dark_Eigenvector.Integral());
        first_index = (4*NthEigenvector-3)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",1.,0.001);
        Chi2Minimizer->SetVariableLimits(first_index,0.5,2.0);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            limit = 0.;
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,-limit,limit);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,-limit,limit);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).real());
        scale = double(Hist_Dark_Eigenvector.Integral());
        first_index = (4*NthEigenvector-2)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",1.,0.001);
        Chi2Minimizer->SetVariableLimits(first_index,0.5,2.0);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            if (mode<NthEigenvector) limit = 0.;
            else limit = 0.5*scale;
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,-limit,limit);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,-limit,limit);
        }

        fill1DHistogram(&Hist_Dark_Eigenvector,eigensolver_dark.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector).imag());
        scale = double(Hist_Dark_Eigenvector.Integral());
        first_index = (4*NthEigenvector-1)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",1.,0.001);
        Chi2Minimizer->SetVariableLimits(first_index,0.5,2.0);
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            limit = 0.;
            Chi2Minimizer->SetVariable(first_index+2*mode-1,"par["+std::to_string(int(first_index+2*mode-1))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode-1,-limit,limit);
            Chi2Minimizer->SetVariable(first_index+2*mode,"par["+std::to_string(int(first_index+2*mode))+"]",0., 0.001);
            Chi2Minimizer->SetVariableLimits(first_index+2*mode,-limit,limit);
        }

        // eigenvalues
        first_index = (4*NthEigenvector-0)*(1+2*n_fourier_modes)+(NthEigenvector-1);
        Chi2Minimizer->SetVariable(first_index, "par["+std::to_string(int(first_index))+"]", eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real(), 0.01);
        //Chi2Minimizer->SetVariableLimits(first_index,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real()*1.0,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real()*1.0);
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
    vector<TH2D> Hist_Redu_MSCLW;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_0_Real;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_0_Imag;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_0_Real;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_0_Imag;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_1_Real;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_1_Imag;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_1_Real;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_1_Imag;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_2_Real;
    vector<TH1D> Hist_FourierCoeff_Eigenvector_2_Imag;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_2_Real;
    vector<TH1D> Hist_FourierCoeff_InvEigenvector_2_Imag;
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
        CurrentEnergy = energy_bins[e];
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        TString filename_data  = "Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_dark  = "Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
        Hist_Bkgd_MSCLW.push_back(TH2D("Hist_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Redu_MSCLW.push_back(TH2D("Hist_Redu_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_FourierCoeff_Eigenvector_0_Real.push_back(TH1D("Hist_FourierCoeff_Eigenvector_0_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_Eigenvector_0_Imag.push_back(TH1D("Hist_FourierCoeff_Eigenvector_0_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_0_Real.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_0_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_0_Imag.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_0_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_Eigenvector_1_Real.push_back(TH1D("Hist_FourierCoeff_Eigenvector_1_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_Eigenvector_1_Imag.push_back(TH1D("Hist_FourierCoeff_Eigenvector_1_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_1_Real.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_1_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_1_Imag.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_1_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_Eigenvector_2_Real.push_back(TH1D("Hist_FourierCoeff_Eigenvector_2_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_Eigenvector_2_Imag.push_back(TH1D("Hist_FourierCoeff_Eigenvector_2_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_2_Real.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_2_Real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
        Hist_FourierCoeff_InvEigenvector_2_Imag.push_back(TH1D("Hist_FourierCoeff_InvEigenvector_2_Imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",6,0,6));
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

        //if (energy_bins[e]<237) continue;
        //if (energy_bins[e]>=282) continue;
        //if (energy_bins[e]<335) continue;
        //if (energy_bins[e]>=398) continue;
        //if (energy_bins[e]<562) continue;
        //if (energy_bins[e]>=794) continue;

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
        NumberOfEigenvectors = 3;
        eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);

        MatrixXcd mtx_eigenval_redu(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        for (int i=0;i<mtx_dark.cols();i++)
        {
            for (int j=0;j<mtx_dark.rows();j++)
            {
                if (i==j && i>=mtx_dark.cols()-NumberOfEigenvectors) 
                {
                    mtx_eigenval_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_redu*eigensolver_data.eigenvectors().inverse();

        fill1DHistogram(&Hist_Data_EigenvectorReal_0.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Data_EigenvectorReal_1.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Data_EigenvectorReal_2.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-3).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_0.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_1.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2).real());
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_2.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3).real());
        //FourierFitEigenvector(&Hist_Data_EigenvectorReal_0.at(e),&Hist_Fit_EigenvectorReal_0.at(e));
        //FourierFitEigenvector(&Hist_Data_EigenvectorReal_1.at(e),&Hist_Fit_EigenvectorReal_1.at(e));
        //FourierFitEigenvector(&Hist_Data_EigenvectorReal_2.at(e),&Hist_Fit_EigenvectorReal_2.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_0.at(e),&Hist_Fit_InvEigenvectorReal_0.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_1.at(e),&Hist_Fit_InvEigenvectorReal_1.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorReal_2.at(e),&Hist_Fit_InvEigenvectorReal_2.at(e));
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
        //FourierFitEigenvector(&Hist_Data_EigenvectorImag_0.at(e),&Hist_Fit_EigenvectorImag_0.at(e));
        //FourierFitEigenvector(&Hist_Data_EigenvectorImag_1.at(e),&Hist_Fit_EigenvectorImag_1.at(e));
        //FourierFitEigenvector(&Hist_Data_EigenvectorImag_2.at(e),&Hist_Fit_EigenvectorImag_2.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_0.at(e),&Hist_Fit_InvEigenvectorImag_0.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_1.at(e),&Hist_Fit_InvEigenvectorImag_1.at(e));
        //FourierFitEigenvector(&Hist_Data_InvEigenvectorImag_2.at(e),&Hist_Fit_InvEigenvectorImag_2.at(e));
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

        NumberOfEigenvectors = 3;
        if (CurrentEnergy>1000.) NumberOfEigenvectors = 2;
        if (CurrentEnergy>2000.) NumberOfEigenvectors = 1;
        ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kVectorBFGS );
        //ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kSteepestDescent );
        //ROOT::Minuit2::Minuit2Minimizer Chi2Minimizer_1st( ROOT::Minuit2::kMigrad );
        //Chi2Minimizer_1st.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
        //Chi2Minimizer_1st.SetMaxIterations(100000); // for GSL
        //Chi2Minimizer_1st.SetTolerance(0.001);
        Chi2Minimizer_1st.SetMaxFunctionCalls(200); // for Minuit/Minuit2
        Chi2Minimizer_1st.SetMaxIterations(200); // for GSL
        Chi2Minimizer_1st.SetTolerance(0.001);

        SetInitialEigenvectors();
        SmoothEigenvectors(mtx_eigenvector_init, mtx_eigenvector_inv_init);

        //ROOT::Math::Functor Chi2Func_1st(&Chi2Function,4*NumberOfEigenvectors*2+NumberOfEigenvectors);
        //Chi2Minimizer_1st.SetFunction(Chi2Func_1st);
        //SetInitialVariables(&Chi2Minimizer_1st, mtx_dark, mtx_data);
        //Chi2Minimizer_1st.Minimize();
        //const double *par = Chi2Minimizer_1st.X();
        //ParametrizeEigenvectors(par);
        //mtx_data_bkgd = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

        n_fourier_modes = NumberOfEigenvectors;
        ROOT::Math::Functor Chi2Func_1st(&FourierChi2Function,4*NumberOfEigenvectors*(1+2*n_fourier_modes)+NumberOfEigenvectors);        
        Chi2Minimizer_1st.SetFunction(Chi2Func_1st);
        FourierSetInitialVariables(&Chi2Minimizer_1st, mtx_dark, mtx_data);
        Chi2Minimizer_1st.Minimize();
        const double *par = Chi2Minimizer_1st.X();
        FourierParametrizeEigenvectors(par);

        FillFourierParameters(par,1,&Hist_FourierCoeff_Eigenvector_0_Real.at(e),&Hist_FourierCoeff_Eigenvector_0_Imag.at(e),&Hist_FourierCoeff_InvEigenvector_0_Real.at(e),&Hist_FourierCoeff_InvEigenvector_0_Imag.at(e));
        FillFourierParameters(par,2,&Hist_FourierCoeff_Eigenvector_1_Real.at(e),&Hist_FourierCoeff_Eigenvector_1_Imag.at(e),&Hist_FourierCoeff_InvEigenvector_1_Real.at(e),&Hist_FourierCoeff_InvEigenvector_1_Imag.at(e));
        FillFourierParameters(par,3,&Hist_FourierCoeff_Eigenvector_2_Real.at(e),&Hist_FourierCoeff_Eigenvector_2_Imag.at(e),&Hist_FourierCoeff_InvEigenvector_2_Real.at(e),&Hist_FourierCoeff_InvEigenvector_2_Imag.at(e));
        mtx_data_bkgd = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

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
        fill2DHistogram(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        fill2DHistogram(&Hist_Redu_MSCLW.at(e),mtx_data_redu);
        //RestoreEdge(Hist_Dark,&Hist_Bkgd_MSCLW.at(e),-0.5,-0.5);
        //NormalizeHist2D(Hist_Data,&Hist_Bkgd_MSCLW.at(e));
        std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral() << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-1):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-1) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-2):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-2) << std::endl;
        std::cout << "eigensolver_data.eigenvalues()(mtx_data.cols()-3):" << std::endl;
        std::cout << eigensolver_data.eigenvalues()(mtx_data.cols()-3) << std::endl;
        std::cout << "mtx_eigenvalue(mtx_data.cols()-1):" << std::endl;
        std::cout << mtx_eigenvalue.row(mtx_data.cols()-1)(mtx_data.cols()-1) << std::endl;
        std::cout << "mtx_eigenvalue(mtx_data.cols()-2):" << std::endl;
        std::cout << mtx_eigenvalue.row(mtx_data.cols()-2)(mtx_data.cols()-2) << std::endl;
        std::cout << "mtx_eigenvalue(mtx_data.cols()-3):" << std::endl;
        std::cout << mtx_eigenvalue.row(mtx_data.cols()-3)(mtx_data.cols()-3) << std::endl;
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
        Hist_Redu_MSCLW.at(e).Write();
        Hist_FourierCoeff_Eigenvector_0_Real.at(e).Write();
        Hist_FourierCoeff_Eigenvector_0_Imag.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_0_Real.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_0_Imag.at(e).Write();
        Hist_FourierCoeff_Eigenvector_1_Real.at(e).Write();
        Hist_FourierCoeff_Eigenvector_1_Imag.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_1_Real.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_1_Imag.at(e).Write();
        Hist_FourierCoeff_Eigenvector_2_Real.at(e).Write();
        Hist_FourierCoeff_Eigenvector_2_Imag.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_2_Real.at(e).Write();
        Hist_FourierCoeff_InvEigenvector_2_Imag.at(e).Write();
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
