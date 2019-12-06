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

int NthIteration = 0;
int NumberOfEigenvectors = 0;
int n_fourier_modes = 6;
int n_taylor_modes = 6;
double CurrentEnergy = 0.;
MatrixXcd mtx_gamma(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data_redu(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark_redu(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_data_err(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_dark_err(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
VectorXcd vtr_fit(N_bins_for_deconv);
VectorXcd vtr_data(N_bins_for_deconv);
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
double init_deriv_at_zero;
bool signalfree_model;

void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,mtx(binx,biny).real());
        }
    }
}
void fill2DHistogramAbs(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,max(0.,mtx(binx,biny).real()));
        }
    }
}
void Smooth2DHistogram(TH2D* hist_origin, TH2D* hist_smooth, int size)
{
    for (int binx=0;binx<hist_origin->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist_origin->GetNbinsY();biny++)
        {
            int binx_low = max(1,binx+1-size);
            int binx_up = min(hist_origin->GetNbinsX(),binx+1+size);
            int biny_low = max(1,biny+1-size);
            int biny_up = min(hist_origin->GetNbinsY(),biny+1+size);
            hist_smooth->SetBinContent(binx+1,biny+1,hist_origin->Integral(binx_low,binx_up,biny_low,biny_up));
        }
    }
}
void fill2DHistogramError(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinError(binx+1,biny+1,mtx(binx,biny).real());
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
MatrixXcd fillMatrixError(TH2D* hist)
{
    MatrixXcd matrix(hist->GetNbinsX(),hist->GetNbinsY());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            matrix(binx,biny) = hist->GetBinError(binx+1,biny+1);
        }
    }
    return matrix;
}
double SmoothEigenvectors(MatrixXcd* mtx, MatrixXcd* mtx_inv)
{
    TH1D hist_ref = TH1D("hist_ref","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    const std::complex<double> If(0.0, 1.0);
    int n_rebin = 4;
    int n_rebin_inv = 4;
    double extension = (MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)*double(n_rebin);
    double extension_inv = (MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)*double(n_rebin_inv);
    double deriv_at_zero = 0.;
    double deriv_at_zero_norm = 0.;
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        TH1D eigenvec_real = TH1D("eigenvec_real","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
        TH1D eigenvec_imag = TH1D("eigenvec_imag","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
        TH1D eigenvec_inv_real = TH1D("eigenvec_inv_real","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
        TH1D eigenvec_inv_imag = TH1D("eigenvec_inv_imag","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
        fill1DHistogram(&eigenvec_real,mtx->col(N_bins_for_deconv-NthEigenvector).real());
        fill1DHistogram(&eigenvec_imag,mtx->col(N_bins_for_deconv-NthEigenvector).imag());
        fill1DHistogram(&eigenvec_inv_real,mtx_inv->row(N_bins_for_deconv-NthEigenvector).real());
        fill1DHistogram(&eigenvec_inv_imag,mtx_inv->row(N_bins_for_deconv-NthEigenvector).imag());
        TH1D eigenvec_real_ext = TH1D("eigenvec_real_ext","",N_bins_for_deconv+n_rebin,MSCL_plot_lower-extension,MSCL_plot_upper);
        TH1D eigenvec_imag_ext = TH1D("eigenvec_imag_ext","",N_bins_for_deconv+n_rebin,MSCL_plot_lower-extension,MSCL_plot_upper);
        TH1D eigenvec_inv_real_ext = TH1D("eigenvec_inv_real_ext","",N_bins_for_deconv+n_rebin_inv,MSCL_plot_lower-extension_inv,MSCL_plot_upper);
        TH1D eigenvec_inv_imag_ext = TH1D("eigenvec_inv_imag_ext","",N_bins_for_deconv+n_rebin_inv,MSCL_plot_lower-extension_inv,MSCL_plot_upper);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            eigenvec_real_ext.SetBinContent(row+1+n_rebin,eigenvec_real.GetBinContent(row+1));
            eigenvec_real_ext.SetBinError(row+1+n_rebin,eigenvec_real.GetBinError(row+1));
            eigenvec_imag_ext.SetBinContent(row+1+n_rebin,eigenvec_imag.GetBinContent(row+1));
            eigenvec_imag_ext.SetBinError(row+1+n_rebin,eigenvec_imag.GetBinError(row+1));
            eigenvec_inv_real_ext.SetBinContent(row+1+n_rebin_inv,eigenvec_inv_real.GetBinContent(row+1));
            eigenvec_inv_real_ext.SetBinError(row+1+n_rebin_inv,eigenvec_inv_real.GetBinError(row+1));
            eigenvec_inv_imag_ext.SetBinContent(row+1+n_rebin_inv,eigenvec_inv_imag.GetBinContent(row+1));
            eigenvec_inv_imag_ext.SetBinError(row+1+n_rebin_inv,eigenvec_inv_imag.GetBinError(row+1));
        }
        eigenvec_real_ext.Rebin(n_rebin);
        eigenvec_imag_ext.Rebin(n_rebin);
        TSpline3 spline_eigenvec_real(&eigenvec_real_ext);
        TSpline3 spline_eigenvec_imag(&eigenvec_imag_ext);
        int col_fix = N_bins_for_deconv-NthEigenvector;
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            double xx = hist_ref.GetBinCenter(row+1);
            mtx->operator()(row,col_fix) = spline_eigenvec_real.Eval(xx)/double(n_rebin); 
            mtx->operator()(row,col_fix) += If*spline_eigenvec_imag.Eval(xx)/double(n_rebin);
        }
        eigenvec_inv_real_ext.Rebin(n_rebin_inv);
        eigenvec_inv_imag_ext.Rebin(n_rebin_inv);
        TSpline3 spline_eigenvec_inv_real(&eigenvec_inv_real_ext);
        TSpline3 spline_eigenvec_inv_imag(&eigenvec_inv_imag_ext);
        int row_fix = N_bins_for_deconv-NthEigenvector;
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            double xx = hist_ref.GetBinCenter(col+1);
            mtx_inv->operator()(row_fix,col) = spline_eigenvec_inv_real.Eval(xx)/double(n_rebin_inv); 
            mtx_inv->operator()(row_fix,col) += If*spline_eigenvec_inv_imag.Eval(xx)/double(n_rebin_inv);
        }
    }
    return 0.;
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

double SignalLogLikelihood(TH2D* hist_data, TH2D* hist_gamma, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(1.);
    int biny_blind = hist_data->GetYaxis()->FindBin(1.);
    double log_likelihood = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx<binx_blind && by<biny_blind)
            {
                double data = hist_data->GetBinContent(bx,by);
                double gamma = hist_gamma->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                model += gamma;
                double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
                double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
                //double weight = exp(-0.5*dx*dx-0.5*dy*dy);
                double weight = 1.;
                model = max(0.001,model);
                double log_factorial = 0.;
                if (data>0) log_factorial = data*(log(data)-1.) + 0.5*log(2*M_PI*data);
                double log_likelihood_this = -2.*(data*log(model)-log_factorial-model);
                if (isnan(log_likelihood_this) || isinf(log_likelihood_this))
                {
                    std::cout << "log_likelihood = nan (or inf) !!!" << std::endl;
                    std::cout << "model = " << model << std::endl;
                    std::cout << "data = " << data << std::endl;
                }
                log_likelihood += weight*log_likelihood_this;
            }
        }
    }
    return log_likelihood;
}

double BlindedLogLikelihood(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(1.);
    int biny_blind = hist_data->GetYaxis()->FindBin(1.);
    double log_likelihood = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>=binx_blind || by>=biny_blind)
            {
                double data = hist_data->GetBinContent(bx,by);
                double dark = hist_dark->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                //model = model*dark;
                double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
                double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
                //double weight = exp(-0.5*dx*dx-0.5*dy*dy);
                double weight = 1.;
                model = max(0.001,model);
                double log_factorial = 0.;
                if (data>0) log_factorial = data*(log(data)-1.) + 0.5*log(2*M_PI*data);
                double log_likelihood_this = -2.*(data*log(model)-log_factorial-model);
                if (isnan(log_likelihood_this) || isinf(log_likelihood_this))
                {
                    std::cout << "log_likelihood = nan (or inf) !!!" << std::endl;
                    std::cout << "model = " << model << std::endl;
                    std::cout << "data = " << data << std::endl;
                }
                log_likelihood += weight*log_likelihood_this;
            }
        }
    }
    return log_likelihood;
}

double ActivationFunction(double x, double ref)
{
    return max((x-ref)/ref,0.);
}

double FirstDerivative(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind*2);
    int biny_upper = hist_data->GetYaxis()->FindBin(MSCW_cut_blind*2);
    int binx_lower = hist_data->GetXaxis()->FindBin(-0.5);
    int biny_lower = hist_data->GetYaxis()->FindBin(-0.5);
    double deriv_model = 0.;
    double deriv_dark = 0.;
    double n_bins = 0.;
    for (int bx=1+1;bx<=hist_data->GetNbinsX()-1;bx++)
    {
        for (int by=1+1;by<=hist_data->GetNbinsY()-1;by++)
        {
            double data = 0.;
            data += pow(hist_data->GetBinContent(bx+1,by)-hist_data->GetBinContent(bx,by),2);
            data += pow(hist_data->GetBinContent(bx,by)-hist_data->GetBinContent(bx-1,by),2);
            data += pow(hist_data->GetBinContent(bx,by+1)-hist_data->GetBinContent(bx,by),2);
            data += pow(hist_data->GetBinContent(bx,by)-hist_data->GetBinContent(bx,by-1),2);
            double dark = 0.;
            dark += pow(hist_dark->GetBinContent(bx+1,by)-hist_dark->GetBinContent(bx,by),2);
            dark += pow(hist_dark->GetBinContent(bx,by)-hist_dark->GetBinContent(bx-1,by),2);
            dark += pow(hist_dark->GetBinContent(bx,by+1)-hist_dark->GetBinContent(bx,by),2);
            dark += pow(hist_dark->GetBinContent(bx,by)-hist_dark->GetBinContent(bx,by-1),2);
            double model = 0.;
            model += pow(hist_model->GetBinContent(bx+1,by)-hist_model->GetBinContent(bx,by),2);
            model += pow(hist_model->GetBinContent(bx,by)-hist_model->GetBinContent(bx-1,by),2);
            model += pow(hist_model->GetBinContent(bx,by+1)-hist_model->GetBinContent(bx,by),2);
            model += pow(hist_model->GetBinContent(bx,by)-hist_model->GetBinContent(bx,by-1),2);
            double weight = 1.;
            double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
            double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
            double width = 0.1;
            if (bx<binx_blind && by<biny_blind)
            {
                deriv_model += model;
                deriv_dark += dark;
                n_bins += 1.;
            }
        }
    }
    return n_bins*ActivationFunction(deriv_model, deriv_dark);
}

double SignalChi2(TH2D* hist_data, TH2D* hist_gamma, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx<binx_blind && by<biny_blind)
            {
                double data = hist_data->GetBinContent(bx,by);
                double gamma = hist_gamma->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                double weight = 1.;
                double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
                double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
                double width = 0.1;
                //double data_err = max(1.0,pow(data,0.5));
                //double model_err = max(1.0,pow(abs(model),0.5));
                //double gamma_err = max(1.0,pow(abs(gamma),0.5));
                //weight = 1./(data_err*data_err+model_err*model_err+gamma_err*gamma_err);
                //double data_err = max(1.0,pow(data,0.5));
                //double model_err = max(1.0,pow(abs(model+gamma),0.5));
                //weight = 1./(data_err*data_err+model_err*model_err);
                double data_err = max(1.0,pow(max(0.,data-gamma),0.5));
                double model_err = max(1.0,pow(max(0.,model),0.5));
                weight = 1./(data_err*data_err+model_err*model_err);
                double chi2_this = weight*pow(data-model-gamma,2);
                if (isnan(chi2_this))
                {
                    if (isnan(weight)) std::cout << "weight==nan!!!" << std::endl;
                    if (isnan(pow(data-model,2))) std::cout << "pow(data-model,2)==nan!!!" << std::endl;
                    std::cout << "model = " << model << std::endl;
                    std::cout << "data = " << data << std::endl;
                    std::cout << "gamma = " << gamma << std::endl;
                    continue;
                }
                chi2 += chi2_this;
            }
        }
    }
    return chi2;
}

double BlindedChi2(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    //int binx_upper = hist_data->GetXaxis()->FindBin(3.);
    //int biny_upper = hist_data->GetXaxis()->FindBin(1.5);
    double chi2 = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>=binx_blind || by>=biny_blind)
            {
                //if (bx>=binx_upper || by>=biny_upper) continue;
                double data = hist_data->GetBinContent(bx,by);
                double dark = hist_dark->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                //model = model*dark;
                double weight = 1.;
                double dx = hist_data->GetXaxis()->GetBinCenter(bx)-(1.);
                double dy = hist_data->GetYaxis()->GetBinCenter(by)-(1.);
                double width = 0.1;
                double data_err = max(1.0,pow(data,0.5));
                double model_err = max(1.0,pow(abs(model),0.5));
                weight = 1./(data_err*data_err+model_err*model_err);
                double chi2_this = weight*pow(data-model,2);
                if (isnan(chi2_this))
                {
                    if (isnan(weight)) std::cout << "weight==nan!!!" << std::endl;
                    if (isnan(pow(data-model,2))) std::cout << "pow(data-model,2)==nan!!!" << std::endl;
                    std::cout << "model = " << model << std::endl;
                    std::cout << "data = " << data << std::endl;
                    continue;
                }
                if (bx>=binx_blind && by>=biny_blind)
                {
                    //continue;
                    chi2 += chi2_this;
                }
                else if (bx>=binx_blind && by<biny_blind)
                {
                    chi2 += chi2_this;
                }
                else
                {
                    chi2 += chi2_this;
                }
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

    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //for (int NthEigenvector=1;NthEigenvector<=N_bins_for_deconv;NthEigenvector++)
    {
        mtx_eigenvector.col(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector);
        mtx_eigenvector_inv.row(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector);
        mtx_eigenvalue(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector);
    }

    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    //for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //{
    //    col_fix = N_bins_for_deconv-NthEigenvector;
    //    first_index = (2*NthEigenvector-2)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    for (int row=0;row<N_bins_for_deconv;row++)
    //    {
    //        double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
    //        double x_norm = par[first_index]*2.0*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
    //        mtx_eigenvector(row,col_fix) += par[first_index+1]/2.;
    //        for (int mode=1;mode<=n_fourier_modes;mode++)
    //        {
    //            mtx_eigenvector(row,col_fix) += par[first_index+1+2*mode-1]*cos(double(mode)*x_norm);
    //            mtx_eigenvector(row,col_fix) += par[first_index+1+2*mode]*sin(double(mode)*x_norm);
    //        }
    //    }

    //    row_fix = N_bins_for_deconv-NthEigenvector;
    //    first_index = (2*NthEigenvector-1)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    for (int col=0;col<N_bins_for_deconv;col++)
    //    {
    //        double x = double(col+1)*(MSCL_plot_upper-MSCL_plot_lower)/double(N_bins_for_deconv)+MSCL_plot_lower;
    //        double x_norm = par[first_index]*1.0*M_PI*(x-MSCL_plot_lower)/(MSCL_plot_upper-MSCL_plot_lower);
    //        mtx_eigenvector_inv(row_fix,col) += par[first_index+1]/2.;
    //        for (int mode=1;mode<=n_fourier_modes;mode++)
    //        {
    //            mtx_eigenvector_inv(row_fix,col) += par[first_index+1+2*mode-1]*cos(double(mode)*x_norm);
    //            mtx_eigenvector_inv(row_fix,col) += par[first_index+1+2*mode]*sin(double(mode)*x_norm);
    //        }
    //    }

    //    // build eigenvalue matrix
    //    first_index = (2*NthEigenvector-0)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
    //    {
    //        mtx_eigenvalue(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvalue) = par[first_index+NthEigenvalue-1];
    //    }
    //}
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector(row,col_fix) += par[first_index+row];
        }

        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv(row_fix,col) += par[first_index+col];
        }

        // build eigenvalue matrix
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            mtx_eigenvalue(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvalue) += par[first_index+NthEigenvalue-1]*par[0];
        }
    }
}


double FourierChi2Function(const double *par)
{

    FourierParametrizeEigenvectors(par);
    //double deriv_at_zero = SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);
    //if (0.8*init_deriv_at_zero>deriv_at_zero) return 1e10;
    //if (1.2*init_deriv_at_zero<deriv_at_zero) return 1e10;
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_model = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;

    TH2D hist_gamma = TH2D("hist_gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_diff = TH2D("hist_diff","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark = TH2D("hist_dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    fill2DHistogramAbs(&hist_gamma,mtx_gamma);
    fill2DHistogramAbs(&hist_data,mtx_data);
    fill2DHistogramAbs(&hist_dark,mtx_dark);
    fill2DHistogramAbs(&hist_model,mtx_model);

    int binx_lower = hist_diff.GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_diff.GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = hist_diff.GetXaxis()->FindBin(MSCL_cut_upper*2)-1;
    int biny_lower = hist_diff.GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_diff.GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = hist_diff.GetYaxis()->FindBin(MSCW_cut_upper*2)-1;
    hist_diff.Add(&hist_data);
    hist_diff.Add(&hist_model,-1.);
    //double gamma_total = hist_diff.Integral(binx_lower,binx_blind,biny_lower,biny_blind);
    double gamma_total = hist_diff.Integral();
    gamma_total = max(0.,gamma_total);
    //double scale = gamma_total/double(hist_gamma.Integral(binx_lower,binx_blind,biny_lower,biny_blind));
    double scale = gamma_total/double(hist_gamma.Integral());
    hist_gamma.Scale(scale);
    if (signalfree_model) hist_gamma.Scale(0);

    double chi2 = BlindedChi2(&hist_data,&hist_dark,&hist_model);
    chi2 += SignalChi2(&hist_data,&hist_gamma,&hist_model);
    //
    //double chi2 = BlindedLogLikelihood(&hist_data,&hist_dark,&hist_model);
    //chi2 += SignalLogLikelihood(&hist_data,&hist_gamma,&hist_model);
    //
    //std::cout << NthIteration << "-th iteration, chi2 = " << chi2 << std::endl; 
    //NthIteration += 1;

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
            mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row);
            //mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row);
        }
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col);
            //mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col);
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

void FourierSet2ndIterationVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, const double *par, int type)
{

    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale = 0.;
    double limit = 0.; 
    
    //for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //{

    //    first_index = (2*NthEigenvector-2)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    limit = 0.1;
    //    if (type!=2) limit = 0.;
    //    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",par[first_index],0.001);
    //    Chi2Minimizer->SetVariableLimits(first_index,par[first_index],par[first_index]);
    //    Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",par[first_index+1],0.001);
    //    Chi2Minimizer->SetVariableLimits(first_index+1,par[first_index+1]*(0.-limit),par[first_index+1]*(0.+limit));
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        limit = 0.1;
    //        if (type!=2) limit = 0.;
    //        if (mode<NthEigenvector) limit = 0.;
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode-1,"par["+std::to_string(int(first_index+1+2*mode-1))+"]",par[first_index+1+2*mode-1],0.001);
    //        Chi2Minimizer->SetVariableLimits(first_index+1+2*mode-1,par[first_index+1+2*mode-1]*(0.-limit),par[first_index+1+2*mode-1]*(0.+limit));
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode,"par["+std::to_string(int(first_index+1+2*mode))+"]",par[first_index+1+2*mode],0.001);
    //        Chi2Minimizer->SetVariableLimits(first_index+1+2*mode,par[first_index+1+2*mode]*(0.-limit),par[first_index+1+2*mode]*(0.+limit));
    //    }

    //    first_index = (2*NthEigenvector-1)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    limit = 0.1;
    //    if (type!=1) limit = 0.;
    //    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",par[first_index],0.001);
    //    Chi2Minimizer->SetVariableLimits(first_index,par[first_index],par[first_index]);
    //    Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",par[first_index+1],0.001);
    //    Chi2Minimizer->SetVariableLimits(first_index+1,par[first_index+1]*(0.-limit),par[first_index+1]*(0.+limit));
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        limit = 0.1;
    //        if (type!=1) limit = 0.;
    //        if (mode<NthEigenvector) limit = 0.;
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode-1,"par["+std::to_string(int(first_index+1+2*mode-1))+"]",par[first_index+1+2*mode-1],0.001);
    //        Chi2Minimizer->SetVariableLimits(first_index+1+2*mode-1,par[first_index+1+2*mode-1]*(0.-limit),par[first_index+1+2*mode-1]*(0.+limit));
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode,"par["+std::to_string(int(first_index+1+2*mode))+"]",par[first_index+1+2*mode],0.001);
    //        Chi2Minimizer->SetVariableLimits(first_index+1+2*mode,par[first_index+1+2*mode]*(0.-limit),par[first_index+1+2*mode]*(0.+limit));
    //    }

    //    // eigenvalues
    //    double input_value = 0;
    //    first_index = (2*NthEigenvector-0)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
    //    {
    //        input_value = par[first_index+NthEigenvalue-1];
    //        limit = 0.1*par[first_index];
    //        Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.001*limit);
    //        Chi2Minimizer->SetVariableLimits(first_index+NthEigenvalue-1,input_value-limit,input_value+limit);
    //    }
    //}

    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        limit = 0.2;
        if (type!=2) limit = 0.;
        first_index = (4*NthEigenvector-4)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",par[first_index+row],0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+row,-limit,limit);
        }
        first_index = (4*NthEigenvector-3)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",par[first_index+row],0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+row,-limit,limit);
        }

        limit = 0.2;
        if (type!=1) limit = 0.;
        first_index = (4*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",par[first_index+col],0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+col,-limit,limit);
        }
        first_index = (4*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",par[first_index+col],0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+col,-limit,limit);
        }

        // eigenvalues
        double input_value = 0;
        first_index = (4*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            input_value = par[first_index+NthEigenvalue-1];
            //limit = 0.02*eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real();
            limit = 0.;
            Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.001*limit);
            Chi2Minimizer->SetVariableLimits(first_index+NthEigenvalue-1,input_value-limit,input_value+limit);
        }
    }

}
double FourierChi2Function1D_0(const double *par)
{
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        vtr_fit(row) = 0.;
    }
    int first_index = 0;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
        double x_norm = par[first_index]*2.0*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
        vtr_fit(row) += par[first_index+1]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            vtr_fit(row) += par[first_index+1+2*mode-1]*cos(double(mode)*x_norm);
            vtr_fit(row) += par[first_index+1+2*mode]*sin(double(mode)*x_norm);
        }
    }
    double chi2 = 0.;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        chi2 += pow(vtr_fit(row).real()-vtr_data(row).real(),2);
    }
    return chi2;
}
double FourierChi2Function1D_1(const double *par)
{
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        vtr_fit(row) = 0.;
    }
    int first_index = 0;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        double x = double(row+1)*(MSCW_plot_upper-MSCW_plot_lower)/double(N_bins_for_deconv)+MSCW_plot_lower;
        double x_norm = par[first_index]*1.0*M_PI*(x-MSCW_plot_lower)/(MSCW_plot_upper-MSCW_plot_lower);
        vtr_fit(row) += par[first_index+1]/2.;
        for (int mode=1;mode<=n_fourier_modes;mode++)
        {
            vtr_fit(row) += par[first_index+1+2*mode-1]*cos(double(mode)*x_norm);
            vtr_fit(row) += par[first_index+1+2*mode]*sin(double(mode)*x_norm);
        }
    }
    double chi2 = 0.;
    for (int row=0;row<N_bins_for_deconv;row++)
    {
        chi2 += pow(vtr_fit(row).real()-vtr_data(row).real(),2);
    }
    return chi2;
}
void FourierSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer)
{

    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);

    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    double scale = 0.;
    double limit = 0.; 

    ROOT::Math::GSLMinimizer Chi2Minimizer1D( ROOT::Math::kVectorBFGS );
    Chi2Minimizer1D.SetMaxFunctionCalls(200); // for Minuit/Minuit2
    Chi2Minimizer1D.SetMaxIterations(200); // for GSL
    Chi2Minimizer1D.SetTolerance(0.001);
    const double *par_1D;

    
    //for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //{

    //    std::cout << NthEigenvector << "-th eigenvector" << std::endl;
    //    vtr_data = eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector);
    //    ROOT::Math::Functor Chi2Func_0(&FourierChi2Function1D_0,(2+2*n_fourier_modes));
    //    Chi2Minimizer1D.SetFunction(Chi2Func_0);
    //    first_index = 0;
    //    Chi2Minimizer1D.SetVariable(first_index,"par_1D["+std::to_string(int(first_index))+"]",1.,0.001);
    //    Chi2Minimizer1D.SetVariable(first_index+1,"par_1D["+std::to_string(int(first_index+1))+"]",1.,0.001);
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        Chi2Minimizer1D.SetVariable(first_index+1+2*mode-1,"par_1D["+std::to_string(int(first_index+1+2*mode-1))+"]",1.,0.001);
    //        Chi2Minimizer1D.SetVariable(first_index+1+2*mode,"par_1D["+std::to_string(int(first_index+1+2*mode))+"]",1.,0.001);
    //    }
    //    Chi2Minimizer1D.Minimize();
    //    par_1D = Chi2Minimizer1D.X();
    //    std::cout << "coeff " << 0 << " = " << par_1D[0] << std::endl;
    //    std::cout << "coeff " << 1 << " = " << par_1D[1] << std::endl;

    //    first_index = (2*NthEigenvector-2)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",par_1D[0],0.001);
    //    Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",par_1D[1],0.001);
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        std::cout << "coeff " << 1+2*mode-1 << " = " << par_1D[1+2*mode-1] << std::endl;
    //        std::cout << "coeff " << 1+2*mode << " = " << par_1D[1+2*mode] << std::endl;
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode-1,"par["+std::to_string(int(first_index+1+2*mode-1))+"]",par_1D[1+2*mode-1],0.001);
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode,"par["+std::to_string(int(first_index+1+2*mode))+"]",par_1D[1+2*mode],0.001);
    //    }

    //    std::cout << NthEigenvector << "-th inverse eigenvector" << std::endl;
    //    vtr_data = eigensolver_data.eigenvectors().inverse().row(mtx_dark.cols()-NthEigenvector);
    //    ROOT::Math::Functor Chi2Func_1(&FourierChi2Function1D_1,(2+2*n_fourier_modes));
    //    Chi2Minimizer1D.SetFunction(Chi2Func_1);
    //    first_index = 0;
    //    Chi2Minimizer1D.SetVariable(first_index,"par_1D["+std::to_string(int(first_index))+"]",1.,0.001);
    //    Chi2Minimizer1D.SetVariable(first_index+1,"par_1D["+std::to_string(int(first_index+1))+"]",1.,0.001);
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        Chi2Minimizer1D.SetVariable(first_index+1+2*mode-1,"par_1D["+std::to_string(int(first_index+1+2*mode-1))+"]",1.,0.001);
    //        Chi2Minimizer1D.SetVariable(first_index+1+2*mode,"par_1D["+std::to_string(int(first_index+1+2*mode))+"]",1.,0.001);
    //    }
    //    Chi2Minimizer1D.Minimize();
    //    par_1D = Chi2Minimizer1D.X();
    //    std::cout << "coeff " << 0 << " = " << par_1D[0] << std::endl;
    //    std::cout << "coeff " << 1 << " = " << par_1D[1] << std::endl;

    //    first_index = (2*NthEigenvector-1)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    Chi2Minimizer->SetVariable(first_index,"par["+std::to_string(int(first_index))+"]",par_1D[0],0.001);
    //    Chi2Minimizer->SetVariable(first_index+1,"par["+std::to_string(int(first_index+1))+"]",par_1D[1],0.001);
    //    for (int mode=1;mode<=n_fourier_modes;mode++)
    //    {
    //        std::cout << "coeff " << 1+2*mode-1 << " = " << par_1D[1+2*mode-1] << std::endl;
    //        std::cout << "coeff " << 1+2*mode << " = " << par_1D[1+2*mode] << std::endl;
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode-1,"par["+std::to_string(int(first_index+1+2*mode-1))+"]",par_1D[1+2*mode-1],0.001);
    //        Chi2Minimizer->SetVariable(first_index+1+2*mode,"par["+std::to_string(int(first_index+1+2*mode))+"]",par_1D[1+2*mode],0.001);
    //    }

    //    // eigenvalues
    //    double input_value = 0.;
    //    first_index = (2*NthEigenvector-0)*(2+2*n_fourier_modes)+NumberOfEigenvectors*(NthEigenvector-1);
    //    for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
    //    {
    //        input_value = eigensolver_data.eigenvalues()(N_bins_for_deconv-NthEigenvector).real();
    //        if (NthEigenvalue!=NthEigenvector) input_value = 0.;
    //        Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.001);
    //    }
    //}

    Chi2Minimizer->SetVariable(0, "par["+std::to_string(int(0))+"]", 1, 0.001);
    Chi2Minimizer->SetVariableLimits(0,0.9,1.1);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        limit = 0.2;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",0.,0.005);
            //if (NthEigenvector>1) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
        }

        limit = 0.2;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",0.,0.005);
            //if (NthEigenvector>1) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
        }

        // eigenvalues
        double input_value = 0.;
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            input_value = 0.;
            if (NthEigenvalue==NthEigenvector && NthEigenvalue==1)
            {
                //limit = 0.2*eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvector).real();
                limit = 0.;
            }
            else
            {
                //limit = 0.1*eigensolver_dark.eigenvalues()(N_bins_for_deconv-1).real();
                //if (NthEigenvalue>NthEigenvector) limit = 0.;
                //limit = 0.;
                limit = 0.05*eigensolver_dark.eigenvalues()(N_bins_for_deconv-1).real();
            }
            Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+NthEigenvalue-1,input_value-limit,input_value+limit);
        }
    }

}
void NormalizeHist2D(TH2D* hist_ref, TH2D* hist_sub)
{
    int binx_lower = hist_ref->GetXaxis()->FindBin(MSCW_cut_lower);
    int binx_upper = hist_ref->GetXaxis()->FindBin(MSCW_cut_upper)-1;
    int biny_lower = hist_ref->GetYaxis()->FindBin(MSCW_cut_upper);
    int biny_upper = hist_ref->GetYaxis()->FindBin(3.*MSCW_cut_upper)-1;
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
void NetflixMethodPrediction(string target_data, double PercentCrab, double tel_elev_lower_input, double tel_elev_upper_input, bool isON, double MSCW_cut_input, double MSCL_cut_input, int rank)
{

    signalfree_model = false;
    //if (target_data=="Proton") signalfree_model = true;

    TH1::SetDefaultSumw2();
    sprintf(target, "%s", target_data.c_str());
    MSCW_cut_blind = MSCW_cut_input;
    MSCL_cut_blind = MSCL_cut_input;
    TelElev_lower = tel_elev_lower_input;
    TelElev_upper = tel_elev_upper_input;
    Theta2_cut_lower = 0.2;
    Theta2_cut_upper = 10.;
    if (isON) Theta2_cut_lower = 0.;
    TString file_tag;
    if (isON) file_tag = "ON";
    else file_tag = "OFF";

    vector<TH2D> Hist_Bkgd_MSCLW;
    vector<TH2D> Hist_GammaRDBM_MSCLW;
    vector<TH2D> Hist_Redu_MSCLW;
    vector<TH1D> Hist_Data_Eigenvalues_real;
    vector<TH1D> Hist_Dark_Eigenvalues_real;
    vector<TH1D> Hist_Data_Eigenvalues_imag;
    vector<TH1D> Hist_Dark_Eigenvalues_imag;
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
    TFile InputDataFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+file_tag+".root");
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "energy = " << energy_bins[e] << std::endl;

        N_bins_for_deconv = N_bins_for_deconv_at_E[e];

        CurrentEnergy = energy_bins[e];
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        TString filename_gamma  = "Hist_GammaData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        if (PercentCrab>0.) filename_gamma  = "Hist_GammaMC_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_data  = "Hist_Data_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TString filename_dark  = "Hist_Dark_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_GammaMC = (TH2D*)InputDataFile.Get(filename_gamma);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
        Hist_Bkgd_MSCLW.push_back(TH2D("Hist_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaRDBM_MSCLW.push_back(TH2D("Hist_GammaRDBM_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Redu_MSCLW.push_back(TH2D("Hist_Redu_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Data_Eigenvalues_real.push_back(TH1D("Hist_Data_Eigenvalues_real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Eigenvalues_real.push_back(TH1D("Hist_Dark_Eigenvalues_real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_Eigenvalues_imag.push_back(TH1D("Hist_Data_Eigenvalues_imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Eigenvalues_imag.push_back(TH1D("Hist_Dark_Eigenvalues_imag_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
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
        //if (energy_bins[e]<398) continue;
        //if (energy_bins[e]>=473) continue;
        //if (energy_bins[e]<473) continue;
        //if (energy_bins[e]>=562) continue;
        //if (energy_bins[e]<562) continue;
        //if (energy_bins[e]>=794) continue;
        //if (energy_bins[e]<1585) continue;
        //if (energy_bins[e]>=2239) continue;
        //if (energy_bins[e]<2239) continue;
        //if (CurrentEnergy<794.) continue;

        MatrixXcd mtx_data_blind(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_dark_blind(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_data_bkgd(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        mtx_gamma = fillMatrix(Hist_GammaMC);
        mtx_data = fillMatrix(Hist_Data);
        mtx_dark = fillMatrix(Hist_Dark);
        mtx_data_err = fillMatrixError(Hist_Data);
        mtx_dark_err = fillMatrixError(Hist_Dark);
        int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCL_cut_lower);
        int binx_blind = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int binx_upper = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind*2)-1;
        int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_lower);
        int biny_blind = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_upper = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind*2)-1;
        mtx_data_blind = MatrixBlind(mtx_data, binx_lower, biny_lower, binx_blind, biny_blind);
        mtx_dark_blind = MatrixBlind(mtx_dark, binx_lower, biny_lower, binx_blind, biny_blind);

        n_fourier_modes = 6;
        n_taylor_modes = 6;

        NumberOfEigenvectors = rank;
        //if (CurrentEnergy>=562.) NumberOfEigenvectors = 2;

        eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);

        MatrixXcd mtx_eigenval_data_redu(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        MatrixXcd mtx_eigenval_dark_redu(Hist_Data->GetNbinsX(),Hist_Data->GetNbinsY());
        for (int i=0;i<mtx_dark.cols();i++)
        {
            for (int j=0;j<mtx_dark.rows();j++)
            {
                if (i==j && i>=mtx_dark.cols()-NumberOfEigenvectors) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                    mtx_eigenval_dark_redu(i,j) = eigensolver_dark.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                    mtx_eigenval_dark_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        mtx_dark_redu = eigensolver_dark.eigenvectors()*mtx_eigenval_dark_redu*eigensolver_dark.eigenvectors().inverse();

        for (int NthEigenvalue=1;NthEigenvalue<=N_bins_for_deconv;NthEigenvalue++)
        {
            Hist_Data_Eigenvalues_real.at(e).SetBinContent(NthEigenvalue,eigensolver_data.eigenvalues()(N_bins_for_deconv-NthEigenvalue).real());
            Hist_Dark_Eigenvalues_real.at(e).SetBinContent(NthEigenvalue,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvalue).real());
            Hist_Data_Eigenvalues_imag.at(e).SetBinContent(NthEigenvalue,eigensolver_data.eigenvalues()(N_bins_for_deconv-NthEigenvalue).imag());
            Hist_Dark_Eigenvalues_imag.at(e).SetBinContent(NthEigenvalue,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvalue).imag());
        }

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

        SetInitialEigenvectors();
        //init_deriv_at_zero = SmoothEigenvectors(&mtx_eigenvector_init, &mtx_eigenvector_inv_init);
        //std::cout << "initial deriv_at_zero = " << init_deriv_at_zero << std::endl;

        n_fourier_modes = 6;
        ROOT::Math::Functor Chi2Func(&FourierChi2Function,1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors*NumberOfEigenvectors); 
        std::cout << "total n paramters = " << 1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors << std::endl;

        ROOT::Math::GSLMinimizer Chi2Minimizer_0th( ROOT::Math::kSteepestDescent );
        Chi2Minimizer_0th.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
        //Chi2Minimizer_0th.SetMaxIterations(10000); // for GSL
        Chi2Minimizer_0th.SetMaxIterations(100); // for GSL
        Chi2Minimizer_0th.SetTolerance(0.001);
        Chi2Minimizer_0th.SetPrintLevel(2);
        Chi2Minimizer_0th.SetFunction(Chi2Func);
        FourierSetInitialVariables(&Chi2Minimizer_0th);
        const double *par_0th = Chi2Minimizer_0th.X();
        NthIteration = 0;
        double init_chi2 = FourierChi2Function(par_0th);
        std::cout << "initial chi2 = " << init_chi2 << std::endl;
        Chi2Minimizer_0th.SetTolerance(0.01*double(N_bins_for_deconv*N_bins_for_deconv));
        Chi2Minimizer_0th.Minimize();
        par_0th = Chi2Minimizer_0th.X();
        std::cout << "final chi2 = " << FourierChi2Function(par_0th) << std::endl;
        FourierParametrizeEigenvectors(par_0th);
        //double final_deriv_at_zero = SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);
        //std::cout << "final deriv_at_zero = " << final_deriv_at_zero << std::endl;

        //ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kVectorBFGS );
        //Chi2Minimizer_1st.SetMaxFunctionCalls(200); // for Minuit/Minuit2
        //Chi2Minimizer_1st.SetMaxIterations(200); // for GSL
        //Chi2Minimizer_1st.SetTolerance(0.001);
        //Chi2Minimizer_1st.SetFunction(Chi2Func);
        //FourierSet2ndIterationVariables(&Chi2Minimizer_1st,par_0th,1);
        //Chi2Minimizer_1st.Minimize();
        //const double *par_1st = Chi2Minimizer_1st.X();
        //FourierParametrizeEigenvectors(par_1st);
        ////SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);

        //ROOT::Math::GSLMinimizer Chi2Minimizer_2nd( ROOT::Math::kVectorBFGS );
        //Chi2Minimizer_2nd.SetMaxFunctionCalls(200); // for Minuit/Minuit2
        //Chi2Minimizer_2nd.SetMaxIterations(200); // for GSL
        //Chi2Minimizer_2nd.SetTolerance(0.001);
        //Chi2Minimizer_2nd.SetFunction(Chi2Func);
        //FourierSet2ndIterationVariables(&Chi2Minimizer_2nd,par_1st,2);
        //Chi2Minimizer_2nd.Minimize();
        //const double *par_2nd = Chi2Minimizer_2nd.X();
        //FourierParametrizeEigenvectors(par_2nd);
        ////SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);

        //ROOT::Math::GSLMinimizer Chi2Minimizer_3rd( ROOT::Math::kVectorBFGS );
        //Chi2Minimizer_3rd.SetMaxFunctionCalls(200); // for Minuit/Minuit2
        //Chi2Minimizer_3rd.SetMaxIterations(200); // for GSL
        //Chi2Minimizer_3rd.SetTolerance(0.001);
        //Chi2Minimizer_3rd.SetFunction(Chi2Func);
        //FourierSet2ndIterationVariables(&Chi2Minimizer_3rd,par_2nd,2);
        //Chi2Minimizer_3rd.Minimize();
        //const double *par_3rd = Chi2Minimizer_3rd.X();
        //FourierParametrizeEigenvectors(par_3rd);
        ////SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);

        mtx_data_bkgd = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;


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
        fill2DHistogramAbs(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        fill2DHistogramAbs(&Hist_GammaRDBM_MSCLW.at(e),mtx_gamma);

        double gamma_total = Hist_Data->Integral(binx_lower,binx_blind,biny_lower,biny_blind)-Hist_Bkgd_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind);
        gamma_total = max(0.,gamma_total);
        double scale_gamma = double(gamma_total)/double(Hist_GammaRDBM_MSCLW.at(e).Integral(binx_lower,binx_blind,biny_lower,biny_blind));
        Hist_GammaRDBM_MSCLW.at(e).Scale(scale_gamma);
        if (signalfree_model) Hist_GammaRDBM_MSCLW.at(e).Scale(0);

        fill2DHistogramAbs(&Hist_Redu_MSCLW.at(e),mtx_data_redu);
        //RestoreEdge(Hist_Dark,&Hist_Bkgd_MSCLW.at(e),-0.5,-0.5);
        //NormalizeHist2D(Hist_Data,&Hist_Bkgd_MSCLW.at(e));
        std::cout << "Hist_Redu_MSCLW.at(e).Integral() = " << Hist_Redu_MSCLW.at(e).Integral(1,20,1,20) << std::endl;
        std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral(1,20,1,20) << std::endl;
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
        std::cout << "mtx_eigenvalue:" << std::endl;
        std::cout << mtx_eigenvalue.block(mtx_data.rows()-NumberOfEigenvectors,mtx_data.cols()-NumberOfEigenvectors,NumberOfEigenvectors,NumberOfEigenvectors) << std::endl;

    }
    InputDataFile.Close();


    TFile OutputFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+file_tag+".root","update");
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_GammaRDBM_MSCLW.at(e).Write();
        Hist_Bkgd_MSCLW.at(e).Write();
        Hist_Redu_MSCLW.at(e).Write();
        Hist_Data_Eigenvalues_real.at(e).Write();
        Hist_Dark_Eigenvalues_real.at(e).Write();
        Hist_Data_Eigenvalues_imag.at(e).Write();
        Hist_Dark_Eigenvalues_imag.at(e).Write();
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
