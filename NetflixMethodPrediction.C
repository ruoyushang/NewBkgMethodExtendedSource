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
MatrixXcd mtx_data_bkgd(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_dark(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_init(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv_vari(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvector_inv(N_bins_for_deconv,N_bins_for_deconv);
MatrixXcd mtx_eigenvalue(N_bins_for_deconv,N_bins_for_deconv);
VectorXcd vtr_fit(N_bins_for_deconv);
VectorXcd vtr_data(N_bins_for_deconv);
VectorXcd vtr_bkgd_eigenvalues(N_bins_for_deconv);
VectorXcd vtr_data_eigenvalues(N_bins_for_deconv);
VectorXcd vtr_dark_eigenvalues(N_bins_for_deconv);
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd;
ComplexEigenSolver<MatrixXcd> eigensolver_bkgd_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_dark;
ComplexEigenSolver<MatrixXcd> eigensolver_data;
ComplexEigenSolver<MatrixXcd> eigensolver_dark_transpose;
ComplexEigenSolver<MatrixXcd> eigensolver_data_transpose;
double init_deriv_at_zero;
double ratio_empty_bins;
int binx_blind_global;
int biny_blind_global;

void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            if (!invert_y) hist->SetBinContent(binx+1,biny+1,mtx(binx,biny).real());
            else hist->SetBinContent(binx+1,hist->GetNbinsY()-biny,mtx(binx,biny).real());
        }
    }
}
std::complex<double> MatrixIntegral(MatrixXcd mtx, int binx_0, int binx_1, int biny_0, int biny_1)
{
    std::complex<double> integral (0.0, 0.0);
    for (int binx=binx_0;binx<=binx_1;binx++)
    {
        for (int biny=biny_0;biny<=biny_1;biny++)
        {
            integral += mtx(binx,biny);
        }
    }
    return integral;
}
void fill2DHistogramAbs(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            if (!invert_y) hist->SetBinContent(binx+1,biny+1,max(0.,mtx(binx,biny).real()));
            else hist->SetBinContent(binx+1,hist->GetNbinsY()-biny,max(0.,mtx(binx,biny).real()));
        }
    }
}
MatrixXcd GetSubEigenvectors(MatrixXcd mtx_input, int region)
{
    int start_row;
    int start_col;
    int size_row;
    int size_col;
    if (region==0 || region==-1)
    {                
        start_row = 0;
        start_col = 0;
        size_row = (binx_blind_global)-start_row;
        size_col = mtx_input.cols();
    }                
    if (region==1 || region==-1)
    {                
        start_row = binx_blind_global;
        start_col = 0;
        size_row = mtx_input.rows()-(binx_blind_global);
        size_col = mtx_input.cols();
    }                
    MatrixXcd mtx_output = mtx_input.block(start_row,start_col,size_row,size_col);
    return mtx_output;
}
MatrixXcd GetSubmatrix(MatrixXcd mtx_input, int region)
{
    // regions:
    //   0 | 1
    //   -----
    //   2 | 3
    int start_row;
    int start_col;
    int size_row;
    int size_col;
    if (region==0)
    {
        start_row = 0;
        start_col = 0;
        size_row = (binx_blind_global)-start_row;
        size_col = (biny_blind_global)-start_col;
    }
    if (region==3)
    {
        start_row = binx_blind_global;
        start_col = biny_blind_global;
        size_row = mtx_input.rows()-(binx_blind_global);
        size_col = mtx_input.cols()-(biny_blind_global);
    }
    if (region==1)
    {
        start_row = 0;
        start_col = biny_blind_global;
        size_row = (binx_blind_global)-start_row;
        size_col = mtx_input.cols()-(biny_blind_global);
    }
    if (region==2)
    {
        start_row = binx_blind_global;
        start_col = 0;
        size_row = mtx_input.rows()-(binx_blind_global);
        size_col = (biny_blind_global)-start_col;
    }
    MatrixXcd mtx_output = mtx_input.block(start_row,start_col,size_row,size_col);
    return mtx_output;
}
MatrixXcd GetAlphaMatrix(MatrixXcd mtx_input, int region)
{
    MatrixXcd mtx_sub = GetSubmatrix(mtx_input, region);
    MatrixXcd mtx_alpha = mtx_sub.transpose()*mtx_sub;
    if (region==1) mtx_alpha = mtx_sub*mtx_sub.transpose();
    return mtx_alpha;
}
MatrixXcd GetUnitMatrix(int rows, int cols)
{
    MatrixXcd mtx_unit(rows,cols);
    for (int row=0;row<rows;row++)
    {
        for (int col=0;col<cols;col++)
        {
            mtx_unit(row,col) = 0.;
            if (col==row) mtx_unit(row,col) = 1.;
        }
    }
    return mtx_unit;
}
MatrixXcd GetProjectionMatrix(MatrixXcd mtx_input, int region)
{
    MatrixXcd mtx_sub = GetSubmatrix(mtx_input, region);
    MatrixXcd mtx_alpha = GetAlphaMatrix(mtx_input,region);

    MatrixXcd mtx_unit = GetUnitMatrix(mtx_sub.rows(),mtx_sub.cols());

    MatrixXcd mtx_output(mtx_sub.rows(),mtx_sub.cols());
    mtx_output = mtx_unit-mtx_sub*mtx_alpha.inverse()*mtx_sub.transpose();
    if (region==1) mtx_output = mtx_unit-mtx_sub.transpose()*mtx_alpha.inverse()*mtx_sub;
    return mtx_output;
}
MatrixXcd GetBetaMatrix(MatrixXcd mtx_input, int idx_1, int idx_2)
{
    MatrixXcd mtx_sub_1 = GetSubmatrix(mtx_input, idx_1);
    MatrixXcd mtx_alpha_1 = GetAlphaMatrix(mtx_input,idx_1);
    MatrixXcd mtx_sub_2 = GetSubmatrix(mtx_input, idx_2);
    MatrixXcd mtx_alpha_2 = GetAlphaMatrix(mtx_input,idx_2);
    MatrixXcd mtx_beta = mtx_sub_2*mtx_alpha_2.inverse().transpose()*mtx_alpha_1.inverse()*mtx_sub_1;
    return mtx_beta;
}
MatrixXcd GetQMatrix(MatrixXcd mtx_input, int idx_1, int idx_2)
{
    MatrixXcd mtx_sub_3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_beta = GetBetaMatrix(mtx_input, 1, 2);
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_unit = GetUnitMatrix(mtx_sub_3.rows(),mtx_sub_3.cols());
    MatrixXcd mtx_Q = (mtx_sub_3-eigensolver_input.eigenvalues()(idx_1)*mtx_unit)*mtx_beta.transpose()*(mtx_sub_3-eigensolver_input.eigenvalues()(idx_2)*mtx_unit);
    return mtx_Q;
}
MatrixXcd GetUnitEigenvalueMatrix(MatrixXcd mtx_input, int entry, int dim)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_lambda(dim,dim);
    for (int row=0;row<dim;row++)
    {
        for (int col=0;col<dim;col++)
        {
            mtx_lambda(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda(row,col) = eigensolver_input.eigenvalues()(entry);
            }
        }
    }
    return mtx_lambda;
}
MatrixXcd GetEigenvalueMatrix(MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    MatrixXcd mtx_lambda(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambda(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda(row,col) = eigensolver_input.eigenvalues()(col);
            }
        }
    }
    return mtx_lambda;
}
MatrixXcd MakeRealEigenvectors(MatrixXcd mtx_input)
{
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            double real = mtx_input(row,col).real();
            double imag = mtx_input(row,col).imag();
            mtx_output(row,col) = pow(real*real+imag*imag,0.5);
        }
    }
    //return mtx_output;
    return mtx_input;
}
MatrixXcd GetLambdaNuMatrix(MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r = eigensolver_input.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_transpose.eigenvectors();
    mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    MatrixXcd mtx_lambdanu(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambdanu(row,col) = 0.;
            if (row==col)
            {
                mtx_lambdanu(row,col) = eigensolver_input.eigenvalues()(col)/mtx_H(row,col);
            }
        }
    }
    return mtx_lambdanu;
}
MatrixXcd GetEigenvalueWeightedVectors(MatrixXcd mtx_eigenvalues, MatrixXcd mtx_input)
{
    MatrixXcd mtx_output(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_output(row,col) = mtx_input(row,col)*mtx_eigenvalues(col,col);
        }
    }
    return mtx_output;
}

MatrixXcd GetWeightedMatrix(MatrixXcd mtx_input, MatrixXcd mtx_weight, int fix_row)
{
    MatrixXcd mtx_output(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_output(row,col) = mtx_input(row,col)*mtx_weight(fix_row,col);
        }
    }
    return mtx_output;
}

void PrintItReversely(MatrixXcd mtx_input)
{
    MatrixXcd mtx_reverse(mtx_input.rows(),mtx_input.cols());
    IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
    //for (int row=0;row<mtx_input.rows();row++)
    //{
    //    for (int col=0;col<mtx_input.cols();col++)
    //    {
    //        mtx_reverse(row,col) = mtx_input(mtx_input.rows()-1-row,mtx_input.cols()-1-col);
    //    }
    //}
    //std::cout << mtx_reverse.format(CleanFmt) << std::endl;
    std::cout << mtx_input.format(CleanFmt) << std::endl;
}
MatrixXcd Convert_hadron_to_gamma_vector_v2(MatrixXcd mtx_input, MatrixXcd mtx_vc, MatrixXcd mtx_eigenvalue, int LeftOrRight)
{
    // LeftOrRight = 1, left
    // LeftOrRight = 2, right
    MatrixXcd mtx_sub = GetSubmatrix(mtx_input, LeftOrRight);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_alpha = GetAlphaMatrix(mtx_input, LeftOrRight);
    MatrixXcd mtx_lambda_vc = GetEigenvalueWeightedVectors(mtx_eigenvalue,mtx_vc);
    if (LeftOrRight==2)
    {
        MatrixXcd mtx_M3_vc = mtx_M3*mtx_vc;
        MatrixXcd mtx_alpha2inv_M2t = mtx_alpha.inverse()*mtx_sub.transpose();
        return -mtx_alpha2inv_M2t*(mtx_M3_vc-mtx_lambda_vc);
    }
    else
    {
        MatrixXcd mtx_M3t_vc = mtx_M3.transpose()*mtx_vc;
        MatrixXcd mtx_alpha1inv_M1 = mtx_alpha.inverse()*mtx_sub;
        return -mtx_alpha1inv_M1*(mtx_M3t_vc-mtx_lambda_vc);
    }

}
MatrixXcd Convert_hadron_to_gamma_vector(MatrixXcd mtx_eigenvector, MatrixXcd mtx_eigenvalue, int LeftOrRight)
{
    // LeftOrRight = 1, left
    // LeftOrRight = 2, right
    MatrixXcd mtx_vg = GetSubEigenvectors(mtx_eigenvector,0);
    MatrixXcd mtx_vc = GetSubEigenvectors(mtx_eigenvector,1);
    MatrixXcd mtx_sub = GetSubmatrix(mtx_data, LeftOrRight);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_data, 3);
    MatrixXcd mtx_alpha = GetAlphaMatrix(mtx_data, LeftOrRight);
    MatrixXcd mtx_lambda_vc = GetEigenvalueWeightedVectors(mtx_eigenvalue,mtx_vc);
    if (LeftOrRight==2)
    {
        MatrixXcd mtx_M3_vc = mtx_M3*mtx_vc;
        MatrixXcd mtx_alpha2inv_M2t = mtx_alpha.inverse()*mtx_sub.transpose();
        //std::cout << "mtx_vg = " << std::endl;
        //PrintItReversely(mtx_vg);
        //std::cout << "mtx_alpha2inv_M2t*(mtx_M3_vc-mtx_lambda_vc) = " << std::endl;
        //PrintItReversely(mtx_alpha2inv_M2t*(mtx_M3_vc-mtx_lambda_vc));
        return -mtx_alpha2inv_M2t*(mtx_M3_vc-mtx_lambda_vc);
    }
    else
    {
        MatrixXcd mtx_M3t_vc = mtx_M3.transpose()*mtx_vc;
        MatrixXcd mtx_alpha1inv_M1 = mtx_alpha.inverse()*mtx_sub;
        //std::cout << "mtx_vg = " << std::endl;
        //PrintItReversely(mtx_vg);
        //std::cout << "mtx_alpha1inv_M1*(mtx_M3t_vc-mtx_lambda_vc) = " << std::endl;
        //PrintItReversely(mtx_alpha1inv_M1*(mtx_M3t_vc-mtx_lambda_vc));
        return -mtx_alpha1inv_M1*(mtx_M3t_vc-mtx_lambda_vc);
    }

}
void VerifyBetaMatrix(MatrixXcd mtx_input)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r = eigensolver_input.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_transpose.eigenvectors();
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_U_rg = GetSubEigenvectors(mtx_U_r,0);
    MatrixXcd mtx_U_lg = GetSubEigenvectors(mtx_U_l,0);
    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_input);
    MatrixXcd mtx_lambdanu_transpose = GetLambdaNuMatrix(mtx_input.transpose());
    std::cout << "mtx_lambdanu = " << std::endl;
    PrintItReversely(mtx_lambdanu);
    std::cout << "mtx_lambdanu_transpose = " << std::endl;
    PrintItReversely(mtx_lambdanu_transpose);
    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_alpha2 = GetAlphaMatrix(mtx_input, 2);
    MatrixXcd mtx_beta_11 = mtx_M1.transpose()*mtx_alpha1.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_12 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_22 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();

    std::cout << "mtx_R = " << std::endl;
    PrintItReversely(mtx_R);
    MatrixXcd mtx_lambdanu_lc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc);
    MatrixXcd mtx_lambdanu_lc_transpose = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc).transpose();
    MatrixXcd mtx_R_predict_term1 = mtx_U_rc.transpose().conjugate()*mtx_U_rc;
    MatrixXcd mtx_R_predict_term2 = (mtx_U_rc.transpose().conjugate()*mtx_M3-mtx_R*mtx_lambdanu_lc_transpose)*mtx_beta_11*(mtx_M3.transpose()*mtx_U_rc-mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate());
    MatrixXcd mtx_R_predict = mtx_R_predict_term1 + mtx_R_predict_term2;
    std::cout << "mtx_R_predict = " << std::endl;
    PrintItReversely(mtx_R_predict);

    std::cout << "mtx_L = " << std::endl;
    PrintItReversely(mtx_L);
    MatrixXcd mtx_lambdanu_rc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc);
    MatrixXcd mtx_lambdanu_rc_transpose = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc).transpose();
    MatrixXcd mtx_L_predict_term1 = mtx_U_lc.transpose().conjugate()*mtx_U_lc;
    MatrixXcd mtx_L_predict_term2 = (mtx_U_lc.transpose().conjugate()*mtx_M3.transpose()-mtx_L*mtx_lambdanu_rc_transpose)*mtx_beta_22*(mtx_M3*mtx_U_lc-mtx_lambdanu_rc.conjugate()*mtx_L.transpose().conjugate());
    MatrixXcd mtx_L_predict = mtx_L_predict_term1 + mtx_L_predict_term2;
    std::cout << "mtx_L_predict = " << std::endl;
    PrintItReversely(mtx_L_predict);

    std::cout << "mtx_H = " << std::endl;
    PrintItReversely(mtx_H);
    MatrixXcd mtx_H_predict_term1 = mtx_U_lc.transpose()*mtx_U_rc;
    MatrixXcd mtx_H_predict_term2 = (mtx_U_lc.transpose()*mtx_M3.transpose()-mtx_L.conjugate()*mtx_lambdanu_rc_transpose.conjugate())*mtx_beta_12*(mtx_M3.transpose()*mtx_U_rc-mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate());
    MatrixXcd mtx_H_predict = mtx_H_predict_term1 + mtx_H_predict_term2;
    std::cout << "mtx_H_predict = " << std::endl;
    PrintItReversely(mtx_H_predict);
}
void VerifyHMatrix(MatrixXcd mtx_input)
{
    MatrixXcd mtx_lambda = GetEigenvalueMatrix(mtx_input);
    MatrixXcd mtx_lambda_transpose = GetEigenvalueMatrix(mtx_input.transpose());
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r = eigensolver_input.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_transpose.eigenvectors();
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    std::cout << "mtx_H = " << std::endl;
    PrintItReversely(mtx_H);
    MatrixXcd mtx_rg_predict = Convert_hadron_to_gamma_vector(mtx_U_r, mtx_lambda, 2);
    MatrixXcd mtx_lg_predict = Convert_hadron_to_gamma_vector(mtx_U_l, mtx_lambda_transpose, 1);
    MatrixXcd mtx_H_from_Q = mtx_lg_predict.transpose()*mtx_rg_predict + mtx_U_lc.transpose()*mtx_U_rc;
    std::cout << "mtx_H_from_Q = " << std::endl;
    PrintItReversely(mtx_H_from_Q);
    MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
    std::cout << "mtx_R = " << std::endl;
    PrintItReversely(mtx_R);
    MatrixXcd mtx_R_from_Q = mtx_rg_predict.transpose().conjugate()*mtx_rg_predict + mtx_U_rc.transpose().conjugate()*mtx_U_rc;
    std::cout << "mtx_R_from_Q = " << std::endl;
    PrintItReversely(mtx_R_from_Q);
    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    std::cout << "mtx_L = " << std::endl;
    PrintItReversely(mtx_L);
    MatrixXcd mtx_L_from_Q = mtx_lg_predict.transpose().conjugate()*mtx_lg_predict + mtx_U_lc.transpose().conjugate()*mtx_U_lc;
    std::cout << "mtx_L_from_Q = " << std::endl;
    PrintItReversely(mtx_L_from_Q);
}
void VerifySmallgMatrix(MatrixXcd mtx_input, int NthEigenvector)
{
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r = eigensolver_input.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_transpose.eigenvectors();
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_U_rg = GetSubEigenvectors(mtx_U_r,0);
    MatrixXcd mtx_U_lg = GetSubEigenvectors(mtx_U_l,0);
    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_lambda = GetUnitEigenvalueMatrix(mtx_input,mtx_input.cols()-NthEigenvector,mtx_M3.rows());
    MatrixXcd mtx_lambda_transpose = GetUnitEigenvalueMatrix(mtx_input.transpose(),mtx_input.cols()-NthEigenvector,mtx_M3.rows());
    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_input);
    MatrixXcd mtx_lambdanu_transpose = GetLambdaNuMatrix(mtx_input.transpose());

    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    MatrixXcd mtx_lM = mtx_U_l.transpose().conjugate()*mtx_input.transpose();
    MatrixXcd mtx_lambdanu_L_r = mtx_L*mtx_lambdanu*mtx_U_r.transpose();
    std::cout << "mtx_lM = " << std::endl;
    PrintItReversely(mtx_lM);
    std::cout << "mtx_lambdanu_L_r = " << std::endl;
    PrintItReversely(mtx_lambdanu_L_r);

    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_Smallg1 = mtx_M2*mtx_alpha1.inverse()*mtx_M1;

    std::cout << "mtx_U_lg = " << std::endl;
    PrintItReversely(mtx_U_lg);
    MatrixXcd mtx_lg_predict = Convert_hadron_to_gamma_vector(mtx_U_l, GetEigenvalueMatrix(mtx_input), 1);
    std::cout << "mtx_lg_predict = " << std::endl;
    PrintItReversely(mtx_lg_predict);

    std::cout << "mtx_lg_predict.transpose().conjugate()*mtx_M2.transpose() + mtx_U_lc.transpose().conjugate()*mtx_M3.transpose() = " << std::endl;
    PrintItReversely(mtx_lg_predict.transpose().conjugate()*mtx_M2.transpose() + mtx_U_lc.transpose().conjugate()*mtx_M3.transpose());
    std::cout << "mtx_lg_predict.transpose().conjugate()*mtx_M2.transpose() = " << std::endl;
    PrintItReversely(mtx_lg_predict.transpose().conjugate()*mtx_M2.transpose());
    std::cout << "mtx_U_lc.col(mtx_input.cols()-NthEigenvector).transpose().conjugate()*((mtx_M3-mtx_lambda_transpose.conjugate())*mtx_Smallg1.transpose()) = " << std::endl;
    PrintItReversely(mtx_U_lc.col(mtx_input.cols()-NthEigenvector).transpose().conjugate()*((mtx_M3-mtx_lambda_transpose.conjugate())*mtx_Smallg1.transpose()));

    MatrixXcd mtx_left = mtx_U_lc.col(mtx_input.cols()-NthEigenvector).transpose().conjugate()*(mtx_M3.transpose()-(mtx_M3-mtx_lambda_transpose.conjugate())*mtx_Smallg1.transpose());
    MatrixXcd mtx_right = mtx_L*mtx_lambdanu*mtx_U_rc.transpose();
    std::cout << "mtx_left = " << std::endl;
    PrintItReversely(mtx_left);
    std::cout << "mtx_right = " << std::endl;
    PrintItReversely(mtx_right);
}
MatrixXcd CutoffEigenvalueMatrix(MatrixXcd mtx_input, int cutoff_entry)
{
    MatrixXcd mtx_output = mtx_input;
    for (int entry=1;entry<=mtx_input.rows();entry++)
    {
        if (entry>cutoff_entry)
        {
            mtx_output(mtx_input.rows()-entry,mtx_input.rows()-entry) = 0.;
        }
    }
    return mtx_output;
}
void VerifyLeftRightVectorCoupling(MatrixXcd mtx_input, int entry_cutoff)
{
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "entry_cutoff = " << entry_cutoff << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r = eigensolver_input.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_input_transpose.eigenvectors();
    mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_U_rg = GetSubEigenvectors(mtx_U_r,0);
    MatrixXcd mtx_U_lg = GetSubEigenvectors(mtx_U_l,0);
    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    MatrixXcd mtx_lambda = GetEigenvalueMatrix(mtx_input);
    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_input);
    mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, entry_cutoff);
    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_alpha2 = GetAlphaMatrix(mtx_input, 2);

    MatrixXcd mtx_Smallg1 = mtx_M2*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_BigG1 = mtx_M3-mtx_Smallg1*mtx_M3.transpose();
    MatrixXcd mtx_G1_rc = mtx_BigG1*mtx_U_rc;
    MatrixXcd mtx_lambda_rc = GetEigenvalueWeightedVectors(mtx_lambda,mtx_U_rc);
    std::cout << "mtx_G1_rc-mtx_lambda_rc = " << std::endl;
    std::cout << (mtx_G1_rc-mtx_lambda_rc).block(0,mtx_input.cols()-3,mtx_M3.rows(),3) << std::endl;
    MatrixXcd mtx_lambdanu_lc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc);
    MatrixXcd mtx_lambdanu_R_lc = mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate();
    std::cout << "-mtx_Smallg1*mtx_lambdanu_R_lc = " << std::endl;
    std::cout << (-mtx_Smallg1*mtx_lambdanu_R_lc).block(0,mtx_input.cols()-3,mtx_M3.rows(),3) << std::endl;

    MatrixXcd mtx_Smallg2 = mtx_M1.transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();
    MatrixXcd mtx_BigG2 = mtx_M3.transpose()-mtx_Smallg2*mtx_M3;
    MatrixXcd mtx_G2_lc = mtx_BigG2*mtx_U_lc;
    MatrixXcd mtx_lambda_lc = GetEigenvalueWeightedVectors(mtx_lambda,mtx_U_lc);
    std::cout << "mtx_G2_lc-mtx_lambda_lc = " << std::endl;
    std::cout << (mtx_G2_lc-mtx_lambda_lc).block(0,mtx_input.cols()-3,mtx_M3.rows(),3) << std::endl;
    MatrixXcd mtx_lambdanu_rc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc);
    MatrixXcd mtx_lambdanu_L_rc = mtx_lambdanu_rc.conjugate()*mtx_L.transpose().conjugate();
    std::cout << "-mtx_Smallg2*mtx_lambdanu_L_rc = " << std::endl;
    std::cout << (-mtx_Smallg2*mtx_lambdanu_L_rc).block(0,mtx_input.cols()-3,mtx_M3.rows(),3) << std::endl;

}
MatrixXcd SolveEigenvectorsIteratively_v2(MatrixXcd mtx_input, MatrixXcd mtx_init, int entry_cutoff)
{
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "entry_cutoff = " << entry_cutoff << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_init;
    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_init);
    ComplexEigenSolver<MatrixXcd> eigensolver_init_transpose;
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_init.transpose());
    MatrixXcd mtx_U_r = eigensolver_init.eigenvectors();
    MatrixXcd mtx_U_l = eigensolver_init_transpose.eigenvectors();
    mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_U_rg = GetSubEigenvectors(mtx_U_r,0);
    MatrixXcd mtx_U_lg = GetSubEigenvectors(mtx_U_l,0);
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
    MatrixXcd mtx_lambda = GetEigenvalueMatrix(mtx_init);
    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_init);
    mtx_lambda = CutoffEigenvalueMatrix(mtx_lambda, entry_cutoff);
    mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, entry_cutoff);

    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_alpha2 = GetAlphaMatrix(mtx_input, 2);
    MatrixXcd mtx_Smallg1 = mtx_M2*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_BigG1 = mtx_M3-mtx_Smallg1*mtx_M3.transpose();
    MatrixXcd mtx_Smallg2 = mtx_M1.transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();
    MatrixXcd mtx_BigG2 = mtx_M3.transpose()-mtx_Smallg2*mtx_M3;
    MatrixXcd mtx_unit = GetUnitMatrix(mtx_M3.rows(),mtx_M3.cols());
    MatrixXcd mtx_beta_11 = mtx_M1.transpose()*mtx_alpha1.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_12 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_22 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();

    MatrixXcd mtx_U_rc_new = MatrixXcd::Zero(mtx_U_rc.rows(),mtx_U_rc.cols());
    MatrixXcd mtx_U_lc_new = MatrixXcd::Zero(mtx_U_lc.rows(),mtx_U_lc.cols());
    for (int entry_0=0;entry_0<entry_cutoff;entry_0++)
    {
        MatrixXcd mtx_G1_lambda = mtx_BigG1-mtx_lambda(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_0-1)*mtx_unit;
        VectorXcd mtx_lambdanu_R_g1_lc = VectorXcd::Zero(mtx_U_lc.rows());
        //for (int entry_1=0;entry_1<mtx_input.rows();entry_1++)
        for (int entry_1=0;entry_1<entry_cutoff;entry_1++)
        {
            mtx_lambdanu_R_g1_lc += -(mtx_lambdanu.conjugate()(mtx_input.rows()-entry_1-1,mtx_input.rows()-entry_1-1)*mtx_R(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_1-1)*mtx_Smallg1*mtx_U_lc.col(mtx_input.rows()-entry_1-1).conjugate()).col(0);
        }
        //std::cout << "mtx_G1_lambda*mtx_U_rc.col(mtx_input.rows()-entry_0-1) = " << std::endl;
        //PrintItReversely(mtx_G1_lambda*mtx_U_rc.col(mtx_input.rows()-entry_0-1));
        //std::cout << "mtx_lambdanu_R_g1_lc = " << std::endl;
        //PrintItReversely(mtx_lambdanu_R_g1_lc);
        //mtx_U_rc_new.col(mtx_input.rows()-entry_0-1) = mtx_G1_lambda.inverse()*mtx_lambdanu_R_g1_lc;
        mtx_U_rc_new.col(mtx_input.rows()-entry_0-1) = mtx_G1_lambda.bdcSvd(ComputeThinU | ComputeThinV).solve(mtx_lambdanu_R_g1_lc);

        MatrixXcd mtx_G2_lambda = mtx_BigG2-mtx_lambda(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_0-1)*mtx_unit;
        //MatrixXcd mtx_lambdanu_L_g2_rc = MatrixXcd::Zero(mtx_U_rc.rows(),1);
        VectorXcd mtx_lambdanu_L_g2_rc = VectorXcd::Zero(mtx_U_rc.rows());
        for (int entry_1=0;entry_1<entry_cutoff;entry_1++)
        {
            mtx_lambdanu_L_g2_rc += -(mtx_lambdanu.conjugate()(mtx_input.rows()-entry_1-1,mtx_input.rows()-entry_1-1)*mtx_L.conjugate()(mtx_input.rows()-entry_1-1,mtx_input.rows()-entry_0-1)*mtx_Smallg2*mtx_U_rc.col(mtx_input.rows()-entry_1-1).conjugate()).col(0);
        }
        //mtx_U_lc_new.col(mtx_input.rows()-entry_0-1) = mtx_G2_lambda.inverse()*mtx_lambdanu_L_g2_rc;
        mtx_U_lc_new.col(mtx_input.rows()-entry_0-1) = mtx_G2_lambda.bdcSvd(ComputeThinU | ComputeThinV).solve(mtx_lambdanu_L_g2_rc);
    }
    mtx_U_rc_new = MakeRealEigenvectors(mtx_U_rc_new);
    mtx_U_lc_new = MakeRealEigenvectors(mtx_U_lc_new);
    std::cout << "mtx_U_rc = " << std::endl;
    std::cout << mtx_U_rc.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_U_rc_new = " << std::endl;
    std::cout << mtx_U_rc_new.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_U_lc = " << std::endl;
    std::cout << mtx_U_lc.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_U_lc_new = " << std::endl;
    std::cout << mtx_U_lc_new.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;

    ////MatrixXcd mtx_lambdanu_lc_new = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc_new);
    ////MatrixXcd mtx_lambdanu_lc_transpose_new = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc_new).transpose();
    //MatrixXcd mtx_lambdanu_lc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc);
    //MatrixXcd mtx_lambdanu_lc_transpose = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc).transpose();
    //MatrixXcd mtx_R_predict_term1 = mtx_U_rc_new.transpose().conjugate()*mtx_U_rc_new;
    ////MatrixXcd mtx_R_predict_term2 = (mtx_U_rc_new.transpose().conjugate()*mtx_M3-mtx_R*mtx_lambdanu_lc_transpose_new)*mtx_beta_11*(mtx_M3.transpose()*mtx_U_rc_new-mtx_lambdanu_lc_new.conjugate()*mtx_R.transpose().conjugate());
    //MatrixXcd mtx_R_predict_term2 = (mtx_U_rc.transpose().conjugate()*mtx_M3-mtx_R*mtx_lambdanu_lc_transpose)*mtx_beta_11*(mtx_M3.transpose()*mtx_U_rc-mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate());
    //MatrixXcd mtx_R_predict = mtx_R_predict_term1 + mtx_R_predict_term2;
    //std::cout << "mtx_R = " << std::endl;
    //PrintItReversely(mtx_R);
    //std::cout << "mtx_R_predict = " << std::endl;
    //PrintItReversely(mtx_R_predict);

    ////MatrixXcd mtx_lambdanu_rc_new = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc_new);
    ////MatrixXcd mtx_lambdanu_rc_transpose_new = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc_new).transpose();
    //MatrixXcd mtx_lambdanu_rc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc);
    //MatrixXcd mtx_lambdanu_rc_transpose = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc).transpose();
    //MatrixXcd mtx_L_predict_term1 = mtx_U_lc_new.transpose().conjugate()*mtx_U_lc_new;
    ////MatrixXcd mtx_L_predict_term2 = (mtx_U_lc_new.transpose().conjugate()*mtx_M3.transpose()-mtx_L*mtx_lambdanu_rc_transpose_new)*mtx_beta_22*(mtx_M3*mtx_U_lc_new-mtx_lambdanu_rc_new.conjugate()*mtx_L.transpose().conjugate());
    //MatrixXcd mtx_L_predict_term2 = (mtx_U_lc.transpose().conjugate()*mtx_M3.transpose()-mtx_L*mtx_lambdanu_rc_transpose)*mtx_beta_22*(mtx_M3*mtx_U_lc-mtx_lambdanu_rc.conjugate()*mtx_L.transpose().conjugate());
    //MatrixXcd mtx_L_predict = mtx_L_predict_term1 + mtx_L_predict_term2;
    //std::cout << "mtx_L = " << std::endl;
    //PrintItReversely(mtx_L);
    //std::cout << "mtx_L_predict = " << std::endl;
    //PrintItReversely(mtx_L_predict);

    //MatrixXcd mtx_H_predict_term1 = mtx_U_lc_new.transpose()*mtx_U_rc_new;
    ////MatrixXcd mtx_H_predict_term2 = (mtx_U_lc_new.transpose()*mtx_M3.transpose()-mtx_L.conjugate()*mtx_lambdanu_rc_transpose_new.conjugate())*mtx_beta_12*(mtx_M3.transpose()*mtx_U_rc_new-mtx_lambdanu_lc_new.conjugate()*mtx_R.transpose().conjugate());
    //MatrixXcd mtx_H_predict_term2 = (mtx_U_lc.transpose()*mtx_M3.transpose()-mtx_L.conjugate()*mtx_lambdanu_rc_transpose.conjugate())*mtx_beta_12*(mtx_M3.transpose()*mtx_U_rc-mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate());
    //MatrixXcd mtx_H_predict = mtx_H_predict_term1 + mtx_H_predict_term2;
    //std::cout << "mtx_H = " << std::endl;
    //PrintItReversely(mtx_H);
    //std::cout << "mtx_H_predict = " << std::endl;
    //PrintItReversely(mtx_H_predict);

    //MatrixXcd mtx_lambda_new = MatrixXcd::Zero(mtx_lambda.rows(),mtx_lambda.cols());
    //for (int row=0;row<mtx_lambda_new.rows();row++)
    //{
    //    for (int col=0;col<mtx_lambda_new.cols();col++)
    //    {
    //        mtx_lambda_new(row,col) = 0.;
    //        if (row==col)
    //        {
    //            mtx_lambda_new(row,col) = mtx_lambdanu(row,col)*mtx_H_predict(row,col);
    //        }
    //    }
    //}
    //std::cout << "mtx_lambda = " << std::endl;
    //PrintItReversely(mtx_lambda);
    //std::cout << "mtx_lambda_new = " << std::endl;
    //PrintItReversely(mtx_lambda_new);

    MatrixXcd mtx_rg_predict = Convert_hadron_to_gamma_vector_v2(mtx_input, mtx_U_rc_new, mtx_lambda, 2);
    MatrixXcd mtx_lg_predict = Convert_hadron_to_gamma_vector_v2(mtx_input, mtx_U_lc_new, mtx_lambda, 1);
    std::cout << "mtx_U_rg = " << std::endl;
    std::cout << mtx_U_rg.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_rg_predict = " << std::endl;
    std::cout << mtx_rg_predict.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_U_lg = " << std::endl;
    std::cout << mtx_U_lg.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    std::cout << "mtx_lg_predict = " << std::endl;
    std::cout << mtx_lg_predict.block(0,mtx_input.cols()-entry_cutoff,mtx_M3.rows(),entry_cutoff) << std::endl;
    MatrixXcd mtx_U_r_new = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    MatrixXcd mtx_U_l_new = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    int start_row;
    int start_col;
    int size_row;
    int size_col;
    start_row = 0;
    start_col = 0;
    size_row = (binx_blind_global)-start_row;
    size_col = mtx_input.cols();
    mtx_U_r_new.block(start_row,start_col,size_row,size_col) = mtx_rg_predict;
    mtx_U_l_new.block(start_row,start_col,size_row,size_col) = mtx_lg_predict;
    start_row = binx_blind_global;
    start_col = 0;
    size_row = mtx_input.rows()-(binx_blind_global);
    size_col = mtx_input.cols();
    mtx_U_r_new.block(start_row,start_col,size_row,size_col) = mtx_U_rc_new;
    mtx_U_l_new.block(start_row,start_col,size_row,size_col) = mtx_U_lc_new;

    MatrixXcd mtx_output = mtx_U_r_new*mtx_lambdanu*mtx_U_l_new.transpose();
    return mtx_output;

}
void SolveEigenvectorsIteratively(MatrixXcd mtx_input, MatrixXcd mtx_init, int entry_cutoff)
{

    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r_input = eigensolver_input.eigenvectors();
    mtx_U_r_input = MakeRealEigenvectors(mtx_U_r_input);
    MatrixXcd mtx_U_l_input = eigensolver_input_transpose.eigenvectors();
    mtx_U_l_input = MakeRealEigenvectors(mtx_U_l_input);
    MatrixXcd mtx_U_rc_input = GetSubEigenvectors(mtx_U_r_input,1);
    MatrixXcd mtx_U_lc_input = GetSubEigenvectors(mtx_U_l_input,1);
    MatrixXcd mtx_R_input = mtx_U_r_input.transpose().conjugate()*mtx_U_r_input;
    MatrixXcd mtx_L_input = mtx_U_l_input.transpose().conjugate()*mtx_U_l_input;
    MatrixXcd mtx_lambda_input = GetEigenvalueMatrix(mtx_input);
    MatrixXcd mtx_lambdanu_input = GetLambdaNuMatrix(mtx_input);

    MatrixXcd mtx_M1_init = GetSubmatrix(mtx_init, 1);
    MatrixXcd mtx_M2_init = GetSubmatrix(mtx_init, 2);
    MatrixXcd mtx_M3_init = GetSubmatrix(mtx_init, 3);
    MatrixXcd mtx_alpha1_init = GetAlphaMatrix(mtx_init, 1);
    MatrixXcd mtx_alpha2_init = GetAlphaMatrix(mtx_init, 2);
    MatrixXcd mtx_Smallg1_init = mtx_M2_init*mtx_alpha1_init.inverse()*mtx_M1_init;
    MatrixXcd mtx_BigG1_init = mtx_M3_init-mtx_Smallg1_init*mtx_M3_init.transpose();
    MatrixXcd mtx_Smallg2_init = mtx_M1_init.transpose()*mtx_alpha2_init.inverse()*mtx_M2_init.transpose();
    MatrixXcd mtx_BigG2_init = mtx_M3_init.transpose()-mtx_Smallg2_init*mtx_M3_init;

    ComplexEigenSolver<MatrixXcd> eigensolver_init;
    eigensolver_init = ComplexEigenSolver<MatrixXcd>(mtx_init);
    ComplexEigenSolver<MatrixXcd> eigensolver_init_transpose;
    eigensolver_init_transpose = ComplexEigenSolver<MatrixXcd>(mtx_init.transpose());
    MatrixXcd mtx_U_r = eigensolver_init.eigenvectors();
    mtx_U_r = MakeRealEigenvectors(mtx_U_r);
    MatrixXcd mtx_U_l = eigensolver_init_transpose.eigenvectors();
    mtx_U_l = MakeRealEigenvectors(mtx_U_l);
    MatrixXcd mtx_U_rc = GetSubEigenvectors(mtx_U_r,1);
    MatrixXcd mtx_U_lc = GetSubEigenvectors(mtx_U_l,1);
    MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
    MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
    MatrixXcd mtx_lambda = GetEigenvalueMatrix(mtx_init);
    MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_init);

    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);
    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_alpha2 = GetAlphaMatrix(mtx_input, 2);
    MatrixXcd mtx_Smallg1 = mtx_M2*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_BigG1 = mtx_M3-mtx_Smallg1*mtx_M3.transpose();
    MatrixXcd mtx_Smallg2 = mtx_M1.transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();
    MatrixXcd mtx_BigG2 = mtx_M3.transpose()-mtx_Smallg2*mtx_M3;
    MatrixXcd mtx_unit = GetUnitMatrix(mtx_M3.rows(),mtx_M3.cols());
    //mtx_Smallg1 = mtx_unit;
    //mtx_Smallg2 = mtx_unit;

    MatrixXcd mtx_Big = MatrixXcd::Zero(mtx_M3.rows()*2*entry_cutoff,mtx_M3.rows()*2*entry_cutoff);
    MatrixXcd mtx_G1_lambda = MatrixXcd::Zero(mtx_M3.rows(),mtx_M3.rows());
    MatrixXcd mtx_G2_lambda = MatrixXcd::Zero(mtx_M3.rows(),mtx_M3.rows());
    MatrixXcd mtx_R_g1 = MatrixXcd::Zero(mtx_M3.rows(),mtx_M3.rows());
    MatrixXcd mtx_L_g2 = MatrixXcd::Zero(mtx_M3.rows(),mtx_M3.rows());
    int first_row;
    int first_col;
    for (int entry_0=0;entry_0<entry_cutoff;entry_0++)
    {
        for (int entry_1=0;entry_1<entry_cutoff;entry_1++)
        {
            mtx_G1_lambda = mtx_BigG1-mtx_lambda(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_1-1)*mtx_unit;
            mtx_G2_lambda = mtx_BigG2-mtx_lambda(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_1-1)*mtx_unit;
            mtx_R_g1 = mtx_lambdanu(mtx_input.rows()-entry_1-1,mtx_input.rows()-entry_1-1)*mtx_R(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_1-1)*mtx_Smallg1;
            mtx_L_g2 = mtx_lambdanu(mtx_input.rows()-entry_1-1,mtx_input.rows()-entry_1-1)*mtx_L(mtx_input.rows()-entry_0-1,mtx_input.rows()-entry_1-1)*mtx_Smallg2;
            if (entry_0==entry_1)
            {
                first_row = entry_0*mtx_M3.rows();
                first_col = entry_1*mtx_M3.rows();
                mtx_Big.block(first_row,first_col,mtx_M3.rows(),mtx_M3.rows()) = mtx_G1_lambda;
                first_row = entry_0*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
                first_col = entry_1*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
                mtx_Big.block(first_row,first_col,mtx_M3.rows(),mtx_M3.rows()) = mtx_G2_lambda;
            }
            first_row = entry_0*mtx_M3.rows();
            first_col = entry_1*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
            mtx_Big.block(first_row,first_col,mtx_M3.rows(),mtx_M3.rows()) = mtx_R_g1;
            first_row = entry_0*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
            first_col = entry_1*mtx_M3.rows();
            mtx_Big.block(first_row,first_col,mtx_M3.rows(),mtx_M3.rows()) = mtx_L_g2;
        }
    }
    //std::cout << "mtx_Big = " << std::endl;
    //PrintItReversely(mtx_Big);

    mtx_lambdanu = CutoffEigenvalueMatrix(mtx_lambdanu, entry_cutoff);
    MatrixXcd mtx_G1_rc = mtx_BigG1*mtx_U_rc;
    MatrixXcd mtx_lambda_rc = GetEigenvalueWeightedVectors(mtx_lambda,mtx_U_rc);
    MatrixXcd mtx_G2_lc = mtx_BigG2*mtx_U_lc;
    MatrixXcd mtx_lambda_lc = GetEigenvalueWeightedVectors(mtx_lambda,mtx_U_lc);
    MatrixXcd mtx_lambdanu_lc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_lc);
    MatrixXcd mtx_lambdanu_R_lc = mtx_lambdanu_lc.conjugate()*mtx_R.transpose().conjugate();
    MatrixXcd mtx_remain_lc = mtx_G1_rc-mtx_lambda_rc+mtx_Smallg1*mtx_lambdanu_R_lc;
    MatrixXcd mtx_lambdanu_rc = GetEigenvalueWeightedVectors(mtx_lambdanu,mtx_U_rc);
    MatrixXcd mtx_lambdanu_L_rc = mtx_lambdanu_rc.conjugate()*mtx_L.transpose().conjugate();
    MatrixXcd mtx_remain_rc = mtx_G2_lc-mtx_lambda_lc+mtx_Smallg2*mtx_lambdanu_L_rc;
    //std::cout << "mtx_remain_lc = " << std::endl;
    //PrintItReversely(mtx_remain_lc);
    //std::cout << "mtx_remain_rc = " << std::endl;
    //PrintItReversely(mtx_remain_rc);

    VectorXcd vtx_remain = VectorXcd::Zero(mtx_M3.rows()*2*entry_cutoff);
    VectorXcd vtx_truth = VectorXcd::Zero(mtx_M3.rows()*2*entry_cutoff);
    for (int entry_0=0;entry_0<entry_cutoff;entry_0++)
    {
        first_row = entry_0*mtx_M3.rows();
        vtx_remain.segment(first_row,mtx_M3.rows()) = mtx_remain_lc.col(mtx_input.rows()-entry_0-1).conjugate();
        vtx_truth.segment(first_row,mtx_M3.rows()) = mtx_U_rc.col(mtx_input.rows()-entry_0-1).conjugate();
        first_row = entry_0*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
        vtx_remain.segment(first_row,mtx_M3.rows()) = mtx_remain_rc.col(mtx_input.rows()-entry_0-1);
        vtx_truth.segment(first_row,mtx_M3.rows()) = mtx_U_lc.col(mtx_input.rows()-entry_0-1);
    }
    //std::cout << "vtx_remain = " << std::endl;
    //std::cout << vtx_remain << std::endl;
    //std::cout << "vtx_truth = " << std::endl;
    //std::cout << vtx_truth << std::endl;
    //std::cout << "mtx_Big*vtx_truth = " << std::endl;
    //std::cout << mtx_Big*vtx_truth << std::endl;

    VectorXcd vtx_solution = mtx_Big.bdcSvd(ComputeThinU | ComputeThinV).solve(vtx_remain);
    //VectorXcd vtx_solution = (mtx_Big.transpose() * mtx_Big).ldlt().solve(mtx_Big.transpose() * vtx_remain);
    //MatrixXcd vtx_solution = mtx_Big.fullPivLu().kernel();
    for (int entry_0=0;entry_0<entry_cutoff;entry_0++)
    {
        std::cout << "vtx_solution rc = " << std::endl;
        first_row = entry_0*mtx_M3.rows();
        std::cout << vtx_solution.segment(first_row,mtx_M3.rows()) << std::endl;
        std::cout << "vtx_solution lc = " << std::endl;
        first_row = entry_0*mtx_M3.rows()+entry_cutoff*mtx_M3.rows();
        std::cout << vtx_solution.segment(first_row,mtx_M3.rows()) << std::endl;
    }
    std::cout << "mtx_U_rc_input = " << std::endl;
    std::cout << mtx_U_rc_input << std::endl;
    std::cout << "mtx_U_lc_input = " << std::endl;
    std::cout << mtx_U_lc_input << std::endl;
    std::cout << "mtx_U_rc (init) = " << std::endl;
    std::cout << mtx_U_rc << std::endl;
    std::cout << "mtx_U_lc (init) = " << std::endl;
    std::cout << mtx_U_lc << std::endl;

}
pair<std::complex<double>,std::complex<double>> GetIterativeMatrix(MatrixXcd mtx_input, int NthEigenvector, std::complex<double> par_R, std::complex<double> par_L)
{
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    ComplexEigenSolver<MatrixXcd> eigensolver_input;
    eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose;
    eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_U_r_input = eigensolver_input.eigenvectors();
    mtx_U_r_input = MakeRealEigenvectors(mtx_U_r_input);
    MatrixXcd mtx_U_l_input = eigensolver_input_transpose.eigenvectors();
    mtx_U_l_input = MakeRealEigenvectors(mtx_U_l_input);
    MatrixXcd mtx_U_rc_input = GetSubEigenvectors(mtx_U_r_input,1);
    MatrixXcd mtx_U_lc_input = GetSubEigenvectors(mtx_U_l_input,1);

    MatrixXcd mtx_M1 = GetSubmatrix(mtx_input, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_input, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_input, 3);

    MatrixXcd mtx_alpha1 = GetAlphaMatrix(mtx_input, 1);
    MatrixXcd mtx_Smallg1 = mtx_M2*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_BigG1 = mtx_M3-mtx_Smallg1*mtx_M3.transpose();

    MatrixXcd mtx_alpha2 = GetAlphaMatrix(mtx_input, 2);
    MatrixXcd mtx_Smallg2 = mtx_M1.transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();
    MatrixXcd mtx_BigG2 = mtx_M3.transpose()-mtx_Smallg2*mtx_M3;
    MatrixXcd mtx_unit = GetUnitMatrix(mtx_M3.rows(),mtx_M3.cols());

    MatrixXcd mtx_result = MatrixXcd::Zero(mtx_M3.rows()+mtx_M3.rows(),mtx_M3.rows()+mtx_M3.rows());
    mtx_result.block(0,0,mtx_M3.rows(),mtx_M3.rows()) = mtx_BigG1;
    mtx_result.block(mtx_M3.rows(),mtx_M3.rows(),mtx_M3.rows(),mtx_M3.rows()) = mtx_BigG2;
    mtx_result.block(0,mtx_M3.rows(),mtx_M3.rows(),mtx_M3.rows()) = par_R*mtx_Smallg1;
    mtx_result.block(mtx_M3.rows(),0,mtx_M3.rows(),mtx_M3.rows()) = par_L*mtx_Smallg2;
    //std::cout << "the big matrix = " << std::endl;
    //PrintItReversely(mtx_result);

    ComplexEigenSolver<MatrixXcd> eigensolver_iterative;
    eigensolver_iterative = ComplexEigenSolver<MatrixXcd>(mtx_result);
    std::cout << "new lambda = " << eigensolver_iterative.eigenvalues()(mtx_result.rows()-1) << std::endl;
    //std::cout << "new eigenvector = " << std::endl;
    //std::cout << eigensolver_iterative.eigenvectors().block(0,mtx_result.cols()-1,mtx_result.rows(),1) << std::endl;
    std::complex<double> lambda_new = eigensolver_iterative.eigenvalues()(mtx_result.rows()-1);
    MatrixXcd mtx_U_rc_new = eigensolver_iterative.eigenvectors().block(0,mtx_result.cols()-1,mtx_M3.rows(),1);
    mtx_U_rc_new = MakeRealEigenvectors(mtx_U_rc_new);
    MatrixXcd mtx_U_lc_new = eigensolver_iterative.eigenvectors().block(mtx_M3.rows(),mtx_result.cols()-1,mtx_M3.rows(),1);
    mtx_U_lc_new = MakeRealEigenvectors(mtx_U_lc_new);
    MatrixXcd mtx_beta_11 = mtx_M1.transpose()*mtx_alpha1.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_12 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha1.inverse()*mtx_M1;
    MatrixXcd mtx_beta_22 = mtx_M2*mtx_alpha2.inverse().transpose()*mtx_alpha2.inverse()*mtx_M2.transpose();
    std::cout << "mtx_U_lc_new = " << std::endl;
    std::cout << mtx_U_lc_new << std::endl;
    std::cout << "mtx_U_rc_new = " << std::endl;
    std::cout << mtx_U_rc_new << std::endl;

    MatrixXcd mtx_R_new = mtx_U_rc_new.transpose().conjugate()*mtx_U_rc_new;
    mtx_R_new += (mtx_U_rc_new.transpose().conjugate()*mtx_M3-par_R*mtx_U_lc_new.transpose())*mtx_beta_11*(mtx_M3.transpose()*mtx_U_rc_new-std::conj(par_R)*mtx_U_lc_new.conjugate());
    std::cout << "mtx_R_new = " << std::endl;
    PrintItReversely(mtx_R_new);

    MatrixXcd mtx_L_new = mtx_U_lc_new.transpose().conjugate()*mtx_U_lc_new;
    mtx_L_new += (mtx_U_lc_new.transpose().conjugate()*mtx_M3.transpose()-par_L*mtx_U_rc_new.transpose())*mtx_beta_22*(mtx_M3*mtx_U_lc_new-std::conj(par_L)*mtx_U_rc_new.conjugate());
    std::cout << "mtx_L_new = " << std::endl;
    PrintItReversely(mtx_L_new);

    MatrixXcd mtx_H_new = mtx_U_lc_new.transpose()*mtx_U_rc_new;
    mtx_H_new += (mtx_U_lc_new.transpose()*mtx_M3.transpose()-std::conj(par_L)*mtx_U_rc_new.transpose().conjugate())*mtx_beta_12*(mtx_M3.transpose()*mtx_U_rc_new-std::conj(par_R)*mtx_U_lc_new.conjugate());
    std::cout << "mtx_H_new = " << std::endl;
    PrintItReversely(mtx_H_new);

    return std::make_pair(mtx_R_new(0,0)*lambda_new/mtx_H_new(0,0),mtx_L_new(0,0)*lambda_new/mtx_H_new(0,0));
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
void fill1DHistogram(TH1D* hist,VectorXcd vtr)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        double real = vtr(binx).real();
        double imag = vtr(binx).imag();
        double ampl = pow(real*real+imag*imag,0.5);
        if (real!=0.)
        {
            ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
        }
        hist->SetBinContent(binx+1,ampl);
        hist->SetBinError(binx+1,0.01*ampl);
    }
}
double NormalizationFactorForEigenvector(VectorXcd vtr, bool invert)
{
    TH1D hist_ampl = TH1D("hist_ampl","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper);
    fill1DHistogram(&hist_ampl,vtr);
    int bin_blind = hist_ampl.GetXaxis()->FindBin(1.);
    double integral_1 = 0.;
    if (!invert) integral_1 = hist_ampl.Integral(bin_blind,hist_ampl.GetNbinsX());
    else integral_1 = hist_ampl.Integral(1,bin_blind);
    return integral_1;
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
            if (!invert_y) matrix(binx,biny) = hist->GetBinContent(binx+1,biny+1);
            else matrix(binx,biny) = hist->GetBinContent(binx+1,hist->GetNbinsY()-biny);
        }
    }
    return matrix;
}
MatrixXcd BuildModelMatrix()
{
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    if (linear) 
    {
        mtx_model = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;
        mtx_model += mtx_eigenvector_vari*mtx_eigenvalue_init*mtx_eigenvector_inv_init;
        mtx_model += mtx_eigenvector_init*mtx_eigenvalue_vari*mtx_eigenvector_inv_init;
        mtx_model += mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_vari;
    }
    else
    {
        mtx_eigenvalue = mtx_eigenvalue_init+mtx_eigenvalue_vari;
        mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
        mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
        mtx_model = mtx_eigenvector*mtx_eigenvalue*mtx_eigenvector_inv;
    }
    return mtx_model;
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

double Factorial_Stirling_Series(double n)
{
    if (n==0) return 1.;
    if (n==1) return 1.;
    return pow(2*M_PI*n,0.5)*pow(n,n)/exp(n)*(1.+1./(12.*n));
}

double BlindedLogLikelihood(TH2D* hist_data, TH2D* hist_dark, TH2D* hist_model)
{
    int binx_blind = hist_data->GetXaxis()->FindBin(MSCL_cut_blind);
    int biny_blind = hist_data->GetYaxis()->FindBin(MSCW_cut_blind);
    int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind+1.);
    int biny_upper = hist_data->GetXaxis()->FindBin(MSCW_cut_blind+1.);
    //int binx_blind = hist_data->GetXaxis()->FindBin(1.);
    //int biny_blind = hist_data->GetYaxis()->FindBin(1.);
    //int binx_upper = hist_data->GetXaxis()->FindBin(3.);
    //int biny_upper = hist_data->GetYaxis()->FindBin(3.);
    //int binx_upper = hist_data->GetXaxis()->FindBin(3.);
    //int biny_upper = hist_data->GetYaxis()->FindBin(2.);
    double likelihood = 1.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>=binx_blind || by>=biny_blind)
            {
                if (bx>=binx_upper || by>=biny_upper) continue;
                double data = hist_data->GetBinContent(bx,by);
                double dark = hist_dark->GetBinContent(bx,by);
                double model = hist_model->GetBinContent(bx,by);
                double likelihood_this = pow(model,data)*exp(-model)/Factorial_Stirling_Series(data);
                if (isnan(likelihood_this))
                {
                    std::cout << "pow(model,data) = " << pow(model,data) << std::endl;
                    std::cout << "exp(-model) = " << exp(-model) << std::endl;
                    std::cout << "Factorial_Stirling_Series(data) = " << Factorial_Stirling_Series(data) << std::endl;
                }
                if (likelihood_this>0.) likelihood = likelihood*likelihood_this;
            }
        }
    }
    return 1./likelihood;
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
    //int binx_blind = hist_data->GetXaxis()->FindBin(1.);
    //int biny_blind = hist_data->GetYaxis()->FindBin(1.);
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
                double data_err = max(1.0,pow(data,0.5));
                double model_err = max(1.0,pow(abs(model),0.5));
                weight = 1./(data_err*data_err+model_err*model_err);
                //weight = 1./(data*data+model*model);
                //weight = 1./(data_err*data_err);
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
    //int binx_upper = hist_data->GetXaxis()->FindBin(MSCL_cut_blind+1.);
    //int biny_upper = hist_data->GetXaxis()->FindBin(MSCW_cut_blind+1.);
    //int binx_blind = hist_data->GetXaxis()->FindBin(1.);
    //int biny_blind = hist_data->GetYaxis()->FindBin(1.);
    //int binx_upper = hist_data->GetXaxis()->FindBin(3.);
    //int biny_upper = hist_data->GetYaxis()->FindBin(3.);
    int binx_upper = hist_data->GetNbinsX();
    int biny_upper = hist_data->GetNbinsY();
    double chi2 = 0.;
    double nbins = 0.;
    for (int bx=1;bx<=hist_data->GetNbinsX();bx++)
    {
        for (int by=1;by<=hist_data->GetNbinsY();by++)
        {
            if (bx>=binx_blind || by>=biny_blind)
            {
                if (bx>=binx_upper || by>=biny_upper) continue;
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
                //weight = 1./(data_err*data_err+model_err*model_err);
                //weight = 1./(data*data+model*model);
                //weight = 1./(data_err*data_err);
                if (data-model<0.) weight = 2.;
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
                nbins += 1.;
            }
        }
    }

    //mtx_data_bkgd = fillMatrix(hist_model);
    //eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);
    //for (int NthEigenvector=2;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //{
    //    double eigenvalue_ratio_bkgd = eigensolver_bkgd.eigenvalues()(mtx_dark.cols()-NthEigenvector).real()-eigensolver_bkgd.eigenvalues()(mtx_dark.cols()-1).real();
    //    double eigenvalue_ratio_dark = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector).real()-eigensolver_dark.eigenvalues()(mtx_dark.cols()-1).real();
    //    double eigenvalue_diff = (eigenvalue_ratio_bkgd-eigenvalue_ratio_dark)*eigensolver_dark.eigenvalues()(mtx_dark.cols()-1).real();
    //    chi2 += pow(eigenvalue_diff,2);
    //}
    
    return chi2;
}

void NetflixParametrizeEigenvectorsEigenBasis(const double *par)
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

    //mtx_eigenvector = mtx_eigenvector_init;
    //mtx_eigenvector_inv = mtx_eigenvector_inv_init;
    //mtx_eigenvalue = mtx_eigenvalue_init;

    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_vari(row,col_fix) = 0.;
            for (int NthBasis=1;NthBasis<=NumberOfEigenvectors;NthBasis++)
            {
                mtx_eigenvector_vari(row,col_fix) += par[first_index+NthBasis-1]*mtx_eigenvector_dark(row,N_bins_for_deconv-NthBasis);
            }
        }

        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_vari(row_fix,col) = 0.;
            for (int NthBasis=1;NthBasis<=NumberOfEigenvectors;NthBasis++)
            {
                mtx_eigenvector_inv_vari(row_fix,col) += par[first_index+NthBasis-1]*mtx_eigenvector_inv_dark(N_bins_for_deconv-NthBasis,col);
            }
        }

        // build eigenvalue matrix
        first_index = 1+(2*NthEigenvector-0)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            mtx_eigenvalue_vari(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvalue) = par[first_index+NthEigenvalue-1];
        }
    }
}

void NetflixParametrizeEigenvectors(const double *par)
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
    //{
    //    mtx_eigenvector.col(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector);
    //    mtx_eigenvector_inv.row(mtx_dark.cols()-NthEigenvector) = mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector);
    //    mtx_eigenvalue(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector);
    //}
    //mtx_eigenvector = mtx_eigenvector_init;
    //mtx_eigenvector_inv = mtx_eigenvector_inv_init;
    //mtx_eigenvalue = mtx_eigenvalue_init;

    const std::complex<double> If(0.0, 1.0);
    // build eigenvector matrix
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;
    
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {
        col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_vari(row,col_fix) = par[first_index+row];
        }

        row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_vari(row_fix,col) = par[first_index+col];
        }

        // build eigenvalue matrix
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            mtx_eigenvalue_vari(N_bins_for_deconv-NthEigenvector,N_bins_for_deconv-NthEigenvalue) = par[first_index+NthEigenvalue-1];
        }
    }
}


double NetflixChi2Function(const double *par)
{

    if (eigenbasis)
    {
        NetflixParametrizeEigenvectorsEigenBasis(par);
    }
    else
    {
        NetflixParametrizeEigenvectors(par);
    }
    //double deriv_at_zero = SmoothEigenvectors(&mtx_eigenvector, &mtx_eigenvector_inv);
    //if (0.8*init_deriv_at_zero>deriv_at_zero) return 1e10;
    //if (1.2*init_deriv_at_zero<deriv_at_zero) return 1e10;
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    mtx_model = BuildModelMatrix();

    TH2D hist_gamma = TH2D("hist_gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_diff = TH2D("hist_diff","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_data = TH2D("hist_data","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_dark = TH2D("hist_dark","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);

    fill2DHistogramAbs(&hist_gamma,mtx_gamma);
    fill2DHistogram(&hist_data,mtx_data);
    fill2DHistogram(&hist_dark,mtx_dark);
    fill2DHistogram(&hist_model,mtx_model);

    int binx_lower = hist_diff.GetXaxis()->FindBin(MSCL_cut_lower);
    int binx_blind = hist_diff.GetXaxis()->FindBin(MSCL_cut_blind)-1;
    int binx_upper = hist_diff.GetXaxis()->FindBin(1.)-1;
    int biny_lower = hist_diff.GetYaxis()->FindBin(MSCW_cut_lower);
    int biny_blind = hist_diff.GetYaxis()->FindBin(MSCW_cut_blind)-1;
    int biny_upper = hist_diff.GetYaxis()->FindBin(1.)-1;
    hist_diff.Add(&hist_data);
    hist_diff.Add(&hist_model,-1.);
    double gamma_total = hist_diff.Integral(binx_lower,binx_upper,biny_lower,biny_upper);
    //double gamma_total = hist_diff.Integral();
    gamma_total = max(0.,gamma_total);
    double scale = gamma_total/double(hist_gamma.Integral(binx_lower,binx_upper,biny_lower,biny_upper));
    //double scale = gamma_total/double(hist_gamma.Integral());
    hist_gamma.Scale(scale);

    double chi2 = 0.;
    chi2 += BlindedChi2(&hist_data,&hist_dark,&hist_model);
    //if (signal_model) chi2 += SignalChi2(&hist_data,&hist_gamma,&hist_model);
    
    //double chi2 = BlindedLogLikelihood(&hist_data,&hist_dark,&hist_model);
    //chi2 += SignalLogLikelihood(&hist_data,&hist_gamma,&hist_model);

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
double RatioEmptyBins(MatrixXcd mtx, int binx_blind, int biny_blind)
{
    double ratio = 0.;
    double n_empty_bins = 0.;
    double n_filled_bins = 0.;
    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            if (row<=binx_blind && col<=biny_blind) continue; // signal region
            if (row>binx_blind && col>biny_blind) continue;
            if (mtx(row,col).real()<=100.)
            {
                n_empty_bins += 1.;
            }
            else
            {
                n_filled_bins += 1.;
            }
        }
    }
    return n_empty_bins/(n_empty_bins+n_filled_bins);
}
void SetInitialEigenvectors(int binx_blind, int biny_blind)
{
    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    const std::complex<double> If(0.0, 1.0);
    int col_fix = 0;
    int row_fix = 0;
    int first_index = 0;

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_init(row,col) = 0.;
            mtx_eigenvector_vari(row,col) = 0.;
            mtx_eigenvector_inv_init(row,col) = 0.;
            mtx_eigenvector_inv_vari(row,col) = 0.;
            mtx_eigenvalue_init(row,col) = 0.;
            mtx_eigenvalue_vari(row,col) = 0.;
            double real = eigensolver_dark.eigenvectors().col(col)(row).real();
            double imag = eigensolver_dark.eigenvectors().col(col)(row).imag();
            double ampl = pow(real*real+imag*imag,0.5);
            if (real!=0.)
            {
                ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
            }
            mtx_eigenvector_dark.col(col)(row) = ampl;
            double real_inv = eigensolver_dark.eigenvectors().inverse().row(row)(col).real();
            double imag_inv = eigensolver_dark.eigenvectors().inverse().row(row)(col).imag();
            double ampl_inv = pow(real_inv*real_inv+imag_inv*imag_inv,0.5);
            if (real_inv!=0.)
            {
                ampl_inv = real_inv/abs(real_inv)*pow(real_inv*real_inv+imag_inv*imag_inv,0.5);
            }
            mtx_eigenvector_inv_dark.row(row)(col) = ampl_inv;
        }
    }

    for (int NthEigenvector=NumberOfEigenvectors+1;NthEigenvector<=N_bins_for_deconv;NthEigenvector++)
    {
        mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row);
        }
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col);
        }
    }

    double first_eigenvalue_dark = eigensolver_dark.eigenvalues()(mtx_dark.cols()-1).real();
    double first_eigenvalue_data = eigensolver_data.eigenvalues()(mtx_dark.cols()-1).real();
    double eigenvalue_ratio = first_eigenvalue_data/first_eigenvalue_dark;
    //if (dark_vector<2) eigenvalue_ratio = 1.;
    std::cout << "eigenvalue_ratio = " << eigenvalue_ratio << std::endl;
    //for (int NthEigenvector=1;NthEigenvector<=N_bins_for_deconv;NthEigenvector++)
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        //mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector);
        if (NthEigenvector>=dark_vector)
        {
            mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector);
        }
        else
        {
            mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_data.eigenvalues()(mtx_dark.cols()-NthEigenvector);
        }

        double sign = 1.;
        if (eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector).dot(eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector)).real()<0.)
        {
            sign = -1.;
        }
        double sign_inv = 1.;
        if (eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector).dot(eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)).real()<0.)
        {
            sign_inv = -1.;
        }
        double norm_data = NormalizationFactorForEigenvector(eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector),false);
        double norm_dark = NormalizationFactorForEigenvector(eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector),false);
        std::cout << NthEigenvector << "-th eigenvector, norm_data = " << norm_data << std::endl;
        std::cout << NthEigenvector << "-th eigenvector, norm_dark = " << norm_dark << std::endl;
        double norm_data_inv = NormalizationFactorForEigenvector(eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector),invert_y);
        double norm_dark_inv = NormalizationFactorForEigenvector(eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector),invert_y);
        std::cout << NthEigenvector << "-th inv. eigenvector, norm_data_inv = " << norm_data_inv << std::endl;
        std::cout << NthEigenvector << "-th inv. eigenvector, norm_dark_inv = " << norm_dark_inv << std::endl;
        double norm_scale = norm_data/norm_dark;
        double norm_scale_inv = norm_data_inv/norm_dark_inv;
        if (abs(norm_scale)>2. || abs(norm_scale_inv)>2.) 
        {
            norm_scale = sign;
            norm_scale_inv = sign_inv;
        }
        if (abs(norm_scale)<0.5 || abs(norm_scale_inv)<0.5) 
        {
            norm_scale = sign;
            norm_scale_inv = sign_inv;
        }
        //if (dark_vector<2)
        //{
        //    norm_scale = 1.;
        //    norm_scale_inv = 1.;
        //}
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            //if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=2 && row<biny_blind)
            //if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=2)
            if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=dark_vector)
            {
                double real = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).real();
                double imag = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).imag();
                double ampl = pow(real*real+imag*imag,0.5);
                if (real!=0.)
                {
                    ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
                }
                mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = ampl;
            }
            else
            {
                double real = eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).real();
                double imag = eigensolver_data.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).imag();
                double ampl = pow(real*real+imag*imag,0.5);
                if (real!=0.)
                {
                    ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
                }
                mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = ampl;
            }
        }
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            //if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=2 && col<binx_blind)
            //if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=2)
            if (NthEigenvector<=NumberOfEigenvectors && NthEigenvector>=dark_vector)
            {
                double real = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).real();
                double imag = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).imag();
                double ampl = pow(real*real+imag*imag,0.5);
                if (real!=0.)
                {
                    ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
                }
                mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = ampl;
            }
            else
            {
                double real = eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).real();
                double imag = eigensolver_data.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).imag();
                double ampl = pow(real*real+imag*imag,0.5);
                if (real!=0.)
                {
                    ampl = real/abs(real)*pow(real*real+imag*imag,0.5);
                }
                mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = ampl;
            }
        }

    }

    //bool isComplexRoot = false;
    //for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    //{
    //    if (!isComplexRoot)
    //    {
    //        double real_this = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector).real();
    //        double imag_this = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector).imag();
    //        double real_next = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector-1).real();
    //        double imag_next = eigensolver_dark.eigenvalues()(mtx_dark.cols()-NthEigenvector-1).imag();
    //        std::cout << "At rank " << NthEigenvector << std::endl;
    //        std::cout << "real_this " << real_this << std::endl;
    //        std::cout << "imag_this " << imag_this << std::endl;
    //        std::cout << "real_next " << real_next << std::endl;
    //        std::cout << "imag_next " << imag_next << std::endl;
    //        if (real_next/real_this>=0.999 && abs(imag_next/imag_this)>=0.999) // complex conjugate
    //        {
    //            std::cout << "Find a complex conjugate at rank " << NthEigenvector << std::endl;
    //            isComplexRoot = true;
    //            if (NthEigenvector==NumberOfEigenvectors)
    //            {
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = 0.;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector-1,mtx_dark.cols()-NthEigenvector-1) = 0.;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector-1,mtx_dark.cols()-NthEigenvector) = 0.;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector-1) = 0.;
    //            }
    //            else
    //            {
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = real_this;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector-1,mtx_dark.cols()-NthEigenvector-1) = real_this;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector-1,mtx_dark.cols()-NthEigenvector) = -1.*imag_this;
    //                mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector-1) = imag_this;
    //                for (int row=0;row<N_bins_for_deconv;row++)
    //                {
    //                    mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector)(row) = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).real();
    //                    mtx_eigenvector_init.col(mtx_dark.cols()-NthEigenvector-1)(row) = eigensolver_dark.eigenvectors().col(mtx_dark.cols()-NthEigenvector)(row).imag();
    //                }
    //                for (int col=0;col<N_bins_for_deconv;col++)
    //                {
    //                    mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector)(col) = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).real();
    //                    mtx_eigenvector_inv_init.row(mtx_dark.cols()-NthEigenvector-1)(col) = eigensolver_dark.eigenvectors().inverse().row(mtx_dark.rows()-NthEigenvector)(col).imag();
    //                }
    //            }
    //        }
    //    }
    //    else
    //    {
    //        isComplexRoot = false;
    //    }
    //}
    

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
double NetflixChi2Function1D_0(const double *par)
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
double NetflixChi2Function1D_1(const double *par)
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
void NetflixSetInitialVariablesEigenBasis(ROOT::Math::GSLMinimizer* Chi2Minimizer, int binx_blind, int biny_blind, int fix_which, int which_to_fit)
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
    

    Chi2Minimizer->SetVariable(0, "par["+std::to_string(int(0))+"]", 1., 0.001);
    Chi2Minimizer->SetVariableLimits(0,1.,1.);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        limit = 0.0;
        first_index = 1+(2*NthEigenvector-2)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int entry=0;entry<NumberOfEigenvectors;entry++)
        {
            Chi2Minimizer->SetVariable(first_index+entry,"par["+std::to_string(int(first_index+entry))+"]",0.,0.001);
            if (fix_which==0) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
            if (NthEigenvector>=cutoff_mode) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
        }

        limit = 0.0;
        first_index = 1+(2*NthEigenvector-1)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int entry=0;entry<NumberOfEigenvectors;entry++)
        {
            Chi2Minimizer->SetVariable(first_index+entry,"par["+std::to_string(int(first_index+entry))+"]",0.,0.001);
            if (fix_which==1) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
            if (NthEigenvector>=cutoff_mode) Chi2Minimizer->SetVariableLimits(first_index+entry,0.,0.);
        }

        // eigenvalues
        double input_value = 0.;
        first_index = 1+(2*NthEigenvector-0)*(NumberOfEigenvectors)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            input_value = 0.;
            limit = 0.;
            Chi2Minimizer->SetVariable(first_index+NthEigenvalue-1, "par["+std::to_string(int(first_index+NthEigenvalue-1))+"]", input_value, 0.01*limit);
            Chi2Minimizer->SetVariableLimits(first_index+NthEigenvalue-1,input_value-limit,input_value+limit);
        }
    }

}
void NetflixSetInitialVariables(ROOT::Math::GSLMinimizer* Chi2Minimizer, int binx_blind, int biny_blind, int fix_which, int which_to_fit)
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
    

    Chi2Minimizer->SetVariable(0, "par["+std::to_string(int(0))+"]", 1., 0.001);
    Chi2Minimizer->SetVariableLimits(0,1.,1.);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        limit = 0.0;
        int col_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-2)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            Chi2Minimizer->SetVariable(first_index+row,"par["+std::to_string(int(first_index+row))+"]",0.,0.01);
            if (fix_which==0) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+row,0.,0.);
            //if (row>=biny_blind) 
            //{
            //    Chi2Minimizer->SetVariableLimits(first_index+row,0.0,0.0);
            //}
            //else
            //{
            //    double vari = mtx_eigenvector_vari(row,col_fix).real();
            //    Chi2Minimizer->SetVariableLimits(first_index+row,-vari,vari);
            //}
            //if (ratio_empty_bins>0.6) 
            //{
            //    Chi2Minimizer->SetVariableLimits(first_index+row,0.0,0.0);
            //}
        }

        limit = 0.0;
        int row_fix = N_bins_for_deconv-NthEigenvector;
        first_index = 1+(2*NthEigenvector-1)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int col=0;col<N_bins_for_deconv;col++)
        {
            Chi2Minimizer->SetVariable(first_index+col,"par["+std::to_string(int(first_index+col))+"]",0.,0.01);
            if (fix_which==1) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
            if (which_to_fit!=NthEigenvector && which_to_fit!=-1) Chi2Minimizer->SetVariableLimits(first_index+col,0.,0.);
            //if (col>=binx_blind) 
            //{
            //    Chi2Minimizer->SetVariableLimits(first_index+col,0.0,0.0);
            //}
            //else
            //{
            //    double vari = mtx_eigenvector_vari(row_fix,col).real();
            //    Chi2Minimizer->SetVariableLimits(first_index+col,-vari,vari);
            //}
            //if (ratio_empty_bins>0.6) 
            //{
            //    Chi2Minimizer->SetVariableLimits(first_index+col,0.0,0.0);
            //}
        }

        // eigenvalues
        double input_value = 0.;
        first_index = 1+(2*NthEigenvector-0)*(N_bins_for_deconv)+NumberOfEigenvectors*(NthEigenvector-1);
        for (int NthEigenvalue=1;NthEigenvalue<=NumberOfEigenvectors;NthEigenvalue++)
        {
            input_value = 0.;
            limit = 0.;
            //if (NthEigenvalue==NthEigenvector && NthEigenvalue==1)
            //{
            //    //limit = 0.;
            //    limit = 0.02*eigensolver_dark.eigenvalues()(N_bins_for_deconv-1).real();
            //}
            //else if (NthEigenvalue!=NthEigenvector)
            //{
            //    limit = 0.;
            //}
            //else
            //{
            //    limit = 0.02*eigensolver_dark.eigenvalues()(N_bins_for_deconv-1).real();
            //    if (mtx_eigenvalue_init(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvalue)==0.)
            //    {
            //        limit = 0.;
            //    }
            //}
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

double NuclearNormChi2Function(const double *par)
{
    
    // build model
    MatrixXcd mtx_model(N_bins_for_deconv,N_bins_for_deconv);
    TH2D hist_model_tmp = TH2D("hist_model_tmp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_gamma_tmp = TH2D("hist_gamma_tmp","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&hist_model_tmp,mtx_data);
    fill2DHistogram(&hist_gamma_tmp,mtx_gamma);

    VectorXcd vtr_ref_eigenvalues = vtr_dark_eigenvalues;

    hist_model_tmp.Add(&hist_gamma_tmp,par[0]);
    mtx_model = fillMatrix(&hist_model_tmp);

    eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_model);
    vtr_bkgd_eigenvalues = eigensolver_bkgd.eigenvalues();

    double chi2 = 0.;
    for (int entry=1;entry<=NumberOfEigenvectors;entry++)
    {
        chi2 += abs(vtr_ref_eigenvalues(N_bins_for_deconv-entry).real());
    }
    //for (int entry=1;entry<=NumberOfEigenvectors;entry++)
    //{
    //    double rel_error_real = 1.; // for small eigenvalues, the uncertainty is 1 order of magnitude.
    //    double rel_error_imag = 1.; // for small eigenvalues, the uncertainty is 1 order of magnitude.
    //    if (entry==1) rel_error_real = 0.005;
    //    if (entry==2) rel_error_real = 0.01;
    //    if (entry==3) rel_error_real = 0.02;
    //    double eigenvalue_ampl = pow(vtr_ref_eigenvalues(N_bins_for_deconv-entry).real()*vtr_ref_eigenvalues(N_bins_for_deconv-entry).real()+vtr_ref_eigenvalues(N_bins_for_deconv-entry).imag()*vtr_ref_eigenvalues(N_bins_for_deconv-entry).imag(),0.5);
    //    double error_real = rel_error_real*eigenvalue_ampl;
    //    double error_imag = rel_error_imag*eigenvalue_ampl;
    //    double chi2_real = vtr_bkgd_eigenvalues(N_bins_for_deconv-entry).real()-vtr_ref_eigenvalues(N_bins_for_deconv-entry).real();
    //    chi2_real = chi2_real*chi2_real/error_real;
    //    double chi2_imag = vtr_bkgd_eigenvalues(N_bins_for_deconv-entry).imag()-vtr_ref_eigenvalues(N_bins_for_deconv-entry).imag();
    //    chi2_imag = chi2_imag*chi2_imag/error_imag;
    //    chi2 += chi2_real+chi2_imag;
    //}

    return chi2;

}
void NuclearNormMinimizationMethod(int binx_blind, int biny_blind)
{

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_data_bkgd(row,col) = 0.;
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }

    eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    vtr_dark_eigenvalues = eigensolver_dark.eigenvalues();
    vtr_data_eigenvalues = eigensolver_data.eigenvalues();

    TH2D hist_model = TH2D("hist_model","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    TH2D hist_gamma = TH2D("hist_gamma","",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper);
    fill2DHistogram(&hist_model,mtx_data);
    fill2DHistogram(&hist_gamma,mtx_gamma);

    std::cout << "total n paramters = " << 1 << std::endl;
    ROOT::Math::Functor Chi2Func(&NuclearNormChi2Function,1); 

    // Choose method upon creation between:
    // kConjugateFR, kConjugatePR, kVectorBFGS,
    // kVectorBFGS2, kSteepestDescent

    //ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kSteepestDescent );
    ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kVectorBFGS2 );
    Chi2Minimizer_1st.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    Chi2Minimizer_1st.SetMaxIterations(100); // for GSL
    //Chi2Minimizer_1st.SetTolerance(0.001);
    Chi2Minimizer_1st.SetTolerance(0.1);
    Chi2Minimizer_1st.SetPrintLevel(0);
    //Chi2Minimizer_1st.SetPrintLevel(2);
    Chi2Minimizer_1st.SetFunction(Chi2Func);

    double initial = 0.;
    double step_size = 0.0001*hist_model.Integral()/hist_gamma.Integral();
    Chi2Minimizer_1st.SetVariable(0,"par["+std::to_string(int(0))+"]", initial, step_size);
    //Chi2Minimizer_1st.SetVariableLimits(0,initial-limit,initial);

    const double *par_1st = Chi2Minimizer_1st.X();
    std::cout << "initial chi2 = " << NuclearNormChi2Function(par_1st) << std::endl;
    Chi2Minimizer_1st.Minimize();
    par_1st = Chi2Minimizer_1st.X();
    std::cout << "final chi2 = " << NuclearNormChi2Function(par_1st) << std::endl;

    hist_model.Add(&hist_gamma,par_1st[0]);
    mtx_data_bkgd = fillMatrix(&hist_model);

    eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);
    for (int NthEigenvector=1;NthEigenvector<=NumberOfEigenvectors;NthEigenvector++)
    {

        if (NthEigenvector<=NumberOfEigenvectors)
        {
            mtx_eigenvalue(mtx_dark.cols()-NthEigenvector,mtx_dark.cols()-NthEigenvector) = eigensolver_bkgd.eigenvalues()(mtx_dark.cols()-NthEigenvector);
        }
    }

}
void SingleTimeMinimization(int fix_which, int which_to_fit)
{

    ROOT::Math::Functor Chi2Func;
    // use this is you use EigenBasis
    if (eigenbasis)
    {
        std::cout << "total n paramters = " << 1+2*NumberOfEigenvectors*(NumberOfEigenvectors)+NumberOfEigenvectors << std::endl;
        Chi2Func = ROOT::Math::Functor(&NetflixChi2Function,1+2*NumberOfEigenvectors*(NumberOfEigenvectors)+NumberOfEigenvectors*NumberOfEigenvectors); 
    }
    else
    {
        std::cout << "total n paramters = " << 1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors << std::endl;
        Chi2Func = ROOT::Math::Functor(&NetflixChi2Function,1+2*NumberOfEigenvectors*(N_bins_for_deconv)+NumberOfEigenvectors*NumberOfEigenvectors); 
    }

    // Choose method upon creation between:
    // kConjugateFR, kConjugatePR, kVectorBFGS,
    // kVectorBFGS2, kSteepestDescent

    //ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kSteepestDescent );
    ROOT::Math::GSLMinimizer Chi2Minimizer_1st( ROOT::Math::kVectorBFGS2 );
    Chi2Minimizer_1st.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    Chi2Minimizer_1st.SetMaxIterations(100); // for GSL
    //Chi2Minimizer_1st.SetTolerance(0.001);
    Chi2Minimizer_1st.SetTolerance(0.1);
    Chi2Minimizer_1st.SetPrintLevel(0);
    Chi2Minimizer_1st.SetFunction(Chi2Func);
    if (eigenbasis)
    {
        NetflixSetInitialVariablesEigenBasis(&Chi2Minimizer_1st,binx_blind_global,biny_blind_global,fix_which,which_to_fit);
    }
    else
    {
        NetflixSetInitialVariables(&Chi2Minimizer_1st,binx_blind_global,biny_blind_global,fix_which,which_to_fit);
    }
    const double *par_1st = Chi2Minimizer_1st.X();
    NthIteration = 0;
    std::cout << "initial chi2 = " << NetflixChi2Function(par_1st) << std::endl;
    Chi2Minimizer_1st.SetTolerance(0.01*double(N_bins_for_deconv*N_bins_for_deconv));
    Chi2Minimizer_1st.Minimize();
    par_1st = Chi2Minimizer_1st.X();
    std::cout << "final chi2 = " << NetflixChi2Function(par_1st) << std::endl;
    if (eigenbasis)
    {
        NetflixParametrizeEigenvectorsEigenBasis(par_1st);
    }
    else
    {
        NetflixParametrizeEigenvectors(par_1st);
    }

    mtx_eigenvalue = mtx_eigenvalue_init+mtx_eigenvalue_vari;
    mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
    mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
    mtx_data_bkgd = BuildModelMatrix();

}
MatrixXcd SpectralDecompositionMethod(MatrixXcd mtx_input, int LeftOrRight, int entry_cutoff, int nth_entry)
{
    // LeftOrRight = 1, left
    // LeftOrRight = 2, right

    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows(),mtx_input.cols());
    ComplexEigenSolver<MatrixXcd> eigensolver_input = ComplexEigenSolver<MatrixXcd>(mtx_input);
    ComplexEigenSolver<MatrixXcd> eigensolver_input_transpose = ComplexEigenSolver<MatrixXcd>(mtx_input.transpose());
    MatrixXcd mtx_r_init = eigensolver_input.eigenvectors();
    MatrixXcd mtx_l_init = eigensolver_input_transpose.eigenvectors();
    mtx_r_init = MakeRealEigenvectors(mtx_r_init);
    mtx_l_init = MakeRealEigenvectors(mtx_l_init);
    MatrixXcd mtx_rc_init = GetSubEigenvectors(mtx_r_init,1);
    MatrixXcd mtx_lc_init = GetSubEigenvectors(mtx_l_init,1);
    MatrixXcd mtx_R = mtx_r_init.transpose().conjugate()*mtx_r_init;
    MatrixXcd mtx_L = mtx_l_init.transpose().conjugate()*mtx_l_init;
    MatrixXcd mtx_H = mtx_l_init.transpose()*mtx_r_init;
    MatrixXcd mtx_lambda(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_lambda(row,col) = 0.;
            if (row==col)
            {
                mtx_lambda(row,col) = eigensolver_input.eigenvalues()(col);
            }
        }
    }
    mtx_lambda = CutoffEigenvalueMatrix(mtx_lambda, entry_cutoff);
    MatrixXcd mtx_S(mtx_input.rows(),mtx_input.cols());
    for (int row=0;row<mtx_input.rows();row++)
    {
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_S(row,col) = 0.;
            if (row==col)
            {
                mtx_S(row,col) = eigensolver_input.eigenvalues()(col)/mtx_H(row,col);
            }
        }
    }
    mtx_S = CutoffEigenvalueMatrix(mtx_S, entry_cutoff);

    MatrixXcd mtx_remain = mtx_input - mtx_r_init*mtx_S*mtx_l_init.transpose();

    MatrixXcd mtx_M1 = GetSubmatrix(mtx_data, 1);
    MatrixXcd mtx_M2 = GetSubmatrix(mtx_data, 2);
    MatrixXcd mtx_M3 = GetSubmatrix(mtx_data, 3);
    MatrixXcd mtx_M1_init = GetSubmatrix(mtx_input-mtx_remain, 1);
    MatrixXcd mtx_M2_init = GetSubmatrix(mtx_input-mtx_remain, 2);
    MatrixXcd mtx_M3_init = GetSubmatrix(mtx_input-mtx_remain, 3);

    VectorXcd vtr_Delta = VectorXcd::Zero((mtx_M2.cols()+mtx_M3.cols())*mtx_M3.rows());
    MatrixXcd mtx_Big = MatrixXcd::Zero((mtx_M2.cols()+mtx_M3.cols())*mtx_M3.rows(),mtx_input.cols());
    MatrixXcd mtx_unit = GetUnitMatrix(mtx_input.rows(),mtx_input.cols());
    MatrixXcd mtx_l_vari = MatrixXcd::Zero(mtx_l_init.rows(),mtx_l_init.cols());
    MatrixXcd mtx_r_vari = MatrixXcd::Zero(mtx_r_init.rows(),mtx_r_init.cols());
    MatrixXcd mtx_l_final = mtx_l_init;
    MatrixXcd mtx_r_final = mtx_r_init;
    // LeftOrRight = 1, left
    // LeftOrRight = 2, right
    if (LeftOrRight==1)
    {
        for (int row=0;row<mtx_M3.rows();row++)
        {
            int first_idx = row*(mtx_M2.cols()+mtx_M3.cols());
            vtr_Delta.segment(first_idx,mtx_M2.cols()) = (mtx_M2-mtx_M2_init).transpose().col(mtx_M2.rows()-row-1);
            vtr_Delta.segment(first_idx+mtx_M2.cols(),mtx_M3.cols()) = (mtx_M3-mtx_M3_init).transpose().col(mtx_M2.rows()-row-1);
        }
        for (int row=0;row<mtx_M3.rows();row++)
        {
            int first_idx = row*(mtx_input.cols());
            mtx_Big.block(first_idx,0,mtx_input.rows(),mtx_input.cols()) = mtx_unit*mtx_S(mtx_input.rows()-nth_entry,mtx_input.cols()-nth_entry)*mtx_rc_init(mtx_M3.rows()-row-1,mtx_input.cols()-nth_entry);
        }
        mtx_l_vari.col(mtx_l_vari.cols()-nth_entry) = mtx_Big.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_Delta);
        mtx_l_final += mtx_l_vari;
        //MatrixXcd mtx_S_l_transpose = GetEigenvalueWeightedVectors(mtx_S,mtx_l_final).transpose();
        //mtx_r_final.col(mtx_r_final.cols()-nth_entry) = ((mtx_R*mtx_S_l_transpose)*mtx_input.inverse()).transpose().conjugate().col(mtx_r_final.cols()-nth_entry);
        //std::cout << "mtx_r_init = " << std::endl;
        //PrintItReversely(mtx_r_init);
        //std::cout << "mtx_r_final = " << std::endl;
        //PrintItReversely(mtx_r_final);
    }
    else
    {
        for (int col=0;col<mtx_M3.cols();col++)
        {
            int first_idx = col*(mtx_M1.rows()+mtx_M3.rows());
            vtr_Delta.segment(first_idx,mtx_M1.rows()) = (mtx_M1-mtx_M1_init).col(mtx_M1.cols()-col-1);
            vtr_Delta.segment(first_idx+mtx_M1.rows(),mtx_M3.rows()) = (mtx_M3-mtx_M3_init).col(mtx_M1.cols()-col-1);
        }
        for (int col=0;col<mtx_M3.cols();col++)
        {
            int first_idx = col*(mtx_input.rows());
            mtx_Big.block(first_idx,0,mtx_input.cols(),mtx_input.rows()) = mtx_unit*mtx_S(mtx_input.rows()-nth_entry,mtx_input.cols()-nth_entry)*mtx_lc_init(mtx_M3.cols()-col-1,mtx_input.rows()-nth_entry);
        }
        mtx_r_vari.col(mtx_r_vari.cols()-nth_entry) = mtx_Big.bdcSvd(ComputeThinU | ComputeThinV).solve(vtr_Delta);
        mtx_r_final += mtx_r_vari;
        //MatrixXcd mtx_S_r_transpose = GetEigenvalueWeightedVectors(mtx_S,mtx_r_final).transpose();
        //mtx_l_final.col(mtx_l_final.cols()-nth_entry) = ((mtx_L*mtx_S_r_transpose)*(mtx_input.transpose()).inverse()).transpose().conjugate().col(mtx_l_final.cols()-nth_entry);
    }
    mtx_output = mtx_r_final*mtx_S*mtx_l_final.transpose() + mtx_remain;


    return mtx_output;
}
void MatrixFactorizationMethod()
{

    for (int col=0;col<N_bins_for_deconv;col++)
    {
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            mtx_data_bkgd(row,col) = 0.;
            mtx_eigenvector(row,col) = 0.;
            mtx_eigenvector_inv(row,col) = 0.;
            mtx_eigenvalue(row,col) = 0.;
        }
    }

    SetInitialEigenvectors(binx_blind_global,biny_blind_global);
    mtx_data_bkgd = mtx_eigenvector_init*mtx_eigenvalue_init*mtx_eigenvector_inv_init;

    //SingleTimeMinimization(-1,-1);
    for (int iteration=0;iteration<20;iteration++)
    {
        std::cout << "iteration = " << iteration << std::endl;
        SingleTimeMinimization(0,1);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
        SingleTimeMinimization(1,1);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
    }
    for (int iteration=0;iteration<20;iteration++)
    {
        std::cout << "iteration = " << iteration << std::endl;
        SingleTimeMinimization(0,2);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
        SingleTimeMinimization(1,2);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
    }
    for (int iteration=0;iteration<10;iteration++)
    {
        std::cout << "iteration = " << iteration << std::endl;
        SingleTimeMinimization(0,1);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
        SingleTimeMinimization(1,1);
        mtx_eigenvalue_init = mtx_eigenvalue;
        mtx_eigenvector_init = mtx_eigenvector;
        mtx_eigenvector_inv_init = mtx_eigenvector_inv;
    }

}

void PrintMatrixDiffRatio(MatrixXcd mtx_a, MatrixXcd mtx_b)
{
    MatrixXcd mtx_a_abs(mtx_a.rows(),mtx_a.cols());
    MatrixXcd mtx_b_abs(mtx_a.rows(),mtx_a.cols());
    MatrixXcd mtx_ratio(mtx_a.rows(),mtx_a.cols());
    IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
    for (int row=0;row<mtx_a.rows();row++)
    {
        for (int col=0;col<mtx_a.cols();col++)
        {
            mtx_a_abs(row,col) = pow(mtx_a(row,col).real()*mtx_a(row,col).real()+mtx_a(row,col).imag()*mtx_a(row,col).imag(),0.5);
            mtx_b_abs(row,col) = pow(mtx_b(row,col).real()*mtx_b(row,col).real()+mtx_b(row,col).imag()*mtx_b(row,col).imag(),0.5);
            mtx_ratio(row,col) = abs(mtx_a_abs(row,col)-mtx_b_abs(row,col))/mtx_a_abs(row,col);
        }
    }
    std::cout << "ref = " << std::endl;
    PrintItReversely(mtx_a_abs);
    std::cout << "diff ratio = " << std::endl;
    PrintItReversely(mtx_ratio);
}

void NetflixMethodPrediction(string target_data, double PercentCrab, double tel_elev_lower_input, double tel_elev_upper_input, bool isON, double MSCW_cut_input, double MSCL_cut_input, int rank)
{

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

    MSCW_plot_upper = gamma_hadron_dim_ratio*(MSCW_cut_input-MSCW_plot_lower)+MSCW_cut_input;
    MSCL_plot_upper = gamma_hadron_dim_ratio*(MSCL_cut_input-MSCL_plot_lower)+MSCL_cut_input;

    vector<TH2D> Hist_Redu_MSCLW;
    vector<TH2D> Hist_Rank0_MSCLW;
    vector<TH2D> Hist_Rank1_MSCLW;
    vector<TH2D> Hist_Rank2_MSCLW;
    vector<TH2D> Hist_Rank3_MSCLW;
    vector<TH2D> Hist_Spectral_Rank0_MSCLW;
    vector<TH2D> Hist_Spectral_Rank1_MSCLW;
    vector<TH2D> Hist_Spectral_Rank2_MSCLW;
    vector<TH2D> Hist_Spectral_Rank3_MSCLW;
    vector<TH2D> Hist_Bkgd_Rank0_MSCLW;
    vector<TH2D> Hist_Bkgd_Rank1_MSCLW;
    vector<TH2D> Hist_Bkgd_Rank2_MSCLW;
    vector<TH2D> Hist_Bkgd_Rank3_MSCLW;
    vector<TH2D> Hist_Bkgd_MSCLW;
    vector<TH2D> Hist_GammaRDBM_MSCLW;
    vector<TH1D> Hist_Data_Eigenvalues_real;
    vector<TH1D> Hist_Dark_Eigenvalues_real;
    vector<TH1D> Hist_Fit_Eigenvalues_real;
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

    TFile InputDataFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+file_tag+".root");

    double NSB;
    TTree* InfoTree = nullptr;
    InfoTree = (TTree*) InputDataFile.Get("InfoTree");
    InfoTree->SetBranchAddress("NSB",&NSB);
    InfoTree->GetEntry(0);

    TString filename_data_incl  = "Hist_Data_MSCLW_incl";
    TH2D* Hist_Data_incl = (TH2D*)InputDataFile.Get(filename_data_incl);
    MatrixXcd mtx_eigenval_data_redu(Hist_Data_incl->GetNbinsX(),Hist_Data_incl->GetNbinsY());
    mtx_data = fillMatrix(Hist_Data_incl);
    eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
    for (int r=0;r<8;r++)
    {
        char rank_idx[50];
        sprintf(rank_idx, "%i", r);
        Hist_Redu_MSCLW.push_back(TH2D("Hist_Redu_MSCLW_Rank"+TString(rank_idx),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j && i>=mtx_data.cols()-r-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else if (i==j && i==mtx_data.cols()-r-2)
                {
                    double real_this = eigensolver_data.eigenvalues()(i).real();
                    double imag_this = eigensolver_data.eigenvalues()(i).imag();
                    double real_next = eigensolver_data.eigenvalues()(i+1).real();
                    double imag_next = eigensolver_data.eigenvalues()(i+1).imag();
                    if (real_this/real_next>=0.999 && abs(imag_this/imag_next)>=0.999) // complex conjugate
                    {
                        mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                    }
                    else
                    {
                        mtx_eigenval_data_redu(i,j) = 0;
                    }
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        fill2DHistogramAbs(&Hist_Redu_MSCLW.at(r),mtx_data_redu);
    }

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
        TString filename_true  = "Hist_TrueBkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up);
        TH2D* Hist_GammaMC = (TH2D*)InputDataFile.Get(filename_gamma);
        TH2D* Hist_Data = (TH2D*)InputDataFile.Get(filename_data);
        TH2D* Hist_Dark = (TH2D*)InputDataFile.Get(filename_dark);
        TH2D* Hist_TrueBkgd = (TH2D*)InputDataFile.Get(filename_true);
        Hist_Rank0_MSCLW.push_back(TH2D("Hist_Rank0_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank1_MSCLW.push_back(TH2D("Hist_Rank1_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank2_MSCLW.push_back(TH2D("Hist_Rank2_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Rank3_MSCLW.push_back(TH2D("Hist_Rank3_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Spectral_Rank0_MSCLW.push_back(TH2D("Hist_Spectral_Rank0_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Spectral_Rank1_MSCLW.push_back(TH2D("Hist_Spectral_Rank1_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Spectral_Rank2_MSCLW.push_back(TH2D("Hist_Spectral_Rank2_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Spectral_Rank3_MSCLW.push_back(TH2D("Hist_Spectral_Rank3_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Bkgd_Rank0_MSCLW.push_back(TH2D("Hist_Bkgd_Rank0_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Bkgd_Rank1_MSCLW.push_back(TH2D("Hist_Bkgd_Rank1_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Bkgd_Rank2_MSCLW.push_back(TH2D("Hist_Bkgd_Rank2_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Bkgd_Rank3_MSCLW.push_back(TH2D("Hist_Bkgd_Rank3_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Bkgd_MSCLW.push_back(TH2D("Hist_Bkgd_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_GammaRDBM_MSCLW.push_back(TH2D("Hist_GammaRDBM_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,MSCL_plot_lower,MSCL_plot_upper,N_bins_for_deconv,MSCW_plot_lower,MSCW_plot_upper));
        Hist_Data_Eigenvalues_real.push_back(TH1D("Hist_Data_Eigenvalues_real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_Eigenvalues_real.push_back(TH1D("Hist_Dark_Eigenvalues_real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_Eigenvalues_real.push_back(TH1D("Hist_Fit_Eigenvalues_real_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_EigenvectorReal_0.push_back(TH1D("Hist_Data_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_EigenvectorReal_0.push_back(TH1D("Hist_Fit_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_EigenvectorReal_0.push_back(TH1D("Hist_Dark_EigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_EigenvectorReal_1.push_back(TH1D("Hist_Data_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_EigenvectorReal_1.push_back(TH1D("Hist_Fit_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_EigenvectorReal_1.push_back(TH1D("Hist_Dark_EigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_EigenvectorReal_2.push_back(TH1D("Hist_Data_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_EigenvectorReal_2.push_back(TH1D("Hist_Fit_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_EigenvectorReal_2.push_back(TH1D("Hist_Dark_EigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_InvEigenvectorReal_0.push_back(TH1D("Hist_Data_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_InvEigenvectorReal_0.push_back(TH1D("Hist_Fit_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_InvEigenvectorReal_0.push_back(TH1D("Hist_Dark_InvEigenvectorReal_0_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_InvEigenvectorReal_1.push_back(TH1D("Hist_Data_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_InvEigenvectorReal_1.push_back(TH1D("Hist_Fit_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_InvEigenvectorReal_1.push_back(TH1D("Hist_Dark_InvEigenvectorReal_1_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Data_InvEigenvectorReal_2.push_back(TH1D("Hist_Data_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Fit_InvEigenvectorReal_2.push_back(TH1D("Hist_Fit_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));
        Hist_Dark_InvEigenvectorReal_2.push_back(TH1D("Hist_Dark_InvEigenvectorReal_2_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",N_bins_for_deconv,0,N_bins_for_deconv));


        //Hist_Data->Scale(1./double(Hist_Data->Integral()));
        //Hist_Dark->Scale(1./double(Hist_Dark->Integral()));
        mtx_gamma = fillMatrix(Hist_GammaMC);
        mtx_data = fillMatrix(Hist_Data);
        mtx_dark = fillMatrix(Hist_Dark);
        eigensolver_data = ComplexEigenSolver<MatrixXcd>(mtx_data);
        eigensolver_data_transpose = ComplexEigenSolver<MatrixXcd>(mtx_data.transpose());
        eigensolver_dark = ComplexEigenSolver<MatrixXcd>(mtx_dark);
        eigensolver_dark_transpose = ComplexEigenSolver<MatrixXcd>(mtx_dark.transpose());

        int binx_lower = Hist_Data->GetXaxis()->FindBin(MSCL_cut_lower);
        binx_blind_global = Hist_Data->GetXaxis()->FindBin(MSCL_cut_blind)-1;
        int binx_upper = Hist_Data->GetXaxis()->FindBin(1.)-1;
        int biny_lower = Hist_Data->GetYaxis()->FindBin(MSCW_cut_lower);
        biny_blind_global = Hist_Data->GetYaxis()->FindBin(MSCW_cut_blind)-1;
        int biny_upper = Hist_Data->GetYaxis()->FindBin(1.)-1;

        n_fourier_modes = 6;
        n_taylor_modes = 6;

        NumberOfEigenvectors = rank;
        //if (CurrentEnergy>=1000.) NumberOfEigenvectors = 2;

        MatrixXcd mtx_U_r = eigensolver_data.eigenvectors();
        MatrixXcd mtx_U_l = eigensolver_data_transpose.eigenvectors();
        MatrixXcd mtx_H = mtx_U_l.transpose()*mtx_U_r;
        MatrixXcd mtx_R = mtx_U_r.transpose().conjugate()*mtx_U_r;
        MatrixXcd mtx_L = mtx_U_l.transpose().conjugate()*mtx_U_l;
        IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
        std::cout << "mtx_H = " << std::endl;
        PrintItReversely(mtx_H);
        std::cout << "mtx_R = " << std::endl;
        PrintItReversely(mtx_R);
        std::cout << "mtx_L = " << std::endl;
        PrintItReversely(mtx_L);
        MatrixXcd mtx_lambda(mtx_data.rows(),mtx_data.cols());
        MatrixXcd mtx_lambda_redu(mtx_data.rows(),mtx_data.cols());
        MatrixXcd mtx_S(N_bins_for_deconv,N_bins_for_deconv);
        MatrixXcd mtx_S_redu(N_bins_for_deconv,N_bins_for_deconv);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            for (int col=0;col<N_bins_for_deconv;col++)
            {
                mtx_S(row,col) = 0.;
                mtx_S_redu(row,col) = 0.;
                mtx_lambda_redu(row,col) = 0.;
                if (row==col)
                {
                    mtx_S(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                    mtx_lambda(row,col) = eigensolver_data.eigenvalues()(col);
                    if (row>=N_bins_for_deconv-1-3)
                    {
                        mtx_S_redu(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                        mtx_lambda_redu(row,col) = eigensolver_data.eigenvalues()(col);
                    }
                }
            }
        }
        std::cout << "mtx_S_redu = " << std::endl;
        PrintItReversely(mtx_S_redu);
        MatrixXcd mtx_data_spectral = mtx_U_r*mtx_S_redu*mtx_U_l.transpose();
        std::cout << "mtx_data - mtx_data_spectral = " << std::endl;
        std::cout << (mtx_data-mtx_data_spectral).format(CleanFmt) << std::endl;
        MatrixXcd mtx_rM = mtx_U_r.transpose().conjugate()*mtx_data;
        std::cout << "relative diff of mtx_rM - mtx_L.conjugate()*mtx_S.conjugate()*mtx_U_r.transpose() = " << std::endl;
        PrintMatrixDiffRatio(mtx_rM, mtx_R.conjugate()*mtx_S_redu.conjugate()*mtx_U_l.transpose());
        MatrixXcd mtx_lM = mtx_U_l.transpose().conjugate()*mtx_data.transpose();
        std::cout << "relative diff of mtx_lM - mtx_R.conjugate()*mtx_S.conjugate()*mtx_U_l.transpose() = " << std::endl;
        PrintMatrixDiffRatio(mtx_lM, mtx_L.conjugate()*mtx_S_redu.conjugate()*mtx_U_r.transpose());

        std::cout << "mtx_gamma = " << std::endl;
        std::cout << mtx_gamma << std::endl;
        std::cout << "mtx_data = " << std::endl;
        std::cout << mtx_data << std::endl;
        std::cout << "mtx_data (region 1) = " << std::endl;
        std::cout << GetSubmatrix(mtx_data, 1) << std::endl;
        std::cout << "mtx_data (region 2) = " << std::endl;
        std::cout << GetSubmatrix(mtx_data, 2) << std::endl;
        std::cout << "mtx_data (region 3) = " << std::endl;
        std::cout << GetSubmatrix(mtx_data, 3) << std::endl;


        MatrixXcd mtx_alpha_1 = GetAlphaMatrix(mtx_data, 1);
        std::cout << "mtx_alpha_1 = " << std::endl;
        std::cout << mtx_alpha_1 << std::endl;
        eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_alpha_1);
        std::cout << "alpha 1 eigenvalues = " << std::endl;
        std::cout << eigensolver_bkgd.eigenvalues() << std::endl;
        MatrixXcd mtx_proj_1 = GetProjectionMatrix(mtx_data, 1);
        std::cout << "mtx_proj_1 = " << std::endl;
        std::cout << mtx_proj_1 << std::endl;

        MatrixXcd mtx_alpha_2 = GetAlphaMatrix(mtx_data, 2);
        std::cout << "mtx_alpha_2 = " << std::endl;
        std::cout << mtx_alpha_2 << std::endl;
        eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_alpha_2);
        std::cout << "alpha 2 eigenvalues = " << std::endl;
        std::cout << eigensolver_bkgd.eigenvalues() << std::endl;
        MatrixXcd mtx_proj_2 = GetProjectionMatrix(mtx_data, 2);
        std::cout << "mtx_proj_2 = " << std::endl;
        std::cout << mtx_proj_2 << std::endl;

        MatrixXcd mtx_rg_predict = Convert_hadron_to_gamma_vector(mtx_U_r, mtx_lambda, 2);
        MatrixXcd mtx_lg_predict = Convert_hadron_to_gamma_vector(mtx_U_l, mtx_lambda, 1);
        //VerifyHMatrix(mtx_data);
        //VerifySmallgMatrix(mtx_data,1);
        //VerifyBetaMatrix(mtx_data);
        //VerifyLeftRightVectorCoupling(mtx_data,mtx_data.rows());
        //VerifyLeftRightVectorCoupling(mtx_data,3);
        //VerifyLeftRightVectorCoupling(mtx_data,1);
        //SolveEigenvectorsIteratively(mtx_data, mtx_dark, 3);
        //MatrixXcd mtx_predict = SolveEigenvectorsIteratively_v2(mtx_data,mtx_data,1);
        //mtx_predict = SolveEigenvectorsIteratively_v2(mtx_data,mtx_data,3);
        //mtx_predict = SolveEigenvectorsIteratively_v2(mtx_data,mtx_data,6);
        
        //MatrixXcd mtx_lambdanu = GetLambdaNuMatrix(mtx_data);
        //std::complex<double> par_R = mtx_lambdanu(mtx_data.rows()-1,mtx_data.cols()-1);
        //std::complex<double> par_L = mtx_lambdanu(mtx_data.rows()-1,mtx_data.cols()-1);
        ////par_R += 0.1*par_R;
        ////par_L += 0.1*par_L;
        //std::cout << "inital par_R = " << par_R << std::endl;
        //std::cout << "inital par_L = " << par_L << std::endl;
        //pair<std::complex<double>,std::complex<double>> pars = GetIterativeMatrix(mtx_data,1,par_R,par_L);
        //std::cout << "new par_R = " << pars.first << std::endl;
        //std::cout << "new par_L = " << pars.second << std::endl;
        //for (int iteration=1;iteration<=100;iteration++)
        //{
        //    pars = GetIterativeMatrix(mtx_data,1,pars.first,pars.second);
        //    std::cout << "new par_R = " << pars.first << std::endl;
        //    std::cout << "new par_L = " << pars.second << std::endl;
        //}

        double count_gamma_like = Hist_Data->Integral(binx_lower,binx_blind_global,biny_lower,biny_blind_global);
        double count_total = Hist_Data->Integral();
        //if (file_tag == "OFF")
        //{
        //    char textfile[50] = "";
        //    sprintf(textfile, "../eigenvector_%s_E%i_Elev%i.txt", target, int(energy_bins[e]), int(TelElev_lower));
        //    std::ofstream myfile;
        //    myfile.open (textfile);
        //    sprintf(textfile, "../eigenvector_%s_E%i_Elev%i_coma.txt", target, int(energy_bins[e]), int(TelElev_lower));
        //    std::ofstream myfile2;
        //    myfile2.open (textfile);
        //    myfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        //    myfile << "NSB:" << NSB << std::endl;
        //    myfile << "total counts:" << count_total << std::endl;
        //    myfile << "gamma-like counts:" << count_gamma_like << std::endl;
        //    myfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        //    myfile << "eigenvalues:" << std::endl;
        //    myfile << eigensolver_data.eigenvalues() << std::endl << std::endl;
        //    bool is_complex_pair = false;
        //    for (int r=0;r<N_bins_for_deconv;r++)
        //    {
        //        is_complex_pair = false;
        //        MatrixXcd mtx_eigenval_data_temp(Hist_Data_incl->GetNbinsX(),Hist_Data_incl->GetNbinsY());
        //        for (int i=0;i<mtx_data.cols();i++)
        //        {
        //            for (int j=0;j<mtx_data.rows();j++)
        //            {
        //                if (i==j && i==mtx_data.cols()-r-1) 
        //                {
        //                    mtx_eigenval_data_temp(i,j) = eigensolver_data.eigenvalues()(i);
        //                }
        //                //else if (i==j && i==mtx_data.cols()-r-2)
        //                //{
        //                //    double real_this = eigensolver_data.eigenvalues()(i).real();
        //                //    double imag_this = eigensolver_data.eigenvalues()(i).imag();
        //                //    double real_next = eigensolver_data.eigenvalues()(i+1).real();
        //                //    double imag_next = eigensolver_data.eigenvalues()(i+1).imag();
        //                //    if (real_this/real_next>=0.999 && abs(imag_this/imag_next)>=0.999) // complex conjugate
        //                //    {
        //                //        mtx_eigenval_data_temp(i,j) = eigensolver_data.eigenvalues()(i);
        //                //        is_complex_pair = true;
        //                //    }
        //                //    else
        //                //    {
        //                //        mtx_eigenval_data_temp(i,j) = 0;
        //                //    }
        //                //}
        //                else
        //                {
        //                    mtx_eigenval_data_temp(i,j) = 0;
        //                }
        //            }
        //        }
        //        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_temp*eigensolver_data.eigenvectors().inverse();
        //        myfile << r << "th largest eigenvalue, contribution in total and in gamma-ray region:" << std::endl;
        //        myfile << "eigenvalue: " << eigensolver_data.eigenvalues()(mtx_data.cols()-r-1) << std::endl;
        //        myfile << "total region: " << MatrixIntegral(mtx_data_redu,0,N_bins_for_deconv-1,0,N_bins_for_deconv-1) << std::endl;
        //        myfile << "signal region: " << MatrixIntegral(mtx_data_redu,binx_lower-1,binx_blind-1,biny_lower-1,biny_blind-1) << std::endl;
        //        myfile2 << eigensolver_data.eigenvalues()(mtx_data.cols()-r-1).real() << "," << eigensolver_data.eigenvalues()(mtx_data.cols()-r-1).imag() << "," << MatrixIntegral(mtx_data_redu,0,N_bins_for_deconv-1,0,N_bins_for_deconv-1).real() << "," << MatrixIntegral(mtx_data_redu,0,N_bins_for_deconv-1,0,N_bins_for_deconv-1).imag() << "," << MatrixIntegral(mtx_data_redu,binx_lower-1,binx_blind-1,biny_lower-1,biny_blind-1).real() << "," << MatrixIntegral(mtx_data_redu,binx_lower-1,binx_blind-1,biny_lower-1,biny_blind-1).imag() << std::endl;
        //    }
        //  //myfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        //  //myfile << "eigenvectors:" << std::endl;
        //  //myfile << eigensolver_data.eigenvectors().block(0,mtx_data.cols()-6,mtx_data.rows(),6) << std::endl << std::endl;
        //  //myfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        //  //myfile << "inverse eigenvectors:" << std::endl;
        //  //myfile << eigensolver_data.eigenvectors().inverse().block(mtx_data.rows()-6,0,6,mtx_data.cols()) << std::endl << std::endl;
        //  //myfile << "U^{-1}MU: "<< std::endl;
        //  //myfile << eigensolver_data.eigenvectors().inverse()*mtx_data*eigensolver_data.eigenvectors() << std::endl << std::endl;
        //  myfile.close();
        //  myfile2.close();
        //}

        int binx_cut = Hist_Data->GetXaxis()->FindBin(1.)-1;
        int biny_cut = Hist_Data->GetYaxis()->FindBin(1.)-1;
        ratio_empty_bins = RatioEmptyBins(mtx_data,binx_cut,biny_cut);
        std::cout << "ratio_empty_bins = " << ratio_empty_bins << std::endl;

        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j && i==mtx_data.cols()-0-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        fill2DHistogram(&Hist_Rank0_MSCLW.at(e),mtx_data_redu);
        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j && i==mtx_data.cols()-1-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        fill2DHistogram(&Hist_Rank1_MSCLW.at(e),mtx_data_redu);
        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j && i==mtx_data.cols()-2-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        fill2DHistogram(&Hist_Rank2_MSCLW.at(e),mtx_data_redu);
        for (int i=0;i<mtx_data.cols();i++)
        {
            for (int j=0;j<mtx_data.rows();j++)
            {
                if (i==j && i==mtx_data.cols()-3-1) 
                {
                    mtx_eigenval_data_redu(i,j) = eigensolver_data.eigenvalues()(i);
                }
                else
                {
                    mtx_eigenval_data_redu(i,j) = 0;
                }
            }
        }
        mtx_data_redu = eigensolver_data.eigenvectors()*mtx_eigenval_data_redu*eigensolver_data.eigenvectors().inverse();
        fill2DHistogram(&Hist_Rank3_MSCLW.at(e),mtx_data_redu);

        for (int row=0;row<N_bins_for_deconv;row++)
        {
            for (int col=0;col<N_bins_for_deconv;col++)
            {
                mtx_S_redu(row,col) = 0.;
                if (row==col)
                {
                    if (row==N_bins_for_deconv-1-0)
                    {
                        mtx_S_redu(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                    }
                }
            }
        }
        mtx_data_spectral = mtx_U_r*mtx_S_redu*mtx_U_l.transpose();
        fill2DHistogram(&Hist_Spectral_Rank0_MSCLW.at(e),mtx_data_spectral);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            for (int col=0;col<N_bins_for_deconv;col++)
            {
                mtx_S_redu(row,col) = 0.;
                if (row==col)
                {
                    if (row==N_bins_for_deconv-1-1)
                    {
                        mtx_S_redu(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                    }
                }
            }
        }
        mtx_data_spectral = mtx_U_r*mtx_S_redu*mtx_U_l.transpose();
        fill2DHistogram(&Hist_Spectral_Rank1_MSCLW.at(e),mtx_data_spectral);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            for (int col=0;col<N_bins_for_deconv;col++)
            {
                mtx_S_redu(row,col) = 0.;
                if (row==col)
                {
                    if (row==N_bins_for_deconv-1-2)
                    {
                        mtx_S_redu(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                    }
                }
            }
        }
        mtx_data_spectral = mtx_U_r*mtx_S_redu*mtx_U_l.transpose();
        fill2DHistogram(&Hist_Spectral_Rank2_MSCLW.at(e),mtx_data_spectral);
        for (int row=0;row<N_bins_for_deconv;row++)
        {
            for (int col=0;col<N_bins_for_deconv;col++)
            {
                mtx_S_redu(row,col) = 0.;
                if (row==col)
                {
                    if (row==N_bins_for_deconv-1-3)
                    {
                        mtx_S_redu(row,col) = eigensolver_data.eigenvalues()(col)/mtx_H(row,col);
                    }
                }
            }
        }
        mtx_data_spectral = mtx_U_r*mtx_S_redu*mtx_U_l.transpose();
        fill2DHistogram(&Hist_Spectral_Rank3_MSCLW.at(e),mtx_data_spectral);

        for (int NthEigenvalue=1;NthEigenvalue<=N_bins_for_deconv;NthEigenvalue++)
        {
            Hist_Data_Eigenvalues_real.at(e).SetBinContent(NthEigenvalue,eigensolver_data.eigenvalues()(N_bins_for_deconv-NthEigenvalue).real());
            Hist_Dark_Eigenvalues_real.at(e).SetBinContent(NthEigenvalue,eigensolver_dark.eigenvalues()(N_bins_for_deconv-NthEigenvalue).real());
        }

        SetInitialEigenvectors(binx_blind_global,biny_blind_global);
        fill1DHistogram(&Hist_Data_EigenvectorReal_0.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Data_EigenvectorReal_1.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Data_EigenvectorReal_2.at(e),eigensolver_data.eigenvectors().col(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_0.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_1.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Data_InvEigenvectorReal_2.at(e),eigensolver_data.eigenvectors().inverse().row(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Dark_EigenvectorReal_0.at(e),mtx_eigenvector_init.col(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Dark_EigenvectorReal_1.at(e),mtx_eigenvector_init.col(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Dark_EigenvectorReal_2.at(e),mtx_eigenvector_init.col(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_0.at(e),mtx_eigenvector_inv_init.row(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_1.at(e),mtx_eigenvector_inv_init.row(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Dark_InvEigenvectorReal_2.at(e),mtx_eigenvector_inv_init.row(mtx_data.cols()-3));

        int col_fix = 0;
        int row_fix = 0;
        int first_index = 0;

        mtx_data_bkgd = SpectralDecompositionMethod(mtx_dark, 1, 3, 1);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 2, 3, 1);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 1, 3, 2);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 2, 3, 2);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 1, 3, 1);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 2, 3, 1);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 1, 3, 2);
        mtx_data_bkgd = SpectralDecompositionMethod(mtx_data_bkgd, 2, 3, 2);

        //MatrixFactorizationMethod();
        //NuclearNormMinimizationMethod(binx_blind_global,biny_blind_global);
        fill2DHistogramAbs(&Hist_Bkgd_MSCLW.at(e),mtx_data_bkgd);
        std::cout << "Hist_Bkgd_MSCLW.at(e).Integral() = " << Hist_Bkgd_MSCLW.at(e).Integral(1,15,1,15) << std::endl;
        eigensolver_bkgd = ComplexEigenSolver<MatrixXcd>(mtx_data_bkgd);

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenval_data_redu(row,col) = 0.;
            }
        }
        mtx_eigenval_data_redu(N_bins_for_deconv-1,N_bins_for_deconv-1) = mtx_eigenvalue_init(N_bins_for_deconv-1,N_bins_for_deconv-1);
        mtx_eigenval_data_redu(N_bins_for_deconv-1,N_bins_for_deconv-1) += mtx_eigenvalue_vari(N_bins_for_deconv-1,N_bins_for_deconv-1);
        if (linear)
        {
            mtx_data_redu = mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_vari*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_vari;
        }
        else
        {
            mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
            mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
            mtx_data_redu = mtx_eigenvector*mtx_eigenval_data_redu*mtx_eigenvector_inv;
        }
        fill2DHistogramAbs(&Hist_Bkgd_Rank0_MSCLW.at(e),mtx_data_redu);

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenval_data_redu(row,col) = 0.;
            }
        }
        mtx_eigenval_data_redu(N_bins_for_deconv-2,N_bins_for_deconv-2) = mtx_eigenvalue_init(N_bins_for_deconv-2,N_bins_for_deconv-2);
        mtx_eigenval_data_redu(N_bins_for_deconv-2,N_bins_for_deconv-2) += mtx_eigenvalue_vari(N_bins_for_deconv-2,N_bins_for_deconv-2);
        if (linear)
        {
            mtx_data_redu = mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_vari*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_vari;
        }
        else
        {
            mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
            mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
            mtx_data_redu = mtx_eigenvector*mtx_eigenval_data_redu*mtx_eigenvector_inv;
        }
        fill2DHistogramAbs(&Hist_Bkgd_Rank1_MSCLW.at(e),mtx_data_redu);

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenval_data_redu(row,col) = 0.;
            }
        }
        mtx_eigenval_data_redu(N_bins_for_deconv-3,N_bins_for_deconv-3) = mtx_eigenvalue_init(N_bins_for_deconv-3,N_bins_for_deconv-3);
        mtx_eigenval_data_redu(N_bins_for_deconv-3,N_bins_for_deconv-3) += mtx_eigenvalue_vari(N_bins_for_deconv-3,N_bins_for_deconv-3);
        if (linear)
        {
            mtx_data_redu = mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_vari*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_vari;
        }
        else
        {
            mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
            mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
            mtx_data_redu = mtx_eigenvector*mtx_eigenval_data_redu*mtx_eigenvector_inv;
        }
        fill2DHistogramAbs(&Hist_Bkgd_Rank2_MSCLW.at(e),mtx_data_redu);

        for (int col=0;col<N_bins_for_deconv;col++)
        {
            for (int row=0;row<N_bins_for_deconv;row++)
            {
                mtx_eigenval_data_redu(row,col) = 0.;
            }
        }
        mtx_eigenval_data_redu(N_bins_for_deconv-4,N_bins_for_deconv-4) = mtx_eigenvalue_init(N_bins_for_deconv-4,N_bins_for_deconv-4);
        mtx_eigenval_data_redu(N_bins_for_deconv-4,N_bins_for_deconv-4) += mtx_eigenvalue_vari(N_bins_for_deconv-4,N_bins_for_deconv-4);
        if (linear)
        {
            mtx_data_redu = mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_vari*mtx_eigenval_data_redu*mtx_eigenvector_inv_init;
            mtx_data_redu += mtx_eigenvector_init*mtx_eigenval_data_redu*mtx_eigenvector_inv_vari;
        }
        else
        {
            mtx_eigenvector = mtx_eigenvector_init+mtx_eigenvector_vari;
            mtx_eigenvector_inv = mtx_eigenvector_inv_init+mtx_eigenvector_inv_vari;
            mtx_data_redu = mtx_eigenvector*mtx_eigenval_data_redu*mtx_eigenvector_inv;
        }
        fill2DHistogramAbs(&Hist_Bkgd_Rank3_MSCLW.at(e),mtx_data_redu);

        for (int NthEigenvalue=1;NthEigenvalue<=N_bins_for_deconv;NthEigenvalue++)
        {
            Hist_Fit_Eigenvalues_real.at(e).SetBinContent(NthEigenvalue,mtx_eigenvalue_init(N_bins_for_deconv-NthEigenvalue,N_bins_for_deconv-NthEigenvalue).real()+NthEigenvalue,mtx_eigenvalue_vari(N_bins_for_deconv-NthEigenvalue,N_bins_for_deconv-NthEigenvalue).real());
        }

        fill1DHistogram(&Hist_Fit_EigenvectorReal_0.at(e),mtx_eigenvector.col(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Fit_EigenvectorReal_1.at(e),mtx_eigenvector.col(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Fit_EigenvectorReal_2.at(e),mtx_eigenvector.col(mtx_data.cols()-3));
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_0.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-1));
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_1.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-2));
        fill1DHistogram(&Hist_Fit_InvEigenvectorReal_2.at(e),mtx_eigenvector_inv.row(mtx_data.cols()-3));
        fill2DHistogramAbs(&Hist_GammaRDBM_MSCLW.at(e),mtx_gamma);

        double gamma_total = Hist_Data->Integral(binx_lower,binx_upper,biny_lower,biny_upper)-Hist_Bkgd_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper);
        gamma_total = max(0.,gamma_total);
        double scale_gamma = double(gamma_total)/double(Hist_GammaRDBM_MSCLW.at(e).Integral(binx_lower,binx_upper,biny_lower,biny_upper));
        Hist_GammaRDBM_MSCLW.at(e).Scale(scale_gamma);

        std::cout << "Hist_TrueBkgd->Integral() = " << Hist_TrueBkgd->Integral(1,15,1,15) << std::endl;
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
        std::cout << "mtx_eigenvalue_init:" << std::endl;
        std::cout << mtx_eigenvalue_init.block(mtx_data.rows()-NumberOfEigenvectors,mtx_data.cols()-NumberOfEigenvectors,NumberOfEigenvectors,NumberOfEigenvectors) << std::endl;
        std::cout << "mtx_eigenvalue:" << std::endl;
        std::cout << mtx_eigenvalue.block(mtx_data.rows()-NumberOfEigenvectors,mtx_data.cols()-NumberOfEigenvectors,NumberOfEigenvectors,NumberOfEigenvectors) << std::endl;
        std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1):" << std::endl;
        std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-1) << std::endl;
        std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2):" << std::endl;
        std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-2) << std::endl;
        std::cout << "eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3):" << std::endl;
        std::cout << eigensolver_bkgd.eigenvalues()(mtx_data.cols()-3) << std::endl;

    }
    InputDataFile.Close();


    TFile OutputFile("../Netflix_"+TString(target_data)+"_Crab"+std::to_string(int(PercentCrab))+"_TelElev"+std::to_string(int(TelElev_lower))+"to"+std::to_string(int(TelElev_upper))+"_"+file_tag+".root","update");
    for (int r=0;r<8;r++)
    {
        Hist_Redu_MSCLW.at(r).Write();
    }
    for (int e=0;e<N_energy_bins;e++)
    {
        Hist_GammaRDBM_MSCLW.at(e).Write();
        Hist_Rank0_MSCLW.at(e).Write();
        Hist_Rank1_MSCLW.at(e).Write();
        Hist_Rank2_MSCLW.at(e).Write();
        Hist_Rank3_MSCLW.at(e).Write();
        Hist_Spectral_Rank0_MSCLW.at(e).Write();
        Hist_Spectral_Rank1_MSCLW.at(e).Write();
        Hist_Spectral_Rank2_MSCLW.at(e).Write();
        Hist_Spectral_Rank3_MSCLW.at(e).Write();
        Hist_Bkgd_Rank0_MSCLW.at(e).Write();
        Hist_Bkgd_Rank1_MSCLW.at(e).Write();
        Hist_Bkgd_Rank2_MSCLW.at(e).Write();
        Hist_Bkgd_Rank3_MSCLW.at(e).Write();
        Hist_Bkgd_MSCLW.at(e).Write();
        Hist_Data_Eigenvalues_real.at(e).Write();
        Hist_Dark_Eigenvalues_real.at(e).Write();
        Hist_Fit_Eigenvalues_real.at(e).Write();
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
    }
    OutputFile.Close();

}
