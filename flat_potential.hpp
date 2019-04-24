#ifndef _FLAT_POTENTIAL_HPP
#define _FLAT_POTENTIAL_HPP

#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "boost/math/special_functions/erf.hpp"

/*
 ****************************************************************
 * 2D FLAT POTENTIAL 
 * 
 *  H = [ -cos(theta), sin(theta) * exp(i * phi) ]
 *      [  sin(theta) * exp(i * phi), cos(theta) ]
 *
 *  theta = pi / 2 * (erf(B * x) + 1)
 *  phi = W * y
 *
 ****************************************************************
 *  API INFO
 ****************************************************************
 *
 *  - void output_potenial_params() 
 *      
 *      THIS FUNCTION OUTPUTS POTENTIAL PARAMETRS
 *
 *  - vector< complex<double> > cal_H(const double x)
 *
 *      THIS FUNCTION CALCULATES HAMILTONIAN FOR A GIVEN POSITION
 *
 *  - vector< complex<double> > cal_nablaH(const double x)
 *
 *      THIS FUNCTION CALCULATES HAMILTONIAN GRADIANT FOR A GIVEN POSITION
 *
 *  - void cal_info_nume(const double x, vector<double>& eva, vector< complex<double> >& dc, vector< complex<double> >& F, vector< complex<double> >& lastevt)
 *
 *      THIS FUNCTION CALCULATES EIGENVALUES, EIGENVECTORS, FORCE MATRIX AND DERIVATIVE COUPLINGS FOR A GIVEN POSITION
 *  
 ****************************************************************
 */

namespace {
    using std::vector;
    using std::complex;

    double param_A = 0.02;
    double param_B = 3.0;
    double param_W = 0.5;

    void output_potential_params() {
        ioer::info("# Potential parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " W = ", param_W);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 3, 
                "set_potenial_params: potential paramter vector size must be >= 3");
        param_A = params[0];
        param_B = params[1];
        param_W = params[2];
    }

    double cal_theta(const vector<double>& r) {
        const double x = r[0];
        return 0.5 * M_PI * (boost::math::erf(param_B * x) + 1.0);
    }

    vector<double> cal_nablatheta(const vector<double>& r) {
        const double x = r[0];
        vector<double> nablatheta(r.size(), 0.0);
        nablatheta[0] = std::sqrt(M_PI) * param_B * std::exp(-param_B * param_B * x * x);
        return nablatheta;
    }

    double cal_phi(const vector<double>& r) {
        const double y = r[1];
        return param_W * y;
    }

    vector<double> cal_nablaphi(const vector<double>& r) {
        vector<double> nablaphi(r.size(), 0.0);
        nablaphi[1] = param_W;
        return nablaphi;
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const double theta = cal_theta(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const double CC = cos(theta);
        const double SS = sin(theta);

        vector< complex<double> > H(4, 0.0);
        H[0+0*2] = -CC;
        H[1+1*2] = CC;
        H[0+1*2] = SS * eip;
        H[1+0*2] = conj(H[0+1*2]);
        H *= param_A;

        return H;
    }

    vector< vector< complex<double> > > cal_nablaH(const vector<double>& r) {
        const int ndim = r.size();
        const int edim = 2;

        const double theta = cal_theta(r);
        const double phi = cal_phi(r);
        const vector<double> nablatheta = cal_nablatheta(r);
        const vector<double> nablaphi = cal_nablaphi(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * phi);
        const double CC = cos(theta);
        const double SS = sin(theta);

        vector< vector< complex<double> > > nablaH(ndim);

        for (int ix(0); ix < ndim; ++ix) {
            vector< complex<double> >& nablaH_ix = nablaH[ix];
            nablaH_ix.resize(edim * edim);

            nablaH_ix[0+0*edim] = SS * nablatheta[ix];
            nablaH_ix[1+1*edim] = -SS * nablatheta[ix];
            nablaH_ix[0+1*edim] = eip * (CC * nablatheta[ix] + matrixop::IMAGIZ * SS * nablaphi[ix]);
            nablaH_ix[1+0*edim] = conj(nablaH_ix[0+1*edim]);

            nablaH_ix *= param_A;
        }

        return nablaH;
    }

    void cal_info_nume(const vector<double>& r,
            vector<double>& eva, 
            vector< vector< complex<double> > >& dc,
            vector< vector< complex<double> > >& F,
            vector< complex<double> >& lastevt)
    {
        /*
         * input:   position vector r
         *
         * output:  eigenvalues
         *          (phase aligned) derivative coupling matrix
         *          force matrix
         *
         * in/out:  last step eigenvectors lastevt (can be a empty vector)
         *              on exit, lastevt is replaced by current step evt
         *
         * -- calculation is performed in the numerical way
         */

        const vector< complex<double> > H = cal_H(r);
        const int ndim = r.size();
        const int edim = static_cast<int>(std::sqrt(static_cast<double>(H.size())));
        vector< complex<double> > evt;
        matrixop::hdiag(cal_H(r), eva, evt);

        // correct phase
        if (not lastevt.empty()) {
            auto tmp = matrixop::matCmat(lastevt, evt, edim);
            for (int j = 0; j < edim; ++j) {
                complex<double> eip = tmp[j+j*edim] / abs(tmp[j+j*edim]);
                for (int k = 0; k < edim; ++k) {
                    evt[k+j*edim] /= eip;
                }
            }
        }

        // F, dc
        const vector< vector< complex<double> > > nablaH = cal_nablaH(r);
        dc.resize(ndim);
        F.resize(ndim);

        for (int ix = 0; ix < ndim; ++ix) {
            dc[ix] = matrixop::matCmatmat(evt, nablaH[ix], evt, edim, edim);
            F[ix].assign(edim * edim, 0.0);

            for (int j = 0; j < edim; ++j) {
                for (int k = 0; k < edim; ++k) {
                    F[ix][j+k*edim] = -dc[ix][j+k*edim];
                    if (j == k) {
                        dc[ix][j+k*edim] = 0.0;
                    }
                    else {
                        dc[ix][j+k*edim] /= (eva[k] - eva[j]);
                    }
                }
            }
        }

        // save evt to lastevt
        lastevt = std::move(evt);
    }
};

#endif
