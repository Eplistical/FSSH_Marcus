#ifndef _SPIN_BOSON_2D_POTENTIAL_HPP
#define _SPIN_BOSON_2D_POTENTIAL_HPP


#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/ioer.hpp"


/*
 ****************************************************************
 * 2D SPIN BOSON POTENTIAL 
 *
 *  H00 = 0.5 * mx * wx**2 * x**2 + 0.5 * my * wy**2 * y**2
 *  H11 = 0.5 * mx * wx**2 * (x-gx)**2 + 0.5 * my * wy**2 * (g-gy)**2 + dG0
 *  H01 = V * exp(i * W * y)
 *  H10 = conj(H01)
 *
 ****************************************************************
 *  API INFO
 ****************************************************************
 *
 *  - void output_potenial_params() 
 *      
 *      THIS FUNCTION OUTPUTS POTENTIAL PARAMETRS
 *
 *  - vector< complex<double> > cal_H(const vector<double>& r)
 *
 *      THIS FUNCTION CALCULATES HAMILTONIAN FOR A GIVEN POSITION
 *
 *  - vector< vector< complex<double> > > cal_nablaH(const vector<double>& r)
 *
 *      THIS FUNCTION CALCULATES HAMILTONIAN GRADIANT FOR A GIVEN POSITION
 *
 *  - void cal_info_nume(const vector<double>& r, vector<double>& eva, vector< vector< complex<double> > >& dc, vector< vector< complex<double> > >& F, vector< complex<double> >& lastevt)
 *
 *      THIS FUNCTION CALCULATES EIGENVALUES, EIGENVECTORS, FORCE MATRIX AND DERIVATIVE COUPLINGS FOR A GIVEN POSITION
 *  
 ****************************************************************
 */

namespace {
    using std::vector;
    using std::complex;

    vector<double> param_mass { 1.0, 1.0 };
    vector<double> param_omega { 4.375e-5, 4.375e-5 };
    vector<double> param_g { 4997.3, 2524.7 };
    double param_dG0 = -0.018;
    double param_V = 2.5e-5;
    double param_W = 0.0;

    void output_potential_params() {
        ioer::info("# 2D Spin Boson Potential parameters: ", 
                    " mass = ", param_mass,
                    " omega = ", param_omega,
                    " g = ", param_g,
                    " dG0 = ", param_dG0,
                    " V = ", param_V,
                    " W = ", param_W,
                    ""
                    );
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 9, 
                "set_potenial_params: potential paramter vector size must be >= 9");
        param_mass[0] = params[0];
        param_mass[1] = params[1];
        param_omega[0] = params[2];
        param_omega[1] = params[3];
        param_g[0] = params[4];
        param_g[1] = params[5];
        param_dG0 = params[6];
        param_V = params[7];
        param_W = params[8];
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const double mx = param_mass[0];
        const double my = param_mass[1];
        const double wx = param_omega[0];
        const double wy = param_omega[1];
        const double gx = param_g[0];
        const double gy = param_g[1];
        const double dG0 = param_dG0;
        const double V = param_V;
        const double W = param_W;

        vector< complex<double> > H(4, 0.0);
        H[0+0*2] = 0.5 * mx * pow(wx * x, 2) + 0.5 * my * pow(wy * y, 2);
        H[1+1*2] = 0.5 * mx * pow(wx * (x - gx), 2) + 0.5 * my * pow(wy * (y - gy), 2) + dG0;
        H[0+1*2] = V * exp(matrixop::IMAGIZ * W * y);
        H[1+0*2] = conj(H[0+1*2]);

        return H;
    }

    vector< vector< complex<double> > > cal_nablaH(const vector<double>& r) {
        const int ndim = 2;

        const double x = r[0];
        const double y = r[1];
        const double mx = param_mass[0];
        const double my = param_mass[1];
        const double wx = param_omega[0];
        const double wy = param_omega[1];
        const double gx = param_g[0];
        const double gy = param_g[1];
        const double dG0 = param_dG0;
        const double V = param_V;
        const double W = param_W;

        // initialize nablaH
        vector< vector< complex<double> > > nablaH(ndim);
        for (int ix(0); ix < ndim; ++ix) {
            nablaH[ix].resize(4, 0.0);
        }

        vector< complex<double> >& Hx = nablaH[0];
        vector< complex<double> >& Hy = nablaH[1];

        Hx.resize(4);
        Hx[0+0*2] = mx * wx * wx * x;
        Hx[1+1*2] = mx * wx * wx * (x - gx);
        Hx[0+1*2] = 0.0;
        Hx[1+0*2] = conj(Hx[0+1*2]);

        Hy.resize(4);
        Hy[0+0*2] = my * wy * wy * y;
        Hy[1+1*2] = my * wy * wy * (y - gy);
        Hy[0+1*2] = matrixop::IMAGIZ * W * V * exp(matrixop::IMAGIZ * W * y);
        Hy[1+0*2] = conj(Hy[0+1*2]);

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

#endif // _SPIN_BOSON_2D_POTENTIAL_HPP
