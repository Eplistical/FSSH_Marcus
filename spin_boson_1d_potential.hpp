#ifndef _SPIN_BOSON_1D_POTENTIAL_HPP
#define _SPIN_BOSON_1D_POTENTIAL_HPP


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
 * 1D SPIN BOSON POTENTIAL 
 * 
 *  H = [ 0.5 * mass * omega**2 * x**2                      V  ]
 *      [ V             0.5 * mass * omega**2 * (x-g)**2 + dG0 ]
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

    double param_mass = 1.0;
    double param_omega = 4.375e-5;
    double param_g = 4997.3;
    double param_dG0 = -0.018;
    double param_V = 2.5e-5;

    void output_potential_params() {
        ioer::info("# 1D Spin Boson Potential Parameters: ", 
                    " mass  = ", param_mass,
                    " omega = ", param_omega,
                    " g     = ", param_g,
                    " dG0   = ", param_dG0,
                    " V     = ", param_V,
                    ""
                    );
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 5, 
                "set_potenial_params: potential paramter vector size must be >= 5");
        param_mass = params[0];
        param_omega = params[1];
        param_g = params[2];
        param_dG0 = params[3];
        param_V = params[4];
    }


    vector< complex<double> > cal_H(const double x) {
        const double m = param_mass;
        const double w = param_omega;
        const double g = param_g;
        const double V = param_V;

        vector< complex<double> > H(4, 0.0);
        H[0+0*2] = 0.5 * m * pow(w * x, 2);
        H[1+1*2] = 0.5 * m * pow(w * (x - g), 2);
        H[0+1*2] = V;
        H[1+0*2] = conj(H[0+1*2]);

        return H;
    }

    vector< complex<double> > cal_nablaH(const double x) {
        const double m = param_mass;
        const double w = param_omega;
        const double g = param_g;
        const double V = param_V;

        vector< complex<double> > nablaH(4);
        vector< complex<double> >& Hx = nablaH;

        Hx.resize(4);
        Hx[0+0*2] = m * w * w * x;
        Hx[1+1*2] = m * w * w * (x - g);
        Hx[0+1*2] = 0.0;
        Hx[1+0*2] = conj(Hx[0+1*2]);

        return nablaH;
    }

    void cal_info_nume(const double x,
            vector<double>& eva, 
            vector< complex<double> >& dc,
            vector< complex<double> >& F,
            vector< complex<double> >& lastevt)
    {
        const int edim = 2;

        vector< complex<double> > evt;
        matrixop::hdiag(cal_H(x), eva, evt);
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
        dc = matrixop::matCmatmat(evt, cal_nablaH(x), evt, edim, edim);
        F.assign(edim * edim, 0.0);

        for (int j = 0; j < edim; ++j) {
            for (int k = 0; k < edim; ++k) {
                F[j+k*edim] = -dc[j+k*edim];
                if (j == k) {
                    dc[j+k*edim] = 0.0;
                }
                else {
                    dc[j+k*edim] /= (eva[k] - eva[j]);
                }
            }
        }

        // save evt to lastevt
        lastevt = std::move(evt);
    }
};

#endif 
