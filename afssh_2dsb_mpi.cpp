#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include <utility>
#include "misc/fmtstring.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/matrixop.hpp"
#include "misc/MPIer.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"
#include "spin_boson_2d_potential.hpp"

// A-FSSH WITH LANDRY'S ALGORITHM
//  LANDRY, SUBOTNIK, J. Chem. Phys. 137, 22A513 (2012)

// -- GLOBAL -- //


using namespace std;
namespace po = boost::program_options;
using state_t = vector< complex<double> >;

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};


// const paramters
const complex<double> zI(0.0, 1.0);
int ndim = 2;
int edim = 2;

const vector<double> mass { 1.0, 1.0};
const vector<double> omega { 4.375e-5, 0.0};
const vector<double> g { 4997.3, 0.0};
const double kT = 9.5e-4;
const vector<double> fric_gamma = { 1.5e-4, 0.0 };

const vector<double> sigma_r { 0.0, 0.0 };
const vector<double> sigma_p = sqrt(kT * mass);

// tunable parameters
double V = 2.5e-5;
double dG0 = -0.018;

vector<double> init_r { 0.0, 0.0};
vector<double> init_p { 0.0, 0.0};
vector<double> init_s { 1.0, 0.0 };

double dt = 1.0;
int Nstep = 10000;
int output_step = 100;
int Ntraj = 2000;
bool enable_hop = true;
bool enable_deco = true;
bool enable_prev = true;
int seed = 0;

// cache
vector<double> eva;
vector< vector< complex<double> > > dc;
vector< vector< complex<double> > > F;
vector< complex<double> > lastevt;


// -- FUNCTIONS -- //


bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("V", po::value<double>(&V), "coupling magnitude V")
        ("dG0", po::value<double>(&dG0), "driving force")
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init_r vector")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init_p vector")
        ("init_s", po::value<decltype(init_s)>(&init_s)->multitoken(), "init_s vector")
        ("Ntraj", po::value<decltype(Ntraj)>(&Ntraj), "# traj")
        ("Nstep", po::value<decltype(Nstep)>(&Nstep), "# step")
        ("output_step", po::value<decltype(output_step)>(&output_step), "# step for output")
        ("dt", po::value<decltype(dt)>(&dt), "single time step")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_deco", po::value<bool>(&enable_deco), "enable decoherence")
        ("enable_prev", po::value<bool>(&enable_prev), "enable momentum reversal")
        ("seed", po::value<decltype(seed)>(&seed), "random seed")
        ;
    po::variables_map vm; 
    //po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

void init_state(state_t& state, const vector<double>& init_r, const vector<double>& init_p, 
        const vector<double>& mass, const vector<double> init_s) 
{
    // state = r(ndim), p(ndim), c(edim), s(1), rmom(ndim * edim * edim), pmom(ndim * edim * edim)
    state.resize(ndim * 2 + edim + 1 + ndim * edim * edim * 2, matrixop::ZEROZ);

    // nuclear DoF (r, p)
    for (int i(0); i < ndim; ++i) {
        state[i].real(randomer::normal(init_r[i], sigma_r[i]));
        state[ndim + i].real(randomer::normal(init_p[i], sigma_p[i]));
    }

    // electronic DoF (c, s)
    vector<double> s_normalized = init_s / sum(init_s);
    for (int i(0); i < edim; ++i) {
        state[ndim * 2 + i] = sqrt(s_normalized[i]);
    }
    state[ndim * 2 + edim].real(randomer::discrete(s_normalized.begin(), s_normalized.end()));

    // moment DoF (rmom, pmom) --> 0.0
}

void integrator(state_t& state, const double dt) {
    // extract info
    vector< complex<double> > r(state.begin(), state.begin() + ndim);
    vector< complex<double> > p(state.begin() + ndim, state.begin() + ndim * 2);
    vector< complex<double> > c(state.begin() + ndim * 2, state.begin() + ndim * 2 + edim);
    int s = static_cast<int>(state[ndim * 2 + edim].real());
    vector< complex<double> > rmom(state.begin() + ndim * 2 + edim + 1, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim);
    vector< complex<double> > pmom(state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim * 2);

    // ele part, RK4
    
    auto el_rk4_func = [&dt, &p]
        (const vector< complex<double> >& c, vector< complex<double> >& k) 
        {
            k = -zI * c * eva;
            for (int i(0); i < ndim; ++i) {
                k -= p[i] / mass[i] * matrixop::matvec(dc[i], c);
            }
            k *= dt;
        };

    vector< complex<double> > k1, k2, k3, k4;
    el_rk4_func(c, k1);
    el_rk4_func(c + 0.5 * k1, k2);
    el_rk4_func(c + 0.5 * k2, k3);
    el_rk4_func(c + k3, k4);

    c += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    copy(c.begin(), c.end(), state.begin() + ndim * 2);

    // moments part, RK4

    if (enable_deco) {
        vector< complex<double> > Emat(edim * edim);
        for (int j(0); j < edim; ++j) {
            Emat[j+j*edim] = eva[j];
        }

        vector< complex<double> > sigma = outer_product(c, conj(c));
        vector< vector< complex<double> > > dFsigma(ndim);

        vector< complex<double> > dF;
        for (int i(0); i < ndim; ++i) {
            dF = F[i] - F[i][s+s*edim] * matrixop::eye(edim);
            dFsigma[i] = matrixop::anticommutator(dF, sigma);
        }

        auto deco_rk4_func = [&Emat, &p, &dFsigma, &dt]
            (
             const int i, const int s,
             const vector< complex<double> >& rmom, 
             const vector< complex<double> >& pmom, 
             vector< complex<double> >& l, 
             vector< complex<double> >& m 
            ) 
            {
                vector< complex<double> > rtmp(edim * edim, matrixop::ZEROZ);
                vector< complex<double> > ptmp(edim * edim, matrixop::ZEROZ);
                for (int j(0); j < ndim; ++j) {
                    rtmp += p[j] / mass[j] * matrixop::commutator(dc[j], rmom);
                    ptmp += p[j] / mass[j] * matrixop::commutator(dc[j], pmom);
                } 
                vector< complex<double> > TR = -zI * matrixop::commutator(Emat, rmom) + pmom / mass[i] - rtmp;
                vector< complex<double> > TP = -zI * matrixop::commutator(Emat, pmom) + 0.5 * dFsigma[i] - ptmp;

                l = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
                m = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));
            };

        vector< complex<double> > l1, l2, l3, l4; // rmom
        vector< complex<double> > m1, m2, m3, m4; // pmom
        vector< complex<double> > irmom(edim * edim);
        vector< complex<double> > ipmom(edim * edim);
        for (int i(0); i < ndim; ++i) {
            copy(rmom.begin() + i * edim * edim, rmom.begin() + (i + 1) * edim * edim, irmom.begin());
            copy(pmom.begin() + i * edim * edim, pmom.begin() + (i + 1) * edim * edim, ipmom.begin());

            deco_rk4_func(i, s, irmom, ipmom, l1, m1);
            deco_rk4_func(i, s, irmom + 0.5 * l1, ipmom + 0.5 * m1, l2, m2);
            deco_rk4_func(i, s, irmom + 0.5 * l2, ipmom + 0.5 * m2, l3, m3);
            deco_rk4_func(i, s, irmom + l3, ipmom + m3, l4, m4);

            irmom += l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;
            ipmom += m1 / 6.0 + m2 / 3.0 + m3 / 3.0 + m4 / 6.0;

            copy(irmom.begin(), irmom.end(), rmom.begin() + i * edim * edim);
            copy(ipmom.begin(), ipmom.end(), pmom.begin() + i * edim * edim);
        }
        copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1);
        copy(pmom.begin(), pmom.end(), state.begin() + ndim * 2 + edim + 1 + edim * edim * ndim);
    }

    // nuclear part, VV

    vector<double> random_force = randomer::vnormal(ndim) * sqrt(2.0 * kT / dt * fric_gamma);
    vector< complex<double> > force(ndim);
    for (int i(0); i < ndim; ++i) {
        force[i] = F[i][s + s * edim].real() - fric_gamma[i] * p[i] + random_force[i];
    }
    p += 0.5 * dt * force;

    r += dt * p / mass;

    cal_info_nume(real(r), eva, dc, F, lastevt);

    for (int i(0); i < ndim; ++i) {
        force[i] = F[i][s + s * edim].real() - fric_gamma[i] * p[i] + random_force[i];
    }
    p += 0.5 * dt * force;

    copy(r.begin(), r.end(), state.begin());
    copy(p.begin(), p.end(), state.begin() + ndim);
}

int hopper(state_t& state) 
{
    // extract info
    vector< complex<double> > r(state.begin(), state.begin() + ndim);
    vector< complex<double> > p(state.begin() + ndim, state.begin() + ndim * 2);
    vector< complex<double> > c(state.begin() + ndim * 2, state.begin() + ndim * 2 + edim);
    int s = static_cast<int>(state[ndim * 2 + edim].real());
    vector< complex<double> > rmom(state.begin() + ndim * 2 + edim + 1, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim);
    vector< complex<double> > pmom(state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim * 2);

    const int from = s;
    const int to = 1 - s;

    // calc hop prob
    complex<double> vd = 0.0;
    for (int i(0); i < ndim; ++i) {
        vd += p[i] / mass[i] * dc[i][to+from*edim];
    }
    double g = -2 * dt * (c[from] * conj(c[to]) * vd).real() / (c[from] * conj(c[from])).real();
    if (g >= 1.0) {
        ioer::info(misc::fmtstring("# hopper: WARNING -- g = %.3f is larger than 1", g));
    }
    double dE = eva[to] - eva[from];

    // random number
    if (randomer::rand() < g) {
        // momentum-rescaling direction: (x-direction)
        vector<double> n(ndim, 0.0);

        /*
        // Mehod #2
        const vector<double> dcR { dc[0][from+to*edim].real(), dc[1][from+to*edim].real()  };
        const vector<double> dcI { dc[0][from+to*edim].imag(), dc[1][from+to*edim].imag()  };
        const double diff_norm2 = norm2(dcR) - norm2(dcI);
        const double twice_eta0 = std::atan(-2 * (dcR[0] * dcI[0] + dcR[1] * dcI[1]) / diff_norm2);
        double eta;
        if (cos(twice_eta0) * diff_norm2 > 0.0) {
            eta = 0.5 * twice_eta0;
        }
        else {
            eta = 0.5 * twice_eta0 + 0.5 * M_PI;

        }
        const complex<double> eieta = exp(zI * eta);
        n[0] = (eieta * dc[0][from+to*edim]).real();
        n[1] = (eieta * dc[1][from+to*edim]).real();
        */
        n[0] = dc[0][from+to*edim].real();
        n[1] = dc[1][from+to*edim].real();

        // hop
        if (norm(n) > 1e-40) {
            vector<double> pn = component(real(p), n);
            double pn_norm = norm(pn); 
            double tmp = pn_norm * pn_norm - 2 * mass[0] * dE; // masses along all dimension should be identical
            if (tmp > 0.0) {
                // hop accepted
                double pn_norm_new = sqrt(tmp);
                p += (pn_norm_new - pn_norm) / pn_norm * pn;

                // replace p & s
                copy(p.begin(), p.end(), state.begin() + ndim);
                state[2 * ndim + edim] = to;

                // moments
                rmom.assign(edim * edim * ndim, matrixop::ZEROZ);
                for (int k(0); k < edim; ++k) {
                    if (k == to or abs(c[k]) < 1e-8) {
                        for (int i(0); i < ndim; ++i) {
                            pmom[k+k*edim+i*edim*edim] = matrixop::ZEROZ;
                        }
                    }
                    else {
                        double tmp = pow(pn_norm_new, 2) + 2.0 * mass[0] * (eva[to] - eva[k]);
                        if (tmp <= 0.0) {
                            for (int i(0); i < ndim; ++i) {
                                pmom[k+k*edim+i*edim*edim] = matrixop::ZEROZ;
                            }
                        }
                        else {
                            vector<double> dpkk = (sqrt(tmp) - pn_norm_new) * pow(abs(c[k]), 2) / pn_norm * pn;
                            for (int i(0); i < ndim; ++i) {
                                pmom[k+k*edim+i*edim*edim] = dpkk[i];
                            }
                        }
                    }
                }
                copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1);
                copy(pmom.begin(), pmom.end(), state.begin() + ndim * 2 + edim + 1 + edim * edim * ndim);

                return (from < to) ? HOP_UP : HOP_DN;
            }
            else {
                // hop frustrated
                if (enable_prev) {
                    double F0d01 = 0.0;
                    double F1d01 = 0.0;
                    double pd01 = 0.0;
                    for (int i(0); i < ndim; ++i) {
                        F0d01 += (F[i][0+0*edim] * dc[i][0+1*edim]).real();
                        F1d01 += (F[i][1+1*edim] * dc[i][0+1*edim]).real();
                        pd01 += (p[i] * dc[i][0+1*edim]).real();
                        if (F0d01 * F1d01 < 0.0 and pd01 * F1d01 < 0.0)
                        {
                            // momentum reversal: along d direction 
                            // TODO: WHAT IF dc IS COMPLEX?
                            p = p - pn * 2.0;
                            copy(p.begin(), p.end(), state.begin() + ndim);
                        }
                    }
                }
                return HOP_FR;
            }
        }
    }
    // hop rejected
    return HOP_RJ;
}

void decoherencer(state_t& state, const double xi = 1.0) {
    // extract info
    vector< complex<double> > r(state.begin(), state.begin() + ndim);
    vector< complex<double> > p(state.begin() + ndim, state.begin() + ndim * 2);
    vector< complex<double> > c(state.begin() + ndim * 2, state.begin() + ndim * 2 + edim);
    int s = static_cast<int>(state[ndim * 2 + edim].real());
    vector< complex<double> > rmom(state.begin() + ndim * 2 + edim + 1, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim);
    vector< complex<double> > pmom(state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim, state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim * 2);

    // collapse & reset
    for (int n(0); n < edim; ++n) {
        if (n != s) {
            double Pcollapse = 0.0;
            double Preset = 0.0;
            complex<double> Fsn_rmom = 0.0;
            for (int i(0); i < ndim; ++i) {
                Pcollapse += 0.5 * (F[i][n+n*edim] - F[i][s+s*edim]).real() * rmom[n+n*edim+i*edim*edim].real();
                Fsn_rmom += F[i][s+n*edim] * rmom[n+n*edim+i*edim*edim];
            }
            Preset = -Pcollapse;
            Pcollapse -= 2.0 * xi * abs(Fsn_rmom);

            if (randomer::rand() < Pcollapse) {
                // Landry 2012 Eqn (45) -- TYPO --
                c[n] = 0.0;
                c[s] = c[s] / abs(c[s]) * sqrt(pow(abs(c[s]), 2) + pow(abs(c[n]), 2));
                copy(c.begin(), c.end(), state.begin() + ndim * 2);

                for (int j(0); j < edim; ++j) {
                    for (int i(0); i < ndim; ++i) {
                        rmom[j+n*edim+i*edim*edim] = 0.0;
                        pmom[j+n*edim+i*edim*edim] = 0.0;
                        rmom[n+j*edim+i*edim*edim] = 0.0;
                        pmom[n+j*edim+i*edim*edim] = 0.0;
                    }
                } 
                copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1);
                copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim);
            }
            if (randomer::rand() < Preset) {
                for (int j(0); j < edim; ++j) {
                    for (int i(0); i < ndim; ++i) {
                        rmom[j+n*edim+i*edim*edim] = 0.0;
                        pmom[j+n*edim+i*edim*edim] = 0.0;
                        rmom[n+j*edim+i*edim*edim] = 0.0;
                        pmom[n+j*edim+i*edim*edim] = 0.0;
                    }
                } 
                copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1);
                copy(rmom.begin(), rmom.end(), state.begin() + ndim * 2 + edim + 1 + ndim * edim * edim);
            }
        }
    }
}

bool check_end(const state_t& state) {
    return false;
}

void afssh_nd_mpi() {
    // potential para
    set_potenial_params( vector<double> { mass[0], mass[1], omega[0], omega[1], g[0], g[1], dG0, V });

    // assign jobs
    const vector<int> my_jobs = MPIer::assign_job(Ntraj);
    const int my_Ntraj = my_jobs.size();

    // initialize
    vector<state_t> state(my_Ntraj);
    for (int itraj(0); itraj < my_Ntraj; ++itraj) {
        init_state(state[itraj], init_r, init_p, mass, init_s);
    }

    // statistics
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    vector<double> hop_count_summary(10, 0.0);

    double n0d = 0.0, n1d = 0.0;
    double KE = 0.0, PE = 0.0;

    // recorders
    int Nrec = Nstep / output_step;
    map<string, vector<double> > record {
        { "n0d", vector<double>(Nrec, 0.0) },
        { "n1d", vector<double>(Nrec, 0.0) },
        { "KE", vector<double>(Nrec, 0.0) },
        { "PE", vector<double>(Nrec, 0.0) },
    };


    //  ----  MAIN LOOP  ---- //


    vector< vector< complex<double> > > lastevt_save(my_Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < my_Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                vector< complex<double> > r(state[itraj].begin(), state[itraj].begin() + ndim);
                cal_info_nume(real(r), eva, dc, F, lastevt);
                // hopper
                if (enable_hop) {
                    int hopflag = hopper(state[itraj]);
                    switch (hopflag) {
                        case HOP_UP : { hopup += 1.0; hop_count[itraj] += 1.0; break; }
                        case HOP_DN : { hopdn += 1.0; hop_count[itraj] += 1.0; break; }
                        case HOP_FR : { hopfr += 1.0; break; }
                        case HOP_RJ : { hoprj += 1.0; break; }
                        default : break;
                    }
                }
                // decoherence
                if (enable_deco) {
                    decoherencer(state[itraj]);
                }
                // integrate t -> t + dt
                integrator(state[itraj], dt);
                // save lastevt
                lastevt_save[itraj] = move(lastevt);
            }
        }

        if (istep % output_step == 0) {
            // record data
            n0d = n1d = 0.0;
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0d, &n1d, &KE, &PE](const state_t& st) { 
                        // extract info
                        vector< complex<double> > r(st.begin(), st.begin() + ndim);
                        vector< complex<double> > p(st.begin() + ndim, st.begin() + ndim * 2);
                        vector< complex<double> > c(st.begin() + ndim * 2, st.begin() + ndim * 2 + edim);
                        int s = static_cast<int>(st[ndim * 2 + edim].real());

                        vector<double> evatmp;
                        vector< complex<double> > U;
                        matrixop::hdiag(cal_H(real(r)), evatmp, U);

                        n0d += pow(abs(U[0+s*2]), 2) + 2 * (U[0+0*2] * c[0] * conj(c[1]) * conj(U[0+1*2])).real();
                        n1d += pow(abs(U[1+s*2]), 2) + 2 * (U[1+0*2] * c[0] * conj(c[1]) * conj(U[1+1*2])).real();

                        // energy
                        KE += 0.5 * sum(pow(abs(p), 2) / mass);
                        PE += evatmp[s]; 
                    });

            const int irec = istep / output_step;
            record["n0d"][irec] = n0d;
            record["n1d"][irec] = n1d;
            record["KE"][irec] = KE;
            record["PE"][irec] = PE;
        }
    }
    MPIer::barrier();


    // ----  COLLECT DATA  ---- //
    

    for (int r(1); r < MPIer::size; ++r) {
        if (MPIer::rank == r) {
            MPIer::send(0, record["n0d"]);
            MPIer::send(0, record["n1d"]);
            MPIer::send(0, record["KE"]);
            MPIer::send(0, record["PE"]);
        }
        else if (MPIer::master) {
            vector<double> vbuf;
            MPIer::recv(r, vbuf); record["n0d"] += vbuf;
            MPIer::recv(r, vbuf); record["n1d"] += vbuf;
            MPIer::recv(r, vbuf); record["KE"] += vbuf;
            MPIer::recv(r, vbuf); record["PE"] += vbuf;
        }
        MPIer::barrier();
    }
    MPIer::barrier();


    // ----  PROCESS & COLLECT HOP STATISTICS DATA  ---- //


    for_each(hop_count.begin(), hop_count.end(),
            [&hop_count_summary](double x) { 
                if (hop_count_summary.size() < static_cast<int>(x) + 1) {
                    hop_count_summary.resize(static_cast<int>(x) + 1);
                }
                hop_count_summary[static_cast<int>(x)] += 1.0; 
            });
    MPIer::barrier();

    for (int r = 1; r < MPIer::size; ++r) {
        if (MPIer::rank == r) {
            MPIer::send(0, hopup, hopdn, hopfr, hoprj, hop_count_summary);
        }
        else if (MPIer::master) {
            double buf;
            vector<double> vbuf;

            MPIer::recv(r, buf); hopup += buf;
            MPIer::recv(r, buf); hopdn += buf;
            MPIer::recv(r, buf); hopfr += buf;
            MPIer::recv(r, buf); hoprj += buf;

            MPIer::recv(r, vbuf); 

            if (hop_count_summary.size() < vbuf.size()) {
                hop_count_summary.resize(vbuf.size());
                hop_count_summary += vbuf;
            }
            else {
                for (int i(0); i < vbuf.size(); ++i) {
                    hop_count_summary[i] += vbuf[i];
                }
            }
        }
        MPIer::barrier();
    }
    MPIer::barrier();


    // ----  OUTPUT  ---- //


    if (MPIer::master) {
        // Output parameters
        output_potential_params();
        ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, 
                    " mass = ", mass, " kT = ", kT, 
                    " omega = ", omega, " g = ", g, " dG0 = ", dG0, 
                    " fric_gamma = ", fric_gamma, 
                    " V = ", V, 
                    " init_r = ", init_r, " init_p = ", init_p, 
                    " sigma_r = ", sigma_r, " sigma_p = ", sigma_p, 
                    " init_s = ", init_s,
                    " enable_hop = ", enable_hop,
                    " enable_deco = ", enable_deco,
                    " enable_prev = ", enable_prev,
                    " seed = ", seed
                );
        ioer::info("# ln(V) = ", log(V));
        ioer::tabout('#', "t", "n0d", "n1d", "KE", "PE", "Etot");

        for (int irec = 0; irec < Nrec; ++irec) {
            ioer::tabout(
                    "", irec * output_step * dt, 
                    record["n0d"][irec] / Ntraj, 
                    record["n1d"][irec] / Ntraj, 
                    record["KE"][irec] / Ntraj, 
                    record["PE"][irec] / Ntraj, 
                    (record["KE"][irec] + record["PE"][irec]) / Ntraj
                    );
        }

        // Output hop info
        ioer::info("# hop statistics: ");
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
    }
    MPIer::barrier();
}

int check() {
    vector<double> r(ndim, 0.0);
    int N = 500;
    vector<double> xarr = linspace(-10000.0, 15000.0, N);
    vector<double> yarr = linspace(-8.0, 8.0, N);

    for (int ix(0); ix < N; ++ix) {
        for (int iy(0); iy < N; ++iy) {
            double x = xarr[ix];
            double y = yarr[iy];
            r[0] = x;
            r[1] = y;

            cal_info_nume(r, eva, dc, F, lastevt);
            ioer::tabout(
                    x, y, 
                    eva[0], eva[1], 
                    F[0][0+0*2].real(), F[1][0+0*2].real(),
                    F[0][1+1*2].real(), F[1][1+1*2].real(),
                    real(dc[0][0+1*2]), imag(dc[0][0+1*2]), 
                    real(dc[1][0+1*2]), imag(dc[1][0+1*2])
                    );
        }
    }
    return 0;
}


// -- MAIN -- //


int main(int argc, char** argv) {
    MPIer::setup();
    if (argc < 2) {
        if (MPIer::master) 
            ioer::info("use --help for detailed info");
    }
    else if (string(argv[1]) == "check") {
        if (MPIer::master) 
            check();
    }
    else {
        if (argparse(argc, argv) == true) {
            randomer::seed(MPIer::assign_random_seed(seed));
            if (MPIer::master) timer::tic();
            afssh_nd_mpi();
            if (MPIer::master) ioer::info("# ", timer::toc());
        }
    }
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
