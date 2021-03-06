#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
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
#include "spin_boson_1d_potential.hpp"

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
const int ndim = 1;
const int edim = 2;

const double mass = 1.0;
const double kT = 9.5e-4;
const double omega = 4.375e-5;
const double Er = 0.0239;
const double M = sqrt(0.5 * Er * mass * omega * omega);
const double g = 2.0 * M / mass / omega / omega;
const double fric_gamma = 1.5e-4; // 0.0024; 

const double sigma_x = 0.0;
const double sigma_p = sqrt(kT * mass);

// tunable parameters
double V = 2.5e-5;
double dG0 = -0.018;

double init_x = 0.0;
double init_p = 0.0;
double init_s = 0.0;

double dt = 1.0;
int Nstep = 10000;
int output_step = 100;
int Ntraj = 5000;
bool enable_hop = true;
bool enable_deco = true;
bool enable_prev = true;
int seed = 0;

// cache
vector<double> eva(edim);
vector< complex<double> > lastevt;
vector< complex<double> > F(edim * edim);
vector< complex<double> > dc(edim * edim);


// -- FUNCTIONS -- //


bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("V", po::value<double>(&V), "coupling magnitude V")
        ("dG0", po::value<double>(&dG0), "driving force")
        ("init_x", po::value<double>(&init_x), "init x")
        ("init_p", po::value<double>(&init_p), "potential para init_p")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_deco", po::value<bool>(&enable_deco), "enable decoherence")
        ("enable_prev", po::value<bool>(&enable_prev), "enable momentum reversal")
        ("seed", po::value<int>(&seed), "random seed. Default 0")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

void init_state(state_t& state) {
    // state = r(1), p(1), c(2), s(1), rmom(4), pmom(4)
    state.resize(13, matrixop::ZEROZ);

    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_p, sigma_p)); 
    state[2].real(sqrt(1.0 - init_s));
    state[3].real(sqrt(init_s));
    state[4].real((randomer::rand() < init_s) ? 1.0 : 0.0);
}

bool check_end(const state_t& state) {
    return false;
}

void integrator(state_t& state, const double dt) {
    // extract info
    // state = r(1), p(1), c(2), s(1), rmom(4), pmom(4)
    double x = state[0].real();
    double p = state[1].real();
    vector< complex<double> > c(state.begin() + 2, state.begin() + 4);
    int s = static_cast<int>(state[4].real());
    vector< complex<double> > rmom(state.begin() + 5, state.begin() + 9);
    vector< complex<double> > pmom(state.begin() + 9, state.begin() + 13);

    // ele part, RK4

    auto el_rk4_func = [&dt, &p]
        ( const vector< complex<double> >& c, vector< complex<double> >& k ) 
        {
            k = dt * (-zI * c * eva - p / mass * matrixop::matvec(dc, c));
        };

    vector< complex<double> > k1, k2, k3, k4; // c

    el_rk4_func(c, k1);
    el_rk4_func(c + 0.5 * k1, k2);
    el_rk4_func(c + 0.5 * k2, k3);
    el_rk4_func(c + k3, k4);

    c += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    copy(c.begin(), c.end(), state.begin() + 2);

    // moments part, RK4

    if (enable_deco) {
        vector< complex<double> > dF, sigma, dFsigma;

        vector< complex<double> > Emat(edim * edim);
        for (int i(0); i < edim; ++i) {
            Emat[i+i*edim] = eva[i];
        }

        dF = F - F[s+s*edim] * matrixop::eye(edim);
        sigma.resize(edim * edim);
        for (int j(0); j < edim; ++j) {
            for (int k(0); k < edim; ++k) {
                sigma[j+k*edim] = c[j] * conj(c[k]);
            }
        }
        dFsigma = matrixop::anticommutator(dF, sigma);

        auto deco_rk4_func = [&dt, &Emat, &p, &dFsigma]
            (
                const int s,
                const vector< complex<double> >& rmom, 
                const vector< complex<double> >& pmom, 
                vector< complex<double> >& l, 
                vector< complex<double> >& m 
            ) 
            {
                vector< complex<double> > TR = -zI * matrixop::commutator(Emat, rmom) + pmom / mass -  p / mass * matrixop::commutator(dc, rmom);
                vector< complex<double> > TP = -zI * matrixop::commutator(Emat, pmom) + 0.5 * dFsigma - p / mass * matrixop::commutator(dc, pmom);
                l = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
                m = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));
            };

        vector< complex<double> > l1, l2, l3, l4; // rmom
        vector< complex<double> > m1, m2, m3, m4; // pmom

        deco_rk4_func(s, rmom, pmom, l1, m1);
        deco_rk4_func(s, rmom + 0.5 * l1, pmom + 0.5 * m1, l2, m2);
        deco_rk4_func(s, rmom + 0.5 * l2, pmom + 0.5 * m2, l3, m3);
        deco_rk4_func(s, rmom + l3, pmom + m3, l4, m4);

        rmom += l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;
        pmom += m1 / 6.0 + m2 / 3.0 + m3 / 3.0 + m4 / 6.0;

        copy(rmom.begin(), rmom.end(), state.begin() + 5);
        copy(pmom.begin(), pmom.end(), state.begin() + 9);
    }

    // nuclear part, VV

    const double random_force = randomer::normal(0.0, sqrt(2 * fric_gamma * kT / dt));
    p += 0.5 * dt * (F[s+s*edim].real() - fric_gamma * p + random_force);
    x += dt * p / mass;
    cal_info_nume(x, eva, dc, F, lastevt);
    p += 0.5 * dt * (F[s+s*edim].real() - fric_gamma * p + random_force);

    state[0].real(x);
    state[1].real(p);
}

int hopper(state_t& state) {
    // extract info
    // state = r(1), p(1), c(2), s(1), rmom(4), pmom(4)
    double x = state[0].real();
    double p = state[1].real();
    vector< complex<double> > c(state.begin() + 2, state.begin() + 4);
    int s = static_cast<int>(state[4].real());
    vector< complex<double> > rmom(state.begin() + 5, state.begin() + 9);
    vector< complex<double> > pmom(state.begin() + 9, state.begin() + 13);

    // calc hop prob
    double g = -2 * dt * (c[s] * conj(c[1-s]) * (p / mass * dc[1-s+s*2])).real() / (c[s] * conj(c[s])).real();
    if (g >= 1.0) {
        ioer::info(misc::fmtstring("# hopper: WARNING -- g = %.3f is larger than 1", g));
    }
    double dE = eva[1-s] - eva[s];
    // hop
    if (randomer::rand() < g) {
        double tmp = p * p - 2 * mass * dE;
        if (tmp > 0.0) {
            // hop accepted
            double pnew = sqrt(tmp);
            state[1] = pnew * (p < 0.0 ? -1 : 1); 
            state[4].real(1.0 - s); 

            // moments
            if (enable_deco) {
                rmom.assign(edim * edim, 0.0);

                for (int k(0); k < edim; ++k) {
                    if (k == s or abs(c[k]) < 1e-8) {
                        pmom[k+k*edim] = 0.0;
                    }
                    else {
                        double ptmp = pnew * pnew + 2 * mass * (eva[s] - eva[k]);
                        if (ptmp < 0.0) {
                            pmom[k+k*edim] = 0.0;
                        }
                        else {
                            pmom[k+k*edim] = pow(abs(c[k]), 2) * (sqrt(ptmp) - pnew) * 
                                (p > 0.0 ? 1.0 : -1.0);
                        }
                    }
                }

                copy(rmom.begin(), rmom.end(), state.begin() + 5);
                copy(pmom.begin(), pmom.end(), state.begin() + 9);
            }

            return (s == 0) ? HOP_UP : HOP_DN;
        }
        else {
            // hop frustrated
            if (enable_prev) {
                double F0d01 = (F[0+0*edim] * dc[0+1*edim]).real();
                double F1d01 = (F[1+1*edim] * dc[0+1*edim]).real();
                double vd01 = p / mass * dc[0+1*2].real();
                if (F0d01 * F1d01 < 0.0 and vd01 * F1d01 < 0.0)
                {
                    // momentum reversal
                    state[1] *= -1.0;
                }
            }
            return HOP_FR; 
        }
    }   
    // hop rejected
    return HOP_RJ;
}

void decoherencer(state_t& state, const double xi = 1.0) {
    // extract info
    // state = r(1), p(1), c(2), s(1), rmom(4), pmom(4)
    double x = state[0].real();
    double p = state[1].real();
    vector< complex<double> > c(state.begin() + 2, state.begin() + 4);
    int s = static_cast<int>(state[4].real());
    vector< complex<double> > rmom(state.begin() + 5, state.begin() + 9);
    vector< complex<double> > pmom(state.begin() + 9, state.begin() + 13);

    // collapse & reset
    for (int n(0); n < edim; ++n) {
        if (n != s) {
            const double Pcollapse = dt * (0.5 * ((F[n+n*edim] - F[s+s*edim]).real() * rmom[n+n*edim].real()) - 2 * xi * abs(F[s+n*edim] * rmom[n+n*edim]));
            const double Preset = dt * (-0.5 * ((F[n+n*edim] - F[s+s*edim]).real() * rmom[n+n*edim].real()));

            if (randomer::rand() < Pcollapse) {
                // Landry 2012 Eqn (45) -- TYPO --
                c[n] = 0.0;
                c[s] = c[s] / abs(c[s]) * sqrt(pow(abs(c[s]), 2) + pow(abs(c[n]), 2));
                copy(c.begin(), c.end(), state.begin() + 2);

                for (int j(0); j < edim; ++j) {
                    rmom[j+n*edim] = 0.0;
                    pmom[j+n*edim] = 0.0;
                    rmom[n+j*edim] = 0.0;
                    pmom[n+j*edim] = 0.0;
                } 
                copy(rmom.begin(), rmom.end(), state.begin() + 5);
                copy(pmom.begin(), pmom.end(), state.begin() + 9);
            }
            if (randomer::rand() < Preset) {
                for (int j(0); j < edim; ++j) {
                    rmom[j+n*edim] = 0.0;
                    pmom[j+n*edim] = 0.0;
                    rmom[n+j*edim] = 0.0;
                    pmom[n+j*edim] = 0.0;
                } 
                copy(rmom.begin(), rmom.end(), state.begin() + 5);
                copy(pmom.begin(), pmom.end(), state.begin() + 9);
            }
        }
    }
}


void afssh_1d_mpi() {
    // potential para
    set_potenial_params(vector<double> { mass, omega, g, dG0, V });

    // assign jobs
    const vector<int> my_jobs = MPIer::assign_job(Ntraj);
    const int my_Ntraj = my_jobs.size();

    // initialize
    vector<state_t> state(my_Ntraj);
    for (auto& st : state) {
        init_state(st);
    }

    // statistics
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    vector<double> hop_count_summary(10, 0.0);

    double n0d = 0.0, n1d = 0.0;
    double KE = 0.0, PE = 0.0;

    // recorder
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
                cal_info_nume(state[itraj][0].real(), eva, dc, F, lastevt);
                // hopper
                if (enable_hop) {
                    int hopflag = hopper(state[itraj]);
                    switch (hopflag) {
                        case HOP_UP : hopup += 1.0; hop_count[itraj] += 1.0; break;
                        case HOP_DN : hopdn += 1.0; hop_count[itraj] += 1.0; break;
                        case HOP_FR : hopfr += 1.0; break;
                        case HOP_RJ : hoprj += 1.0; break;
                        default : break;
                    }
                }
                // decoherence
                if (enable_deco) {
                    decoherencer(state[itraj]);
                }
                // integrater t -> t + dt
                integrator(state[itraj], dt);
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }

        if (istep % output_step == 0) {
            // record data
            n0d = n1d = 0.0;
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0d, &n1d, &KE, &PE](const state_t& st) { 
                        double x = st[0].real();
                        double p = st[1].real();
                        vector< complex<double> > c(st.begin() + 2, st.begin() + 4);
                        int s = static_cast<int>(st[4].real());

                        vector<double> evatmp;
                        vector< complex<double> > U;
                        matrixop::hdiag(cal_H(x), evatmp, U);

                        n0d += pow(abs(U[0+s*2]), 2) + 2 * (U[0+0*2] * c[0] * conj(c[1]) * conj(U[0+1*2])).real();
                        n1d += pow(abs(U[1+s*2]), 2) + 2 * (U[1+0*2] * c[0] * conj(c[1]) * conj(U[1+1*2])).real();

                        // energy
                        KE += p * p / mass / 2;
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
        // para & header
        output_potential_params();
        ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step,
                    " mass = ", mass, " kT = ", kT, 
                    " omega = ", omega, " g = ", g,  " Er = ", Er, " dG0 = ", dG0, 
                    " fric_gamma = ", fric_gamma, 
                    " V = ", V, 
                    " init_x = ", init_x, " init_p = ", init_p, 
                    " sigma_x = ", sigma_x, " sigma_p = ", sigma_p, 
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

        // output hop info
        ioer::info("# hop statistics: ");
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
    }
    MPIer::barrier();
}

void check() {
    // check surf
    for (double x = -10000; x < 15000; x += 10) {
        cal_info_nume(x, eva, dc, F, lastevt);
        ioer::tabout(x, eva, F[0+0*edim].real(), F[1+1*edim].real(), dc[0+1*edim].real(), dc[0+1*edim].imag());
    }
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
            afssh_1d_mpi();
            if (MPIer::master) ioer::info("# ", timer::toc());
        }
    }
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
