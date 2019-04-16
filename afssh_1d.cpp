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
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"
#include "spin_boson_1d_potential.hpp"


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
const double omega = 3.5e-4;
const double Er = 0.0239;
const double M = sqrt(0.5 * Er * mass * omega * omega);
const double g = 2.0 * M / mass / omega / omega;
const double dG0 = -0.015;
const double fric_gamma = 0.0024; 

const double sigma_x = 0.0;
//const double sigma_x = sqrt(kT / mass / omega / omega); 
const double sigma_px = sqrt(kT * mass);

// tunable parameters
double V = 1.49e-5;
double init_x = 0.0;
double init_p = 0.0;
double init_s = 0.0;

double dt = 2.5;
int Nstep = 10000;
int output_step = 100;
int Ntraj = 5000;
bool enable_hop = true;
int seed = 0;

// cache
vector<double> eva(edim);
vector< complex<double> > lastevt(edim * edim);
vector< complex<double> > F(edim * edim);
vector< complex<double> > dc(edim * edim);


// -- FUNCTIONS -- //


bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("V", po::value<double>(&V), "coupling magnitude V")
        ("init_x", po::value<double>(&init_x), "init x")
        ("init_p", po::value<double>(&init_p), "potential para init_p")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
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
    state[1].real(randomer::normal(init_p, sigma_px)); 
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

    double v = p / mass;
    vector< complex<double> > Emat(edim * edim);
    for (int i(0); i < edim; ++i) {
        Emat[i+i*edim] = eva[i];
    }

    vector< complex<double> > k1, k2, k3, k4; // c

    k1 = dt * (-zI * c * eva - v * matrixop::matvec(dc, c));
    k2 = dt * (-zI * (c + 0.5 * k1) * eva - v * matrixop::matvec(dc, c + 0.5 * k1));
    k3 = dt * (-zI * (c + 0.5 * k2) * eva - v * matrixop::matvec(dc, c + 0.5 * k2));
    k4 = dt * (-zI * (c + k3) * eva - v * matrixop::matvec(dc, c + k3));

    c += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    copy(c.begin(), c.end(), state.begin() + 2);

    // moments part, RK4

    vector< complex<double> > l1, l2, l3, l4; // rmom
    vector< complex<double> > m1, m2, m3, m4; // pmom
    vector< complex<double> > TR, TP, dF, sigma, dFsigma;

    dF = F - F[s+s*edim] * matrixop::eye(edim);
    sigma.resize(edim * edim);
    for (int j(0); j < edim; ++j) {
        for (int k(0); k < edim; ++k) {
            sigma[j+k*edim] = c[j] * conj(c[k]);
        }
    }
    dFsigma = matrixop::anticommutator(dF, sigma);
    

    TR = -zI * matrixop::commutator(Emat, rmom) - v * matrixop::commutator(dc, rmom) + pmom / mass;
    TP = -zI * matrixop::commutator(Emat, pmom) - v * matrixop::commutator(dc, pmom) + 0.5 * dFsigma;
    l1 = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
    m1 = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));

    TR = -zI * matrixop::commutator(Emat, rmom + 0.5 * l1) - v * matrixop::commutator(dc, rmom + 0.5 * l1) + (pmom + 0.5 * m1) / mass;
    TP = -zI * matrixop::commutator(Emat, pmom + 0.5 * m1) - v * matrixop::commutator(dc, pmom + 0.5 * m1) + 0.5 * dFsigma;
    l2 = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
    m2 = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));

    TR = -zI * matrixop::commutator(Emat, rmom + 0.5 * l2) - v * matrixop::commutator(dc, rmom + 0.5 * l2) + (pmom + 0.5 * m2) / mass;
    TP = -zI * matrixop::commutator(Emat, pmom + 0.5 * m2) - v * matrixop::commutator(dc, pmom + 0.5 * m2) + 0.5 * dFsigma;
    l3 = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
    m3 = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));

    TR = -zI * matrixop::commutator(Emat, rmom + l3) - v * matrixop::commutator(dc, rmom + l3) + (pmom + m3) / mass;
    TP = -zI * matrixop::commutator(Emat, pmom + m3) - v * matrixop::commutator(dc, pmom + m3) + 0.5 * dFsigma;
    l4 = dt * (TR - TR[s+s*edim] * matrixop::eye(edim));
    m4 = dt * (TP - TP[s+s*edim] * matrixop::eye(edim));

    rmom += l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;
    pmom += m1 / 6.0 + m2 / 3.0 + m3 / 3.0 + m4 / 6.0;

    copy(rmom.begin(), rmom.end(), state.begin() + 5);
    copy(pmom.begin(), pmom.end(), state.begin() + 9);


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
    double dE = eva[1-s] - eva[s];
    // hop
    if (randomer::rand() < g) {
        double tmp = p * p - 2 * mass * dE;
        if (tmp > 0.0) {
            // hop accepted
            double pnew = sqrt(tmp);
            state[1] = pnew * (p < 0.0 ? -1 : 1); 
            state[6].real(1.0 - s); 
            return (s == 0) ? HOP_UP : HOP_DN;
        }
        else {
            // hop frustrated
            double F0d01 = (F[0+0*edim] * dc[0+1*edim]).real();
            double F1d01 = (F[1+1*edim] * dc[0+1*edim]).real();
            double vd01 = p / mass * dc[0+1*2].real();
            if (F0d01 * F1d01 < 0.0 and vd01 * F1d01 < 0.0)
            {
                // momentum reversal
                state[1] *= -1.0;
            }
            return HOP_FR; 
        }
    }   
    // hop rejected
    return HOP_RJ;
}

void afssh() {
    // initialize
    vector<state_t> state(Ntraj);
    for (auto& st : state) {
        init_state(st);
    }
    // statistics
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<int> hop_rec(Ntraj, 0);
    double n0d = 0.0, n1d = 0.0;
    double KE = 0.0, PE = 0.0;
    // cache
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    // main loop
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info_nume(state[itraj][0].real(), eva, dc, F, lastevt);
                // hopper
                if (enable_hop) {
                    int hopflag = hopper(state[itraj]);
                    switch (hopflag) {
                        case HOP_UP : hopup += 1.0; hop_rec[itraj] += 1; break;
                        case HOP_DN : hopdn += 1.0; hop_rec[itraj] += 1; break;
                        case HOP_FR : hopfr += 1.0; break;
                        case HOP_RJ : hoprj += 1.0; break;
                        default : break;
                    }
                }
                // propagate
                integrator(state[itraj], dt);
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
        }

    }
    // output hop info
    ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
    vector<int> hop_stat(50,0);
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        hop_stat[hop_rec[itraj]] += 1;
    }
    ioer::info("# hop stat: ", hop_stat);
}

void check() {

}


// -- MAIN -- //


int main(int argc, char** argv) {
    if (argc < 2) {
        ioer::info("use --help for detailed info");
        return -1;
    }
    else if (string(argv[1]) == "check") {
        check();
        return 0;
    }
    else {
        if (argparse(argc, argv) == false) {
            return 0;
        }
        randomer::seed(seed);
        timer::tic();
        afssh();
        ioer::info("# ", timer::toc());
    }
    return 0;
}
