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

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;
using state_t = vector< complex<double> >;
const complex<double> zI(0.0, 1.0);

const double mass = 1.0;
const double kT = 9.5e-4;
const double omega = 3.5e-4;
const double Er = 0.0239;
const double M = sqrt(0.5 * Er * mass * omega * omega);
const double g = 2.0 * M / mass / omega / omega;
const double dG0 = -0.015;
const double fric_gamma = 0.0024; // 2 * mass * omega;

double V = 1.49e-5;
double init_x = 0.0;
double init_px = 0.0;
double init_s = 0.0;

const double sigma_x = 0.0;
//const double sigma_x = sqrt(kT / mass / omega / omega); 
const double sigma_px = sqrt(kT * mass);

int Nstep = 40000;
double dt = 2.5;
int output_step = 100;
int Ntraj = 5000;
bool enable_hop = true;

vector< complex<double> > lastevt;
vector<double> eva(2);
vector<double> Fx(2);
vector< complex<double> > dcx(4);
double random_force = 0.0;
string integrator = "vv";
int seed = 0;

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("V", po::value<double>(&V), "coupling magnitude V")
        ("init_x", po::value<double>(&init_x), "init x")
        ("init_px", po::value<double>(&init_px), "potential para init_px")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("integrator", po::value<string>(&integrator), "integrator to use. vv or lgv")
        ("seed", po::value<int>(&seed), "random seed. Default 0")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    misc::crasher::confirm(integrator == "vv" or integrator == "lgv", "Invalid integrator");
    return true;
}

vector< complex<double> > cal_H(const double x) {
    vector< complex<double> > H(4);
    H[0+0*2] = 0.5 * mass * omega * omega * x * x;
    H[1+1*2] = 0.5 * mass * omega * omega * (x - g) * (x - g) + dG0;
    H[1+0*2] = V;
    H[0+1*2] = conj(H[1+0*2]);
    return H;
}

vector< complex<double> > cal_nablaH(const double x) {
    vector< complex<double> > nablaH(4);
    nablaH[0+0*2] = mass * omega * omega * x;
    nablaH[1+1*2] = mass * omega * omega * (x - g);
    nablaH[1+0*2] = 0.0;
    nablaH[0+1*2] = 0.0;
    return nablaH;
}

void init_state(state_t& state) {
    state.resize(5, matrixop::ZEROZ);
    // state = x, vx, c0, c1, s
    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_px, sigma_px) / mass); 
    state[2].real(sqrt(1.0 - init_s));
    state[3].real(sqrt(init_s));
    state[4].real((randomer::rand() < init_s) ? 1.0 : 0.0);
}

void cal_info_nume(const double x)
{
    // nume
    vector< complex<double> > evt;
    matrixop::hdiag(cal_H(x), eva, evt);
    // correct phase
    if (not lastevt.empty()) {
        auto tmp = matrixop::matCmat(lastevt, evt, 2);
        for (int j = 0; j < 2; ++j) {
            complex<double> eip = tmp[j+j*2] / abs(tmp[j+j*2]);
            for (int k = 0; k < 2; ++k) {
                evt[k+j*2] /= eip;
            }
        }
    }
    // F, dc
    dcx = matrixop::matCmatmat(evt, cal_nablaH(x), evt, 2, 2);
    Fx.assign(2, 0.0);
    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            if (j == k) {
                Fx[j] = -dcx[j+k*2].real();
                dcx[j+k*2] = 0.0;
            }
            else {
                dcx[j+k*2] /= (eva[k] - eva[j]);
            }
        }
    }
    lastevt = move(evt);
}

void cal_info(const double x)
{
    cal_info_nume(x);
}

void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, v, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    int s = static_cast<int>(state[4].real());

    // cal_info(x); // <- make program slower but prevent total energy to drift
    // state_dot
    state_dot[0] = vx;
    state_dot[1] = (Fx[s] - fric_gamma * vx + random_force) / mass;
    state_dot[2] = -zI * c[0] * eva[0] - c[1] * (vx * dcx[0+1*2]) - c[0] * (vx * dcx[0+0*2]);
    state_dot[3] = -zI * c[1] * eva[1] - c[0] * (vx * dcx[1+0*2]) - c[1] * (vx * dcx[1+1*2]);
    state_dot[4] = matrixop::ZEROZ;
}

void integrator_Langevin(state_t& state, const double dt) {
    static bool firsttime = true;
    if (firsttime) {
        ioer::info("# INTEGRATOR: LGV");
        firsttime = false;
    }
    // state = x, vx, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    int s = static_cast<int>(state[4].real());

    // ele part, RK4

    complex<double> k1,k2,k3,k4;
    complex<double> l1,l2,l3,l4;
    k1 = dt * (-zI * c[0] * eva[0] - c[1] * (vx * dcx[0+1*2]));
    l1 = dt * (-zI * c[1] * eva[1] - c[0] * (vx * dcx[1+0*2]));
    k2 = dt * (-zI * (c[0] + 0.5 * k1) * eva[0] - (c[1] + 0.5 * l1) * (vx * dcx[0+1*2]));
    l2 = dt * (-zI * (c[1] + 0.5 * l1) * eva[1] - (c[0] + 0.5 * k1) * (vx * dcx[1+0*2]));
    k3 = dt * (-zI * (c[0] + 0.5 * k2) * eva[0] - (c[1] + 0.5 * l2) * (vx * dcx[0+1*2]));
    l3 = dt * (-zI * (c[1] + 0.5 * l2) * eva[1] - (c[0] + 0.5 * k2) * (vx * dcx[1+0*2]));
    k4 = dt * (-zI * (c[0] + k3) * eva[0] - (c[1] + l3) * (vx * dcx[0+1*2]));
    l4 = dt * (-zI * (c[1] + l3) * eva[1] - (c[0] + k3) * (vx * dcx[1+0*2]));

    state[2] += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    state[3] += l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;


    // nuclear part, Langevin integrator
    const vector<double> randnum = randomer::vnormal(2);
    const double sigma = sqrt(2 * kT * fric_gamma) / mass;

    const double a = Fx[s] / mass;
    const double C = 0.5 * dt * dt * (a - fric_gamma / mass * vx) + 
                    sigma * pow(dt, 1.5) * (0.5 * randnum[0] + 0.5 / sqrt(3) * randnum[1]);
    const double x_new = x + dt * vx + C;

    cal_info(x_new);

    const double a_new = Fx[s] / mass;
    const double vx_new = vx + 0.5 * dt * (a_new + a) - dt * fric_gamma / mass * vx + sigma * sqrt(dt) * randnum[0] - fric_gamma / mass * C;

    state[0].real(x_new);
    state[1].real(vx_new);
}

void integrator_VV(state_t& state, const double dt) {
    static bool firsttime = true;
    if (firsttime) {
        ioer::info("# INTEGRATOR: VV");
        firsttime = false;
    }
    // state = x, vx, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    int s = static_cast<int>(state[4].real());

    // ele part, RK4

    complex<double> k1,k2,k3,k4;
    complex<double> l1,l2,l3,l4;
    k1 = dt * (-zI * c[0] * eva[0] - c[1] * (vx * dcx[0+1*2]));
    l1 = dt * (-zI * c[1] * eva[1] - c[0] * (vx * dcx[1+0*2]));
    k2 = dt * (-zI * (c[0] + 0.5 * k1) * eva[0] - (c[1] + 0.5 * l1) * (vx * dcx[0+1*2]));
    l2 = dt * (-zI * (c[1] + 0.5 * l1) * eva[1] - (c[0] + 0.5 * k1) * (vx * dcx[1+0*2]));
    k3 = dt * (-zI * (c[0] + 0.5 * k2) * eva[0] - (c[1] + 0.5 * l2) * (vx * dcx[0+1*2]));
    l3 = dt * (-zI * (c[1] + 0.5 * l2) * eva[1] - (c[0] + 0.5 * k2) * (vx * dcx[1+0*2]));
    k4 = dt * (-zI * (c[0] + k3) * eva[0] - (c[1] + l3) * (vx * dcx[0+1*2]));
    l4 = dt * (-zI * (c[1] + l3) * eva[1] - (c[0] + k3) * (vx * dcx[1+0*2]));

    state[2] += k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    state[3] += l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;


    // nuclear part, VV
    random_force = randomer::normal(0.0, sqrt(2 * fric_gamma * kT / dt));
    vx += 0.5 * dt * (Fx[s] - fric_gamma * vx + random_force) / mass;
    x += dt * vx;
    cal_info(x);
    vx += 0.5 * dt * (Fx[s] - fric_gamma * vx + random_force) / mass;

    state[0].real(x);
    state[1].real(vx);
}

int hopper(state_t& state) {
    // state = x, vx, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    int s = static_cast<int>(state[4].real());
    // calc hop prob
    double g = -2 * dt * (c[s] * conj(c[1-s]) * (vx * dcx[1-s+s*2])).real() / (c[s] * conj(c[s])).real();
    double dE = eva[1-s] - eva[s];
    // hop
    if (randomer::rand() < g) {
        double tmp = vx * vx - 2 * dE / mass;
        if (tmp > 0.0) {
            // hop accepted, modify p and s
            double vnew = sqrt(tmp);
            state[1] = vnew * (vx < 0.0 ? -1 : 1); 
            state[4].real(1.0 - s); 
            return (s == 0) ? HOP_UP : HOP_DN;
        }
        else {
            // hop frustrated, momentum reversal
            double F0d01 = Fx[0] * dcx[0+1*2].real();
            double F1d01 = Fx[1] * dcx[0+1*2].real();
            double vxd01 = vx * dcx[0+1*2].real();
            if (F0d01 * F1d01 < 0.0 and vxd01 * F1d01 < 0.0)
            {
                state[1].real(-vx);
            }
            return HOP_FR; 
        }
    }   
    // hop rejected
    return HOP_RJ;
}

bool check_end(const state_t& state) {
    return false;
}

void fssh() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        init_state(state[itraj]);
    }
    // main loop
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    // statistics
    double n0d = 0.0, n1d = 0.0;
    double KE = 0.0, PE = 0.0;
    vector<int> hop_rec(Ntraj, 0);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info(state[itraj][0].real());
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
                /*
                random_force = randomer::normal(0.0, sqrt(2 * fric_gamma * kT / dt));
                rk4.do_step(sys, state[itraj], istep * dt, dt);
                */
                if (integrator == "vv") {
                    integrator_VV(state[itraj], dt);
                }
                else if (integrator == "lgv") {
                    integrator_Langevin(state[itraj], dt);
                }
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis

            n0d = n1d = 0.0;
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0d, &n1d, &KE, &PE](const state_t& st) { 
                        double x = st[0].real();
                        double vx = st[1].real();
                        vector< complex<double> > c { st[2], st[3] };
                        int s = static_cast<int>(st[4].real());

                        // diabatic population -- Method 3 from Landry 2013 communication 
                        vector<double> eva;
                        vector< complex<double> > U;
                        matrixop::hdiag(cal_H(x), eva, U);

                        n0d += pow(abs(U[0+s*2]), 2) + 2 * (U[0+0*2] * c[0] * conj(c[1]) * conj(U[0+1*2])).real();
                        n1d += pow(abs(U[1+s*2]), 2) + 2 * (U[1+0*2] * c[0] * conj(c[1]) * conj(U[1+1*2])).real();

                        // energy
                        KE += 0.5 * mass * vx * vx;
                        PE += eva[s]; 
                    });

            // output
            if (istep == 0) {
                // para & header
                ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step,
                            " mass = ", mass, " kT = ", kT, 
                            " omega = ", omega, " g = ", g,  " Er = ", Er, " dG0 = ", dG0, 
                            " fric_gamma = ", fric_gamma, 
                            " V = ", V, 
                            " init_x = ", init_x, " init_px = ", init_px, 
                            " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                            " init_s = ", init_s,
                            " seed = ", seed
                        );
                ioer::tabout('#', "t", "n0d", "n1d", "KE", "PE", "Etot");
            }
            n0d /= Ntraj;
            n1d /= Ntraj;
            KE /= Ntraj;
            PE /= Ntraj;
            ioer::tabout("", istep * dt, n0d, n1d, KE, PE, KE + PE);

            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                break;
            }
        }

    }
    // hop info
    ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
    vector<int> hop_stat(50,0);
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        hop_stat[hop_rec[itraj]] += 1;
    }
    ioer::info("# hop stat: ", hop_stat);
}

void check() {
    // check surf
    for (double x = -1000; x < 1500; x += 10.0) {
        cal_info(x);
        ioer::tabout(x, eva, Fx, abs(dcx) / 10);
    }

    // check init
    state_t st;
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        init_state(st);
        ioer::tabout(real(st));
    }
}

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
        fssh();
        ioer::info("# ", timer::toc());
    }
    return 0;
}
