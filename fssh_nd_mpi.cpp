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
#include "marcus_potential.hpp"

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};

using namespace std;
namespace po = boost::program_options;
using boost::math::erf;
using state_t = vector< complex<double> >;

int ndim = 2;
int edim = 2;

double kT = 0.01;

vector<double> mass { 2000.0, 2000.0};
vector<double> init_r { 0.0, 0.0};
vector<double> init_p { 0.0, 0.0};
vector<double> sigma_r { 0.0, 0.0 };
vector<double> sigma_p = sqrt(kT * mass);
vector<double> init_s { 1.0, 0.0 };
vector<double> potential_params;
vector<double> fric_gamma = 2.0 * mass * param_omega;


double dt = 0.1;
int Nstep = 10000;
int output_step = 100;
int Ntraj = 2000;
int seed = 0;

bool testmode = false;

vector<double> eva;
vector< vector< complex<double> > > dc;
vector< vector< complex<double> > > F;
vector< complex<double> > lastevt;

bool argparse(int argc, char** argv) 
{
    /*
     * parse input arguments
     */
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("Ntraj", po::value<decltype(Ntraj)>(&Ntraj), "# traj")
        ("Nstep", po::value<decltype(Nstep)>(&Nstep), "# step")
        ("output_step", po::value<decltype(output_step)>(&output_step), "# step for output")
        ("dt", po::value<decltype(dt)>(&dt), "single time step")
        ("kT", po::value<decltype(kT)>(&kT), "temperature kT")
        ("mass", po::value<decltype(mass)>(&mass)->multitoken(), "mass vector")
        ("init_r", po::value<decltype(init_r)>(&init_r)->multitoken(), "init_r vector")
        ("init_p", po::value<decltype(init_p)>(&init_p)->multitoken(), "init_p vector")
        ("sigma_r", po::value<decltype(sigma_r)>(&sigma_r)->multitoken(), "sigma_r vector")
        ("sigma_p", po::value<decltype(sigma_p)>(&sigma_p)->multitoken(), "sigma_p vector")
        ("init_s", po::value<decltype(init_s)>(&init_s)->multitoken(), "init_s vector")
        ("potential_params", po::value<decltype(potential_params)>(&potential_params)->multitoken(), "potential_params vector")
        ("testmode", po::value<decltype(testmode)>(&testmode), "run test mode")
        ("seed", po::value<decltype(seed)>(&seed), "random seed")
        ;
    po::variables_map vm; 
    //po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    ndim = mass.size();
    edim = init_s.size();

    if (not potential_params.empty()) {
        set_potenial_params(potential_params);
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

void init_state(state_t& state, const vector<double>& init_r, const vector<double>& init_p, 
        const vector<double>& mass, const vector<double> init_s) 
{
    /*
     * initialize state information for a single trajectory
     */

    // check info
    const int ndim = mass.size();
    const int edim = init_s.size();
    misc::crasher::confirm(init_r.size() == ndim and init_p.size() == ndim,
                                "init_state: mass, init_r, init_p must have identical sizes");
    misc::crasher::confirm(std::all_of(init_s.begin(), init_s.end(), [](double si) { return si >= 0.0; }), 
                                "init_state: all init_s elements must be non-negative" );

    // state -> (r, p, c, s)
    const int Lstate = ndim * 2 + edim + 1; 
    state.resize(Lstate, matrixop::ZEROZ);

    // init nuclear DoF (r, p)
    for (int i(0); i < ndim; ++i) {
        state[i].real(randomer::normal(init_r[i], sigma_r[i]));
        state[ndim + i].real(randomer::normal(init_p[i], sigma_p[i]));
    }

    // init electronic DoF (c, s)
    vector<double> s_normalized = init_s / sum(init_s);
    for (int i(0); i < edim; ++i) {
        state[ndim * 2 + i] = sqrt(s_normalized[i]);
    }
    state[ndim * 2 + edim].real(randomer::discrete(s_normalized.begin(), s_normalized.end()));

}

void VV_integrator(state_t& state, const vector<double>& mass, const double t, const double dt) 
{
    /*
     *  VV integrator
     */

    // check
    const int ndim = mass.size();
    const int edim = state.size() - 2 * ndim - 1;

    // extract information
    vector<double> r(ndim);
    vector<double> p(ndim);
    vector< complex<double> > c(edim);
    int s;
    for (int i(0); i < ndim; ++i) {
        r[i] = state[i].real();
        p[i] = state[ndim + i].real();
    }
    for (int i(0); i < edim; ++i) {
        c[i] = state[ndim * 2 + i];
    }
    s = static_cast<int>(state[ndim * 2 + edim].real());

    const double half_dt = 0.5 * dt;

    // electronic part -- RK4

    cal_info_nume(r, eva, dc, F, lastevt);

    vector< complex<double> > rk4_mat(edim * edim, 0.0);
    for (int j(0); j < edim; ++j) {
        for (int k(0); k < edim; ++k) {
            for (int i(0); i < ndim; ++i) {
                rk4_mat[j+k*edim] -= p[i] * dc[i][j+k*edim] / mass[i];
            }
        }
    }
    for (int j(0); j < edim; ++j) {
        rk4_mat[j+j*edim] -= matrixop::IMAGIZ * eva[j];
    }

    vector< complex<double> > k1, k2, k3, k4;
    k1 = dt * matrixop::matmat(rk4_mat, c, edim);
    k2 = dt * matrixop::matmat(rk4_mat, c + 0.5 * k1, edim);
    k3 = dt * matrixop::matmat(rk4_mat, c + 0.5 * k2, edim);
    k4 = dt * matrixop::matmat(rk4_mat, c + k3, edim);
    c += 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    // nuclear part -- VV integrator
    vector<double> force(ndim);
    vector<double> random_force(ndim);
    for (int i(0); i < ndim; ++i) {
        random_force[i] = randomer::normal(0.0, sqrt(2.0 * fric_gamma[i] * kT / dt));
    }

    complex<double> vdotdc(0.0, 0.0);
    for (int k(0); k < ndim; ++k) {
        vdotdc += p[k] / mass[k] * dc[k][(1-s)+s*edim];
    }

    for (int i(0); i < ndim; ++i) {
        // adiabatic force
        force[i] = F[i][s + s * edim].real();
        // Berry force
        force[i] += 2 * (dc[i][s+(1-s)*edim] * vdotdc).imag();
        // friction & random force
        force[i] += -fric_gamma[i] * p[i] + random_force[i];
    }
    p += half_dt * force;

    r += dt * p / mass;

    cal_info_nume(r, eva, dc, F, lastevt);

    vdotdc = complex<double>(0.0, 0.0);
    for (int k(0); k < ndim; ++k) {
        vdotdc += p[k] / mass[k] * dc[k][(1-s)+s*edim];
    }

    for (int i(0); i < ndim; ++i) {
        // adiabatic force
        force[i] = F[i][s + s * edim].real();
        // Berry force
        force[i] += 2 * (dc[i][s+(1-s)*edim] * vdotdc).imag();
        // friction & random force
        force[i] += -fric_gamma[i] * p[i] + random_force[i];
    }
    p += half_dt * force;

    // update 
    for (int i(0); i < ndim; ++i) {
        state[i].real(r[i]);
        state[ndim + i].real(p[i]);
    }
    for (int i(0); i < edim; ++i) {
        state[ndim * 2 + i] = c[i];
    }
}

int hopper_2state(state_t& state, const vector<double>& mass) 
{
    /*
     * 2-state hopper 
     */

    // check
    const int ndim = mass.size();
    const int edim = state.size() - 2 * ndim - 1;
    misc::crasher::confirm(edim == 2, "hopper: this hopper can only be apllied to edim = 2");
    for (int i(0); i < ndim; ++i) {
        misc::crasher::confirm(mass[i] == mass[0], "hopper: this hopper can work w/ identical mass on each dim");
    }

    // extract information
    vector<double> r(ndim);
    vector<double> p(ndim);
    vector< complex<double> > c(edim);
    int s;
    for (int i(0); i < ndim; ++i) {
        r[i] = state[i].real();
        p[i] = state[ndim + i].real();
    }
    for (int i(0); i < edim; ++i) {
        c[i] = state[ndim * 2 + i];
    }
    s = static_cast<int>(state[ndim * 2 + edim].real());

    const int from = s;
    const int to = 1 - s;


    // calc hop prob
    complex<double> vd = 0.0;
    for (int i(0); i < ndim; ++i) {
        vd += p[i] / mass[i] * dc[i][to+from*edim];
    }
    double g = -2 * dt * (c[from] * conj(c[to]) * vd).real() / (c[from] * conj(c[from])).real();
    double dE = eva[to] - eva[from];


    // random number
    if (randomer::rand() < g) {
        // momentum-rescaling direction: (x-direction)
        vector<double> n(ndim, 0.0);

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
        const complex<double> eieta = exp(matrixop::IMAGIZ * eta);
        n[0] = (eieta * dc[0][from+to*edim]).real();
        n[1] = (eieta * dc[1][from+to*edim]).real();

        // hop
        if (norm(n) > 1e-40) {
            vector<double> pn = component(p, n);
            double pn_norm = norm(pn); 
            double tmp = pn_norm * pn_norm - 2 * mass[0] * dE; // masses along all dimension should be identical
            if (tmp > 0.0) {
                // hop accepted
                double pn_norm_new = sqrt(tmp);
                p += (pn_norm_new - pn_norm) / pn_norm * pn;

                // replace p & s
                for (int i(0); i < ndim; ++i) {
                    state[ndim + i].real(p[i]);
                }
                state[2 * ndim + edim] = to;
                return (from < to) ? HOP_UP : HOP_DN;
            }
            else {
                // hop frustrated
                return HOP_FR;
            }
        }
    }
    // hop rejected
    return HOP_RJ;
}


struct observer {
    /*
     * class for recording observables
     */

    public:
        int m_Nrec, m_irec;
        map< string, vector<double> > m_data_arr;
        const vector<string> m_keys { 
            "n0d", "n1d", 
            "KE", "PE"
        };

    public:
        observer(const int Nrec) 
            : m_Nrec(Nrec), m_irec(0)
        {
            for (const string& key : m_keys) {
                m_data_arr.insert( std::make_pair(key, vector<double>(m_Nrec, 0.0)) );
            }
        }

        ~observer() = default;

    public:
        void add_record(const vector<state_t>& states, const vector<double>& mass) { 
            /*
             * record a row of data
             */

            // check 
            const int ndim = mass.size();
            const int edim = states[0].size() - 2 * ndim - 1;

            misc::crasher::confirm(ndim == 2, "observer::add_record: ndim must be 2");
            misc::crasher::confirm(edim == 2, "observer::add_record: edim must be 2");
            misc::crasher::confirm(m_irec < m_Nrec, "observer::add_record: capacity reached");

            // population & energy
            double n0d = 0.0, n1d = 0.0;
            double KE = 0.0, PE = 0.0;

            for_each(states.begin(), states.end(), 
                    [&n0d, &n1d, &KE, &PE,
                     &ndim, &edim, &mass] (const state_t& st) { 
                        // extract info
                        vector<double> r(ndim);
                        vector<double> p(ndim);
                        vector< complex<double> > c(edim);
                        int s;

                        for (int i(0); i < ndim; ++i) {
                            r[i] = st[i].real();
                            p[i] = st[ndim + i].real();
                        }
                        for (int i(0); i < edim; ++i) {
                            c[i] = st[ndim * 2 + i];
                        }
                        s = static_cast<int>(st[ndim * 2 + edim].real());


                        // calc eva & evt
                        vector<double> eva;
                        vector< complex<double> > U;
                        matrixop::hdiag(cal_H(r), eva, U);

                        // population
                        n0d += pow(abs(U[0+s*2]), 2) + 2 * (U[0+0*2] * c[0] * conj(c[1]) * conj(U[0+1*2])).real();
                        n1d += pow(abs(U[1+s*2]), 2) + 2 * (U[1+0*2] * c[0] * conj(c[1]) * conj(U[1+1*2])).real();

                        // energy
                        for (int i(0); i < ndim; ++i) {
                            KE += 0.5 / mass[i] * p[i] * p[i];
                        }
                        PE += eva[s];
                    });

            m_data_arr["n0d"].at(m_irec) += n0d;
            m_data_arr["n1d"].at(m_irec) += n1d;
            m_data_arr["KE"].at(m_irec) = KE;
            m_data_arr["PE"].at(m_irec) = PE;

            // increase m_irec
            m_irec += 1;
        }

        void fill_unrecorded(const string& mode) { 
            /*
             * fill unrecorded data with the last record
             */
            if (m_irec > 0) {
                double val;
                for (const string& key : m_keys) {
                    if (mode == "zero") {
                        val = 0.0;
                    }
                    else if (mode == "last") {
                        val = m_data_arr[key][m_irec - 1];
                    }

                    fill(m_data_arr[key].begin() + m_irec, m_data_arr[key].end(), val);
                }
            }

            m_irec = m_Nrec;
        }

        void join_record() { 
            /*
             * collect all data to the master process
             */
            for (int r(1); r < MPIer::size; ++r) {
                if (MPIer::rank == r) {
                    for (const string& key : m_keys) {
                        MPIer::send(0, m_data_arr[key]);
                    }
                }
                else if (MPIer::master) {
                    vector<double> vbuf;
                    for (const string& key : m_keys) {
                        MPIer::recv(r, vbuf);
                        m_data_arr[key] += vbuf;
                    }
                }
                MPIer::barrier();
            }
            MPIer::barrier();
        }

        vector<string> get_keys() const {
            /*
             * get all keys
             */
            return m_keys;
        }

        double get_record(const string& key, int irec) const { 
            /*
             * get a piece of data 
             */
            auto it = m_data_arr.find(key);
            misc::crasher::confirm(it != m_data_arr.end(), "observer::get: the key does not exist.");
            return it->second.at(irec);
        }

        map<string, double> get_record(int irec) const { 
            /*
             * get a row of data 
             */
            map<string, double> row;
            for (const string& key : m_keys) {
                row.insert( make_pair(key, m_data_arr.at(key).at(irec)) );
            }
            return row;
        }
};

bool check_end(const state_t& state) {
    return false;
}

template <typename ...Params>
void integrator(Params&&... params)
{
    VV_integrator(std::forward<Params>(params)...);
}

template <typename ...Params>
int hopper(Params&&... params)
{
    return hopper_2state(std::forward<Params>(params)...);
}

void fssh_nd_mpi() {
    /*
     * n-dimensional parallel FSSH
     */

    // assign jobs
    const vector<int> my_jobs = MPIer::assign_job(Ntraj);
    const int my_Ntraj = my_jobs.size();


    // initialize trajectories
    vector<state_t> state(my_Ntraj);
    for (int itraj(0); itraj < my_Ntraj; ++itraj) {
        init_state(state[itraj], init_r, init_p, mass, init_s);
    }


    // hop statistics
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    vector<double> hop_count_summary(50, 0.0);

    // recorders
    int Nrec = Nstep / output_step;
    observer obs(Nrec);


    //  ----  MAIN LOOP  ---- //


    // last evt save
    vector< vector< complex<double> > > lastevt_save(my_Ntraj);

    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < my_Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // integrate t -> t + dt
                integrator(state[itraj], mass, istep * dt, dt);
                // hopper
                int hopflag = hopper(state[itraj], mass);
                switch (hopflag) {
                    case HOP_UP : { hopup += 1.0; hop_count[itraj] += 1.0; break; }
                    case HOP_DN : { hopdn += 1.0; hop_count[itraj] += 1.0; break; }
                    case HOP_FR : { hopfr += 1.0; break; }
                    case HOP_RJ : { hoprj += 1.0; break; }
                    default : break;
                }
                // save lastevt
                lastevt_save[itraj] = move(lastevt);
            }
        }

        if (istep % output_step == 0) {
            // data record
            obs.add_record(state, mass);

            // check end
            const bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                // fill the rest
                obs.fill_unrecorded("last");
                break;
            }
        }
    }
    MPIer::barrier();


    // ----  COLLECT DATA  ---- //
    

    obs.join_record();
    MPIer::barrier();


    // ----  PROCESS & COLLECT HOP STATISTICS DATA  ---- //


    for_each(hop_count.begin(), hop_count.end(),
            [&hop_count_summary](double x) { hop_count_summary[static_cast<int>(x)] += 1.0; });
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
            MPIer::recv(r, vbuf); hop_count_summary += vbuf;
        }
        MPIer::barrier();
    }
    MPIer::barrier();


    // ----  OUTPUT  ---- //


    if (MPIer::master) {
        // Output parameters
        output_potential_params();
        ioer::info("# FSSH parameters: ", 
                " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, 
                " mass = ", mass, 
                " init_r = ", init_r, " init_p = ", init_p, 
                " sigma_r = ", sigma_r, " sigma_p = ", sigma_p, 
                " init_s = ", init_s,
                " fric_gamma = ", fric_gamma,
                " output_step = ", output_step
                );
        // Output header
        ioer::tabout(
                "#", "t", 
                "n0d", "n1d", 
                "KE", "PE",
                "Etot");
        // Output data
        map<string, double> row;
        for (int irec = 0; irec < Nrec; ++irec) {
            row = obs.get_record(irec);

            // average over all traj
            for (const string& key : obs.get_keys()) {
                row[key] /= Ntraj;
            }

            ioer::tabout(
                    "", irec * output_step * dt, 
                    row["n0d"], row["n1d"],
                    row["KE"], row["PE"],
                    row["KE"] + row["PE"]
                    );
        }

        // Output hop info
        ioer::info("# hop statistics: ");
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
    }
    MPIer::barrier();
}

int test() {
    vector<double> r(3, 0.0);
    int N = 500;
    vector<double> xarr = linspace(-10.0, 30.0, N);
    vector<double> yarr = linspace(-8.0, 8.0, N);

    for (int ix(0); ix < N; ++ix) {
        for (int iy(0); iy < N; ++iy) {
            double x = xarr[ix];
            double y = yarr[iy];
            r[0] = x;
            r[1] = y;

            cal_info_nume(r, eva, dc, F, lastevt);
            ioer::tabout(x, y, 
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

int main(int argc, char** argv) {

    MPIer::setup();
    if (argc < 2) {
        if (MPIer::master) ioer::info("use --help for detailed info");
    }
    else {
        // parse args
        if (argparse(argc, argv) == false) {
            return 0;
        }

        // run 
        randomer::seed(MPIer::assign_random_seed(seed));
        if (MPIer::master) timer::tic();

        if (testmode) {
            if (MPIer::master) ioer::info("# TEST MODE ON");
            test();
        }
        else {
            fssh_nd_mpi();
        }

        if (MPIer::master) ioer::info("# ", timer::toc());
    }
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
