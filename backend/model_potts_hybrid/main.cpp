// Browser runner for the original Ryom frontal-posterior dynamics.
// Each posterior memory mu is associated reciprocally with frontal memory mu.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "functions.h"
#include "pnet.h"
#include "rand_gen.h"

void set_params(Potts_params &params, int N) {
    params.N = N;
    params.Cm = 50;
    params.p = 50;
    params.S = 7; params.a = 0.25; params.U = 0.1;
    params.T1 = 20.; params.T2 = 200.; params.T3A = 10.;
    params.T3B = 100000.; params.beta = 11.;
    params.w = 1.1; params.gammaA = 0.5;
}

int main(int argc, char *argv[]) {
    // Separate numerical stdout from the inherited progress diagnostics.
    std::ostream output(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc < 8 || argc > 10) throw std::invalid_argument("Expected wP wF SP SF tauP tauF lambda [steps=5000] [N=256]");
        double w_p = std::stod(argv[1]), w_f = std::stod(argv[2]);
        int S_p = std::stoi(argv[3]), S_f = std::stoi(argv[4]);
        double tau_p = std::stod(argv[5]), tau_f = std::stod(argv[6]);
        double lambda = std::stod(argv[7]);
        int steps = argc > 8 ? std::stoi(argv[8]) : 5000;
        int N = argc > 9 ? std::stoi(argv[9]) : 256;
        if (!std::isfinite(w_p) || !std::isfinite(w_f) || w_p < 0. || w_p > 2. || w_f < 0. || w_f > 2. ||
            S_p < 3 || S_p > 11 || S_f < 3 || S_f > 11 ||
            !std::isfinite(tau_p) || !std::isfinite(tau_f) || tau_p < 100. || tau_p > 800. || tau_f < 100. || tau_f > 1600. ||
            !std::isfinite(lambda) || lambda < 0. || lambda > 1. || steps < 1 || steps > 5000 || N <= 50 || N > 500)
            throw std::invalid_argument("Parameters outside supported bounds");
        Potts_params posterior, frontal;
        set_params(posterior, N); set_params(frontal, N);
        posterior.w = w_p; posterior.S = S_p; posterior.T2 = tau_p;
        frontal.w = w_f; frontal.S = S_f; frontal.T2 = tau_f;
        int p = posterior.p;
        std::vector<int> xi_p(N*p), xi_f(N*p);
        // This inherited helper seeds srand48(seed+1). Arguments 0 and 1
        // reproduce the original pattern files generated with seeds 1 and 2.
        make_random_memory(N, p, S_p, posterior.a, 0, xi_p.data(), 0);
        make_random_memory(N, p, S_f, frontal.a, 0, xi_f.data(), 1);
        PNet p_net(&posterior, 1000, 1);
        p_net.make_J_assoc(p, xi_p.data(), posterior.a, lambda);
        PNet f_net(&frontal, 1000, 1001);
        f_net.make_J_assoc(p, xi_f.data(), frontal.a, lambda);
        std::vector<int> pairs(p*p, 0);
        for (int mu = 0; mu < p; ++mu) pairs[mu*p+mu] = 1;
        std::vector<double> Jf2p(static_cast<size_t>(N)*N*S_p*S_f);
        std::vector<double> Jp2f(Jf2p.size());
        Network_runner manager;
        if (lambda == 1.) {
            manager.make_null_connection(&p_net, &f_net, Jf2p.data());
            manager.make_null_connection(&f_net, &p_net, Jp2f.data());
        } else {
            manager.make_Hebb_connection(&p_net, &f_net, Jf2p.data(), xi_p.data(), xi_f.data(), lambda, pairs.data());
            manager.make_Hebb_connection(&f_net, &p_net, Jp2f.data(), xi_f.data(), xi_p.data(), lambda, pairs.data());
        }
        manager.Runs = steps;
        srand48(1); rlxd_init(2, 1);
        // The inherited runner emits posterior then frontal at each snapshot.
        std::ostringstream sequence;
        output << std::fixed << std::setprecision(4);
        manager.run_two_nets(&p_net, &f_net, Jf2p.data(), Jp2f.data(), xi_p.data(), xi_f.data(), sequence, output, 0, 1);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Simulation error: " << error.what() << '\n';
        return 2;
    }
}
