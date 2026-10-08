// Browser runner for simcode_pilot's frozen 49-memory topology study.
// Uses Pottsviz's streaming architecture and the study's active equations.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "functions.h"
#include "pnet.h"
#include "rand_gen.h"

std::vector<double> load_topology(const char *path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Missing topology CSV");
    std::vector<double> weights;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::string cell;
        std::vector<double> values;
        double degree = 0;
        while (std::getline(row, cell, ',')) {
            double value = std::stod(cell);
            if (value != 0. && value != 1.) throw std::invalid_argument("Topology must be binary");
            values.push_back(value); degree += value;
        }
        if (values.size() != 49 || degree == 0) throw std::invalid_argument("Invalid topology row");
        for (double value : values) weights.push_back(value / degree);
    }
    if (weights.size() != 49*49) throw std::invalid_argument("Expected 49 topology rows");
    return weights; // rows: frontal target; columns: posterior source
}

int main(int argc, char *argv[]) {
    std::ostream output(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc != 20) throw std::invalid_argument("Expected topology wP wF SP SF tauP tauF lambdaP lambdaF U beta a T1 density patternSeed runtimeSeed cue N steps");
        double wP = std::stod(argv[2]), wF = std::stod(argv[3]);
        int SP = std::stoi(argv[4]), SF = std::stoi(argv[5]);
        double tauP = std::stod(argv[6]), tauF = std::stod(argv[7]);
        double lambdaP = std::stod(argv[8]), lambdaF = std::stod(argv[9]);
        double U = std::stod(argv[10]), beta = std::stod(argv[11]), a = std::stod(argv[12]);
        double T1 = std::stod(argv[13]), density = std::stod(argv[14]);
        int patternSeed = std::stoi(argv[15]), runtimeSeed = std::stoi(argv[16]);
        int cue = std::stoi(argv[17]), N = std::stoi(argv[18]), steps = std::stoi(argv[19]);
        for (double value : {wP, wF, tauP, tauF, lambdaP, lambdaF, U, beta, a, T1, density})
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite parameter");
        if (wP < .6 || wP > 2 || wF < .6 || wF > 2 || SP < 3 || SP > 11 || SF < 3 || SF > 11 ||
            tauP < 100 || tauP > 800 || tauF < 100 || tauF > 1600 || tauF <= tauP ||
            lambdaP < 0 || lambdaP > 1 || lambdaF < 0 || lambdaF > 1 ||
            U < 0 || U > .6 || beta < 1 || beta > 21 || a < .1 || a > .4 ||
            T1 < 5 || T1 > 35 || density < .05 || density > .25 ||
            patternSeed < 1 || patternSeed > 2147482000 || runtimeSeed < 1 ||
            cue < 0 || cue >= 49 || N < 50 || N > 500 || steps < 1 || steps > 5000)
            throw std::invalid_argument("Parameters outside supported bounds");
        Potts_params posterior{}, frontal{};
        for (Potts_params *params : {&posterior, &frontal}) {
            params->N = N; params->Cm = static_cast<int>(N*density); params->p = 49;
            params->a = a; params->U = U; params->beta = beta; params->T1 = T1;
            // Inhibition dynamics are disabled in the current pilot study.
            params->T3A = 10.; params->T3B = 100000.; params->gammaA = .5;
        }
        posterior.w = wP; posterior.S = SP; posterior.T2 = tauP;
        frontal.w = wF; frontal.S = SF; frontal.T2 = tauF;
        std::vector<int> xiP(N*49), xiF(N*49);
        // The helper seeds srand48(seed+1), matching the pilot pattern generator.
        make_random_memory(N, 49, SP, a, 0, xiP.data(), patternSeed-1);
        make_random_memory(N, 49, SF, a, 0, xiF.data(), patternSeed);
        PNet pNet(&posterior, 1000, patternSeed);
        pNet.make_J_assoc(49, xiP.data(), a, lambdaP);
        PNet fNet(&frontal, 1000, patternSeed+1000);
        fNet.make_J_assoc(49, xiF.data(), a, lambdaF);
        auto p2f = load_topology(argv[1]);
        std::vector<double> f2p(p2f.size());
        for (int p = 0; p < 49; ++p)
            for (int f = 0; f < 49; ++f) f2p[p*49+f] = p2f[f*49+p];
        std::vector<double> Jf2p(static_cast<size_t>(N)*N*SP*SF), Jp2f(Jf2p.size());
        Network_runner manager;
        if (lambdaP == 1.) manager.make_null_connection(&pNet, &fNet, Jf2p.data());
        else manager.make_Hebb_connection(&pNet, &fNet, Jf2p.data(), xiP.data(), xiF.data(), lambdaP, f2p.data());
        if (lambdaF == 1.) manager.make_null_connection(&fNet, &pNet, Jp2f.data());
        else manager.make_Hebb_connection(&fNet, &pNet, Jp2f.data(), xiF.data(), xiP.data(), lambdaF, p2f.data());
        manager.Runs = steps;
        srand48(runtimeSeed); rlxd_init(2, runtimeSeed);
        std::ostringstream sequence;
        output << std::fixed << std::setprecision(4);
        manager.run_two_nets(&pNet, &fNet, Jf2p.data(), Jp2f.data(), xiP.data(), xiF.data(), sequence, output, cue, 1);
    } catch (const std::exception &error) {
        std::cerr << "Simulation error: " << error.what() << '\n';
        return 2;
    }
}
