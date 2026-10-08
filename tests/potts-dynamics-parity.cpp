// Compile against each backend to compare initialization and unit dynamics.
#include <iomanip>
#include <iostream>
#include "pnet.h"

class Probe : public PNet {
public:
    using PNet::PNet;
    double input = .7;
    void dump(std::ostream &output) {
        output << s0[0] << ' ' << theta0A[0] << ' ' << theta0B[0];
        for (int k = 0; k < params.S; ++k)
            output << ' ' << s[k] << ' ' << theta[k] << ' ' << r[k];
        output << '\n';
    }
protected:
    void compute_field(int i) override {
        for (int k = 0; k < params.S; ++k)
            h[i*params.S+k] = (k == 0 ? input : -.15*k);
    }
};

int main() {
    std::ostream output(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    output << std::setprecision(17);
    for (int states : {3, 7, 11}) for (double tau2 : {100., 400.}) {
        Potts_params params{};
        params.N = 50; params.Cm = 7; params.p = 49;
        params.S = states; params.a = .25; params.U = .1;
        params.w = 1.1; params.beta = 11.; params.T1 = 20.;
        params.T2 = tau2; params.T3A = 10.; params.T3B = 100000.;
        params.gammaA = .5;
        Probe net(&params, 10, 1);
        net.initialise();
        net.dump(output);
        for (int step = 0; step < 300; ++step) {
            // Drive both sustained activity and its decay.
            if (step == 100) net.input = -.4;
            net.update_unit(0, 0, nullptr, nullptr);
            net.dump(output);
        }
    }
}
