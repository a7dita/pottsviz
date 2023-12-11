
/******* Kwang Il Ryom, 16/10/2023
* A latching network of continous Potts units

* I use flattend array: xi[mu][i] = xi[mu*N+i],                                             *
*   J[i][j][k][l] = J[i*params.Cm*S*S + j*S + k*S*params.Cm + l]                *
* Make 3 files, named "xxx_cue0",  "xxx_cue1", "xxx_cue2"

************************************************/
#include <iostream>
#include <fstream>
#include "global.h"
#include <stdlib.h>
#include <cassert>
//********************************************************************
// Put global variables here, though discouraged to do so:

bool show_sim; // whether or not to see simulations on screen

// ------------------------------------------------------------------------
void set_params(Potts_params *params){
    /*
    set default values for paramters
    */
    (*params).N = 500;     // number of units in a net
    (*params).Cm = 75;     // mean connections per neuron
    (*params).p = 100;      // number of patterns
    (*params).S = 7;        // number of Potts states
    (*params).a = 0.25;    // sparsity
    (*params).U = 0.1;      // global threshold
    (*params).T1 = 20.;     // integration of incoming fields
    (*params).T2 = 200.;  // adaptation time scale
    (*params).T3A = 10.;    //fast inhibition
    (*params).T3B = 100000.; //slow inhibition
    (*params).beta = 11.;   // inverse temperature
    (*params).w = 1.1;      // self-reinforcement term
    (*params).gammaA = 0.5; // fast inhibition
}
// save things to a file to test, only used for testing version.
void save_data(double * vec, int n, std::ofstream & outfile){
    int i;
    for(i=0;i<n;i++){
        outfile << *(vec+i) << "\n";
    }
}
// ------------------------------------------------------------------------
int main(int argc, char *argv[]){

    
    int rand_seed = 1990; // random seed for srand48() and rlxd_init()
    
    int S;
    double w;
    double tau2;
    S = atoi(argv[1]);
    w = atof(argv[2]);
    tau2 = atof(argv[3]);
    
    assert(S>1);
    assert(w>=0.0);
    assert(tau2>=50. && tau2<=1000.);
    
    int Runs = 5000; // when to stop simulations
        
    //setting parameters for the network
    Potts_params params; 
    set_params(&params); // default values, see above
    // now overriding default values for our needs
    params.S = S;
    params.w = w;
    params.T2 = tau2;
    
    // making a network
    int *xi;
    xi = new int[params.N*params.p];
    make_random_patterns(params, rand_seed, xi);
    //test_patterns(params, xi);
    PNet net('A', params); // the first argument should be char, any character
    net. construct_CJ(xi); // connecting Potts units
    net. SetupTables(rand_seed); // making a random sequence for updating order of units
    
    /*
    // making a folder for saving files
    char outdir[0x100];
    switch(version){
        case 0:
            snprintf(outdir, sizeof(outdir), "mkdir -p data_w%.1f_gA%.1f",w, gammaA);
            break;
        default:
            std::cout<< "Invalid value: 0<= version <=3." << std::endl;
            exit(8);
    }    
    system(outdir);
    */
    // Running network
    
    Network_runner runner; // a simple class object for running the network
    runner.nCue = 1;
    runner.dt = 5;
    runner.run_net(&net, 0, Runs, rand_seed, xi);
    
    return 0;
}
