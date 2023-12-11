
/********************************************************************************
* Study latching of f-net and p-net when the two networks interact. 
* One frontal pattern is connected with Z posterior patterns.
*   * PILOT VERSION--> only 10 of f-patterns are connected with posterior patterns
* Updating order is inter-mingled between two subnets
* J = (1+lambda)*J_{intra} + (1-lambda)*J_{inter}
*   * lambda = 1.0 ---> isolated networks
*   * betwee: (0, 1) ---> weaker strength for inter-connections
*   * lambda = 0.0 ---> fully interacting.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <cassert>
#include <math.h>
#include <omp.h>

#ifndef FUNCTION_H
#include "functions.h"
#endif /* FUNCTION_H */

#ifndef RANDGEN_H
#include "rand_gen.h"
#endif /* RANDGEN_H */

#ifndef PNET_H
#include "pnet.h"
#endif /* PNET_H */ 
/****************** GLOBAL VARIABLES ********************************
* I NEVER use global variables in my code!
* If you need to use one, include "extern double x" in another file.
********************************************************************/
void set_params(Potts_params *params){
    /*
    Set default values for paramters
    */
    // (*params).N = 50;     // number of units in a net
    (*params).N = 500;     // number of units in a net
    (*params).Cm = (int)(*params).N * 0.15;     // mean connections per neuron
    (*params).S = 7;        // number of Potts states
    params->a = 0.25;
    params->p = (int)(*params).N * 0.20; // for non-associative memory connections
    (*params).U = 0.1;      // global threshold
    (*params).T1 = 20.;     // integration of incoming fields
    (*params).T2 = 200.;  // adaptation time scale
    (*params).T3A = 10.;    //fast inhibition
    (*params).T3B = 100000.; //slow inhibition
    (*params).beta = 11.;   // inverse temperature
    (*params).w = 1.1;      // self-reinforcement term
    (*params).gammaA = 0.5; // fast inhibition
}

//------------------------------------------------------------------------
int main(int argc, char *argv[]){

    omp_set_num_threads(40);

    double w_p = atof(argv[1]);
    double w_f = atof(argv[2]);
    int S_p = atoi(argv[3]);
    int S_f = atoi(argv[4]);
    double tau_p = atof(argv[5]);
    double tau_f = atof(argv[6]);
    double lambda_p = atof(argv[7]); // coeff. for intra-connections
    double lambda_f = atof(argv[8]);
    int rand_seed = atoi(argv[9]); // dynamic seed

    int Z_flag = atoi(argv[10]); // methods of connecting the above
    int Z = atoi(argv[11]); // memories to connect between pre and post net

    int N = 500;
    int p = N * 0.20;
    int pat_seed = 1;
    assert(w_p>=0.0);
    assert(w_f>=0.0);
    assert(S_p>1 && S_f>1);
    assert(tau_p>0 && tau_f>0);
    assert(lambda_p>=-1. && lambda_p<=1.);
    assert(lambda_f>=-1. && lambda_f<=1.);
    assert(rand_seed>0);
    
    
    Potts_params p_params, f_params;
    set_params(&p_params);
    set_params(&f_params);
    p_params.w = w_p;
    p_params.S = S_p;
    p_params.p = p;
    p_params.T2 = tau_p;

    f_params.w = w_f;
    f_params.S = S_f;
    f_params.p = p;
    f_params.T2 = tau_f;
    
    int N_p = p_params.N;
    int N_f = f_params.N;
    
    // reading patterns
    int *xi_p = new int[N_p*p];
    if(1){
        std::ifstream file;
        std::ifstream & file_ref = file;
        char fname[0x100];
        snprintf(fname, sizeof(fname), "backend/model_potts_hybrid/patterns/pat_N%d_p%d_S%d_a%.2f_seed%d",
            N_p, p, S_p, p_params.a, pat_seed+0);
        file.open(fname);
        read_patterns(N_p, p, xi_p, file_ref);
        file.close();
    }
    int *xi_f = new int[N_f*p];
    if(1){
        std::ifstream file;
        std::ifstream & file_ref = file;
        char fname[0x100];
        snprintf(fname, sizeof(fname), "backend/model_potts_hybrid/patterns/pat_N%d_p%d_S%d_a%.2f_seed%d",
            N_f, p, S_f, f_params.a, pat_seed+1);
        file.open(fname);
        read_patterns(N_f, p, xi_f, file_ref);
        file.close();
    }

    // make networks
    PNet p_net(&p_params, 1000, pat_seed); 
    p_net.make_J_assoc(p, xi_p, p_params.a, lambda_p);

    PNet f_net(&f_params, 1000, pat_seed+1000);
    f_net.make_J_assoc(p, xi_f, f_params.a, lambda_f); 

    // #################################################
    Network_runner manager;
    double *Jf2p = new double[N_p*N_f*S_p*S_f];
    double *Jp2f = new double[N_f*N_p*S_f*S_p];

    int table_f2p[p*p];
    int table_p2f[p*p];

    // making table of connections from posterior to frontal - i_frontal * j_posterior
    make_match_internet(table_p2f, Z_flag, Z, p);

    // equivalent to transposing the unflattened 2D matrix - creating the symmetric connection
    for(int i=0;i<p;++i){
        for(int j=0;j<p;++j){
            table_f2p[i*p+j] = table_p2f[j*p+i];
        }
    }

    if(lambda_p==1.0){ // NOTE lambd = 1.0 means no connection. lambda = 0.0 means same strength.
        manager.make_null_connection(&p_net, &f_net, Jf2p);
    }
    else{
        manager.make_Hebb_connection(&p_net, &f_net, Jf2p, xi_p, xi_f, lambda_p, table_f2p);
    }
    if(lambda_f==1.0){
        manager.make_null_connection(&f_net, &p_net, Jp2f);
    }
    else{
        manager.make_Hebb_connection(&f_net, &p_net, Jp2f, xi_f, xi_p, lambda_f, table_p2f);
    }
    manager.Runs = 5000;
    // #################################################
    if(1){
        char make_outdir[0x100];
        snprintf(make_outdir, sizeof(make_outdir), 
        "mkdir -p backend/data/Zf%d_Z%d/lambP%.1f_lambF%.1f_%dS%d_%.2fw%.2f_%.1fT%.1f_p%d_seed%d",
        Z_flag, Z, lambda_p, lambda_f,S_p,S_f,w_p,w_f,tau_p,tau_f,p, pat_seed);
        system(make_outdir);
    }

    srand48(rand_seed);
    rlxd_init(2, rand_seed);

    // #pragma omp parallel for
    for(int cue=0;cue<1;++cue){

        // std::cout << "thread id = " << omp_get_thread_num() << "\t" << "cue = " << cue << "\n" << std::endl;

        char fname1[0x100];
        snprintf(fname1, sizeof(fname1),
            "backend/data/Zf%d_Z%d/lambP%.1f_lambF%.1f_%dS%d_%.2fw%.2f_%.1fT%.1f_p%d_seed%d/seq_cue%d_seed%d",
            Z_flag, Z, lambda_p,lambda_f,S_p,S_f,w_p,w_f,tau_p,tau_f,p,pat_seed,cue,rand_seed);
        std::ofstream file1(fname1);
        file1 << std::fixed << std::setprecision(3);

        // if(cue<2 && rand_seed<3){
        if(cue<p){
            char fname2[0x100];
            snprintf(fname2, sizeof(fname2), 
                "backend/data/Zf%d_Z%d/lambP%.1f_lambF%.1f_%dS%d_%.2fw%.2f_%.1fT%.1f_p%d_seed%d/mall_cue%d_seed%d",
                Z_flag, Z, lambda_p,lambda_f,S_p,S_f,w_p,w_f,tau_p,tau_f,p,pat_seed,cue,rand_seed);
            std::ofstream file2(fname2);
            file2 << std::fixed << std::setprecision(4);

            manager.run_two_nets(&p_net, &f_net, Jf2p, Jp2f, xi_p, xi_f, file1, file2, cue, 1);
            file2.close();
        }
        else{
            manager.run_two_nets(&p_net, &f_net, Jf2p, Jp2f, xi_p, xi_f, file1, file1, cue, 0);
        }
        file1.close();
    }
    
    delete[] xi_p;
    delete[] xi_f;
    delete[] Jf2p;
    delete[] Jp2f;
    return 0;

}
