/********************************************************
*   Member functions of PNet class   
* Turn on DEBUG to see some details
********************************************************/
//#define DEBUG

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
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

// ------------------------------------------------------------------------
PNet::PNet(Potts_params *params, int NumSet, int seed){        
        /*
        Constructor of a Potts network, fully connected
        */
        std::cout << "PNet: parameterized constructor." << std::endl;
        

        this->params = *params;
        this->NumSet = NumSet;
        assert( params->Cm < params->N ); // if Cm==N, then use fully-connected version
        // getting memory for variables
        int N = params->N;
        int S = params->S;
        s = new double[N*S];
        s0 = new double[N];
        theta = new double[N*S];
        theta0A = new double[N];
        theta0B = new double[N];
        h = new double[N*S];
        r = new double[N*S];
        J = new double[N*params->Cm*S*S];
        weight = new double[S+1];
        C = new int[N*params->Cm];
        
        mean_activity = params->a;
        
        Permut = new int[N*NumSet];
        this->SetupTables(seed); // seed srand48()
        this->make_J_null();
        this->make_dilution_matrix();
        
}
// ------------------------------------------------------------------------
PNet::~PNet(){
    /* desctructor */
    std::cout << "PNet destructor." << std::endl;
    delete[] s;
    delete[] s0;
    delete[] theta;
    delete[] theta0A;
    delete[] theta0B;
    delete[] h;
    delete[] J;
    delete[] r;
    delete[] Permut;
    delete[] weight;
    delete[] C;
}
// ------------------------------------------------------------------------
void PNet::initialise(){
    /*
    Initialise every variable with equilibrium values
    */
    
    // LATER REWRITE to compute sigma0 in general
    double sigma0 = 0.0265; // The default value only works for beta=11.0, S=7, U=0.1.
    
    int i, k;
    int N = params.N;
    int S = params.S;
    
    for(i=0;i<N;i++){
        for(k=0;k<S;k++)
            *(s+i*S+k) = sigma0;  
        *(s0+i) = 1.-S*sigma0;
        theta0A[i] = params.gammaA*(1.-s0[i]);
        theta0B[i] = (1.0- params.gammaA)*(1. - s0[i]);
    }

    //compute h, r, theta*/
    for(i=0;i<N;i++){
        this->compute_field(i);
        for(k=0;k<S;k++){
            *(r+i*S+k) = *(h+i*S+k); //values in the stationary state
            *(theta+i*S+k) = *(s+i*S+k); //values in the stationary state
        }
    }
    std::cout << "After initialising network into an equilibrium state." << std::endl;   
}
// ------------------------------------------------------------------------
void PNet::initialise(int cue, const int *xi){
    /*
    Initialise with one pattern, given by cue.
    */
    
    // LATER REWRITE to compute sigma0 in general
    double sigma0 = 0.999; // pattern-wanted state is given this value
    int state;
    int i, k;
    int N = params.N;
    int S = params.S;
    
    for(i=0;i<N;++i){
        for(k=0;k<S;++k)
            *(s+i*S+k) = 0.0;  
        state = *(xi+cue*N+i);
        if(state!=S){
            *(s+i*S+state) = sigma0;
            *(s0+i) = 1.0 - sigma0;
        }
        else{
            *(s0+i) = 1.0;
        }
        theta0A[i] = 0.0;//params.gammaA*(1.-s0[i]);
        theta0B[i] = 0.0;//(1.0- params.gammaA)*(1. - s0[i]);
    }
        
    //compute h, r, theta
    for(i=0;i<N;++i){
        this->compute_field(i);
        for(k=0;k<S;++k){
            r[i*S+k] = h[i*S+k]; //values in the stationary state
            theta[i*S+k] = 0.0;//s[i*S+k]; //values in the stationary state
        }
    }
    
    std::cout << "After initialising with a pattern: cue = " << cue <<std::endl;
    
}
// ------------------------------------------------------------------------
void PNet::make_dilution_matrix(){
    /*
    Make the dilution matrix, based on drand48()--->srand48()
    */
    int i, x;
    int array[params.Cm];
    std::cout << "Random dilution: " << params.Cm << std::endl;
    for(i=0;i<params.N;++i){
        make_Cij(i, 0, params.N, params.Cm, array);
        for(x=0;x<params.Cm;++x)
            C[i*params.Cm+x] = array[x];
    }
}
// ------------------------------------------------------------------------
void PNet::make_J_null(){
    /*
    Make all connections zero
    */
    int i,x,k,l, ind;
    int N = params.N;
    int Cm = params.Cm;
    int S = params.S;
    for(i=0;i<N;++i){
        for(k=0;k<S;++k){
            for(x=0;x<Cm;++x){
                for(l=0;l<S;++l){
                    ind = i*Cm*S*S + k*Cm*S + x*S + l;
                    *(J+ind) = 0.0;
                }
            }
        }
    }
}
// ------------------------------------------------------------------------
void PNet::make_J_assoc(int size, const int * eta, double b, double lambda) {
    /* 
    Connect units in a Hebbian way:
    * size = number of rows in eta, that is, number of memories
    * eta = flat array of size*N
    * b = sparsity of each memory representation
    * lambda = relative coefficient to multiply
    */
    
    int i,x,k,l,mu;
    int ind;
    int S = params.S;
    int Cm = params.Cm;
    double cof1 = b/(double)S;
    double cof2 = 2.0*(double)Cm*b*(1.-cof1);
    double factor1, factor2;
    double temp = 0.0;
    double factor = 1.0 + lambda;
    for(i=0;i<params.N;++i){
        for(k=0;k<S;++k){
            for(x=0;x<Cm;++x){
                for(l=0;l<S;++l){
                    ind = i*Cm*S*S + k*Cm*S + x*S + l;
                    temp = 0.0;
                    for(mu=0; mu<size;++mu){
                        factor1 = (double)(eta[mu*params.N+i]==k)-cof1;
                        factor2 = (double)(eta[mu*params.N+C[i*Cm+x]]==l)-cof1;
                        temp += factor1*factor2;
                    }
                    *(J+ind) = factor*temp/cof2;
                }
            }
        }
    }
    std::cout << "After making Hebbian connections..." << std::endl;
}

// ------------------------------------------------------------------------
void PNet::update_unit(int i, int flag, double *f1, double *f2){
    /*
    Update one unit i of the member variable s[i], given an extern field
    Parameters:
    * f1, f2 = pointers of external fields
    * flag = 0 ---> don't use any external field
    * flag = 1 ---> use only the first external field
    * flag = 2 ---> use both external fields
    */
    int S = params.S;
    int k;
    //self-excitation
    double self = params.w*(1.0-s0[i])/(double)S; 
    this->compute_field(i);
    for(k=0;k<S;k++){
        *(h+i*S+k) += params.w*s[i*S+k]-self; //self-term
        if(flag>0){ // add external field
            *(h+i*S+k) += *(f1+i*S+k);
            if(flag>1){
                *(h+i*S+k) += *(f2+i*S+k);
            } 
        }
        
        theta[i*S+k] += (s[i*S+k]-theta[i*S+k]) / params.T2; 
        r[i*S+k] += (h[i*S+k]-theta[i*S+k]-r[i*S+k]) / params.T1;
    }
    theta0A[i] += (params.gammaA*(1.-s0[i])-theta0A[i])/params.T3A;
    theta0B[i] += ((1.-params.gammaA)*(1.-s0[i])-theta0B[i])/params.T3B;
    
    //Computing s[i,k], avoiding overflow i
    double rmax = -1000.0;
    double r0 = theta0A[i] + theta0B[i] + params.U;
    for(k=0;k<S;k++){
        if(r[i*S+k]>rmax)
            rmax = r[i*S+k];
    }
    if(r0>rmax)
        rmax = r0;
    double Z_ = 0.;
    for(k=0;k<S;k++){
        this->weight[k] = exp(params.beta*(r[i*S+k] - rmax));
        Z_ = Z_ + this->weight[k];
    }
    this->weight[S] = exp(params.beta*(r0 - rmax));
    Z_ = Z_ + this->weight[S];
    double invZ = 1./Z_;
    for(k=0;k<S;k++){
        s[i*S+k] = invZ*this->weight[k];
    }
    s0[i] = invZ*this->weight[S];
}
//-----------------------------------------------------------------
void PNet::compute_field(int i){
    /*
    Compute the internal field
    */
    int x, k, l, index;
    int Cm = params.Cm;
    int cof1 = Cm*params.S*params.S;
    int cof2 = Cm*params.S;
    for(k=0;k<params.S;k++){
        h[i*params.S+k] = 0.;
        for(x=0;x<Cm;x++){
            for(l=0;l<params.S;l++){
                index = i*cof1 + k*cof2 + x*params.S + l;
                h[i*params.S+k] += (*(J+index)) * s[C[i*Cm+x]*params.S+l];
            }
        }
    }
}

// ------------------------------------------------------------------------
void PNet::update_network(int n, int flag, double *f1, double *f2){
    /*
    Update the member variable s[i], given two external fields
    * n = the index of sequence for updating order
    * flag: an integer of how many external field are used, 0, 1 and 2
    * f1 = a pointer of the first external field
    * f2 = a pointer of the second external field
    * DO NOT regulate the mean activity now
    */
    int i, index;
    assert(flag<3 && flag>=0);
    //this->compute_mean_activity();
    //double prev_s0 = 0.0;
    for(i=0;i<params.N;++i){
        index = *(Permut+n*params.N+i);
        this->update_unit(index, flag, f1, f2);
    }
}

//-----------------------------------------------------------------
void PNet::SetupTables(int seed){
    /*
    Make a random sequence of unit indices for updating order
    */
    int i, k;
    int N = params.N;
    int vec[N];
    srand48(seed); // "shuffle" uses drand48()
    
    for(i=0;i<N;++i)
        vec[i] = i;
    for(k=0; k<NumSet;++k){
        shuffle(vec, N);
        for(i=0;i<N;++i)
            *(Permut+k*N+i) = *(vec+i);
    }
}

// ------------------------------------------------------------------------
void PNet::compute_mean_activity(){
    /*
    * Compute the  mean activity of the net.
    * Set the member variable, mean_activity
    */
    double sum = 0.0;
    for(int i=0;i<params.N;++i){
        sum += (1.0-s0[i]);
    }
    mean_activity = sum/(double)params.N;
}
// ------------------------------------------------------------------------
void PNet::compute_m(const int *xi, double *m){
    int  i, k, mu;
    double ma, maa;
    int xi_i; // one state wanted by the current pattern.
    double cof1 = params.a/(double)params.S;
    double cof2 = (double)params.N*(params.a)*(1.-cof1);

    assert(params.p>0);
#pragma omp parallel for private(i, k, mu, ma, maa, xi_i)
    for(mu=0; mu<params.p; mu++){
        maa=0.;
        for(i=0; i<params.N; i++){
            ma=0.;
            xi_i = *(xi+mu*params.N+i);
            for(k=0;k<params.S;k++){
                ma = ma + ((double)(xi_i==k)-cof1)*s[i*params.S+k];  
            }
            maa=maa+ma;
        }
        *(m+mu) = maa/cof2;
    }
}
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// * A CLASS THAT RUNS THE NETWORK AS NEEDED
// ------------------------------------------------------------------------
Network_runner::Network_runner(){
    
    this->nCue = 1;
    this->Runs = 1000;
}
Network_runner::~Network_runner(){
    std::cout << "Simulation is finished." << std::endl;
}
// ------------------------------------------------------------------------
void Network_runner::make_null_connection(PNet *post_net, PNet *pre_net, 
    double * Jmat){
    /*
    Fill the matrix's elements by zeros
    ----------
    Paramters: None
    */
    int N = post_net->params.N;
    int S1 = post_net->params.S;
    int S2 = pre_net->params.S;
    assert(N==pre_net->params.N);
    int i, j, k, l, index;
    
    for(i=0; i<N; ++i){
        for(k=0; k<S1; ++k){
            for(j=0; j<N; ++j){
                for(l=0; l<S2; ++l){
                    index = i*N*S1*S2 + k*N*S2 + j*S2 + l;
                    *(Jmat+index) = 0.0;
                }
            }
        }
    }
}
// ------------------------------------------------------------------------
void Network_runner::make_Hebb_connection(PNet *post_net, PNet *pre_net, double *Jmat,
    const int *xi_post, const int *xi_pre, double lambda, int *match_table){
    /*
    Fill the matrix's elements of Jmat by Hebbian learning rule
    ----------
    Paramters:
        * lambda -> coefficient for intra-connections
        * match_table -> p_post by p_pre matrix of matching patterns between two sub-networks
    */
    
    std::cout << "Connecting two networks..." << std::endl;
    int N1 = post_net->params.N;
    int N2 = pre_net->params.N;
    int S1 = post_net->params.S;
    int S2 = pre_net->params.S;
    int Cm = post_net->params.Cm;
    int p = post_net->params.p;
    assert(p==pre_net->params.p);
    assert(Cm<=N2);
    double factor = 1.0 - lambda; // overall factor for all connections
    double denom = 0.5/(double)Cm;

    double a1 = post_net->params.a;
    double a2 = pre_net->params.a;
    double a = sqrt(a1*a2);

#pragma omp parallel for
    for(int i=0; i<N1; ++i){
        int j, k, l, index;
        int mu, nu;

        double temp = 0.0;
        double den1 = sqrt(1.0-a1/S1);
        double den2 = sqrt(1.0-a2/S2);

        for(k=0; k<S1; ++k){
            for(j=0; j<N2; ++j){
                for(l=0; l<S2; ++l){
                    index = i*N2*S1*S2 + k*N2*S2 + j*S2 + l;
                    temp = 0.0;
                    for(mu=0;mu<p;++mu){
                        for(nu=0;nu<p;++nu){
                            if(match_table[mu*p+nu]!=0){
                                double cof1 = (double)(xi_post[mu*N1+i]==k)-a1/S1;
                                double cof2 = (double)(xi_pre[nu*N2+j]==l)-a2/S2;
                                temp += cof1*cof2; 
                            }
                        }  
                    }
                    temp = temp/a/den1/den2;
                    *(Jmat+index) = factor*temp*denom;
                }
            }
        }
    }
    // dilute the connections
    if(Cm<N2){
        
        int size = N2 - Cm; // number of elements to make zero
                            // FIXME but why to dilute it like this? will not it be biased?
        std::cout << "Diluting inter-net connections: " << size << std::endl;
#pragma omp parallel for
        for(int i=0; i<N1; ++i){
            int array[size];
            make_Cij(-1, 0, N2, size, array);
            for(int k=0; k<S1; ++k){
                for(int j=0; j<size; ++j){
                    for(int l=0; l<S2; ++l){
                        int index = i*N2*S1*S2 + k*N2*S2 + array[j]*S2 + l; // FIXME why array[j]?
                        *(Jmat+index) = 0.0;
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------------
void Network_runner::run_two_nets(PNet *p_net, PNet *f_net, double *Jf2p, double *Jp2f,
        const int *xi_p, const int *xi_f, 
        std::ofstream & buf1, std::ofstream &buf2, int cue, int save_all){
    /*
    Run two network with both f->p  and p->f connections
    Parameters:
    * p_net, f_net: pointers of network objects
    * Jf2p, Jp2f: pointers of connections between two networks
    * xi_p, xi_f: pointers of memory patterns of two networks
    * buf1, buf2: buffers of saving data into the file
    * cue: pattern index for running one simulation
    * save_all: if save_all==1, then save all overlaps to buf2.
    */
    int N1 = p_net->params.N;
    int N2 = f_net->params.N;
    int S1 = p_net->params.S;
    int S2 = f_net->params.S;
    
    int p = p_net->params.p;
    assert(p==f_net->params.p);

    double field_f2p[N1*S1];
    double field_cue2p[N1*S1];
    double field_p2f[N2*S2];
    double field_cue2f[N2*S2];
    int array4updating[N1+N2]; // to save which sub-net to update first
    //update both networks
    double ovlp[p];
    int max_index, rand_num_p, rand_num_f;
    int i, t, index, counter1, counter2;
    for(i=0;i<N1+N2;++i){
        if(i<N1){
            *(array4updating+i) = 1;
        }
        else{
            *(array4updating+i) = 2;
        }
    }

    p_net->initialise();
    f_net->initialise();
    
    std::cout << "Running networks..." << std::endl;
    for(t=0;t<this->Runs+1;++t){
        

        shuffle(array4updating, N1+N2);

        this->make_field_decay(f_net, field_cue2f, cue, xi_f, t);
        this->make_field_decay(p_net, field_cue2p, cue, xi_p, t);

        this->compute_field(p_net, f_net, Jf2p, field_f2p);
        this->compute_field(f_net, p_net, Jp2f, field_p2f);
        rand_num_p = (int)(p_net->NumSet*drand48());
        rand_num_f = (int)(f_net->NumSet*drand48());
        
        counter1 = 0;
        counter2 = 0;
        for(i=0;i<N1+N2;++i){
            if( *(array4updating+i)==1 ){
                index = *(p_net->Permut+rand_num_p*N1+counter1);
                p_net->update_unit(index, 2, field_f2p, field_cue2p);
                ++counter1;
            }
            else{
                index = *(f_net->Permut+rand_num_f*N2+counter2);
                f_net->update_unit(index, 2, field_p2f, field_cue2f);
                ++counter2;
            }
        }
        assert(counter1==N1);
        assert(counter2==N2);
        
        // test the overlaps
        if(t%dt==0){

            p_net->compute_m(xi_p, ovlp);
            buf1 << t <<  "\t";
            if(save_all==1){
                buf2 << t << "\t";
                save_doubles(p, ovlp, buf2);
            }
            save_overlaps(2, p_net->params.p, ovlp, buf1);
            f_net->compute_m(xi_f, ovlp);
            if(save_all==1){
                buf2 << t << "\t";
                save_doubles(p, ovlp, buf2);
            }
            max_index = save_overlaps2(2, f_net->params.p, ovlp, buf1);
            
            buf1 << std::endl;
            if(ovlp[max_index]<0.01 && t>200){
                t = Runs;
            }
        }

        if(t%100==0){
            std::cout << "step = " << t << std::endl;
        }
    }
}
// ------------------------------------------------------------------------
void Network_runner::make_field_decay(PNet *net_ptr, double *field, int cue, 
    const int *xi, int t){
    /*
    Make a vector of field that decays exponentially over time.
    ------- Parameters
    * cue = a pattern index to cue
    * t = effective time, field is zero if t==0 or t>9*tau
    * xi = pointer for the patterns    
    */

    int N = net_ptr->params.N;
    int S = net_ptr->params.S;
    
    int i, k;
    if(t==0 || t>9*tau){
        for(i=0;i<N;++i){
            for(k=0;k<S;++k){
                *(field+i*S+k) = 0.0;
            }
        }
    }
    else{
        double INconst = g*exp(-(double)t/tau); 
        for(i=0;i<N;++i){
            for(k=0;k<S;++k){
                *(field+i*S+k) = 0.0;
            }
            if( *(xi+cue*N+i) != S ){
                field[i*S+*(xi+cue*N+i)] = INconst;
            }
        }
    }
}
// ------------------------------------------------------------------------
void Network_runner::compute_field(PNet *post_net, PNet *pre_net, const double *J
    ,double *field){
    /*
    compute field
    */
    int N = pre_net->params.N;
    assert(N==post_net->params.N);
    int S1 = post_net->params.S;
    int S2 = pre_net->params.S;

#pragma omp parallel for
    for(int i=0;i<N;++i){
        int j, k, l, index;
        for(k=0;k<S1;k++){
            field[i*S1+k] = 0.;
            for(j=0;j<N;j++){
                for(l=0;l<S2;l++){
                    index = i*N*S1*S2 + k*N*S2 + j*S2 + l;
                    field[i*S1+k] += (*(J+index)) * pre_net->s[j*S2+l];
                }
            }
        }
    }
}
// ------------------------------------------------------------------------
/*

*/
