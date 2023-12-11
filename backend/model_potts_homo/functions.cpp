
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "global.h"
#include "rand_gen.h"
#include <cassert>
// global variable defined in "main.cpp"
extern bool show_sim;
//**************************************************************************************************************************
//**** Functions that are used only in this file, and
// thus without any prototypes in "global.h"
//**************************************************************************************************************************

// ------------------------------------------------------------------------
// define a structure
typedef struct {
    int index;
    double field;
    int state;
} mystruct;

// ------------------------------------------------------------------------
// comparing function for qsort()
int cmpf (const void *x, const void *y)
{
  double xx = (*(mystruct*)x).field;
  double yy = (*(mystruct*)y).field;
  /*sort by decreasing order*/
  if (xx > yy) return -1;
  if (xx < yy) return  1;
  return 0;
}

// ------------------------------------------------------------------------
// shuffling function
void shuffle(int *array, int size)
{
    double rand_num;
    if (size > 1)
    {
        int i, j, t;
        double rand_val;
        for (i = 0; i < size - 1; i++)
        {
          rand_num = drand48();
          rand_val = rand_num*(double)(size-i-2) + 1.;
          j = i + (int)rand_val;
          t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

// ------------------------------------------------------------------------
// making a random vector
void make_rand_vector(double *T, int n, double param, int seed)
{
    double vector[n];
    int i;
    
    rlxd_init(1, seed);
    int version = 1;
    switch(version){
        case 1:
            wexpo_dble(vector, n, param);
            for(i=0;i<n;i++) *(T+i) = vector[i];
            break;
        case 2:
	        gauss_dble(vector, n);
	        for(i=0;i<n;i++) *(T+i) = param*vector[i];
	        break;
	    case 3:
	        ranlxd(vector,n);
	        for(i=0;i<n;i++) *(T+i) = param*vector[i];
	        break;
	    default:{
	        std::cout<<"Invalid value: version should be 1, 2 or 3.\n";
	        exit(8);
	    }
    }
}
// ------------------------------------------------------------------------
// Pick $n$ random indices between [nMin, nMax), avoiding $i$
void make_Cij(int i, int nMin, int nMax, int n, int *result)
{
    if(n>nMax-nMin) {
        printf("Error in make_cij\n");
        
        exit(0);
    }
    
    int j;
    int array[nMax-nMin];
    int counter = 0;
    for(j=nMin;j<nMax;j++){
        if (j!=i){
            array[counter]=j;
            counter++;
        }
    }
    shuffle(array, counter);
    for(j=0;j<n;j++) result[j]=array[j];
}

// ------------------------------------------------------------------------
// making random patterns by using qsort()
void make_random_patterns(Potts_params params, int seed, int *xi){
    int i, mu;
    srand48(seed);
    int N = params.N;
    int S = params.S;
    int p = params.p;
    int Na = (int) ((double)N*params.a);
    int rand_state;
    int counter = 0;
    for(mu=0;mu<p;mu++){   
        mystruct X[N];
        for(i=0;i<N;i++){
            //first initialize each pattern to a null (quiet) state
            *(xi+mu*N+i) = S;
            //xi[mu*N+i] = S;
            	
            X[i].index = i;
            X[i].field = drand48();
            X[i].state = (int)((double)S*drand48());
        }
        
        //sort the patterns in order of decreasing field keeping the indices
        qsort(X, N, sizeof(mystruct), cmpf);  

        // set N*a largest of them to become active
        counter = 0;
        for(i=0;i<Na;i++){
            rand_state = X[i].state; 
            //while (rand_state==S) rand_state = (int)((double)S*drand48());
            *(xi+mu*N+X[i].index) = rand_state; 
            //xi[mu*N+X[i].index] = rand_state;
            if (rand_state!=S) counter = counter + 1;
        }
          
    }
    std::cout << "After making random patterns. Na = " << counter << std::endl;
}
//Another version, by using shuffle
void make_random_patterns_V2(Potts_params params, int seed, int *xi){
    int i, mu;
    srand48(seed);
    int N = params.N;
    int S = params.S;
    int p = params.p;
    int Na = (int) ((double)N*params.a);
    int X[N];
    for(i=0;i<N;i++){
        X[i] = i;
    }
    int counter;
    int rand_state;
    
    
    for(mu=0;mu<p;mu++){   
        shuffle(X, N);  
        // set N*a largest of them to become active
        counter = 0;
        for(i=0;i<N;i++){
            rand_state = S;
            if(i<Na){ 
                while (rand_state==S){
                    rand_state = (int)((double)S*drand48());
                }
            }
            
            *(xi+mu*N + X[i]) = rand_state;
            if (rand_state!=S) counter = counter + 1;
        }
        //std::cout << counter << "\t";   
    }
    
    
    std::cout << std::endl << "After making random patterns. Na = " << counter << std::endl;
    
}
// ------------------------------------------------------------------------
void test_patterns(Potts_params params, const int *xi){
    int mu,i;
    
    int N = params.N;
    int p = params.p;
    int S = params.S;
    
    mu = 0;
    int Na;
    for(mu=0;mu<p;mu++){  
        Na = 0;
        for(i=0;i<N;i++)
            if(xi[mu*N+i]<S) Na = Na + 1;
        std::cout << "mu = " << mu << ", Na= " << Na << std::endl;
    }
    
}

//**************************************************************************************************************************
//   Member functions of PNet class   
// their prototypes are in global.h
// *************************************************************************************************************************
// ------------------------------------------------------------------------
Network_runner::Network_runner(){
    
    this->dt = 10; // default
    this->nCue = 1; // default
}
Network_runner::~Network_runner(){
    std::cout << "Simulation is finished." << std::endl;
}
// Run networks
void Network_runner::run_net(PNet *net, int start_cue, int Runs, int seed, const int *xi)
{
    
    srand48(seed);
    Potts_params params = net->params;
    int cue = 0;
    int t; // integer time index = number of updates of the whole net.
    int i; // dummy variables
    int rand_num; // which one of random sequences to use
    int N = params.N;
    
    
    std::ofstream mall;
    std::ofstream & mall_ref = mall; // reference
    
    char buffer[0x100];
    
    //std::cout<< "w = " << params.w << "\t" << "gammaA = " << params.gammaA << std::endl;
    
    for(cue=start_cue; cue<this->nCue+start_cue;cue++){
        
        net->initialise();
        
        
        snprintf(buffer, sizeof(buffer), "./backend/data/mall_S%d_w%.2f_gA%.1f_T%.1f_cue%d", params.S, params.w, params.gammaA, params.T2, cue);
        mall.open(buffer, std::ios::out);
        mall << std::fixed << std::setprecision(5) << std::endl;

        // START DYNAMICS 
        std::cout << "----- dynamics started with cue = "<< cue << " -----" << std::endl;
        
        
        int each_unit; // pick one unit randomly
        for(t=0;t<Runs+1;t++){
            rand_num = (int)(net->NumSet*drand48());
            for(i=0;i<N;i++){
                each_unit = net->Permut[rand_num][i];
                net->update_state(each_unit, t, *(xi+cue*N+each_unit));
            }
            // snapshot the network
            if (t%dt==0) {
                net->snapshot(t, xi, mall_ref);
                if (net->max_ovlp < 0.001 && t>20){
                    t = Runs;
                    //std::cout << t << "\t" << net->max_ovlp << std::endl;
                }
            }
            
        }
        mall.close();
        //std::cout << "---- dynamics finished with cue = "<< cue << " ----" << std::endl;
    }
}


//**************************************************************************************************************************
//   Member functions of PNet class   
// their prototypes are in global.h
// *************************************************************************************************************************

// constructor
PNet::PNet(char name, Potts_params params){
        int i;
        // setting default parameter values
        NumSet = 50; // usually fine with this.
        
        this->name = name;
        // getting memory for variables
        this->params = params;
        s = new double*[params.N];
        for(i=0; i<params.N; i++) s[i] = new double[params.S+1];
        h = new double[params.S];
        
        r = new double*[params.N];
        for(i=0; i<params.N; i++) r[i] = new double[params.S];
        theta = new double*[params.N];
        for(i=0; i<params.N; i++) theta[i] = new double[params.S];
        C = new int*[params.N];
        for(i=0; i<params.N; i++) C[i] = new int[params.Cm];
        theta0A = new double[params.N];
        theta0B = new double[params.N];
        J = new double[params.N*params.Cm*params.S*params.S];
        Permut = new int*[NumSet];
        for(i=0; i<NumSet; i++) Permut[i] = new int[params.N];

        
}
// ------------------------------------------------------------------------
// desctructor
PNet::~PNet(){
        int i;
        for(i=0; i<this->params.N; i++) delete[] s[i];
        delete[] s;
        for(i=0; i<this->params.N; i++) delete[] r[i];
        delete[] r;
        for(i=0; i<this->params.N; i++) delete[] theta[i];
        delete[] theta;
        for(i=0; i<this->params.N; i++) delete[] C[i];
        delete[] C;
        for(i=0; i<NumSet; i++) delete[] Permut[i];
        delete[] Permut;
        
        delete[] theta0A;
        delete[] theta0B;
        delete[] J;
        delete[] h;
}


// ------------------------------------------------------------------------
// initialise with one pattern
void PNet::initialise(int mu, const int *xi){
     int i, l, k, x;
    int N = this->params.N;
    int S = params.S;
    int index;
    for(i=0;i<params.N;i++) { 
        for(k=0;k<S+1;k++) s[i][k]=0.0;
        s[i][xi[mu*params.N+i]] = 1.0;
        theta0A[i] = params.gammaA*(1.-s[i][S]);
        theta0B[i] = (1.0- params.gammaA)*(1.-s[i][S]);  
    }
    
    for(i=0;i<N;i++){
        for(k=0;k<S;k++){
            h[k] = 0.;
            for(x=0;x<params.Cm;x++){
                for(l=0;l<S;l++){
                    index = i*params.Cm*S*S + k*params.Cm*S + x*S + l;
                    h[k] +=  (*(J+index)) * s[C[i][x]][l];
                }
            }
            r[i][k] = h[k]; //values in the stationary state
            theta[i][k] = s[i][k]; //values in the stationary state
        }
    }
    std::cout << "After initialising network into mu = " << mu << std::endl;
    max_mu = params.p + 1;
    max2_mu = params.p + 2;
    max_ovlp = -10.0;
    max2_ovlp = -20.0;
    R = 0.0;
}
// ------------------------------------------------------------------------
// initialise with equilibrium states
void PNet::initialise(){
    // later write to compute sigma0 in general
    double sigma0 = 0.0265; // The default value only works for beta=11.0, S=7, U=0.1.
    
    int index;
    int i, l, k, x;
    int N = this->params.N;
    int S = params.S;
    
    for(i=0;i<N;i++){
        for(k=0;k<S;k++) s[i][k] = sigma0;  
        s[i][S] = 1.-S*s[i][0];
        theta0A[i] = params.gammaA*(1.-s[i][S]);
        theta0B[i] = (1.0- params.gammaA)*(1.-s[i][S]);
    }

    //compute h, r, theta*/
    for(i=0;i<N;i++){
        for(k=0;k<S;k++){
            h[k] = 0.;
            for(x=0;x<params.Cm;x++){
                for(l=0;l<S;l++){  
                    index = i*params.Cm*S*S + k*params.Cm*S + x*S + l;
                    h[k] += (*(J+index)) * s[C[i][x]][l];
                }
            }
            r[i][k] = h[k]; //values in the stationary state
            theta[i][k] = s[i][k]; //values in the stationary state
        }
    }
    std::cout << "After initialising network into an equilibrium state." << std::endl;
    max_mu = params.p + 1;
    max2_mu = params.p + 2;
    max_ovlp = -10.0;
    max2_ovlp = -20.0;
    R = 0.0;
}

// ------------------------------------------------------------------------
void PNet::construct_CJ(const int *xi) 
{
    srand48(1990);
    int temp1[this->params.Cm];
    int i,j, k,l,mu;
    int ind;
    double cof1 = params.a/(double)params.S;
    double cof2 = (double)params.Cm*params.a*(1.-cof1);
    int S = params.S;
    for(i=0; i< params.N; i++) 
    {
        make_Cij(i, 0, params.N, params.Cm, temp1);
        for(j=0;j<params.Cm;j++) C[i][j] = temp1[j];
    }
    std::cout << "Random dilution, after making C[i][j].\n";
    
    for(i=0; i<params.N; i++){
        for(k=0; k<S; k++){
            for(j=0; j<params.Cm; j++){
                for(l=0; l<S; l++){
                    ind = i*params.Cm*S*S + k*params.Cm*S + j*S + l;
                    *(J+ind) = 0.0;
                    for(mu=0; mu<params.p; mu++){
                        *(J+ind) = *(J+ind) + ((double)(*(xi+mu*params.N+i)==k)-cof1)*((double)(*(xi+mu*params.N+C[i][j])==l)-cof1);
                    }
                    *(J+ind) = (*(J+ind)) / cof2;
                }
            }
        }
    }
    std::cout << "After computing J[i][j][k][l]." << std::endl;
}



// ------------------------------------------------------------------------
// i is unit number, n is the time (#.updates of net), 
// xi_i is state of the unit i wanted by cued pattern.
void PNet::update_state(int i, int n, int xi_i)
{
    int x, k, l, kmax;
    double Z, rmax, self, INconst;
    
    int S = params.S;
    int Cm = params.Cm;
    
    //double field[S];
    
    //self-excitation
    self = 0.; 
    for(l=0;l<S;l++) self = self + s[i][l]; 
    self=(params.w/(double)S)*self; 

    /*initial field*/
    double tau = 10.;
    double g = 5.;
    INconst = (double)(n>1)*g*exp(-((n-1)/tau)); 
    int index;
    for(k=0;k<S;k++){
        h[k] = 0.;
        for(x=0;x<Cm;x++){
            for(l=0;l<S;l++){
                index = i*Cm*S*S + k*Cm*S + x*S + l;
                h[k] += (*(J+index)) * s[C[i][x]][l];
            }
        }
        h[k] += params.w*s[i][k]-self;
        if(xi_i==k) h[k] += INconst;
        theta[i][k] += (s[i][k]-theta[i][k]) / params.T2; 
        r[i][k] += (h[k]-theta[i][k]-r[i][k]) / params.T1;
    }

    theta0A[i] += (params.gammaA*(1.-s[i][S])-theta0A[i])/params.T3A;
    theta0B[i] += ((1.-params.gammaA)*(1.-s[i][S])-theta0B[i])/params.T3B;
    
    //Computing s[i][k], avoiding overflow in the exponential
    rmax = -1000.0;
    kmax = -100;
    for(k=0;k<S;k++){
        if(r[i][k]>rmax){
            rmax = r[i][k];
            kmax = k;
        }
    }

    Z=0.;
    for(k=0;k<S;k++){
        Z = Z + exp(params.beta*(r[i][k] - rmax));
    }
    Z = Z + exp(params.beta*(theta0A[i] + theta0B[i] + params.U - rmax));

    double invZ;
    invZ = 1./Z;
    for(k=0;k<S;k++){
        s[i][k] = invZ*exp(params.beta*(r[i][k]-rmax));
    }
    s[i][S] = invZ*exp(params.beta*(theta0A[i] + theta0B[i] + params.U - rmax));
}


// ------------------------------------------------------------------------
// snapshot the network state
void PNet::snapshot(const int t, const int *xi){

    int p = params.p;
    int mu;
    int  Mumax, Mumax2;
	double Mmax, Mmax2;
    double m[p];
    
    this->compute_m(xi, m);
    
    Mmax = -1.;
    Mmax2 = -2.;
    Mumax = p+1;
    Mumax2 = p+2;
    
    //find pattern that has maximal overlap with net
    this -> R = 0.; // order paramter in Amit's paper.
    for(mu=0;mu<p;mu++){
        R = R + m[mu]*m[mu]; 
        if(m[mu]>Mmax){
            Mmax2 = Mmax;
            Mmax = m[mu];
            Mumax2 = Mumax;
            Mumax = mu;
        } 
        else if(m[mu]>Mmax2){
            Mmax2 = m[mu];
            Mumax2 = mu;
        }
    }
    R = R - Mmax*Mmax - Mmax2*Mmax2;
    R = R/(double)p;
    if(show_sim!=0) std::cout << "t = " << t << ", Mumax=" << Mumax << ", Mmax=" << Mmax << "\t" << ", Mmax2=" << Mmax2 <<std::endl;	
    max_mu = Mumax;
    max2_mu = Mumax2;
    max_ovlp = Mmax;
    max2_ovlp = Mmax2;
    	
}

void PNet::snapshot(const int t, const int *xi, std::ofstream &buffer){

    int p = params.p;
    int mu;
    int  Mumax, Mumax2;
	double Mmax, Mmax2;
    double m[p];
    
    this->compute_m(xi, m);
    
    Mmax = -1.;
    Mmax2 = -2.;
    Mumax = p+1;
    Mumax2 = p+2;
    
    //find pattern that has maximal overlap with net
    this -> R = 0.; // order paramter in Amit's paper.
    
    buffer<<t<<"\t";
    for(mu=0;mu<p;mu++){
        buffer<<m[mu]<<"\t";
        R = R + m[mu]*m[mu]; 
        if(m[mu]>Mmax){
            Mmax2 = Mmax;
            Mmax = m[mu];
            Mumax2 = Mumax;
            Mumax = mu;
        } 
        else if(m[mu]>Mmax2){
            Mmax2 = m[mu];
            Mumax2 = mu;
        }
    }
    buffer<<"\n";
    
    R = R - max_ovlp*max_ovlp - max2_ovlp*max2_ovlp;
    R = R/(double)p;
    if (show_sim!=0) std::cout << "t = " << t << ", Mumax=" << Mumax << ", Mmax=" << Mmax << "\t" << ", Mmax2=" << Mmax2 <<std::endl;	
    max_mu = Mumax;
    max2_mu = Mumax2;
    max_ovlp = Mmax;
    max2_ovlp = Mmax2;
    	
}
// ------------------------------------------------------------------------
void PNet::compute_m(const int *xi, double *m){
    int  i, k, mu;
    double ma, maa;
    int xi_i; // one state wanted by the current pattern.
    double cof1 = params.a/(double)params.S;
    double cof2 = (double)params.N*(params.a)*(1.-cof1);
    for(mu=0; mu<params.p; mu++){
        maa=0.;
        for(i=0; i<params.N; i++){
            ma=0.;
            xi_i = *(xi+mu*params.N+i);
            for(k=0;k<params.S;k++){
                ma = ma + ((double)(xi_i==k)-cof1)*s[i][k];  
            }
            maa=maa+ma;
        }
        m[mu] = maa/cof2;
    }
}

// ------------------------------------------------------------------------
// Produce random sequences for updating order of units. Without checking if two sequences are identical.
void PNet::SetupTables(int seed){
    int item, jtem, info;
    int factor, k;  
    double rand_num;
    int N = params.N;
    
    srand48(seed);
    
    for(k=0; k<NumSet; k++){
        item = 0; 
        while(item<N){
            rand_num = drand48();
            info = (int)(N*rand_num);
            factor = 0;
            // Below we check for repeated unit index only for current sequence, ignoring previous sequences.
            while(factor == 0){
                factor = 1;
                for(jtem=0; jtem<item; jtem++) {
                    if(Permut[k][jtem] == info){
                        info = (info + 1)-(int)((info+1)/N)*N; // just equal to (info+1) mod N
                        jtem = item;
                        factor = 0;
                    }
                }
            }
            Permut[k][item] = info;
            item++;
        }
    }
    std::cout<< "After SetupTables." << std::endl;   
}







