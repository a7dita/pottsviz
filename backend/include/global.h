
#include <fstream>
#include <iostream>
/* 
************ Making a random vector for time constants *******
flag=0: gauss-->mean=0, param= std
flag=1: uniform-->range=[0,param)
flag=2: exponential-->PDF=exp(-x/param)/params, mean=param, var=param^2
*/

void extern make_rand_vector(double *T, int n, double param, int seed);

/******************************************************************

Definitions of class and struct.
    -Class PNet_timescale is inherited from PNet, only the function "run_dynamics" is overridden.

******************************************************************/

// ------------------------------------------------------------------------
struct Potts_params{
    int N; // number of units in a net
    int Cm; // mean connections per neuron
    int p; // number of patterns
    int S;  // number of Potts states
    double a;    // sparsity
    double U;     // global threshold
    double T1;   // integration of incoming fields
    double T2;  // adaptation
    double T3A; //fast inhibition
    double T3B; //slow inhibition
    double beta; // inverse temperature
    double w; // self-reinforcement term
    double gammaA; // fast inhibition
};
// ------------------------------------------------------------------------

//void set_params(Potts_params *params);
void make_random_patterns(Potts_params params, int seed, int *xi);
void make_random_patterns_V2(Potts_params params, int seed, int *xi);
void test_patterns(Potts_params params, const int *xi);
// ------------------------------------------------------------------------


class PNet{
    public:
        PNet(char name, Potts_params params);
        ~PNet();
        
        void initialise(); //initialising with static states
        void initialise(int mu, const int* xi); //initialise with one pattern
        void construct_CJ(const int *xi); // connecting units
        void update_state
            (int i, 
            int n, 
            int xi_i); 
            // xi_i is state at unit i wanted by cued pattern
        void run_dynamics(
            int cue, 
            int nUpdates, 
            int seed, 
            const int *xi, 
            const double *T2); // T2 is the table for tau2
        void SetupTables(int seed); // make random seqeunces for updating orders
        void snapshot(const int t, const int *xi); // snapshot the network state
        void snapshot(const int t, const int *xi, std::ofstream &buffer); // return also the overlaps, buffer is a refernce
        void compute_m(const int *xi, double *m); // compute overlaps
        char name; // name of the network
        double max_ovlp; // maximum overlap
        double max2_ovlp; // second overlap
        int max_mu; // index of a pattern with maximum overlap
        int max2_mu; // index of the second largest pattern
        double R; // average of squared overlaps, except the largest two.
    friend void compute_m(PNet *net_ptr, const int *xi, double *m); // compute overlaps for all patterns
    friend double compute_m(PNet *net_ptr, const int *xi, int mu); // compute overlap for one pattern
    friend class Network_runner;
    private:
        Potts_params params; // network parameters
        int NumSet; // number of random sequences for updating orders of units
        double **s; // spin variables
        double *h;	// field, now just for one unit.
        double **r;	// membraine potential
        double *theta0A; //fast GABA A inhibition
        double *theta0B; //slow GABA B inhibition
        double **theta;  //threshold
        double *J; // connectivity, flattened
        int **C; // connectivity
        int **Permut;; // updating order of units.
    
};
/*************************************************************************/
class Network_runner{
    public:
        Network_runner();
        ~Network_runner();
        void run_net(
            PNet *net,
            int start_cue,
            int Runs,
            int seed, 
            const int *xi
        );
        int dt;
        int nCue;
    private:
        
        // Nothing yet
};

