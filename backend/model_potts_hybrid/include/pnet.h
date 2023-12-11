/*
* Prototypes of class PNet is defined here
*/
#ifndef PNET_H
#define PNET_H

#ifndef FUNCTION_H
#include "functions.h"
#endif /* FUNCTION_H */
// ------------------------------------------------------------------------
class PNet{
    public:
        PNet(Potts_params *params, int NumSet, int seed);
        virtual ~PNet();
        void initialise();
        void initialise(int cue, const int *xi);

        void make_dilution_matrix();
        void make_J_null();
        void make_J_assoc(int size, const int *eta, double b, double lambda);
        
        void SetupTables(int seed);
        
        void update_unit(int i, int flag, double *f1, double *f2);
        void update_network(int n, int flag, double *f1, double *f2);
        
        void compute_m(const int *xi, double *m);
        void compute_mean_activity();
        double mean_activity;
        
        friend class Network_runner;
    
    protected: // functions that are listed here will be overriden.
        void virtual compute_field(int i);
        
        Potts_params params; // network parameters
        int NumSet; // number of random sequences for updating orders of units 
        double *s; // Potts spin variables
        double *s0; // quiet variables
        double *theta; // state-specific thresholds 
        double *theta0A; //fast GABA A inhibition
        double *theta0B; //slow GABA B inhibition
        double *h; // local field
        double *J; // connectivity, flattened
        double *r;	// membraine potential
        int *Permut; // updating order of units.
        int *C; // dilution matrix
        double *weight; // for computing activatin function
        

};
// ------------------------------------------------------------------------
class Network_runner{
    public:
        Network_runner();
        ~Network_runner();

        void make_null_connection(PNet *post_net, PNet *pre_net, double * Jmat);
        void make_Hebb_connection(PNet *post_net, PNet *pre_net, double *Jmat,
            const int *xi_post, const int *xi_pre, double lambda, int *match_table);

        void run_two_nets(PNet *p_net, PNet *f_net, double *Jf2p, double *Jp2f,
            const int *xi_p, const int *xi_f, 
            std::ofstream & buf1, std::ofstream &buf2, int cue, int save_all);
        void compute_field(PNet *post_net, PNet *pre_net, const double *J
            , double *field);
        void make_field_decay(PNet *net_ptr, double *field, int cue, const int *xi
            ,int t);
        int Runs; // number of updates
        int nCue;
    private:
        // INSERT PRIVATE MEMEBERS HERE
        double g = 5.0;
        double tau = 40.0;
        int dt = 10; // snapshot interval
};
#endif /* PNET_H */
