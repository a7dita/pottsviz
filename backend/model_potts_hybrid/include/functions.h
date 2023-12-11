#include <iostream>
#include <fstream>
/*
* Functions used in other files are listed below
*/
#ifndef FUNCTION_H
#define FUNCTION_H
// ------------------------------------------------------------------------
typedef struct {
    /* define a structure*/
    int index;
    double field;
    int state;
} mystruct;
// ------------------------------------------------------------------------
struct Potts_params{
    int N; // number of units in a net
    int Cm; // mean connections per neuron
    int S;  // number of Potts states
    int p; // number of patterns: for non-associative memories, p=-1
    double a; // mean activity to maintain
    double U;  // global threshold
    double T1; // integration of incoming fields
    double T2; // adaptation
    double T3A; //fast inhibition
    double T3B; //slow inhibition
    double beta; // inverse temperature
    double w; // self-reinforcement term
    double gammaA; // fast inhibition
};

#ifndef FUNCTION_C
extern int cmpf_descend (const void *x, const void *y);
extern int cmpf_ascend (const void *x, const void *y);
extern int cmp_dbl_descend (const void *x, const void *y);
extern void shuffle(int * array, int size);
extern int random_sample(double vec[], int size);
extern void make_Cij(int i, int nMin, int nMax, int n, int *result);

extern void save_patterns(int N, int p, const int * xi, std::ofstream & outfile);
extern void read_patterns(int N, int p, int *xi, std::ifstream & infile);
extern void read_doubles(int N, int p, double * vec, std::ifstream & file_ref);
extern void save_doubles(int size, double * vec, std::ofstream & file_ref);
extern void make_random_memory(int N, int p, int S, double b, int index, int * xi, int seed);
extern void make_flag(int p, int Z, const int *table, int K, int *flag);
extern void get_units_of_max_field(int N, int S, const double *h,int Na, int *index_array, int *state_array);
extern int get_max_index(int size, double * vec);
extern void display_overlaps(int n_to_show, int size, double * vec);
extern void save_overlaps(int n_to_show, int size, double * vec, std::ofstream &buf);
extern int save_overlaps2(int n_to_show, int size, double * vec, std::ofstream &buf);
extern void make_match_internet(int *match_table, int flag, int Z, int p);

#endif /* FUNCTION_C */

#endif /* FUNCTION_H*/

