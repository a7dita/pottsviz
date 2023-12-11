/******************************************************************
* Kwang Il Ryom, March 2nd 2023
* Functions that do not use the object "pnet_class".
*******************************************************************/
#define FUNCTION_C

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "functions.h"

#ifndef RANDGEN_H
#include "rand_gen.h"
#endif /* RANDGEN_H */
// ------------------------------------------------------------------------
/* global variable defined in other files:
* I do not use any global variables
* example: extern int Z;
*/

// ------------------------------------------------------------------------
int cmpf_descend (const void *x, const void *y)
{
    /* comparing function for qsort(), decreasing order
    * use "mystruct" above
    */
    double xx = (*(mystruct*)x).field;
    double yy = (*(mystruct*)y).field;
    if (xx > yy) return -1;
    if (xx < yy) return  1;
    return 0;
}
// ------------------------------------------------------------------------
int cmpf_ascend (const void *x, const void *y)
{
    /* comparing function for qsort(), increasing order
       * use "mystruct" above
    */
    double xx = (*(mystruct*)x).field;
    double yy = (*(mystruct*)y).field;
    if (xx > yy) return 1;
    if (xx < yy) return  -1;
    return 0;
}
// ------------------------------------------------------------------------
int cmp_dbl_descend (const void *x, const void *y){
    /* comparing function for qsort(), for double, decreasing order */
    double xx = *(double*)x;
    double yy = *(double*)y;
    if(xx>yy) return -1;
    if(xx<yy) return 1;
    return 0;
}
// ------------------------------------------------------------------------
void shuffle(int * array, int size)
{
    /* 
    * shuffle the 1D "array" of "size" 
    * if the "size" is smaller than the actual length of array, the last elements are not touched.
    * if the "size" should NOT be larger than actual size of array.
    */
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
int random_sample(double vec[], int size){
    /*
        Given a weight vector, randomly sample one index among {0,1,...,size-1}
    */
    // check if the elements are non-negative
    int i;
    for(i=0;i<size;++i){
        if(vec[i]<0.0)
            return -1;
    }
    double sum = 0.0;
    double rnum;
    for(i=0;i<size;++i)
        sum += vec[i];
    
    rnum = sum*drand48(); // a number between [0, sum)
    for(i=0;i<size;++i){
        if(rnum < vec[i]) 
            return i;
        else
            rnum -= vec[i];
    }
    return -2; // dummy
}

// ------------------------------------------------------------------------
void make_random_memory(int N, int p, int S, double b, int index, int * xi, int seed){
    /* M random patterns by using qsort() and mystruct
    * Parameters:
        p = number of memories 
        index = starting index for storing memories in "xi"
        seed = random seed
        N = number of neurons
        S = number of Potts states
        b = sparsity of patterns
    * Return:
        xi = flat array of integers, p*N in size
    */
    int i, mu;
    srand48(seed+1); 
    int Na = (int) (N*b);
    for(mu=index;mu<p+index;mu++){   
        mystruct X[N];
        for(i=0;i<N;i++){
            //first initialize each pattern to a null (quiet) state
            *(xi+mu*N+i) = S;
            X[i].index = i;
            X[i].field = drand48();
            X[i].state = (int)((double)S*drand48());
        }
        
        //sort the patterns in order of decreasing field keeping the indices
        qsort(X, N, sizeof(mystruct), cmpf_descend);  

        // set N*a largest of them to become active
        for(i=0;i<Na;i++)
            *(xi+mu*N+X[i].index) = X[i].state;  
    }
    std::cout << "---- After making random patterns." << std::endl;
}
// ------------------------------------------------------------------------
void save_patterns(int N, int p, const int *xi, std::ofstream & outfile){
    /* saving patterns to a file 
    * N = number of Potts units
    * p = number of patterns
    */
    int mu,i;    
    for(mu=0;mu<p;mu++){  
        for(i=0;i<N;i++)
            outfile << *(xi+mu*N+i) << "\t";
        outfile << std::endl;
    }   
}
// ------------------------------------------------------------------------
void read_patterns(int N, int p, int *xi, std::ifstream & infile){
    /* reading patterns from a file 
    * N = number of Potts units
    * p = number of patterns
    */
    int mu,i;
    int number;
    for(mu=0;mu<p;mu++){  
        for(i=0;i<N;i++){
            infile >> number;
            *(xi+mu*N+i) = number;
        }
    }   
}
// ------------------------------------------------------------------------
void read_doubles(int N, int p, double * vec, std::ifstream & file_ref){
    /* reading double variables from a file 
    * N = number of putative columns
    * p = number of putative rows
    */
    int mu,i;
    double x;
    for(mu=0;mu<p;mu++){  
        for(i=0;i<N;i++){
            file_ref >> x;
            *(vec+mu*N+i) = x;
        }
    }   
}
// ------------------------------------------------------------------------
void save_doubles(int size, double * vec, std::ofstream & file_ref){
    /*
    Save double variables to file
    */
    int mu;
    for(mu=0;mu<size;mu++){  
        file_ref << *(vec+mu) << "\t";
    }
    file_ref << std::endl;
}
//---------------------------------------------------------------------------------------
void get_units_of_max_field(int N, int S, const double *h,
                            int Na, int *index_array, int *state_array){
    /*
    Get some units and their states that receive the largest fields
    N = number of units
    S = number of Potts states
    h = flat array of N*S fields
    Na = number of units to pick up
    Return: index_array, state_array
    */
    int i, k;
    mystruct X[N];
    int kmax = -1;
    double hmax = -100000.0;
    for(i=0;i<N;i++){
        kmax = -1;
        hmax = -100000.0;
        for(k=0;k<S;++k){
            if( *(h+i*S+k)>hmax ){
                hmax = *(h+i*S+k);
                kmax = k;
            }
        }
        X[i].index = i;
        X[i].field = hmax;
        X[i].state = kmax;
    }
     //sort the units in order of decreasing field keeping the indices
    qsort(X, N, sizeof(mystruct), cmpf_descend); 
    
    for(i=0;i<Na;++i){
        *(index_array+i) = X[i].index;
        *(state_array+i) = X[i].state;
    }
}
//---------------------------------------------------------------------------------------
void make_Cij(int i, int nMin, int nMax, int n, int *result){
    /*
    Pick $n$ random indices between [nMin, nMax), avoiding $i$
    */
    assert(n>0);
    assert(n+nMin<=nMax);
    int j;
    int array[nMax-nMin];
    int counter = 0;
    for(j=nMin;j<nMax;j++){
        if (j!=i){
            array[counter]=j;
            ++counter;
        }
    }
    shuffle(array, counter);
    for(j=0;j<n;j++) 
        result[j]=array[j];
}
//--------------------------------------------------------------------------------
int get_max_index(int size, double * vec){
    /*
    return the index (first instance) of a maximum value
    */
    double max_val = -100000.0;
    int index = -1;
    for(int i=0;i<size;++i){
        if( *(vec+i)>max_val ){
            index = i;
            max_val = *(vec+i);
        }
    }
    return index;
}
//---------------------------------------------------------------------------
void display_overlaps(int n_to_show, int size, double * vec){
    /*
    Show maximum overlaps on the screen
    --------
    Parameters:
    n_to_show = how many overlaps to show
    size = length of the vector
    vec = the overlap values
    */
    mystruct X[size];
    for(int mu=0;mu<size;++mu){
        X[mu].index = mu;
        X[mu].field = vec[mu];
        X[mu].state = 0;
    }
    qsort(X, size, sizeof(mystruct), cmpf_descend);
    //std::cout << "t = " << t << ", ";
    for(int k=0;k<n_to_show;++k){
        std::cout<<"("<<X[k].index << ", " << X[k].field << ") ";
    }
    
    std::cout << std::endl;
}

//---------------------------------------------------------------------------
void save_overlaps(int n_to_show, int size, double * vec, std::ofstream &buf){
    /*
    save maximum overlaps to the "buf"
    --------
    Parameters:
    n_to_show = how many overlaps to show
    size = length of the vector
    vec = the overlap values
    */
    mystruct X[size];
    for(int mu=0;mu<size;++mu){
        X[mu].index = mu;
        X[mu].field = vec[mu];
        X[mu].state = 0;
    }
    qsort(X, size, sizeof(mystruct), cmpf_descend);
    for(int k=0;k<n_to_show;++k){
        buf << X[k].index << "\t" << X[k].field << "\t";
    }
    
    //buf << std::endl;
}
//----------------------------------------------------------------------------
int save_overlaps2(int n_to_show, int size, double * vec, std::ofstream &buf){
    /*
    save maximum overlaps to the "buf" and return maximum-overlap index
    --------
    Parameters:
    n_to_show = how many overlaps to show
    size = length of the vector
    vec = the overlap values
    */
    mystruct X[size];
    for(int mu=0;mu<size;++mu){
        X[mu].index = mu;
        X[mu].field = vec[mu];
        X[mu].state = 0;
    }
    qsort(X, size, sizeof(mystruct), cmpf_descend);
    for(int k=0;k<n_to_show;++k){
        buf << X[k].index << "\t" << X[k].field << "\t";
    }
    
    return X[0].index;
}
//---------------------------------------------------------------------------
void make_match_internet(int *match_table, int flag, int Z, int p){

    assert(Z*Z<=p);
    assert(flag < 5);

    // initialize the table with zeros
    for(int i=0;i<p;++i){
        for(int j=0;j<p;++j){
            match_table[i*p+j] = 0;
        }
    }
    // flag 0 - connect memorie one to one.
    if (flag == 0) {
        for(int i=0;i<p;++i){
            int j = i;
            match_table[i*p+j] = 1;
        }
    }
    // flag 1 - connecting Z patterns of pre network with Z non-overlapping patterns of post network
    if (flag == 1) {
        int pZ = 10;
        assert(pZ*Z<p);
        for(int i=0;i<pZ;++i){
            for(int j=0;j<Z;++j){
                match_table[i*p + i*Z + j] = 1;
            }
        }
    }
    // flag 2 - ring
    // connecting Z*Z patterns of pre network with Z*Z patterns of post network - each having Z unique connections
    if (flag == 2) {
        for (int i = 0; i < Z*Z; ++i) {
            for (int j = i; j < i+Z; ++j) {
                match_table[i*p + (j % (Z*Z))] = 1;
            }
        }
    }
    // flag 3 - cluster
    // connecting Z*Z patterns of pre network with Z*Z patterns of post network - each having Z unique connections
    if (flag == 3) {
        for (int i0 = 0; i0 < Z; ++i0) {
            for (int i1 = 0; i1 < Z; ++i1) {
                for (int j0 = 0; j0 < Z; ++j0) {
                    for (int j1 = 0; j1 < Z; ++j1) {
                        int i = i0 + Z * i1;
                        int j = j0 + Z * j1;
                        if (i0 == (j0 + j1 * i1) % Z) {
                            match_table[i*p + j] = 1;
                        }
                    }
                }
            }
        }
    }
    // flag 4 - ring
    // connecting Z*Z patterns of pre network with Z*Z patterns of post network - each having Z unique connections
    if (flag == 4) {
        int rep = (int)((p/2)/ (Z*Z));
        for (int r = 0; r < rep; ++r) {
            for (int i = 0; i < Z*Z; ++i) {
                for (int j = i; j < i+Z; ++j) {
                    match_table[(Z*Z * r + i) * p + (Z*Z *r +(j % (Z*Z)))] = 1;
                }
            }
        }
    }

}
