
/********************************************************************************
* Make patterns and save
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <cassert>
#include <math.h>

// ------------------------------------------------------------------------
typedef struct {
    /* define a structure*/
    int index;
    double field;
    int state;
} mystruct;
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
    srand48(seed); 
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

//------------------------------------------------------------------------
int main(int argc, char *argv[]){
    
    int N = atoi(argv[1]);
    int p = atoi(argv[2]);
    int S = atoi(argv[3]);
    double a = atof(argv[4]);
    int rand_seed = atoi(argv[5]);

    assert(p>1);
    assert(N>10);
    assert(S>1);
    assert(a>0.0 && a<=1.0);
    assert(rand_seed>0);
    
    
    int *xi = new int[N*p];
    make_random_memory(N, p, S, a, 0, xi, rand_seed);
    
    std::ofstream file;
    std::ofstream & file_ref = file;
    char fname[0x100];
    snprintf(fname, sizeof(fname), "./pat_N%d_p%d_S%d_a%.2f_seed%d", 
        N, p, S, a, rand_seed);
    file.open(fname, std::ios_base::out);
    save_patterns(N, p, xi, file_ref);
    file.close();
    
    delete[] xi;
    return 0;

}



