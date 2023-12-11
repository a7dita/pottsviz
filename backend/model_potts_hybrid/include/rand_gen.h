/* For random numer generators */
#ifndef RANDGEN_H
#define RANDGEN_H

#ifndef RANLXD_C
    // Seed the generator: level is 1 or 2, seed is between 1 and 2^31-1
    void extern rlxd_init(int level, int seed); 

    // uniform distribution from [0, 1)
    void extern ranlxd(double r[], int n); 

     // Normal distribution, N(0,1). Calls ranlxd inside.
    void extern gauss_dble(double rd[],int n);

    // Exponential dist,  PDF=exp(-x/alpha)/alpha,
    // mean = alpha, var = alpha^2, Calls ranlxd inside.
    void extern wexpo_dble (double r[], int n, double alpha);
#endif /* RANLXD_C */

#endif /* RANDGEN_H */
