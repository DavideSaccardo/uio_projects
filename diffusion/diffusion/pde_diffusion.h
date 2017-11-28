#ifndef PDE_DIFFUSION_H
#define PDE_DIFFUSION_H


class pde_diffusion
{
public:
    int N;
    double alpha;
    double r_cn, r_i, r_e, r_2d;
    void Crank_nicolson(double* u);
    pde_diffusion(int n, double a);
    void Implicit(double *u);
    void Explicit(double *u);
    void two_dimension(double **u, double tolerance, int cutoff);
};

#endif // PDE_DIFFUSION_H
