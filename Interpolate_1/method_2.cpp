#include "window.h"

/*
 * Piecemeal interpolation by cubic Hermite polynomials with the definition
 * of missing boundary conditions using an additional node in the boundary
 * nodes.
 */

void find_koef_2(int n, double *f, double *fd, double *a, double delta)
{
    /*********************************
     * массив а имеет длину 4n
     * массив f имеет длину n+1
     * массив fd имеет длину n+1
     * delta - приращение по х
     *********************************/

    for( int i = 0; i < n; i ++)
    {
        a[i*4]   = f [i];
        a[i*4+1] = fd[i];
        a[i*4+2] = (((f[i+1]-f[i])/delta) - fd[i])/delta;
        a[i*4+3] = (fd[i] + fd[i+1] - 2*((f[i+1]-f[i])/delta))/(delta*delta);
//        printf("a_1%d = %e\n", i, a[i*4]);
//        printf("a_2%d = %e\n", i, a[i*4+1]);
//        printf("a_3%d = %e\n", i, a[i*4+2]);
//        printf("a_4%d = %e\n", i, a[i*4+3]);
    }
}

double find_value_2(double y, double *x, double *a, double delta, double le)
{
    /*********************************
     * массив а имеет длину 4n
     * массив x имеет длину n+1
     * y - точка в которой ищем значение
     * delta - приращение по х
     *********************************/

    int i = (y-le)/delta;

    return a[i*4] + a[i*4+1]*(y-x[i]) + a[i*4+2]*(y-x[i])*(y-x[i]) + a[i*4+3]*(y-x[i])*(y-x[i])*(y-x[i+1]);
}
