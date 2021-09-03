#include "window.h"
#include "head.h"

double find_error(double *x, double (*f) (double, double), double hx, double hy,
                  int cnx, int dnx, int cny, int dny, int my_num, int p,
                  double *buf, /*double x_0, double y_0,*/
                  pthread_barrier_t *barrier, double f_max, int disturb)
{
    int n = 2*cnx+dnx, m = 2*cny+dny;
    int N = (n+1)*(m+1) - (dnx-1)*(dny-1);
    int k1 = (my_num*N)/p, k2 = ((my_num+1)*N)/p;
    int i = 0, j = 0;
    int e1 = 0, e2 = 0, e3 = 0;
    double err = 0, er1 = 0, er2 = 0, er3 = 0, er4 = 0;
    for(int k = k1; k < k2; k++)
    {
        k2ij(n, k, i, j, cnx, dnx, cny, dny);
    /*inside domain*/
        if(i > 0 && i < n && j > 0 && j < m
                &&((i < cnx || i > cnx+dnx)
                   || (i >= cnx && i <= dnx+cnx && (j < cny || j > cny+dny)) ))
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }
    /*left vertical lines of domain*/
        else if( (i == 0 && j > 0 && j < m)
                || (i == cnx+dnx && j > cny && j < cny+dny) )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }
    /*right vertical lines of domain*/
        else if( (i == n && j > 0 && j < m)
                || (i == cnx && j > cny && j < cny+dny) )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
        }
    /*lower horizontal lines of domain*/
        else if( (j == 0 && i > 0 && i < n)
                || (j == cny+dny && i > cnx && i < cnx+dnx) )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }
    /*upper horizontal lines of domain*/
        else if( (j == m && i > 0 && i < n)
                || (j == cny && i > cnx && i < cnx+dnx) )
        {
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }
    /*angles of domain*/
        else if( i == cnx+dnx && j == cny )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == cnx && j == cny+dny )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == cnx && j == cny )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == cnx+dnx && j == cny+dny )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == 0 && j == 0)
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j) - f_max*disturb);

            err = er1 > err ? er1 : err;
            err = er2 > err ? er2 : err;
            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == n && j == m )
        {
            er4 = fabs( x[k]          - f (hx*i,       hy*j));
        }

        else if( i == n && j == 0 )
        {
            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);

            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er1 > err ? er1 : err;
            err = er4 > err ? er4 : err;
        }

        else if( i == 0 && j == m )
        {
            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));
            er4 = fabs( x[k]          - f (hx*i,       hy*j));

            err = er3 > err ? er3 : err;
            err = er4 > err ? er4 : err;
        }
    }

    buf[my_num] = err;
    pthread_barrier_wait(barrier);

    for(int k = 0; k < p; k ++)
    {
        err = buf[k] > err ? buf[k] : err;
    }
    return err;
}














//        k2ij(n, k, i, j, cnx, dnx, cny, dny);
//        if(j < m && i < n && )
//        {
//            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);
//            ij2k(n, m, i+1, j+1, e2, cnx, dnx, cny, dny);
//            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);

//            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));
//            er2 = fabs((x[k]+x[e2])/2 - f (hx*(i+0.5), hy*(j+0.5)));
//            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));

//            err = er1 > err ? er1 : err;
//            err = er2 > err ? er2 : err;
//            err = er3 > err ? er3 : err;
//        }
//        else if(j == m && i < n)
//        {
//            ij2k(n, m, i+1, j,   e3, cnx, dnx, cny, dny);
//
//            er3 = fabs((x[k]+x[e3])/2 - f (hx*(i+0.5), hy*j));

//            err = er3 > err ? er3 : err;
//        }
//        else if(i == n && j < m)
//        {
//            ij2k(n, m, i,   j+1, e1, cnx, dnx, cny, dny);

//            er1 = fabs((x[k]+x[e1])/2 - f (hx*i,       hy*(j+0.5)));

//            err = er1 > err ? er1 : err;
//        }
