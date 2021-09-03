#include "window.h"

//Approximations by Chebyshev polynomials using the least squares method.

void find_koef_1(int kol_nodes, double *y, double *f, double *koef)
{
    double a=0,b=0,c=0,d=0,d_last=0,x1=0,x2=0, ui=0;
    int N = kol_nodes - 1;
    memset(koef, 0, kol_nodes*sizeof (double));
//    for(int i = 0; i < kol_nodes; i++)
//    {
//        printf("alfa[%d] = %lf\n", i, koef[i]);
//    }

    for(int i = 0; i < kol_nodes; i++)
    {
        for(int j = 0; j < N; j ++)
        {
            x1 = (2*y[j]  -y[N]-y[0])/(y[N]-y[0]);
            x2 = (2*y[j+1]-y[N]-y[0])/(y[N]-y[0]);
//            printf("x1 = x[%d] = %lf\n", j, x1);
//            printf("x2 = x[%d] = %lf\n", j+1, x2);
            if(i == 0)
            {
                d_last = d;
                a = -(acos(x2) - acos(x1));
                b = -(sin(acos(x2)) - sin(acos(x1)));
                c = (a*x2 - b)/(x2-x1);
                d = (b-a*x1)/(x2-x1);
            }
            else if(i == 1)
            {
                d_last = d;
                a = -(1./i)*(sin(i*acos(x2))-sin(i*acos(x1)));
                b = -0.5*( acos(x2) + 0.5*sin(2*acos(x2))
                          -acos(x1) - 0.5*sin(2*acos(x1)));
                c = (a*x2 - b)/(x2-x1);
                d = (b-a*x1)/(x2-x1);
            }
            else
            {
                d_last = d;
                a = -(1./i)*(sin(i*acos(x2))-sin(i*acos(x1)));
                b = -0.5*( (1./(i-1))*(sin((i-1)*acos(x2)) - sin((i-1)*acos(x1)))
                          +(1./(i+1))*(sin((i+1)*acos(x2)) - sin((i+1)*acos(x1))) );
                c = (a*x2 - b)/(x2-x1);
                d = (b-a*x1)/(x2-x1);
            }

//            printf("a^%d_%d = %lf\n", i,j,a);
//            printf("b^%d_%d = %lf\n", i,j,b);
//            printf("c^%d_%d = %lf\n", i,j,c);
//            printf("d^%d_%d = %lf\n\n", i,j,d);

            if(j==0)
            {
//                printf("c^%d_%d = %lf\n\n", i,j,c);
                ui = c;
            }
            else
            {
//                printf("c^%d_%d = %lf\n", i,j,c);
//                printf("d^%d_%d = %lf\n\n", i,j-1,d_last);
                ui = c + d_last;
            }
//            printf("u^%d_%d = %lf\n", i, j, ui);
//            printf("f[%d] = %lf\n", j, f[j]);
            koef[i] += ui*f[j];
//            printf("alfa[%d] = %lf\n\n", i, koef[i]);
        }
        ui = d;
//        printf("u^%d_%d = %lf\n", i, N, ui);
//        printf("f[%d] = %lf\n", N, f[N]);
        koef[i] += ui*f[N];
//        printf("alfa[%d] = %lf\n\n", i, koef[i]);
    }

//    printf("alfa[0] = %lf \n",koef[0]);
    koef[0] *= (1./PI/*180*/);
//    printf("alfa[0] = %lf \n",koef[0]);
    for(int i = 1; i < kol_nodes; i ++)
    {
//        printf("alfa[%d] = %lf \n",i,koef[i]);
        koef[i] *= (2./PI/*180*/);
//        printf("alfa[%d] = %lf \n",i,koef[i]);
    }
}

double find_value_1(double y, double a, double b, int kol_nodes, double *koef)
{
    double z = (2*y-b-a)/(b-a);
    double t0 = 1, t1 = z, ti=0;
    double value = 0;
//    int n_little = kol_nodes - 2;
    int N = kol_nodes - 1;
//printf("%lf\n",y);
    if(kol_nodes == 1)
    {
        value = koef[0]*t0;
//        printf("%lf\n",koef[0]);
    }
    else if(kol_nodes >= 2)
    {
        value = koef[0]*t0 + koef[1]*t1;
//        printf("%lf\n",koef[0]);
//        printf("%lf\n",koef[1]);
    }
    for( int i = 2; i < N; i ++)
    {
        ti = (2*z*t1) - t0;
        value += koef[i]*ti;
        t0 = t1;
        t1 = ti;
//        printf("%lf\n",koef[i]);
    }
    return value;
}


//for(int j = 0; j < kol_nodes; j ++)
//{
//    printf("y[%d] = %lf\n", j, y[j]);
//}
//for(int j = 0 ; j < kol_nodes; j ++)
//{
//    printf("%d\n", j);
//    if(j!=N)
//    {
//        x1 = (2*y[j]  -y[N]-y[0])/(y[N]-y[0]);
//        x2 = (2*y[j+1]-y[N]-y[0])/(y[N]-y[0]);
//        printf("x1 = x[%d] = %lf\n", j, x1);
//        printf("x2 = x[%d] = %lf\n", j+1, x2);
//    }
//    for(int i = 0; i < kol_nodes; i++)
//    {
////            if(j!=N)
////            {
//            if(i == 0)
//            {
//                a = -(acos(x2) - acos(x1));
//                b = -(sin(acos(x2)) - sin(acos(x1)));
//                c = (a*x2 - b)/(x2-x1);
//                d = (b-a*x1)/(x2-x1);
//            }
//            else if(i == 1)
//            {
//                d_last = d;
//                a = -(1./i)*(sin(i*acos(x2))-sin(i*acos(x1)));
//                b = -0.5*( acos(x2) + 0.5*sin(2*acos(x2))
//                          -acos(x1) - 0.5*sin(2*acos(x1)));
//                c = (a*x2 - b)/(x2-x1);
//                d = (b-a*x1)/(x2-x1);
//            }
//            else
//            {
//                d_last = d;
//                a = -(1./i)*(sin(i*acos(x2))-sin(i*acos(x1)));
//                b = -0.5*( (1./(i-1))*(sin((i-1)*acos(x2)) - sin((i-1)*acos(x1)))
//                          +(1./(i+1))*(sin((i+1)*acos(x2)) - sin((i+1)*acos(x1))) );
//                c = (a*x2 - b)/(x2-x1);
//                d = (b-a*x1)/(x2-x1);
//            }

//            printf("a^%d_%d = %lf\n", i,j,a);
//            printf("b^%d_%d = %lf\n", i,j,b);
//            printf("c^%d_%d = %lf\n", i,j,c);
//            printf("d^%d_%d = %lf\n", i,j,d);

//            if(j == 0)
//            {
//                ui = c;
//            }
//            else if(j!=N)
//            {
//                ui = c + d_last;
//            }
//            else
//            {
//                ui = d;
//            }
////            }
////            else
////            {
////                ui = d;
////            }
//        printf("u^%d_%d = %lf\n", i, j, ui);
//        printf("f[%d] = %lf\n\n", j, f[j]);
//        koef[i] += ui*f[j];
//    }
//}
