#include "window.h"
#include "head.h"

int sys_solve (double *a, int *ind, double *b, double *x, double *u, double *v,
               double *r /*additional memory*/, double *buf
               /*one for all, length == p*/, double eps, int p, int my_num,
               pthread_barrier_t * barrier)
{
    //x[i] == 0 for i from [0, ind[0]-1 ]
    int max_it = 50, res = 0, it = 0;

//    print_bl(b, ind[0]-1, 1, 1);
    int len = ind[0]-1;
    double norm = 0;
    buf[my_num] = 0;
    for(int i = (my_num*len)/p; i < ((my_num+1)*len)/p; i++)
    {
        if(buf[my_num] < fabs(b[i]))
        {
            buf[my_num] = fabs(b[i]);
        }
    }
    pthread_barrier_wait(barrier);
    for(int i = 0; i < p; i++)
    {
        if(buf[i] > norm)
        {
            norm = buf[i];
        }
    }
    pthread_barrier_wait(barrier);
//    if(my_num == 0)
//    {
//        memset(buf, 0, p*sizeof(double));
//    }
    eps *= norm;

    for(it = 0; it < MAX_IT; it += max_it)
    {
        res = solve_stage(a, ind, b, x, u, v, r /*additional memory*/,
                          buf /*one for all, length == p*/, max_it, eps, p,
                          my_num, barrier);
        if(res >= 0)
        {
            break; //solve finded
        }
    }
    if(it >= MAX_IT)
    {
        return -1;//solve not finded
    }
    return it;
}

int solve_stage(double *a, int *ind, double *b, double *x, double *u, double *v,
                double *r/*additional memory*/,
                double *buf /*one for all, length == p*/, int max_it, double eps,
                int p, int my_num, pthread_barrier_t * barrier)
{
    double c1 = 0, c2 = 0, tau = 0, len = ind[0]-1;
    int it = 0;

    //r = A*x-b
    //r = Ax
    msr_matr_mult_vector(a, ind, x, r, p, my_num);
//    pthread_barrier_wait(barrier);
//    print_bl(r, 1, len, len);
    //r = r - b
    linear_combination(r, b, len, -1, p, my_num);
//    pthread_barrier_wait(barrier);
//    print_bl(r, 1, len, len);

    for(it = 0; it < max_it; it++)
    {
        //v = [M^(-1)]*r, M - diag( A )
        apply_preconditioner(a, r, v, len, p, my_num);
        pthread_barrier_wait(barrier);
//        print_bl(v, 1, len, len);

        //u = A*v
        msr_matr_mult_vector(a, ind, v, u, p, my_num);
        pthread_barrier_wait(barrier);
//        print_bl(u, 1, len, len);

        //c1 = <u,r>
        c1 = scalar_product(u, r, len, buf, p, my_num, barrier);
        pthread_barrier_wait(barrier);
//        cout<<"c1 == "<<c1<<"\n";

        //c2 = <u,u>
        c2 = scalar_product(u, u, len, buf, p, my_num, barrier);
        pthread_barrier_wait(barrier);
//        cout<<"c2 == "<<c2<<"\n";

        //case for exit
        if(fabsl(c1) < eps*eps || fabsl(c2) < eps*eps)
        {
            break;
        }

        tau = c1/c2;
//        cout<<"tau == "<<tau<<"\n";

        //x = x - tau*v
        linear_combination(x, v, len, -tau, p, my_num);
//        pthread_barrier_wait(barrier);
//        print_bl(x, 1, len, len);

        //r = r - tau*u
        linear_combination(r, u, len, -tau, p, my_num);
        pthread_barrier_wait(barrier);
//        print_bl(r, 1, len, len);
    }
    if(it >= max_it)
    {
        return -1;
    }
    return it;
}

void msr_matr_mult_vector(double *a, int *ind, double *x, double *v, int p,
                          int my_num)
{
    int size = ind[0]-1;
    int i1 = (my_num*size)/p, i2 = ((my_num+1)*size)/p;
    int len = 0, start_elem = 0;
    double sum = 0;
    for(int i = i1; i < i2; i++)
    {
//        cout<<"i == "<<i<<"\n";
        sum = a[i]*x[i];
        len = ind[i+1]-ind[i];
        start_elem = ind[i];
        for(int j = 0; j < len; j++)
        {
           sum += a[start_elem+j]*x[ind[start_elem+j]];
        }

//        cout<<a[i]<<" ";
//        for(int j = 0; j < len; j++)
//        {
//           cout<<a[start_elem+j]<<" ";
//        }
//        cout<<"\n";


//        cout<<i<<" ";
//        for(int j = 0; j < len; j++)
//        {
//           cout<<ind[start_elem+j]<<" ";
//        }
//        cout<<"\n";

//        cout<<x[i]<<" ";
//        for(int j = 0; j < len; j++)
//        {
//           cout<<x[ind[start_elem+j]]<<" ";
//        }
//        cout<<"\n";
        v[i] = sum;
    }
}

void linear_combination(double *v, double *b, double len, double koef, int p,
                        int my_num)
{
    int i1 = (my_num*len)/p, i2 = ((my_num+1)*len)/p;
    for(int i = i1; i < i2; i++)
    {
        v[i] = v[i] + koef*b[i];
    }
}

void apply_preconditioner(double *a, double *v, double *u, int size, int p,
                          int my_num)
{
    int i1 = (my_num*size)/p, i2 = ((my_num+1)*size)/p;
    for(int i = i1; i < i2; i++)
    {
        u[i] = v[i] / a[i];
    }
}

double scalar_product(double *u, double *v, int len, double *buf, int p,
                      int my_num, pthread_barrier_t *barrier)
{
    double sum = 0;
    int i1 = (my_num*len)/p, i2 = ((my_num+1)*len)/p;
    for(int i = i1; i < i2; i++)
    {
        sum += u[i]*v[i];
    }
    buf[my_num] = sum;
    pthread_barrier_wait(barrier);
    sum = 0;
    for(int i = 0; i < p; i++)
    {
        sum += buf[i];
    }
    return sum;
}
