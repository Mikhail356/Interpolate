#include "window.h"
#include "head.h"

void *solve (void *input)
{
    Args *arg = (Args*) input;
    Storage *ss = arg->store;
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET( get_nprocs() - arg->my_num - 1,&cpu);
    pthread_setaffinity_np(pthread_self(),sizeof(cpu),&cpu);
//    int k = ss->k, f_change = ss->f_change;
    int cnx = 0, dnx = 0, cny = 0, dny = 0, n = 0, m = 0/*, nx = 0, ny = 0*/,
            N = 0, *ind = nullptr, p = 0, kol_tyk = 0, key_push = 0;
    double *a = nullptr, *b = nullptr, *x_new = nullptr, *u = nullptr,
            *r = nullptr, *v = nullptr, *buf = nullptr, len_x = 0, len_y = 0,
            residual = 0, eps = 0, f_max = 0;
    double (*f) (double, double) = nullptr;
//    a = ss->a;
//    ind = ss->ind;
    pthread_mutex_lock(ss->mutex_all);
    buf = ss->buf;
    p = ss->p;
    pthread_mutex_unlock(ss->mutex_all);

    while(true)
    {
        pthread_barrier_wait(arg->store->barrier_all);
        if(ss->status == -1)
        {
            break;
        }
        if(arg->my_num == 0)
        {
            ss->status = 1;
        }
//        pthread_mutex_lock(ss->mutex_all);
        residual = 0;
        cnx = ss->c*ss->nx;
        dnx = ss->d*ss->nx;
        cny = ss->c*ss->ny;
        dny = ss->d*ss->ny;
        len_x = ss->len_x;
        len_y = ss->len_y;
        eps = ss->eps;
        f = ss->f;
        kol_tyk = ss->disturb;
        f_max = ss->f_max;
        n = 2*cnx+dnx;
        m = 2*cny+dny;
        N = (n+1)*(m+1) - (dnx-1)*(dny-1);
        key_push = ss->key_push;
//        pthread_mutex_unlock(ss->mutex_all);
        pthread_barrier_wait(arg->barrier);
        if(arg->my_num == 0)
        {
            if(key_push!=1)
            {
                if(ss->a != nullptr)
                {
                    delete [] ss->a;
                    ss->a = nullptr;
                }
                if(ss->ind != nullptr)
                {
                    delete [] ss->ind;
                    ss->ind = nullptr;
                }
                allocate_MSR_matrix(n, m, &ss->a, &ss->ind, cnx,
                                    dnx, cny, dny);
            }
            ss->x_new = new double [N];
            ss->b = new double [N];
            ss->u = new double [N];
            ss->v = new double [N];
            ss->r = new double [N];
            memset(ss->x_new, 0, N * sizeof(double));
//            ss->key_push = 0;
        }
        pthread_barrier_wait(arg->barrier);
        x_new = ss->x_new;
        ind = ss->ind;
        a = ss->a;
        b = ss->b;
        u = ss->u;
        r = ss->r;
        v = ss->v;
//        buf = ss->buf;
//        pthread_mutex_unlock(ss->mutex_all);
//        pthread_barrier_wait(arg->barrier);
        if(key_push != 1)
        {
            build_MSR_matrix(n, m, a, ind, p /*# thread*/,
                             arg->my_num /*number of thread*/,
                             cnx, dnx, cny, dny, arg);
        }
        pthread_barrier_wait(arg->barrier);
        find_mm(f, f_max, len_x/n, len_y/m, n, m);
//        cout<<f_max<<endl;
        build_vector(n, m, p /*# thread*/, arg->my_num /*number of thread*/,
                     cnx, dnx, cny, dny, b, f, len_x/n, len_y/m, f_max,
                     kol_tyk);
//        pthread_barrier_wait(arg->barrier);
        sys_solve (a, ind, b, x_new, u, v, r /*additional memory*/, buf,
                   eps, p, arg->my_num, arg->barrier);
//        pthread_barrier_wait(arg->barrier);
        residual = find_error(x_new, f, len_x/n, len_y/m, cnx, dnx, cny,
                              dny, arg->my_num, p, buf, arg->barrier, f_max,
                              kol_tyk);
//        pthread_barrier_wait(arg->barrier);
        pthread_mutex_lock(ss->mutex_all);
        if(arg->my_num == 0)
        {
//            if(key_push != 1)
//            {
//                if(ss->a != nullptr)
//                {
//                    delete [] ss->a;
//                    ss->a = nullptr;
//                }
//                if(ss->ind != nullptr)
//                {
//                    delete [] ss->ind;
//                    ss->ind = nullptr;
//                }
//            }
            delete [] b; delete [] u; delete [] v; delete [] r;
            ss->residual = residual;
            ss->N_new = N;
            ss->x_new = x_new;
            ss->status = 0;
            ss->f_change = 1;
        }
        ss->key_push = 0;
        pthread_mutex_unlock(ss->mutex_all);
//        pthread_barrier_wait(arg->barrier);
    }
    if(arg->my_num == 0)
    {
        delete [] ss->x_new;
        delete [] ss->a;
        delete [] ss->ind;
    }
    return 0;
}

int check_file(char *name, int *c, int *d, double *len_x, double *len_y)
{
    FILE *file;
    file = fopen("%r", name);
    if(fscanf(file,"%d",c)
            ||fscanf(file,"%d",d)
            ||fscanf(file,"%lf",len_x)
            ||fscanf(file,"%lf",len_y))
    {
        return -1;
    }
    return 0;
}

void find_mm(double (*f)(double, double), double &f_max,
             double hx, double hy, int n, int m)
{
    double max1 = 0, max2 = 0;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            max1 = f(i*hx, j*hy);
            if(max1 > max2)
            {
                max2 = max1;
            }
        }
    }
    f_max = max2;
}
