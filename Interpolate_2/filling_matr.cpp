#include "head.h"

/*********************
 * matr indexing from down to up
 * and from left to right
 * ^ >
 * j - num line
 * i - num column
 *********************/


/****************************
 * mapping ij to k and back
 * k - number of knot
 * j,i coordinate of k knot
 * n - # segments on axis x == 2*cnx+dnx
 * m - # segments on axis y == 2*cny+dny
 ****************************/
void ij2k ( int n, int /*m*/, int i, int j, int &k,
            int &cnx, int &/*dnx*/, int &cny, int &dny)
{
    if(j <= cny)
    {
        k = j*(n+1)+i;
    }
    else if(j < cny+dny)
    {
        if(i <= cnx)
        {
            k = (n+1)*(cny+1) + (j-cny-1)*2*(cnx+1) + i;
        }
        else
        {
            k = (n+1)*(cny+1) + (j-cny-1)*2*(cnx+1) + i - (cnx-1);
        }
    }
    else
    {
        k = (n+1)*(cny+1) + 2*(dny-1)*(cnx+1) + (j-cny-dny)*(n+1) + i;
    }

//    printf("%d %d\nk == %d\n\n", j, i, k);
}

void k2ij ( int n, int k, int &i, int &j,
            int &cnx, int &dnx, int &cny, int &dny)
{
    if( k < (cny+1)*(n+1) )
    {
        j = k / (n+1); /*number of row down 2 up*/
        i = k % (n+1); /*number of column left 2 right*/
    }
    else if( k < ((cny+1)*(n+1)) + (2*(dny-1)*(cnx+1)) )
    {
        j = ((k - (cny+1)*(n+1)) / (2*(cnx+1))) + cny+1;
        i = (k - (cny+1)*(n+1)) % (2*(cnx+1));
        if(i > cnx)
        {
            i += (dnx-1);
        }
    }
    else
    {
        j = ((k - (n+1)*(cny+1) - 2*(dny-1)*(cnx+1)) / (n+1)) + cny + dny;
        i = (k - (n+1)*(cny+1) - 2*(dny-1)*(cnx+1)) % (n+1);
    }
//    printf("k == %d\n %d %d\n\n", k, j, i);
}

/******************************************
 * # knot of basic function that lies in domain
 * exect knot on diagonal
 * n - # segments on axis x
 * m - # segments on axis y
 * return # knot without central knot of
 *   ___
 *  /| /|
 * /_|/_|
 * | /| /
 * |/_|/
 * in domain
 *  ____________________________
 * |                            |
 * |     __________________     |
 * |    |                  |    |
 * |    |                  |    |
 * |    |__________________|    |
 * |                            |
 * |____________________________|
 * ****************************************/
int get_num_offdiag(int n, int m, int &cnx, int &dnx, int &cny, int &dny, int k)
{
    int i=0, j=0;
    k2ij(n, k, i, j, cnx, dnx, cny, dny);

    /*inside domain*/
    if(i > 0 && i < n && j > 0 && j < m
            &&((i < cnx || i > cnx+dnx)
               || (i >= cnx && i <= dnx+cnx && (j < cny || j > cny+dny)) ))
    {
        return 6;
    }

    /*verticle lines*/
    if( ((i == 0 || i == n) && j > 0 && j < m)
            ||((i == cnx || i == cnx+dnx) && j > cny && j < cny+dny ) )
    {
        return 4;
    }

    /*horisontal lines*/
    if( ((j == 0 || j == m) && i > 0 && i < n)
            ||((j == cny || j == cny+dny) && i > cnx && i < cnx+dnx) )
    {
        return 4;
    }

    /*upper left and lower right inside angles*/
    if( (i == cnx && j == cny+dny) || (i == cnx+dnx && j == cny) )
    {
        return 6;
    }

    /*lower left and upper right inside angles*/
    if( (i == cnx && j == cny) || (i == cnx+dnx && j == cny+dny))
    {
        return  5;
    }

    /*lower left and upper right outside angles*/
    if( (i == 0 && j == 0) || (i == n && j == m) )
    {
        return 3;
    }

    /*upper left and lower right outside angles*/
    if( (i == 0 && j == m) || (i == n && j == 0) )
    {
        return 2;
    }

    //never go here
    abort();
    return -1000;
}

/***************************
 * # non 0 element in msr matrix
 * n - # segments on axis x
 * m - # segments on axis y
 * *************************/
int get_non_zeros(int n, int m, int &cnx, int &dnx, int &cny, int &dny)
{
    int K = (n+1)*(m+1) - (dnx-1)*(dny-1); //common # of knot
    int nz = 0; //non zero

    for(int k = 0; k < K; k ++)
    {
        nz += get_num_offdiag(n, m, cnx, dnx, cny, dny, k);
    }
    return nz;
    //and then the total size of the MSR matrix = K+1+nz
}

int get_offdiag_elem(int n, int m, int k, double *a_diag, double *a,
                     int *ind/*not more than 6*/, int &cnx, int &dnx,
                     int &cny, int &dny) /*of msr matrix format*/
{
    int i=0, j=0;
    k2ij(n, k, i, j, cnx, dnx, cny, dny);
/*inside domain*/
    if(i > 0 && i < n && j > 0 && j < m
            &&((i < cnx || i > cnx+dnx)
               || (i >= cnx && i <= dnx+cnx && (j < cny || j > cny+dny)) ))
    {
        *a_diag = 0.5;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i,   j-1, ind[3], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[4], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[5], cnx, dnx, cny, dny);
        for( int l = 0; l < 6; l++)
        {
            a[l]=1./12;
        }
        return 6;
    }
/*left vertical lines of domain*/
    if( (i == 0 && j > 0 && j < m)
            || (i == cnx+dnx && j > cny && j < cny+dny) )
    {
        *a_diag = 0.25;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i,   j-1, ind[3], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./24;
        return 4;
    }
/*right vertical lines of domain*/
    if( (i == n && j > 0 && j < m)
            || (i == cnx && j > cny && j < cny+dny) )
    {
        *a_diag = 0.25;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i,   j-1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[3], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./24;
        a[2] = 1./12;
        a[3] = 1./12;
        return 4;
    }
/*lower horizontal lines of domain*/
    if( (j == 0 && i > 0 && i < n)
            || (j == cny+dny && i > cnx && i < cnx+dnx) )
    {
        *a_diag = 0.25;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[3], cnx, dnx, cny, dny);
        a[0] = 1./12;
        a[1] = 1./12;
        a[2] = 1./24;
        a[3] = 1./24;
        return 4;
    }
/*upper horizontal lines of domain*/
    if( (j == m && i > 0 && i < n)
            || (j == cny && i > cnx && i < cnx+dnx) )
    {
        *a_diag = 0.25;
        ij2k(n, m, i+1, j,   ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i  , j-1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[3], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./24;
        return 4;
    }
/*angles of domain*/
    if( i == cnx+dnx && j == cny )
    {
        *a_diag = 5./12;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i  , j-1, ind[3], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[4], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j  , ind[5], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./12;
        a[4] = 1./12;
        a[5] = 1./24;
        return 6;
    }

    if( i == cnx && j == cny+dny )
    {
        *a_diag = 5./12;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i  , j-1, ind[3], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[4], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j  , ind[5], cnx, dnx, cny, dny);
        a[0] = 1./12;
        a[1] = 1./12;
        a[2] = 1./24;
        a[3] = 1./24;
        a[4] = 1./12;
        a[5] = 1./12;
        return 6;
    }

    if( i == cnx && j == cny )
    {
        *a_diag = 1./3;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i,   j-1, ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[3], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[4], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./24;
        a[2] = 1./12;
        a[3] = 1./12;
        a[4] = 1./12;
        return 5;
    }

    if( i == cnx+dnx && j == cny+dny )
    {
        *a_diag = 1./3;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        ij2k(n, m, i  , j-1, ind[3], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[4], cnx, dnx, cny, dny);
        a[0] = 1./12;
        a[1] = 1./12;
        a[2] = 1./12;
        a[3] = 1./24;
        a[4] = 1./24;
        return 5;
    }

    if( i == 0 && j == 0)
    {
        *a_diag = 1./6;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j+1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i+1, j,   ind[2], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./24;
        return 3;
    }

    if( i == n && j == m )
    {
        *a_diag = 1./6;
        ij2k(n, m, i,   j-1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j-1, ind[1], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[2], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./12;
        a[2] = 1./24;
        return 3;
    }

    if( i == n && j == 0 )
    {
        *a_diag = 1./12;
        ij2k(n, m, i,   j+1, ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i-1, j,   ind[1], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./24;
        return 2;
    }

    if( i == 0 && j == m )
    {
        *a_diag = 1./12;
        ij2k(n, m, i+1, j,   ind[0], cnx, dnx, cny, dny);
        ij2k(n, m, i,   j-1, ind[1], cnx, dnx, cny, dny);
        a[0] = 1./24;
        a[1] = 1./24;
        return 2;
    }

    abort();
    return  -1000;
}

int allocate_MSR_matrix( int n, int m, double **p_a, int **p_ind, int &cnx,
                         int &dnx, int &cny, int &dny )
{
    double *a;
    int *ind;
    int N = (n+1)*(m+1) - (dny-1)*(dnx-1);
    int nz = get_non_zeros(n, m, cnx, dnx, cny, dny);
    int len = N+1+nz;
    int l=0, sum=0;

    a = new(std::nothrow) double [len];
    if(a == nullptr)
    {
        return -1;
    }
    ind = new(std::nothrow) int [len];
    if(ind == nullptr)
    {
        delete [] a;
        return  -2;
    }
//    memset(a, 0, len*sizeof (double));
//    memset(ind, 0, len*sizeof (int));

    ind[0] = N+1;
    for(int k = 1; k <= N; k ++)
    {
        l = get_num_offdiag(n, m, cnx, dnx, cny, dny, k-1);
        ind[k] = ind[k-1] + l;
        sum+=l;
    }

    if (nz != sum)
    {
        printf("crashed in allocate_MSR_matrix 1\n");
        abort();
    }
    if (ind[N] != len )
    {
        printf("crashed in allocate_MSR_matrix 2\n");
        abort();
    }

    *p_a = a;
    *p_ind = ind;

//    delete [] a;
//    delete [] ind;
    return 0;
}

/*from this place to do paral*/

void build_MSR_matrix( int n, int m, double *a, int *ind, int p /*# thread*/,
                       int k /*number of thread*/, int &cnx, int &dnx,
                       int &cny, int &dny, Args *arg)
{
    int k1=0, k2=0, sum=0, s=0;
    int N = (n+1)*(m+1) - (dny-1)*(dnx-1);
    k1 = (k * N)/p;
    k2 = ((k+1)*N)/p;


//    int i = 0, j = 0, kkk = 0; // delete it
    for( int l = k1; l < k2; l ++)
    {
        s = get_offdiag_elem(n, m, l, (a+l)/*addres diagonal elem*/,
                             (a+ind[l])/*begin row*/,
                             (ind + ind[l])/*index of not diagonal elem*/,
                             cnx, dnx, cny, dny);

//        k2ij(n, l, i, j, cnx, dnx,cny, dny);
//        ij2k(n, m, i, j, kkk, cnx, dnx, cny, dny);
//        printf("num == %d   %d %d %d\nind  ", l, i, j, kkk);

//        for(int iii = 0; iii < s; iii++)
//        {
//            printf("%d ", ind[ind[0]+sum+iii]);
//        }
//        printf("\nk2ij ");
//        for(int iii = 0; iii < s; iii++)
//        {
//            k2ij(n, ind[ind[0]+sum+iii], i, j, cnx, dnx,cny, dny);
//            printf("%d  %d|  ",i ,j );
//        }
//        printf("\nij2k ");
//        for(int iii = 0; iii < s; iii++)
//        {
//            k2ij(n, ind[ind[0]+sum+iii], i, j, cnx, dnx,cny, dny);
//            ij2k(n, m, i, j, kkk, cnx, dnx, cny, dny);
//            printf("%d ",kkk );
//        }
//        printf("\n\n");

        sum+=s;
    }

    reduce_sum (&sum, arg);/*sum from all threads should be = nz*/

    if(N+1+sum != ind[N])
    {
        printf("crash in build_MSR_matrix\n");
        abort();
    }
}

void reduce_sum (int *sum, Args *arg)
{
    arg->save_int = *sum;
    pthread_barrier_wait( arg->barrier );

//    pthread_mutex_lock( arg->mutex );
//    cout<<arg->my_num<<" "<<arg->save_int<<"\n";
//    pthread_mutex_unlock( arg->mutex );
//    pthread_barrier_wait( arg->barrier );
    if(arg->my_num == 0)
    {
        int tmp = 0;
        for(int l = 0; l < arg->p; l++)
        {
            tmp += arg[l].save_int;
        }

        for(int l = 0; l < arg->p; l++)
        {
            arg[l].save_int = tmp;
        }
    }
    pthread_barrier_wait( arg->barrier );
    *sum = arg->save_int;
}

void print_msr(double *a, int *ind)
{
    int size = ind[0]-1;
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum = 0;
        sum+=a[i];
        for(int j = ind[i]; j < ind[i+1]; j++)
        {
            sum+=a[j];
        }
        printf("\t%lf  ", sum);



        printf("%d  %lf ", i, a[i]);
        for(int j = ind[i]; j < ind[i+1]; j++)
        {
            printf("%lf ", a[j]);
        }
        printf("\n");
    }
}
