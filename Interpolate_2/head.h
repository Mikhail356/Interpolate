#ifndef HEAD_H
#define HEAD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
//#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>

#define ERR_INPUT            -1
#define ERR_OPEN             -2
#define ERR_READ             -3
#define ERR_UNKNOWN_ERROR    -4
#define NOT_ENOUGH_MEMORY    -5
#define BAD_PTHREAD          -6
#define MUTEX_INIT           -7
#define MUTEX_DESTROY        -8
#define MAX_IT 100
#define EPS 1e-16
#define PRINT_SIZE 1000000

using namespace std;

#define LOG(...) std::cerr<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"

struct Storage;
struct Point;
struct Triangle;


class Args
{
    public:
        int my_num = 0;//num of the thread
        int error = 0;//for error in process
        int save_int = 0;//for different operations
        int p = 0;
        double save_double = 0;//for different operations
        Storage *store = nullptr;
        pthread_barrier_t* barrier;
//        pthread_mutex_t* mutex;
};

struct Point
{
    double x = 0;
    double y = 0;
};

struct Triangle
{
    Point p [3];
};

struct Storage
{
    int status = 0;
    double len_x = 0;
    double len_y = 0;
    double eps = 0;//precision by recurent method solve system
    double residual = 0;
    int p = 0;//# threads
    int k = 0;//number function for ploting
    double *a = nullptr;//matrix in MSR format
    int *ind = nullptr;//index in MSR format matrix *a
    double *b = nullptr;//vector in MSR format in right part
//    double *x = nullptr;//vector for painting
    double *x_new = nullptr;//vector of solve: Ax == b
    double *u = nullptr;//add memory for counting system
    double *v = nullptr;//add memory for counting system
    double *r = nullptr;//add memory for counting system
    double (*f) (double, double) = nullptr;//current function
    double f_max = 0;
    double f_min = 0;
    int disturb = 0;//disturbance current function
    double *buf = nullptr;//for geting result by all threads
//    int n = 0;//common # segments on axis x == (2*c+d)*nx
//    int m = 0;//common # segments on axis y == (2*c+d)*ny
    int nx = 0;//input parametr, see domain
    int ny = 0;//input parametr, see domain
    int c = 0;//full # of knot = (2*c+d+1)*nx*ny
    int d = 0;//cutout # of knot = (d-1)*nx*(d-1)*ny
//    int N = 0;//common count of knot
    int N_new = 0;
    int f_change = 0;//indicator that f were changed
    int key_push = 0;

    pthread_barrier_t* barrier_all;
    pthread_mutex_t* mutex_all;
};

/*filling_matr.cpp*/

void build_MSR_matrix(int n, int m, double *a, int *ind, int p /*# thread*/,
                      int k /*number of thread*/, int &cnx, int &dnx,
                      int &cny, int &dny, Args *arg);
void reduce_sum (int *sum, Args *arg);
int allocate_MSR_matrix(int n, int m, double **p_a, int **p_ind, int &cnx,
                        int &dnx, int &cny, int &dny );
int get_offdiag_elem(int n, int m, int k, double *a_diag, double *a,
                     int *ind/*not more than 6*/, int &cnx, int &dnx,
                     int &cny, int &dny);
int get_non_zeros(int n, int m, int &cnx, int &dnx, int &cny, int &dny);
int get_num_offdiag(int n, int m, int &cnx, int &dnx, int &cny, int &dny, int k);
void k2ij ( int n, int k, int &i, int &j,
            int &cnx, int &dnx, int &cny, int &dny);
void ij2k ( int n, int /*m*/, int i, int j, int &k,
            int &cnx, int &/*dnx*/, int &cny, int &dny);
void print_msr(double *a, int *ind);


/*filling_vector.cpp*/


int get_part_vector(int n, int m, int k, double *b, double *w,
                    double (*f) (double, double), int &cnx, int &dnx,
                    int &cny, int &dny, double hx, double hy, double x,
                    double y, double f_max, int kol_tyk);
void build_vector(int n, int m, int p /*# thread*/,
                  int k /*number of thread*/, int &cnx, int &dnx,
                  int &cny, int &dny, double *b, double (*f) (double, double),
                  double hx, double hy, double f_max, int kol_tyk);


/*help_kosole.cpp*/

int check_file(char *name, int *c, int *d, double *len_x, double *len_y);
void *solve (void *input);
void find_mm(double (*f)(double, double), double &f_max,
             double hx, double hy, int n, int m);

/*system_solve.cpp*/

int sys_solve (double *a, int *ind, double *b, double *x, double *u, double *v,
               double *r /*additional memory*/, double *buf
               /*one for all, length == p*/, double eps, int p, int my_num,
               pthread_barrier_t * barrier);
int solve_stage(double *a, int *ind, double *b, double *x, double *u, double *v,
                double *r/*additional memory*/,
                double *buf /*one for all, length == p*/, int max_it, double eps,
                int p, int my_num, pthread_barrier_t * barrier);
void msr_matr_mult_vector(double *a, int *ind, double *x, double *v, int p,
                          int my_num);
void linear_combination(double *v, double *b, double len, double koef, int p,
                        int my_num);
void apply_preconditioner(double *a, double *v, double *u, int size, int p,
                          int my_num);
double scalar_product(double *u, double *v, int len, double *buf, int p,
                      int my_num, pthread_barrier_t *barrier);

/*find_error.cpp*/
double find_error(double *x, double (*f) (double, double), double hx, double hy,
                  int cnx, int dnx, int cny, int dny, int my_num, int p,
                  double *buf, /*double x_0, double y_0,*/
                  pthread_barrier_t *barrier, double f_max, int disturb);

double f_0 (double , double);
double f_1 (double x, double);
double f_2 (double , double y);
double f_3 (double x, double y);
double f_4 (double x, double y);
double f_5 (double x, double y);
double f_6 (double x, double y);
double f_7 (double x, double y);

int compare1(const void* a, const void* b);
int compare2(const void* a, const void* b);
int compare3(const void* a, const void* b, void *eye);


template <typename T>
void print_bl (T *a, int column_len, int row_len, int real_len)
{
    for(int i = 0; i < min(column_len, PRINT_SIZE); i++)
    {
        for(int j = 0; j < min(row_len, PRINT_SIZE); j++)
        {
            cout<<a[i*real_len+j]<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";
}

#endif // HEAD_H

