#include "window.h"
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

int get_part_vector(int n, int m, int k, double *b, double *w,
                    double (*f) (double, double), int &cnx, int &dnx,
                    int &cny, int &dny, double hx, double hy, double x,
                    double y, double f_max, int kol_tyk)
{
    int i=0, j=0;
    k2ij(n, k, i, j, cnx, dnx, cny, dny);
    /*
     * j - num line
     * i - num column
    */
/*inside domain*/
    if(i > 0 && i < n && j > 0 && j < m
            &&((i < cnx || i > cnx+dnx)
               || (i >= cnx && i <= dnx+cnx && (j < cny || j > cny+dny)) ))
    {
        *b =      w[0]  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]  * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[2]  * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y)
                + w[3]  * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]  * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]  * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]  * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]  * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[8]  * f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y)
                + w[9]  * f((i+1)*hx           + x, (j+1)*hy           + y)
                + w[10] * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[11] * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12] * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13] * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14] * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15] * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16] * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17] * f((i-1)*hx           + x, (j  )*hy           + y)
                + w[18] * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }
/*left vertical lines of domain*/
    if( (i == 0 && j > 0 && j < m)
            || (i == cnx+dnx && j > cny && j < cny+dny) )
    {
        *b =      w[0]/2.  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]/2.  * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[2]     * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y)
                + w[3]     * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]/2.  * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[7]/2.  * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[8]     * f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y)
                + w[9]     * f((i+1)*hx           + x, (j+1)*hy           + y)
                + w[10]    * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[11]    * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]    * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]/2. * f((i  )*hx           + x, (j-1)*hy           + y);

        *b = *b / 192.;
        return 0;
    }
/*right vertical lines of domain*/
    if( (i == n && j > 0 && j < m)
            || (i == cnx && j > cny && j < cny+dny) )
    {
        *b =      w[0]/2.  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]/2.  * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[4]/2.  * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]     * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]     * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]/2.  * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[13]/2. * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14]    * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15]    * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16]    * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]    * f((i-1)*hx           + x, (j  )*hy           + y)
                + w[18]    * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }
/*lower horizontal lines of domain*/
    if( (j == 0 && i > 0 && i < n)
            || (j == cny+dny && i > cnx && i < cnx+dnx) )
    {
        *b =      (w[0]/2.  * f((i  )*hx           + x, (j  )*hy           + y))
                + (w[1]     * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y))
                + (w[2]     * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y))
                + (w[3]/2.  * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y))
                + (w[6]/2.  * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y))
                + (w[7]     * f((i  )*hx           + x, (j+1)*hy           + y))
                + (w[8]     * f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y))
                + (w[9]     * f((i+1)*hx           + x, (j+1)*hy           + y))
                + (w[10]    * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y))
                + (w[11]/2. * f((i+1)*hx           + x, (j  )*hy           + y))
                + (w[17]/2. * f((i-1)*hx           + x, (j  )*hy           + y))
                + (w[18]    * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y));

        *b = *b / 192.;
        return 0;
    }
/*upper horizontal lines of domain*/
    if( (j == m && i > 0 && i < n)
            || (j == cny && i > cnx && i < cnx+dnx) )
    {
        *b =      w[0]/2.  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[3]/2.  * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]     * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]     * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]/2.  * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[11]/2. * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]    * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]    * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14]    * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15]    * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16]    * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]/2. * f((i-1)*hx           + x, (j  )*hy           + y);

        *b = *b / 192.;
        return 0;
    }
/*angles of domain*/
    if( i == cnx+dnx && j == cny )
    {
        *b =      w[0]*(5./6.) * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]/2.      * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[2]         * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y)
                + w[3]         * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]         * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]         * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]/2.      * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]/2.      * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[8]         * f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y)
                + w[9]         * f((i+1)*hx           + x, (j+1)*hy           + y)
                + w[10]        * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[11]        * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]        * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]        * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14]        * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15]        * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16]        * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]/2.     * f((i-1)*hx           + x, (j  )*hy           + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == cnx && j == cny+dny )
    {
        *b =      w[0]*(5./6.) * f((i  )*hx           + x, (j  )*hy          + y)
                + w[1]         * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[2]         * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y)
                + w[3]/2.      * f((i  )*hx + (hx/2.) + x, (j  )*hy          + y)
                + w[4]/2.      * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]         * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]         * f((i  )*hx - (hx/2.) + x, (j  )*hy          + y)
                + w[7]         * f((i  )*hx           + x, (j+1)*hy          + y)
                + w[8]         * f((i  )*hx + (hx/2.) + x, (j+1)*hy          + y)
                + w[9]         * f((i+1)*hx           + x, (j+1)*hy          + y)
                + w[10]        * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[11]/2.     * f((i+1)*hx           + x, (j  )*hy          + y)
                + w[13]/2.     * f((i  )*hx           + x, (j-1)*hy          + y)
                + w[14]        * f((i  )*hx - (hx/2.) + x, (j-1)*hy          + y)
                + w[15]        * f((i-1)*hx           + x, (j-1)*hy          + y)
                + w[16]        * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]        * f((i-1)*hx           + x, (j  )*hy          + y)
                + w[18]        * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == cnx && j == cny )
    {
        *b =      w[0]*(2./3.) * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]/2.      * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[3]/2.      * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]         * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]         * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]         * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]/2.      * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[11]/2.     * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]        * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]        * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14]        * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15]        * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16]        * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]        * f((i-1)*hx           + x, (j  )*hy           + y)
                + w[18]        * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == cnx+dnx && j == cny+dny )
    {
        *b =      w[0]*(2./3.)  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]          * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[2]          * f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y)
                + w[3]          * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]/2.       * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[6]/2.       * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]          * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[8]          * f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y)
                + w[9]          * f((i+1)*hx           + x, (j+1)*hy           + y)
                + w[10]         * f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[11]         * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]         * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]/2.      * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[17]/2.      * f((i-1)*hx           + x, (j  )*hy           + y)
                + w[18]         * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == 0 && j == 0)/*-----------------------where add noise-------------------------*/
    {
        *b =      w[0]/3.       *( f((i  )*hx           + x, (j  )*hy           + y) + f_max*kol_tyk*0.1)
                + w[1]/2.       *( f((i  )*hx           + x, (j  )*hy + (hy/2.) + y) + f_max*kol_tyk*0.1)
                + w[2]          *( f((i  )*hx + (hx/2.) + x, (j  )*hy + (hy/2.) + y) + f_max*kol_tyk*0.1)
                + w[3]/2.       *( f((i  )*hx + (hx/2.) + x, (j  )*hy           + y) + f_max*kol_tyk*0.1)
                + w[7]/2.       *( f((i  )*hx           + x, (j+1)*hy           + y) + f_max*kol_tyk*0.1)
                + w[8]          *( f((i  )*hx + (hx/2.) + x, (j+1)*hy           + y) + f_max*kol_tyk*0.1)
                + w[9]          *( f((i+1)*hx           + x, (j+1)*hy           + y) + f_max*kol_tyk*0.1)
                + w[10]         *( f((i+1)*hx           + x, (j  )*hy + (hy/2.) + y) + f_max*kol_tyk*0.1)
                + w[11]/2.      *( f((i+1)*hx           + x, (j  )*hy           + y) + f_max*kol_tyk*0.1);

        *b = *b / 192.;
        return 0;
    }

    if( i == n && j == m )
    {
        *b =      w[0]/3.        * f((i  )*hx           + x, (j  )*hy           + y)
                + w[4]/2.        * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[5]           * f((i  )*hx - (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[6]/2.        * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[13]/2.       * f((i  )*hx           + x, (j-1)*hy           + y)
                + w[14]          * f((i  )*hx - (hx/2.) + x, (j-1)*hy           + y)
                + w[15]          * f((i-1)*hx           + x, (j-1)*hy           + y)
                + w[16]          * f((i-1)*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[17]/2.       * f((i-1)*hx           + x, (j  )*hy           + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == n && j == 0 )
    {
        *b =      w[0]/6.  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[1]/2.  * f((i  )*hx           + x, (j  )*hy + (hy/2.) + y)
                + w[6]/2.  * f((i  )*hx - (hx/2.) + x, (j  )*hy           + y)
                + w[7]/2.  * f((i  )*hx           + x, (j+1)*hy           + y)
                + w[17]/2. * f((i-1)*hx           + x, (j  )*hy           + y)
                + w[18]    * f((i  )*hx - (hx/2.) + x, (j  )*hy + (hy/2.) + y);

        *b = *b / 192.;
        return 0;
    }

    if( i == 0 && j == m )
    {
        *b =      w[0]/6.  * f((i  )*hx           + x, (j  )*hy           + y)
                + w[3]/2.  * f((i  )*hx + (hx/2.) + x, (j  )*hy           + y)
                + w[4]/2.  * f((i  )*hx           + x, (j  )*hy - (hy/2.) + y)
                + w[11]/2. * f((i+1)*hx           + x, (j  )*hy           + y)
                + w[12]    * f((i  )*hx + (hx/2.) + x, (j  )*hy - (hy/2.) + y)
                + w[13]/2. * f((i  )*hx           + x, (j-1)*hy           + y);

        *b = *b / 192.;
        return 0;
    }

    abort();
    return  -1000;
}

/*from this place to do paral*/

void build_vector(int n, int m, int p /*# thread*/,
                  int k /*number of thread*/, int &cnx, int &dnx,
                  int &cny, int &dny, double *b, double (*f) (double, double),
                  double hx, double hy, double f_max, int kol_tyk)
{
    int k1=0, k2=0;
    int N = (n+1)*(m+1) - (dny-1)*(dnx-1);
    k1 = (k * N)/p;
    k2 = ((k+1)*N)/p;

    double *w = new double [19];
    w[0] = 36.;
    for(int i = 1; i <= 6; i++)
    {
        w[i] = 20.;
    }
    for(int i = 7; i <= 17; i+=2)
    {
        w[i] = 2.;
    }
    for(int i = 8; i <= 18; i+=2)
    {
        w[i] = 4.;
    }

//    print_bl(w, 19, 1, 1);
    for( int l = k1; l < k2; l ++)
    {
        get_part_vector(n, m, l/*num by current element*/,
                        (b+l)/*addres elem*/,
                        w/*koef in sum*/, f/*needing function*/,
                        cnx, dnx, cny, dny, hx, hy, 0, 0, f_max, kol_tyk);
    }

    delete [] w;
}
