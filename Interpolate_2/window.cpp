#include "window.h"
#include "head.h"

double f_0 (double , double)    { return  1; }
double f_1 (double x, double)   { return  x; }
double f_2 (double , double y)  { return  y; }
double f_3 (double x, double y) { return  (x + y); }
double f_4 (double x, double y) { return  sqrt (x * x + y * y); }
double f_5 (double x, double y) { return  (x * x + y * y); }
double f_6 (double x, double y) { return  exp (x * x - y * y); }
double f_7 (double x, double y) { return  (1./(25 * (x * x + y * y) + 1)); }
Window::Window (QWidget *parent, int &argc, char *argv[])
    : QWidget (parent)
{
    int err = 0;
    if (argc != 7)
    {
//      printf ("Usage: %s [name] nx ny k eps p\n", argv[0]);
        store->nx = 1;
        store->ny = 1;
        store->k = 0;
        store->eps = 1e-12;
        store->p = 4;
        name = "data"/*"/home/misha/Interoplate_2_new/data"*/;
//      return ERR_INPUT;
    }
    else
    {
        if ((store->nx  = atoi (argv [2])) <= 0 ||
            (store->ny  = atoi (argv [3])) <= 0 ||
            (store->k   = atoi (argv [4])) < 0 ||
            (store->eps = atof (argv [5])) <= 1e-16 ||
            (store->p   = atoi (argv [6])) <= 0)
        {
            printf ("Usage: %s [name] nx ny k eps p", argv[0]);
            err = ERR_INPUT;
//            return ERR_INPUT;
        }
        name = argv[1];
    }

//    fstream file;
//    file.open(name.toStdString());
//    if(!file.is_open())
//    {
//        printf("Error open file\n");
////        return ERR_OPEN;
//    }
//    file>>store->c>>store->d>>store->len_x>>store->len_y;
    ifstream file;
    vector<string> lines;
    string line;
    int vector_len = 0;

    file.open(name.toStdString());
    if(!file.is_open())
    {
        printf("Error open file\n");
        store->status = -1;
        err = ERR_OPEN;
    }
    while (getline (file, line))
    {
        line += " ";
        if (line[0] != '#')
        {
            lines.push_back (line);
            vector_len++;
        }
    }
    file.close();

    string data_file;
    for(int i = 0; i < vector_len; i++)
    {
        data_file += lines[i];
    }

    istringstream input (data_file);
    input >> store->c >> store->d >> store->len_x >> store->len_y;


//    file.close();

    if (store->len_x < 1e-16 || store->len_y < 1e-16)
    {
        printf ("Error in input file\n");
        err = ERR_READ;
//        return ERR_READ;
    }

    store->buf = new double [store->p];
    store->N_new = ((2*store->c+store->d)*store->nx + 1)
            *((2*store->c+store->d)*store->ny + 1)
            - ((store->d*store->nx-1)*(store->d*store->ny-1));

//    pthread_barrier_t barrier;
//    pthread_barrier_t barrier_all;
    pthread_barrier_init(&barrier_all, nullptr, store->p+1);
    store->barrier_all = &barrier_all;
    pthread_barrier_init(&barrier, nullptr, store->p);
    pthread_mutex_init(&mutex_all, NULL);
    store->mutex_all = &mutex_all;

    switch(store->k)
    {
    case 0:
        store->f = f_0;
        name = "f = 1";
        break;
    case 1:
        store->f = f_1;
        name = "f = x";
        break;
    case 2:
        store->f = f_2;
        name = "f = y";
        break;
    case 3:
        store->f = f_3;
        name = "f = x + y";
        break;
    case 4:
        store->f = f_4;
        name = "f = sqrt(x^2 + y^2)";
        break;
    case 5:
        store->f = f_5;
        name = "f = x^2 + y^2";
        break;
    case 6:
        store->f = f_6;
        name = "f = exp(x^2-y^2)";
        break;
    case 7:
        store->f = f_7;
        name = "f = 1/(25(x^2 + y^2) + 1)";
        break;
    }


    arg = new Args [store->p];
    thread_id = new pthread_t [store->p];
    f_max = 1;

    if(MAX_C>store->c)
    {
        MAX_D = (int)(store->d * ((double)(MAX_C/store->c)));
    }
    else
    {
        MAX_D = (int)(store->d);
    }
    if(MAX_NX>store->nx)
    {
        MAX_NY = (int)(store->ny * ((double)(MAX_NX/store->nx)));
    }
    else
    {
        MAX_NY = (int)(store->ny);
    }
//    MAX_D = (int)(store->d * ((double)(MAX_C/store->c)));
//    MAX_NY = (int)(store->ny * ((double)(MAX_NX/store->nx)));
    MAX_N = (2*MAX_C+MAX_D)*MAX_NX;
    MAX_M = (2*MAX_C+MAX_D)*MAX_NY;
    N = 1;
    x = new double [N];
    tree = new Triangle [N];
//    memset(x, 0, N*sizeof(double));
//    memset(tree, 0, N*sizeof(Triangle));
    int p = store->p;
    for(int i = 0; i < p; i++)
    {
        arg[i].store = store;
        arg[i].my_num = i;
        arg[i].p = p;
        arg[i].barrier = &barrier;
        pthread_create(&thread_id[i], nullptr, solve, arg+i);
    }
    pthread_barrier_wait(store->barrier_all);
    t.setInterval(100);
    t.setSingleShot(true);
    connect(&t, SIGNAL(timeout()), this, SLOT(update_screen()));
    t.start();
    if(err != 0)
    {
        close();
    }
    update();
}

void
Window :: update_screen()
{
    t.start();
    pthread_mutex_lock(&mutex_all);
    int i = store->f_change;
    pthread_mutex_unlock(&mutex_all);
    if(i == 1)
    {
        update();
    }
}

Window::~Window ()
{
//    pthread_mutex_lock(&mutex_all);
//    store->status = -1;
    int kol_out = store->p;
//    pthread_mutex_unlock(&mutex_all);
    for(int i = 0; i < kol_out; i++)
    {
        pthread_join(thread_id[i],0);
    }
    delete [] tree;
    delete [] thread_id;
    delete [] arg;
    delete [] x;
//    pthread_mutex_lock(&mutex_all);
    delete [] store->buf;
//    delete [] store->x_new;
//    store->x_new = nullptr;
//    pthread_mutex_unlock(&mutex_all);

    pthread_mutex_destroy(&mutex_all);
    pthread_barrier_destroy(&barrier);
    pthread_barrier_destroy(&barrier_all);
}
QSize Window::minimumSizeHint () const
{
    return QSize (100, 100);
}

QSize Window::sizeHint () const
{
    return QSize (1000, 1000);
}

void
Window::keyPressEvent (QKeyEvent * event)
{

  switch (event -> key ())
    {// change function
    case Qt::Key_0:
      pthread_mutex_lock(&mutex_all);
      store->k = (store->k + 1) % func_count;
      store->key_push = 1;
      switch(store->k)
      {
      case 0:
          store->f = f_0;
          name = "f = 1";
          break;
      case 1:
          store->f = f_1;
          name = "f = x";
          break;
      case 2:
          store->f = f_2;
          name = "f = y";
          break;
      case 3:
          store->f = f_3;
          name = "f = x + y";
          break;
      case 4:
          store->f = f_4;
          name = "f = sqrt(x^2 + y^2)";
          break;
      case 5:
          store->f = f_5;
          name = "f = x^2 + y^2";
          break;
      case 6:
          store->f = f_6;
          name = "f = exp(x^2-y^2)";
          break;
      case 7:
          store->f = f_7;
          name = "f = 1/(25(x^2 + y^2) + 1)";
          break;
      }
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
//      setEnabled(true);
      pthread_barrier_wait (store->barrier_all);
      update ();
      break;

      // change render on screen
    case Qt::Key_1:
//      setEnabled(true);
      render_id = (render_id + 1) % render_count;
//      setEnabled(false);
      update ();
      break;

      // division len by 2
    case Qt::Key_2:
      pthread_mutex_lock(&mutex_all);
      store->key_push = 1;
      store->len_x *= 0.5;
      store->len_y *= 0.5;
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait (store->barrier_all);
      update ();
      break;

      // multiple len on 2
    case Qt::Key_3:
      pthread_mutex_lock(&mutex_all);
      store->key_push = 1;
      store->len_x *= 2.;
      store->len_y *= 2.;
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait (store->barrier_all);
      update ();
      break;

      // multiple n on 2
    case Qt::Key_4:
      pthread_mutex_lock(&mutex_all);
      store->nx *= 2;
      store->ny *= 2;
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait(store->barrier_all);
      update ();
      break;

      // division n by 2
    case Qt::Key_5:
      pthread_mutex_lock(&mutex_all);
      if(store->nx > 1 && store->ny > 1)
      {
          store->nx /= 2;
          store->ny /= 2;
      }
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait(store->barrier_all);
      update ();
      break;

      // increase disturbance
    case Qt::Key_6:
      pthread_mutex_lock(&mutex_all);
      store->key_push = 1;
      store->disturb ++;
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait(store->barrier_all);
      update ();
      break;

      // decrease disturbance
    case Qt::Key_7:
      pthread_mutex_lock(&mutex_all);
      store->key_push = 1;
      store->disturb --;
      pthread_mutex_unlock(&mutex_all);
      setEnabled(false);
      pthread_barrier_wait(store->barrier_all);
      update ();
      break;

      // rotate
    case Qt::Key_8:
    case Qt::Key_Left:
      angle = angle - (PI/180.);
      update ();
      break;

      // rotate counterclockewise
    case Qt::Key_9:
    case Qt::Key_Right:
      angle = angle + (PI/180.);
      update ();
      break;

      // close application
    case Qt::Key_Escape:
      close ();

    default:
      break;

    }

}

void
Window::closeEvent(QCloseEvent */*event*/)
{
    pthread_mutex_lock(&mutex_all);
    store->status = -1;
    pthread_mutex_unlock(&mutex_all);
    pthread_barrier_wait(store->barrier_all);
//    int kol_thread = store->p;
//    for(int i = 0; i < kol_thread; i++)
//    {
//        pthread_join
//    }
//    pthread_mutex_lock(&mutex_all);
//    if(store->x_new != nullptr)
//    {
//        delete [] store->x_new;
//        store->x_new = nullptr;
//    }
//    pthread_mutex_unlock(&mutex_all);
}

void
Window::paintEvent (QPaintEvent */*event*/)
{
    QPainter painter (this);
    pthread_mutex_lock(&mutex_all);
    if(store->status == 0)
    {
        cnx = store->c * store->nx;
        cny = store->c * store->ny;
        dnx = store->d * store->nx;
        dny = store->d * store->ny;
        n = (2*store->c+store->d)*store->nx;
        m = (2*store->c+store->d)*store->ny;
        hx = store->len_x / n;
        hy = store->len_y / m;
        x_len = store->len_x;
        y_len = store->len_y;
        store->f_max = f_max;
        f = store->f;
        residual = store->residual;
        disturbance = store->disturb;

        save_int = 2*((min(n, MAX_N)*min(m, MAX_M))
                  - (min(MAX_D*MAX_NX, dnx)*min(MAX_D*MAX_NY, dny)));
//        printf("save_int == %d\n", save_int);
        if(store->f_change == 1)
        {
            setEnabled(true);
            delete [] x;
            delete [] tree;
            x = new double [store->N_new];
            N = store->N_new;
            tree = new Triangle [save_int];
            if(store->x_new != nullptr)
            {
                memcpy(x, store->x_new, N*sizeof(double));
            }
            delete [] store->x_new;
            store->x_new = nullptr;
            store->f_change = 0;
            printf("max{F_min; F_max} = %e\n", f_max);
        }
    }
    pthread_mutex_unlock(&mutex_all);

    pen = QPen (Qt::black, 1, Qt::SolidLine);
    brash.setColor (Qt::white);
    brash.setStyle (Qt::SolidPattern);
    painter.setPen (pen);
    draw_background (painter);
    brash.setColor (Qt::green);

    // angle_in_grad = m_nagle * (2 * PI) / 360

    h_sin = sin(angle);
    h_cos = cos(angle);

    switch (render_id)
    {
    case 0:
        find_min_max (&Window::func_value);
        break;
    case 1:
        find_min_max (&Window::app_value);
        break;
    case 2:
        find_min_max (&Window::err_value);
        break;
    }

    switch (render_id)
    {
    case 0:
        draw_func (painter, &Window::func_value);
        break;
    case 1:
        draw_func (painter, &Window::app_value);
        break;
    case 2:
        draw_func (painter, &Window::err_value);
        break;
    }

    char buf[256];

    painter.setPen (Qt::darkBlue);

    painter.drawText (10, 30, name);

    sprintf (buf, "n = %d; m = %d", n, m);
    painter.drawText (10, 50, buf);

    sprintf (buf, "Residual = %e", residual);
    painter.drawText (10, 70, buf);

    sprintf (buf, "max{F_min; F_max} = %e", f_max);
    painter.drawText (10, 90, buf);

    sprintf (buf, "Len_x = %lf; Len_y = %lf", x_len, y_len);
    painter.drawText (10, 110, buf);

    sprintf (buf, "p = %d", disturbance);
    painter.drawText (10, 130, buf);

    painter.setPen (Qt::darkGray);

    switch (render_id)
    {
    case 0:
        painter.drawText (10, 150, "Function");
        break;
    case 1:
        painter.drawText (10, 150, "Approximation");
        break;
    case 2:
        painter.drawText (10, 150, "Residual");
        break;
    }
//    sprintf (buf, "angle = %lf", angle*180./PI);
//    painter.drawText (10, 170, buf);
}

double
Window :: func_value(double x_, double y_)
{
    if(fabs(x_)<EPS&&fabs(y_)<EPS)
    {
        return f(x_,y_) + disturbance*f_max*0.1;
    }
    return f(x_,y_);
}

double
Window :: app_value(double x_, double y_)
{
    int i = x_/hx;
    int j = y_/hy;
    int k = 0;
    ij2k(n, m, i, j, k, cnx, dnx, cny, dny);
    return x[k];
}

double
Window::err_value (double x_, double y_)
{
  return fabs (func_value (x_, y_) - app_value (x_, y_));
}

void
Window::find_min_max (double (Window::* f) (double, double))
{
    double dx = (double) x_len/min(n,MAX_N);
    double dy = (double) y_len/min(m,MAX_M);

    f_max = 0;
    double x = 0, y = 0, val = 0;

    for (int i = 0; i <= min(n,MAX_N); i++)
    {
        for (int j = 0; j <= min(m,MAX_M); j++)
        {
            x = i * dx;
            y = j * dy;
            val = (this ->* f) (x, y);
            f_max = std::max (f_max, fabs (val));
        }
    }
    f_max = max (f_max, 1e-15);

//    QPoint point;
    double x0 = 0, y0 = 0, x_ = 0,
            y_ = 0, len = 0, tmp1 = 0, tmp2 = 0;
    x_min = x_max = y_min = y_max = 0;
//    x0 = /*x_len*/ /* * 0.5*/0;
//    y0 = /*y_len*//* * 0.5*/0;
    x0 = y0 = 0/*(x_len+y_len) * 0.5*/;
    len = sqrt((x_len*x_len)+(y_len*y_len));
//    len = x_len;

    for (int i = -180; i < 180; i++)
    {
        x_ = x0 + len * sin ((double)i * PI / 180.0);
        y_ = y0 + len * cos ((double)i * PI / 180.0);
//        tmp1 = x_*h_cos - y_*h_sin;
//        tmp2 = x_*h_sin + y_*h_cos;
        tmp1 = x_;
        tmp2 = y_;

        tmp2 += f_max*/*(this->*f)(tmp1, tmp2)*/
                cos(/*angle*/DEF_ANGLE * PI / 180.0);
//        projection(tmp1, tmp2, 0., point);
//        tmp1 = point.x();
        x_min = min(x_min, tmp1);
        x_max = max(x_max, tmp1);

//        projection(tmp1, tmp2, -1., point);
//        tmp2 = point.y();
        y_min = min(y_min, tmp2);

//        projection(tmp1, tmp2, 1., point);
//        tmp2 = point.y();
        y_max = max(y_max, tmp2);
    }
}

void
Window::sort_triangles (Triangle *tr, QPoint eye)
{
    int len = 2 * ((min(MAX_N, n)*min(MAX_M,m)
                    - (min(dnx, MAX_D*MAX_NX)
                       *min(dny,MAX_D*MAX_NY))));
//    for(int i = 0; i < len; i++)
//    {
//        tr[i].p[0].x += eye.x();
//        tr[i].p[1].x += eye.x();
//        tr[i].p[2].x += eye.x();
//        tr[i].p[0].y += eye.y();
//        tr[i].p[1].y += eye.y();
//        tr[i].p[2].y += eye.y();
//    }
    qsort_r (tr, len, sizeof(Triangle), compare3, &eye);
//    for(int i = 0; i < len; i++)
//    {
//        tr[i].p[0].x -= eye.x();
//        tr[i].p[1].x -= eye.x();
//        tr[i].p[2].x -= eye.x();
//        tr[i].p[0].y -= eye.y();
//        tr[i].p[1].y -= eye.y();
//        tr[i].p[2].y -= eye.y();
//    }
}

void
Window::sort_triangles2 (Triangle *tr, Point &eye)
{
    int len = 2 * ((min(MAX_N, n)*min(MAX_M,m)
                    - (min(dnx, MAX_D*MAX_NX)
                       *min(dny,MAX_D*MAX_NY))));
    qsort_r (tr, len, sizeof(Triangle), compare3, &eye);
}
int compare3(const void* a, const void* b, void *e)
{
    Point eye = *(Point*) e;
    Triangle arg1 = *(const Triangle*)a;
    Triangle arg2 = *(const Triangle*)b;
//    double center1 = (arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3
//                    +(arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3;
//    double center2 = (arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3
//                    +(arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3;
//    double center1 = fabs((arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3 - eye.x)
//            + fabs((arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3 - eye.y);
//    double center2 = fabs((arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3 - eye.x)
//            + fabs((arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3 - eye.y);
    double center1 = ((arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3 - eye.x)
            *((arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3 - eye.x)
            +((arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3 - eye.y)
            *((arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3 - eye.y);
    double center2 = ((arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3 - eye.x)
            *((arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3 - eye.x)
            +((arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3 - eye.y)
            *((arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3 - eye.y);
//    center1 = fabs(center1 - eye.x - eye.y);
//    center2 = fabs(center2 - eye.x - eye.y);

    if ( center1 < center2 ) return 1;
    if ( center1 > center2 ) return -1;
    return 0;
}
int compare1(const void* a, const void* b)
{
    Triangle arg1 = *(const Triangle*)a;
    Triangle arg2 = *(const Triangle*)b;
    double center1 = (arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3
                    +(arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3;
    double center2 = (arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3
                    +(arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3;

    if ( center1 < center2 ) return 1;
    if ( center1 > center2 ) return -1;
    return 0;
}

int compare2(const void* a, const void* b)
{
    Triangle arg1 = *(const Triangle*)a;
    Triangle arg2 = *(const Triangle*)b;
    double center1 = (arg1.p[0].x + arg1.p[1].x + arg1.p[2].x)/3
                    +(arg1.p[0].y + arg1.p[1].y + arg1.p[2].y)/3;
    double center2 = (arg2.p[0].x + arg2.p[1].x + arg2.p[2].x)/3
                    +(arg2.p[0].y + arg2.p[1].y + arg2.p[2].y)/3;

    if ( center1 > center2 ) return 1;
    if ( center1 < center2 ) return -1;
    return 0;
}

void
Window::rotate (Triangle *tr, int len)//in R^3
{
    double x0, y0, x1, y1, x2, y2;
    for(int i = 0; i < len; i++)
    {
        x0 = tr[i].p[0].x * h_cos - tr[i].p[0].y * h_sin;
        y0 = tr[i].p[0].x * h_sin + tr[i].p[0].y * h_cos;
        x1 = tr[i].p[1].x * h_cos - tr[i].p[1].y * h_sin;
        y1 = tr[i].p[1].x * h_sin + tr[i].p[1].y * h_cos;
        x2 = tr[i].p[2].x * h_cos - tr[i].p[2].y * h_sin;
        y2 = tr[i].p[2].x * h_sin + tr[i].p[2].y * h_cos;
        tr[i].p[0].x = x0;
        tr[i].p[0].y = y0;
        tr[i].p[0].x = x1;
        tr[i].p[0].y = y1;
        tr[i].p[0].x = x2;
        tr[i].p[0].y = y2;
    }
}

void
Window::corotate (Triangle *tr, int len)//in R^3
{
    double x0, y0, x1, y1, x2, y2;
    double angle_ = 180;
    double h_cos_ = cos(angle_), h_sin_ = sin(angle_);
    for(int i = 0; i < len; i++)
    {
        x0 = tr[i].p[0].x * h_cos_ - tr[i].p[0].y * h_sin_;
        y0 = tr[i].p[0].x * h_sin_ + tr[i].p[0].y * h_cos_;
        x1 = tr[i].p[1].x * h_cos_ - tr[i].p[1].y * h_sin_;
        y1 = tr[i].p[1].x * h_sin_ + tr[i].p[1].y * h_cos_;
        x2 = tr[i].p[2].x * h_cos_ - tr[i].p[2].y * h_sin_;
        y2 = tr[i].p[2].x * h_sin_ + tr[i].p[2].y * h_cos_;
        tr[i].p[0].x = x0;
        tr[i].p[0].y = y0;
        tr[i].p[0].x = x1;
        tr[i].p[0].y = y1;
        tr[i].p[0].x = x2;
        tr[i].p[0].y = y2;
    }
}

void
Window :: projection (double x, double y, double z, Point &p/*QPoint &p*/)
{
//    double X, Y;
//    X = Y = 0.0;

//    static double a = cos (DEF_ANGLE * PI / 180.0);

    // (1, 0, 0) -> (-c, -s)
//    X += -v_cos * x;
//    Y += -v_sin * x;
    // (0, 1, 0) -> (+c, -s)
//    X += v_cos * y;
//    Y += -v_sin * y;
    // (0, 0, 1) -> ( 0, +c)
//    Y += v_cos * z / a;

    int w = width () - 1;
    int h = height () - 1;

//    printf("(x, y) == (%lf, %lf)\n",x, y);
    int x_ = w * ((x - x_min) / (x_max - x_min));

    int y_ = h
            * ((y_max - y - (z*cos(/*angle*/DEF_ANGLE * PI / 180.0)))
                  /(y_max - y_min+f_max));
//    p.setX(x_);
//    p.setY(y_);
    p.x = x_;
    p.y = y_;
}

void
Window :: projection (double x, double y, double z, QPoint &p)
{
//    double X, Y;
//    X = Y = 0.0;

//    static double a = cos (DEF_ANGLE * PI / 180.0);

    // (1, 0, 0) -> (-c, -s)
//    X += -v_cos * x;
//    Y += -v_sin * x;
    // (0, 1, 0) -> (+c, -s)
//    X += v_cos * y;
//    Y += -v_sin * y;
    // (0, 0, 1) -> ( 0, +c)
//    Y += v_cos * z / a;

    int w = width () - 1;
    int h = height () - 1;

    int x_ = w * ((x - x_min) / (x_max - x_min));

    int y_ = h
            * ((y_max - y - (z*cos(DEF_ANGLE * PI / 180.0)))
                  /(y_max - y_min+f_max));
    p.setX(x_);
    p.setY(y_);
}


void
Window::coord_in_window (QPoint &p)
{
    int w = width () - 1;
    int h = height () - 1;

//    int diagonal = sqrt((x_max - x_min)*(x_max - x_min)
//                        + (y_max - y_min)*(y_max - y_min));
    int x_ = w * ((p.x() - x_min) / /*diagonal*/(x_max - x_min));
//            + 0.5*w*(x_max/(x_max-x_min));
    int y_ = h * ((y_max - p.y()) / /*diagonal*/(y_max - y_min));
//            + 0.5*h*(y_max/(y_max-y_min));

    p.setX(x_);
    p.setY(y_);
}

void
Window::draw_background (QPainter &painter)
{
    QPolygon polygon;
    int border = -10;
    polygon.push_back (QPoint (border, border));
    polygon.push_back (QPoint (width () - border, border));
    polygon.push_back (QPoint (width () - border, height () - border));
    polygon.push_back (QPoint (border, height () - border));
    draw_polygon (painter, polygon);
}
void
Window::draw_polygon (QPainter & painter, const QPolygon & polygon)
{
  QPainterPath path;
  path.addPolygon (polygon);
  painter.fillPath (path, brash);/*
                                  *fill path with color in brash
                                  * execluding boundary */
  painter.drawPolygon (polygon);
}

void
Window::draw_func (QPainter & painter,
                   double (Window::* f) (double, double))
{
    double dx = x_len / (min(n, MAX_N));
    double dy = y_len / (min(m, MAX_M));
    double x = 0, y = 0;
    int i = 0, j = 0/*, k = 0*/, min_m = min(m, MAX_M),
            min_n = min(n, MAX_N),
            cnx_ = MAX_C*MAX_NX, dnx_ = MAX_D*MAX_NX,
            cny_ = MAX_C*MAX_NY, dny_ = MAX_D*MAX_NY;
    int min_cnx = min(cnx_, cnx), min_dnx = min(dnx_, dnx),
            min_cny = min(cny_, cny),
            min_dny = min(dny_, dny);
    int len = ((min_m+1) * (min_n+1))-(min_dnx-1)*(min_dny-1),
            shift = 0;

//    printf("(min_n, min_m) == (%d, %d)\n", min_n, min_m);
//    printf("len == %d\n", len);
    for (int t = 0; t < len; t++)
    {
        k2ij(min_n, t, i, j, min_cnx, min_dnx,
             min_cny, min_dny);
//        printf("t = %d\n", t);
//        printf("(i , j) == (%d, %d)\n", i, j);
//        printf("(dx, dy) == (%lf, %lf)\n", dx, dy);
        if((i == min_cnx && j > min_cny && j < min_cny+min_dny)
                ||(j == min_cny && i > min_cnx && i < min_cnx+min_dnx)
                ||( i == min_cnx && j == min_cny )
                ||( i == 0 && j == min_m )
                ||( i == min_n && j == min_m )
                ||( i == min_n && j == 0 )
                ||(j == min_m && i > 0 && i < min_n)
                ||(i == min_n && j > 0 && j < min_m))
        {
            shift++;
        }
        else // CounterClockWise
        {
            x = i * dx;
            y = j * dy;
            tree[2*(t-shift)].p[1].x = x;
            tree[2*(t-shift)].p[1].y = y;
            tree[2*(t-shift)+1].p[2].x = x;
            tree[2*(t-shift)+1].p[2].y = y;

            x = (i+1) * dx;
            y = j * dy;
            tree[2*(t-shift)+1].p[0].x = x;
            tree[2*(t-shift)+1].p[0].y = y;

            x = i * dx;
            y = (j+1) * dy;
            tree[2*(t-shift)].p[0].x = x;
            tree[2*(t-shift)].p[0].y = y;

            x = (i+1) * dx;
            y = (j+1) * dy;
            tree[2*(t-shift)].p[2].x = x;
            tree[2*(t-shift)].p[2].y = y;
            tree[2*(t-shift)+1].p[1].x = x;
            tree[2*(t-shift)+1].p[1].y = y;
        }
    }

    QPolygon triug (3);
    QPoint point;
    QPoint line [3];
    double x_ = 0, y_ = 0;
    double len_len = x_len+y_len;
    projection (x_, y_, 0, line[0]);
    projection (-1*(2*/*x*/len_len)*h_sin, 1*(2*/*y*/len_len)*h_cos, 0, line[1]);
    projection (1*(2*/*x*/len_len)*h_cos, 1*(2*/*y*/len_len)*h_sin, 0, line[2]);
//    printf("line[0] (x,y) == (%d, %d)\n"
//           "line[1] (x,y) == (%d, %d)\n"
//           "line[2] (x,y) == (%d, %d)\n",
//            line[0].x(), line[0].y(),
//            line[1].x(), line[1].y(),
//            line[2].x(), line[2].y());

    double koefx1, koefx2, koefy1, koefy2;
    koefx1 = (line[1].x()-line[0].x());
    koefy1 = (line[1].y()-line[0].y());
    koefx2 = (line[2].x()-line[0].x());
    koefy2 = (line[2].y()-line[0].y());

    if(((int)(angle*180./PI)%360>=90
            && fabs((int)(angle*180./PI)%360)<=359)
            || ((int)(angle*180./PI)%360<=0
            && (int)(angle*180./PI)%360>=-270))
    {
        painter.setPen(Qt::red);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x(),line[0].y()-height());
        painter.setPen(Qt::blue);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x()+koefx1*1/*(width()+height())*/,
                 line[0].y()+koefy1*1/*(width()+height())*/);
        painter.setPen(Qt::gray);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x()+koefx2*1/*(width()+height())*/,
                line[0].y()+koefy2*1/*(width()+height())*/);
    }

    Point glaz;
    x_ = cos(-angle-PI/2);
    y_ = sin(-angle-PI/2);
    glaz.x = x_*(x_len+y_len);
    glaz.y = y_*(x_len+y_len);

    sort_triangles2 (tree, glaz);

    painter.setPen(Qt::darkGreen/*black*/);
    for(int i = 0; i < 2*(len-shift); i++)
    {
        for(int t = 0; t < 3; t++)
        {
            x_ = tree[i].p[t].x * h_cos - tree[i].p[t].y * h_sin;
            y_ = tree[i].p[t].x * h_sin + tree[i].p[t].y * h_cos;
            projection (x_, y_,
                        (this->*f)(tree[i].p[t].x, tree[i].p[t].y)
                        , /*tree[i].p[t]*/point);
            triug.setPoint(t, point);
        }
        draw_polygon(painter, triug);
    }

    if(((int)(angle*180./PI)%360>=0
            && fabs((int)(angle*180./PI)%360)<=90)
            || ((int)(angle*180./PI)%360<=-270
            && (int)(angle*180./PI)%360>=-359))
    {
        painter.setPen(Qt::red);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x(),line[0].y()-height());
        painter.setPen(Qt::blue);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x()+koefx1*1/*(width()+height())*/
                ,line[0].y()+koefy1*1/*(width()+height())*/);
        painter.setPen(Qt::gray);
        painter.drawLine(line[0].x(),line[0].y(),
                line[0].x()+koefx2*1/*(width()+height())*/
                ,line[0].y()+koefy2*1/*(width()+height())*/);
    }
}
