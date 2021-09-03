#include "window.h"

static
double f_0 (double /*x*/)
{
  return 1;
}

static
double f_1 (double x)
{
  return x;
}

static
double f_2 (double x)
{
  return x * x;
}

static
double f_3 (double x)
{
  return x * x * x;
}

static
double f_4 (double x)
{
  return x * x * x * x;
}

static
double f_5 (double x)
{
  return exp(x);
}

static
double f_6 (double x)
{
  return 1./(25*x*x + 1);
}

static
double fd_0 (double /*x*/)
{
  return 0;
}

static
double fd_1 (double /*x*/)
{
  return 1;
}

static
double fd_2 (double x)
{
  return 2*x;
}

static
double fd_3 (double x)
{
  return 3 * x * x;
}

static
double fd_4 (double x)
{
  return 4 * x * x * x;
}

static
double fd_5 (double x)
{
  return exp(x);
}

static
double fd_6 (double x)
{
  return (-100.*x)/((25*x*x + 1)*(25*x*x + 1));
}

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = -1;
  render_id = 0;
  max_f = 0;
  p = 0;

  change_func ();
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  else if (argc!=5 )
    return -1;

  else if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || sscanf (argv[3], "%d", &n) != 1 || n <= 1
      || sscanf (argv[4], "%d", &func_id) != 1  || func_id < 0  || func_id > 6)
    return -2;

  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;    fd = fd_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;    fd = fd_1;
        break;
      case 2:
        f_name = "f (x) = x^2";
        f = f_2;    fd = fd_2;
        break;
      case 3:
        f_name = "f (x) = x^3";
        f = f_3;    fd = fd_3;
        break;
      case 4:
        f_name = "f (x) = x^4";
        f = f_4;    fd = fd_4;
        break;
      case 5:
        f_name = "f (x) = exp(x)";
        f = f_5;    fd = fd_5;
        break;
      case 6:
        f_name = "f (x) = 1./(25*x^2 + 1)";
        f = f_6;    fd = fd_6;
        break;
    }
  update ();
  sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
  printf("%s\n", buf);
  return 0;
}

/// change current function for drawing
void Window::change_func ()
{
  func_id = (func_id + 1) % 7;

  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;    fd = fd_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;    fd = fd_1;
        break;
      case 2:
        f_name = "f (x) = x^2";
        f = f_2;    fd = fd_2;
        break;
      case 3:
        f_name = "f (x) = x^3";
        f = f_3;    fd = fd_3;
        break;
      case 4:
        f_name = "f (x) = x^4";
        f = f_4;    fd = fd_4;
        break;
      case 5:
        f_name = "f (x) = exp(x)";
        f = f_5;    fd = fd_5;
        break;
      case 6:
        f_name = "f (x) = 1./(25*x^2 + 1)";
        f = f_6;    fd = fd_6;
        break;
    }
  p = 0;
  update ();
  sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
  printf("%s\n", buf);
}

void Window::change_render_func ()
{
    render_id = (render_id+1) % 6;
    update();
    sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
    printf("%s\n", buf);
}

void Window::scale_on_x_mult_2 ()
{
    if(fabs(a) < DBL_MAX/2 && fabs(b) < DBL_MAX/2)
    {
        a*=2;  b*=2;
        update();
        sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
        printf("%s\n", buf);
    }
}
void Window::scale_on_x_reduce_2 ()
{
    if(b-a >= 1.e-6)
    {
        a/=2;  b/=2;
        update();
        sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
        printf("%s\n", buf);
    }
}
void Window::n_mult_2 ()
{
    if(n < INT_MAX/2)
    {
        n*=2;
//        n++;
        update();
    }
}
void Window::n_reduce_2 ()
{
    if(n > 2)
    {
        n/=2;
//        n--;
        update();
    }
}
void Window::p_plus ()
{
    if( p < INT_MAX )
    {
        p++;
        update();
        sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
        printf("%s\n", buf);
    }
}
void Window::p_minus ()
{
    if(-p < INT_MAX)
    {
        p--;
        update();
        sprintf(buf, "max{|Fmin|,|Fmax|} = %e", max_f);
        printf("%s\n", buf);
    }
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{
  QPainter painter (this);
  double x1=0, x2=0, y1=0, y2=0;
  double max_y=0, min_y=0;
//  double delta_y = /*1.*/ (b-a)/ QMainWindow().height();
  double delta_x = /*1.*/ (b-a)/ QMainWindow().width();
  double delta_mm =0, delta_ab = 0;
  QPen pen_black(Qt::black, 0, Qt::SolidLine);
  QPen pen_green(Qt::darkGreen, 0, Qt::SolidLine);
  QPen pen_magenta(Qt::darkMagenta , 0, Qt::SolidLine);
  QPen pen_red(Qt::red, 0, Qt::SolidLine);

  painter.setPen (pen_black);

  double delta = (b-a)/(n-1);
  double residual_1 = 0, residual_2 = 0;
  long long wid = QMainWindow().width()+1;
  max_y = min_y = 0;
  delta_ab = b-a;

  double *func = new double[n+1];
  double *fder = new double[n+1];
  double *x = new double[n+1];
  double *res = new double[4*n];

  if(render_id < 3)
  {
      // calculate min and max for current function
      for (long long i = 0; i < wid; i++)
      {
          x1 = a + i*delta_x;
          y1 = f (x1);
          if (y1 < min_y)
            min_y = y1;
          if (y1 > max_y)
            max_y = y1;
      }
      delta_mm = 1.02*(max_y - min_y);

      // draw approximated line for graph
      painter.setPen (pen_black);
      x1 = a;
      y1 = f (x1);
      for(long long i  = 1; i < wid; i++)
      {
          x2 = a + i*delta_x;
          y2 = f (x2);
          painter.drawLine (QPoint ((width()*(x1-a))/delta_ab, (height()*(max_y - y1))/delta_mm),
                            QPoint ((width()*(x2-a))/delta_ab, (height()*(max_y - y2))/delta_mm));
          x1 = x2, y1 = y2;
      }

      // draw axis
      painter.setPen (pen_red);
      painter.drawLine (QPoint (0, (height()*(max_y))/delta_mm),
                        QPoint (width() , (height()*(max_y))/delta_mm));

      painter.drawLine (QPoint ((width()*(  -a))/delta_ab, 0),
                        QPoint ((width()*(  -a))/delta_ab, height() ));

      if( (render_id == 0 || render_id == 2) && n <= 300)
      {
          // draw aproximation from method_1
          painter.setPen (pen_green);
          delta = (b-a)/(n-1);

          x[0] = a;   func[0] = f(a);
          x[n-1] = b; func[n-1] = f(b);
          for(int i = 1; i < n-1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_1(n, x, func, res);

          x1 = a;
          y1 = find_value_1(x1, a, b, n, res);
          for (long long i = 1; i < wid; i++)
          {
              x2 = a + i*delta_x;
              y2 = find_value_1(x2, a, b, n, res);
              painter.drawLine (QPoint ((width()*(x1-a))/delta_ab, (height()*(max_y - y1))/delta_mm),
                                QPoint ((width()*(x2-a))/delta_ab, (height()*(max_y - y2))/delta_mm));
              x1 = x2, y1 = y2;
              if(fabs (f(x1) - y1) >= residual_1)
              {
                  residual_1 = fabs (f(x1) - y1);
              }
          }
      }
      if(render_id == 1 || render_id == 2)
      {
          // draw aproximation from method_2
          painter.setPen (pen_magenta);
          delta = (b-a)/(n-1);
          x[0] = a;  func[0] = f(x[0]);
          fder[0] = fd(a);

          for(int i = 1; i < n+1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
              fder[i] = fd(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_2(n, func, fder, res, delta);

          x1 = a;
          y1 = find_value_2(a, x, res, delta, a);
          residual_2 = fabs(f(x1) - y1);
          for (long long i = 1; i < wid; i++)
          {
              x2 = a + i*delta_x;
              y2 = find_value_2(x2, x, res, delta, a);
              painter.drawLine (QPoint ((width()*(x1-a))/delta_ab, (height()*(max_y - y1))/delta_mm),
                                QPoint ((width()*(x2-a))/delta_ab, (height()*(max_y - y2))/delta_mm));
              x1 = x2, y1 = y2;
              if( fabs( f(x2) - y2 ) >= residual_2)
              {
                  residual_2 = fabs(f(x2) - y2);
              }
          }
      }
  }
  else
  {
      // calculate min and max for current function
      if( (n <= 300) && (render_id == 3 || render_id == 4) )
      {
          //computation residual_1 for aproximation from method_1
          painter.setPen (pen_green);
          delta = (b-a)/(n-1);
          x[0] = a;   func[0] = f(a);
          x[n-1] = b; func[n-1] = f(b);

          for(int i = 1; i < n-1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_1(n, x, func, res);

          for (long long i = 0; i < wid; i++)
          {
              x2 = a + i*delta_x;
              y2 = find_value_1(x2, a, b, n, res);

              x1 = x2, y1 = y2;
              if(fabs (f(x2) - y2) >= residual_1)
              {
                  residual_1 = fabs (f(x2) - y2);
              }
          }
      }
      if(render_id == 3 || render_id == 5)
      {
          //computation residual_2 for aproximation from method_2
          painter.setPen (pen_magenta);
          delta = (b-a)/(n-1);
          x[0] = a;  func[0] = f(x[0]);
          fder[0] = fd(a);

          for(int i = 1; i < n+1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
              fder[i] = fd(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_2(n, func, fder, res, delta);

          for (long long i = 0; i < wid; i++)
          {
              x2 = a + i*delta_x;
              y2 = fabs( find_value_2(x2, x, res, delta, a) - f(x2) );

              x1 = x2, y1 = y2;
              if( y2 >= residual_2)
              {
                  residual_2 = y2;
              }
          }
      }

      max_y = residual_1 > residual_2 ? residual_1: residual_2;
      delta_mm = 1.02*(max_y - min_y);

      if(delta_mm < EPS)
      {
          delta_mm = EPS*1.02;
          max_y = EPS;
      }
      // draw axis
      painter.setPen (pen_red);
      painter.drawLine (QPoint (0      , (height()*(max_y))/delta_mm),
                        QPoint (width(), (height()*(max_y))/delta_mm));

      painter.drawLine (QPoint ((width()*(  -a))/delta_ab, 0),
                        QPoint ((width()*(  -a))/delta_ab, (height()*(max_y - min_y))/delta_mm));

      if( n <= 300 && (render_id == 3 || render_id == 4))
      {
          // draw aproximation from method_1
          painter.setPen (pen_green);
          delta = (b-a)/(n-1);
          x[0] = a;   func[0] = f(a);
          x[n-1] = b; func[n-1] = f(b);

          for(int i = 1; i < n-1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_1(n, x, func, res);

          x1 = a;
          y1 = fabs (find_value_1(x1, a, b, n, res) - f(x1));
          for (long long i = 1; i < wid; i ++)
          {
              x2 = a + i * delta_x;
              y2 = fabs (find_value_1(x2, a, b, n, res) - f(x2));

              painter.drawLine (QPoint ((width()*(x1-a))/delta_ab, (height()*(max_y - y1))/delta_mm),
                                QPoint ((width()*(x2-a))/delta_ab, (height()*(max_y - y2))/delta_mm));
              x1 = x2, y1 = y2;
          }
      }
      if(render_id == 3 || render_id == 5)
      {
          // draw aproximation from method_2
          painter.setPen (pen_magenta);
          delta = (b-a)/(n-1);
          x[0] = a;  func[0] = f(x[0]);
          fder[0] = fd(a);

          for(int i = 1; i < n+1; i ++)
          {
              x[i] = x[i-1] + delta;
              func[i] = f(x[i]);
              fder[i] = fd(x[i]);
          }
          func[n/2] += 0.1*p*max_f;
          find_koef_2(n, func, fder, res, delta);

          x1 = a;
          y1 = fabs(find_value_2(x1, x, res, delta, a) - f(x1));
          for (long long i = 1; i < wid; i++)
          {
              x2 = a + i*delta_x;
              y2 = fabs( find_value_2(x2, x, res, delta, a) - f(x2) );
              painter.drawLine (QPoint ((width()*(x1-a))/delta_ab, (height()*(max_y - y1))/delta_mm),
                                QPoint ((width()*(x2-a))/delta_ab, (height()*(max_y - y2))/delta_mm));

              x1 = x2, y1 = y2;
          }
      }
  }

  delete [] func;
  delete [] fder;
  delete [] x;
  delete [] res;

  //print on screen another data
  if(render_id == 0 || render_id == 2 || render_id == 3 || render_id == 4)
  {
      painter.setPen ("green");
      sprintf(buf, "Residual_1 = %e", residual_1);
      painter.drawText (0, 150, buf);
  }
  if(render_id == 1 || render_id == 2 || render_id == 3 || render_id == 5)
  {
      painter.setPen (pen_magenta);
      sprintf(buf, "Residual_2 = %e", residual_2);
      painter.drawText (0, 170, buf);
  }

  painter.setPen ("blue");

  painter.drawText (90, 50, f_name);

  sprintf(buf, "k = %d", func_id);
  painter.drawText (0, 50, buf);

  sprintf(buf, "a = %e", a);
  painter.drawText (0, 70, buf);

  sprintf(buf, "b = %e", b);
  painter.drawText (0, 90, buf);

  sprintf(buf, "max{|Fmin|,|Fmax|} = %e", fabs(min_y) > fabs(max_y) ? fabs(min_y) : fabs(max_y));
  painter.drawText (0, 110, buf);

  sprintf(buf, "n = %d", n);
  painter.drawText (0, 130, buf);

  sprintf(buf, "p = %d", p);
  painter.drawText (0, 190, buf);

  sprintf(buf, "Render_id = %d", render_id);
  painter.drawText (0, 210, buf);

  max_f = fabs(min_y) > fabs(max_y) ? fabs(min_y) : fabs(max_y);
}
