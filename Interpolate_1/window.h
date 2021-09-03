#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include <QPainter>

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#include <sstream>

#define DEFAULT_A -1
#define DEFAULT_B 1
#define DEFAULT_N 2
#define PI 3.141592653589793
#define EPS 1.e-16

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int render_id;
  const char *f_name;
  double a, b, max_f;
  int n, p;
  double (*f) (double);
  double (*fd) (double);
  char buf[256];

public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  void find_mm_y(double *func, double *fd, double *res, double *x,
                 double delta, double delta_x, double *error1,
                 double *error2);

public slots:
  void change_func ();
  void change_render_func ();
  void scale_on_x_mult_2 ();
  void scale_on_x_reduce_2 ();
  void n_mult_2 ();
  void n_reduce_2 ();
  void p_plus ();
  void p_minus ();

protected:
  void paintEvent (QPaintEvent *event);
};


void find_koef_1(int n, double *x, double *f, double *a);
double find_value_1(double x, double a, double b, int n, double *koef);
void find_koef_2(int n, double *f, double *fd, double *a, double delta);
double find_value_2(double y, double *x, double *a, double delta, double le);

#endif

