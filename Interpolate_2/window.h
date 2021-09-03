#ifndef WINDOW_H
#define WINDOW_H

#include <QApplication>
#include <QPainter>
#include <QDebug>
#include <QWidget>
#include <QtWidgets/QtWidgets>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "head.h"

#define PI 3.141592653589793
#define MAX_C 3  //fundamental web constants
#define MAX_NX 3 //fundamental web constants
#define DEF_ANGLE 30

class Window : public QWidget
{
    Q_OBJECT

public://things for draw now
    double *x = nullptr;//current vector of solve: Ax == b
    int N = 0;
    double angle = 0;
    double h_sin = 0;
    double h_cos = 0;
//    int k = 0;num function for ploting, lie in storage
    int n = 0;
    int m = 0;
    int dny = 0;
    int dnx = 0;
    int cnx = 0;
    int cny = 0;
    double hx = 0;
    double hy = 0;
    int s = 0;
    int func_count = 8;
    int render_count = 3;
    int render_id = 0;
    int save_int = 0;
    Storage store [1];//things for count *x and N
    Args *arg = nullptr;
    double residual = 0;
    int disturbance = 0;
//    Args *arg = nullptr;

    QPen pen;
    QBrush brash ;
    QTimer t;
    pthread_t *thread_id = nullptr;
    pthread_barrier_t barrier;
    pthread_barrier_t barrier_all;
    pthread_mutex_t mutex_all;
    double (*f) (double, double) = nullptr;
    double draw_max = 0;
    double draw_min = 0;
    double x_len = 0;
    double y_len = 0;
    double x_min = 0;/*min coord x in window*/
    double x_max = 0;/*max coord x in window*/
    double y_min = 0;/*min coord y in window*/
    double y_max = 0;/*max coord y in window*/
    double f_max = 0;//max of current function on the picture
    double norm = 0; //for current function that equals store->k
    Triangle *tree = nullptr;
    int MAX_D = 0;
    int MAX_NY = 0;
    int MAX_N = 0;
    int MAX_M = 0;
    double v_cos = cos(DEF_ANGLE * PI / 180.0);
    double v_sin = sin(DEF_ANGLE * PI / 180.0);
    QString name;
    //    QAction *action = nullptr;

public:
  Window (QWidget *parent, int &argc, char *argv[]);
  ~Window();

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

public slots:
  void update_screen();
//  void change_func ();
//  void change_render_func ();
//  void scale_on_x_mult_2 ();
//  void scale_on_x_reduce_2 ();
//  void n_mult_2 ();
//  void n_reduce_2 ();
//  void p_plus ();
//  void p_minus ();

protected:
  void paintEvent (QPaintEvent *event);
  void keyPressEvent (QKeyEvent * event);
  void closeEvent(QCloseEvent *event);

private:
  double func_value (double x, double y);
  double app_value  (double x, double y);
  double err_value  (double x, double y);
  void draw_background (QPainter &painter);
  void draw_polygon (QPainter & painter, const QPolygon & polygon);
  void find_min_max (double (Window::* f) (double, double));
  void draw_func (QPainter & painter, double (Window::* f) (double, double));
  void sort_triangles (Triangle *tr, QPoint eye);
  void sort_triangles2 (Triangle *tr, Point &eye);
  double center (Triangle &tr);
  void rotate (Triangle *tr, int len);
  void corotate (Triangle *tr, int len);
  void projection (double x, double y, double z, QPoint &p);
  void projection (double x, double y, double z, Point &p);
  void coord_in_window (QPoint &p);
};

//int compare(const void* a, const void* b);

#endif // MAINWINDOW_H
