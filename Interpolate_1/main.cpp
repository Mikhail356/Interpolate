#include "window.h"

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action;

  if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!",
                            "Usage: ./basic_graph a b n k "
                            "or ./basic_graph");
      return -1;
    }

  action = tool_bar->addAction ("&0 Change function", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("0"));
  action = tool_bar->addAction ("&1 Change render", graph_area, SLOT (change_render_func () ));
  action->setShortcut (QString ("1"));
  action = tool_bar->addAction ("&2 Scale x*2", graph_area, SLOT (scale_on_x_mult_2 ()));
  action->setShortcut (QString ("2"));
  action = tool_bar->addAction ("&3 Scale x/2", graph_area, SLOT (scale_on_x_reduce_2 ()));
  action->setShortcut (QString ("3"));
  action = tool_bar->addAction ("&4 N*2", graph_area, SLOT (n_mult_2 ()));
  action->setShortcut (QString ("4"));
  action = tool_bar->addAction ("&5 N/2", graph_area, SLOT (n_reduce_2 ()));
  action->setShortcut (QString ("5"));
  action = tool_bar->addAction ("&6 Value_f_in_x[n/2] ++", graph_area, SLOT (p_plus ()));
  action->setShortcut (QString ("6"));
  action = tool_bar->addAction ("&7 Value_f_in_x[n/2] --", graph_area, SLOT (p_minus ()));
  action->setShortcut (QString ("7"));

  action = tool_bar->addAction ("&Exit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Interpolate_1");

  window->show ();
  app.exec ();
  delete window;
  return 0;
}
