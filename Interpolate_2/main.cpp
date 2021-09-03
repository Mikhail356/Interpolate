#include "window.h"
#include "head.h"

int
main (int argc, char * argv [])
{
  QApplication app (argc, argv);
  Window wind(nullptr, argc, argv);
  app.setApplicationName ("Interpolate 2");
  wind.show ();
  return app.exec ();
}

//  QAction *action;

//  action = tool_bar->addAction ("&0 Change function", graph_area, SLOT (change_func ()));
//  action->setShortcut (QString ("0"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&1 Change render", graph_area, SLOT (change_render_func () ));
//  action->setShortcut (QString ("1"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&2 Scale x*2", graph_area, SLOT (scale_on_x_mult_2 ()));
//  action->setShortcut (QString ("2"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&3 Scale x/2", graph_area, SLOT (scale_on_x_reduce_2 ()));
//  action->setShortcut (QString ("3"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&4 N*2", graph_area, SLOT (n_mult_2 ()));
//  action->setShortcut (QString ("4"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&5 N/2", graph_area, SLOT (n_reduce_2 ()));
//  action->setShortcut (QString ("5"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&6 Value_f_in_x[n/2] ++", graph_area, SLOT (p_plus ()));
//  action->setShortcut (QString ("6"));
//  action->setCheckable(true);
//  action = tool_bar->addAction ("&7 Value_f_in_x[n/2] --", graph_area, SLOT (p_minus ()));
//  action->setShortcut (QString ("7"));
//  action->setCheckable(true);

//  action = tool_bar->addAction ("&Exit", window, SLOT (close ()));
//  action->setShortcut (QString ("Ctrl+X"));
//  action->setCheckable(true);

//  window->setMenuBar (tool_bar);
