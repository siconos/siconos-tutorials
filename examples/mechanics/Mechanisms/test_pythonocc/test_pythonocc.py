from OCC.Display.SimpleGui import *
from OCC.Core.BRepPrimAPI import *
display, start_display, add_menu, add_function_to_menu = init_display()
my_box = BRepPrimAPI_MakeBox(10.,20.,30.).Shape()
display.DisplayShape(my_box)
start_display()
