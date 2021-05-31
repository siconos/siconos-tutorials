from OCC.Display.SimpleGui import *
#set_backend('qt')
#set_backend('X')

display, start_display, add_menu, add_function_to_menu = init_display()
def simple_test(event=None):
    display.Test()
def simple_cylinder(event=None):
    from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    s = BRepPrimAPI_MakeCylinder(60, 200)
    display.DisplayShape(s.Shape())

add_menu('simple test')
add_function_to_menu('simple test',simple_test)
add_function_to_menu('simple test',simple_cylinder)

display.View_Iso()
display.FitAll()
# display loop
start_display()
