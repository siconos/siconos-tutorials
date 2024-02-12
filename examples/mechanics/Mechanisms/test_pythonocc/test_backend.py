from OCC.Display.SimpleGui import *
#set_backend('qt')
#set_backend('X')

display, start_display, add_menu, add_function_to_menu = init_display()
def simple_test(event=None):
    display.Test()
def simple_cylinder(event=None):
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
    s = BRepPrimAPI_MakeCylinder(60, 200)
    display.DisplayShape(s.Shape())
def simple_box(event=None):   
    from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
    my_box = BRepPrimAPI_MakeBox(10.0, 200.0, 300.0).Shape()
    display.DisplayShape(my_box, update=True)

    
def step_file(event=None):
     from OCC.Core.STEPControl import STEPControl_Reader
     from OCC.Core.TopExp import TopExp_Explorer
     from OCC.Core.TopAbs import TopAbs_FACE
     from OCC.Core.StepRepr import StepRepr_RepresentationItem

     reader = STEPControl_Reader()
     tr = reader.WS().TransferReader()
     reader.ReadFile('../Trip_Magnetic/CAD/Case/Case_Needle_aasembly.stp')
     reader.ReadFile('toto.step')
     reader.TransferRoots()
     shape = reader.OneShape()
     display.DisplayShape(shape, update=True)


    

add_menu('simple test')
add_function_to_menu('simple test',simple_test)
add_function_to_menu('simple test',simple_cylinder)
add_function_to_menu('simple test',simple_box)
add_function_to_menu('simple test',step_file)





display.View_Iso()
display.FitAll()
# display loop
start_display()
