# This part has to be reviewed
# Qt finding process has changed and our FindQGLViewer is outdated.

set(QTCOMPONENTS Core Xml OpenGL Gui Widgets)
foreach(comp ${QTCOMPONENTS})
endforeach()
find_package(X11 REQUIRED)
find_package(QGLViewer REQUIRED)

# See "modern cmake", we have to use target_link_libraries(EXE QGLViewer ...)
target_link_libraries(QGLSpheres Qt4::QtCore Qt4::QtXml Qt4::QtOpenGL Qt4::QtGui Qt4::QtWidgets)
target_link_libraries(QGLSpheres QGLViewer)

