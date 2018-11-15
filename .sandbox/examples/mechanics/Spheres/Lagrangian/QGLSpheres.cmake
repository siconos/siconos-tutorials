# This part has to be reviewed
# Qt finding process has changed and our FindQGLViewer is outdated.

set(QTCOMPONENTS QtCore QtXml QtOpenGL QtGui QtWidgets)
find_package(Qt4 REQUIRED ${QTCOMPONENTS}) # Qt5??
find_package(X11 REQUIRED)
find_package(QGLViewer REQUIRED)
# find_package(OpenGL REQUIRED)  Check this https://cmake.org/cmake/help/v3.13/module/FindOpenGL.html
# See "modern cmake", we have to use target_link_libraries(EXE QGLViewer ...)
target_link_libraries(QGLSpheres Qt4::QtCore Qt4::QtXml Qt4::QtOpenGL Qt4::QtGui Qt4::QtWidgets)
target_link_libraries(QGLSpheres QGLViewer)

