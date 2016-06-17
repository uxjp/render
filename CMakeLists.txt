project(render)
cmake_minimum_required(VERSION 2.8)


set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)


set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  # ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules
     /usr/lib/CGAL
    )


#flags
    # for the libs  <vector> and  <unordered_set> junior
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#Libraries

    # Qt library
    find_package(Qt5Widgets)
    find_package(Qt5OpenGL)

    find_package(Qt5Xml REQUIRED) # package necessary for QGLVIEWER!!!!!!!!!!!
    find_package(QGLViewer)

    # OpenGL
    find_package(OpenGL)


    # CGAL and its components
    find_package( CGAL QUIET COMPONENTS Core )
    if ( NOT CGAL_FOUND )
      message(STATUS "This project requires the CGAL library, and will not be compiled.")
      return()
    endif()
    # include helper file
       include( ${CGAL_USE_FILE} )

    # GLEW
    find_package(GLEW REQUIRED)
    if (GLEW_FOUND)
       include_directories(${GLEW_INCLUDE_DIRS})
       link_libraries(${GLEW_LIBRARIES})
    endif()


# executable
set (Sources WidgetOGL.cpp main.cpp)
set (Headers WidgetOGL.h )

set (Shaders loadShaders/shader.hpp  loadShaders/shader.cpp
             loadShaders/shaders/obj.vert loadShaders/shaders/obj.frag)

# link of executable with files
add_executable(${PROJECT_NAME} ${Sources}  ${Headers} ${Shaders} )

include_directories (${CMAKE_CURRENT_BINARY_DIR})
qt5_use_modules(${PROJECT_NAME}  Widgets OpenGL )
qt5_use_modules(${PROJECT_NAME} Xml)


target_link_libraries(${PROJECT_NAME}  ${QT_LIBRARIES}
                                       ${OPENGL_LIBRARIES}
                                       ${QGLVIEWER_LIBRARIES}
                                        ${GLEW_LIBRARIES}
                                       )


