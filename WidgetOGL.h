#ifndef WIDGETOGL_H
#define WIDGETOGL_H


// Include standard headers
    #include <stdio.h>
    #include <stdlib.h> // abs
    #include <math.h> // pow
    #include <string>
    #include <vector>
    #include <iostream>
    #include <fstream>
    #include <algorithm>
    #include <sstream>
    using namespace std;



// Include GLM
    #include <glm/glm.hpp>
    #include "glm/gtc/matrix_transform.hpp"
    #include <glm/gtx/rotate_vector.hpp>

// QT
    #include <QtGui>
    #include <QGLWidget>
    #include <QOpenGLWidget>
    #include <QOpenGLFunctions_3_3_Core>


// CGAL for HPR
    #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
    #include <CGAL/Polyhedron_3.h>
    #include <CGAL/convex_hull_3.h>
    #include <CGAL/Timer.h>

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::Vector_3 Vector_3;
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
    typedef Polyhedron::Vertex_iterator Vertex_iterator;

#include <QGLViewer/camera.h>

struct attributes {
    GLfloat coord[3];
 //   GLfloat uv_tex[2];
};



class WidgetOGL: public QGLWidget,protected QOpenGLFunctions_3_3_Core
{
 Q_OBJECT
public:
    WidgetOGL(const QGLFormat &format, QGLWidget *parent = 0);
    ~WidgetOGL();

// VARIABLES//
public:

    CGAL::Timer timer;

    std::vector<glm::vec3> points_original;
    std::vector<Point_3> m_points;
    std::vector<glm::vec3> points;
    std::vector<glm::vec3> color;

    // mouse
    QPoint lastPos;

    //variables of moviment of mouse
    GLfloat angleX=0;
    GLfloat angleY=0;
    GLfloat angleZ=0;
    GLfloat transX=0;
    GLfloat transY=0;
    GLfloat transZ=0;

    // variables of view
    glm::vec3 cameraPos;
    glm::vec3 camera;
    glm::vec3 cible;
    glm::vec3 up;

    // matrix
    glm::mat4 Model ;
    glm::mat4 View ;
    glm::mat4 Projection;

    // adress shader
    GLuint MVPID ;
    GLuint MVID;

    // gldraws
    int numbers_points;

    //link shader
    GLuint program;

    //VBOs
    GLuint  vbo_attributes_vertex;
    GLuint  vbo_attributes_color;

    // window size
    int width ,height;


    // toogles
    bool HPR_toogle= false;
    bool HPR_action =false;
    bool toogle_read_pixels = false;

    // test
    int i=0;
    double R = 10000,boost=1;

// FUNCTIONS //
protected:

    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

    //functions shader
    GLuint create_shader( char *file_path, GLenum type);
    std::string file_read(char *file_path);
    void print_log(GLuint object);
    GLuint link(GLuint ProgramID, GLuint VertexShaderID, GLuint FragmentShaderID);

    // function get data
    std::vector<glm::vec3> get_data_xyz(char* file_path);
    unsigned char *loadBMP_custom(const char * imagepath);

    // get matrix
    glm::mat4 getMatrixModelViewProjection();
    glm::mat4 getMatrixModelView();
    // get camera position
    glm::vec3 getCameraPosition();

    // mouse
    void mousePressEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void keyPressEvent(QKeyEvent *k);

    // HPR
        // simple HPR
        void HPR(std::vector<Point_3> points,glm::vec3 view_point, double R);
            std::vector<Point_3> transformation(std::vector<Point_3> points,Vector_3 view_point,double R);
            Point_3 transform(Point_3 p , Vector_3 view_point , double R );
            std::vector<glm::vec3>  visibility(Polyhedron poly,Vector_3 view_point,double R);

        // HPR optmal
         void HPR_ray_optmal(std::vector<Point_3> points , glm::vec3 view_p );
             Vector_3 view_point_opposite(std::vector<Point_3> points, Vector_3 view_point );
             double calcul_R_min(std::vector<Point_3> points,Vector_3 view_point);
             double calcul_R_optmal(std::vector<Point_3> points, Vector_3 view_point,
                                               Vector_3 view_point_opp, double Rmin);
             double function_disjoint(std::vector<Point_3> points,Vector_3 view_point,
                                                  Vector_3 view_point_opp, double R);
                int intersection(Polyhedron poly1, Vector_3 view_point ,
                             Polyhedron poly2, Vector_3 view_point_opp, double R);
                    bool compare_point(Point_3 p1,Point_3 p2);
             double grad_function_disjoint(double F_R_current, double F_R_before,
                                                  double R_current, double R_before);

    // Reading pixels
    void read_pixels();

};

#endif // WIDGETOGL_H
