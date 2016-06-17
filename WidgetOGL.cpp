// Include GLEW
//#include <GL/glew.h>

#include <QKeyEvent>

#include "WidgetOGL.h"

#define PI 3.14159265

#define BUFFER_OFFSET(i)   ((char *)NULL + (i))

// construtor
WidgetOGL::WidgetOGL(const QGLFormat& format, QGLWidget *parent)
    :QGLWidget(format,parent)
{

 width = 256;
 height = 256;

 camera = glm::vec3(0,0,3);
 cible = glm::vec3(0,0,0);
 up = glm::vec3(0,1,0);

 Model = glm::mat4(1.0f);
 View = glm::lookAt(camera, // Camera is at (4,3,3), in World Space
                    cible, // and looks at the origin
                    up);  // Head is up (set to 0,-1,0 to look upside-down)

 Projection = glm::perspective(45.0f,
                              (float)width/ (float)height,
                              0.01f,
                              100.0f);



}
// Desconstrutor
WidgetOGL::~WidgetOGL()
{
    // clean
     glDeleteBuffers(1, &vbo_attributes_vertex);
     // clean
      glDeleteBuffers(1, &vbo_attributes_color);


}


//modulated functions related to shader----------
GLuint WidgetOGL::create_shader( char* file_path, GLenum type)
{
 std::string shaderCode = file_read(file_path);
  if (shaderCode.empty()) {
    fprintf(stderr, "Erro abrindo: %s: ", file_path); perror("");
    return 0;
  }

//  printf ("%s", shaderCode.c_str());
  printf("Compiling : %s ...\n"  ,file_path);


  GLuint sourceId = glCreateShader(type);
  const GLchar *source = shaderCode.c_str();
  glShaderSource(sourceId, 1, &source , NULL);
   // free( (void*) source);
  glCompileShader(sourceId);


  GLint compile_ok = GL_FALSE;
  glGetShaderiv(sourceId, GL_COMPILE_STATUS, &compile_ok);
      if (compile_ok == GL_FALSE) {
        fprintf(stderr, "%s:", file_path);

        print_log(sourceId);// function log created

        glDeleteShader(sourceId);
        return 0;
      }

  return sourceId;


}
std::string WidgetOGL::file_read( char* file_path)
{
    std::string shaderCode;
        std::ifstream shaderStream(file_path, std::ios::in);
        if(shaderStream.is_open()){
            std::string Line = "";
            while(getline(shaderStream, Line))
                shaderCode += "\n" + Line;
            shaderStream.close();
        }else{
            printf("Impossible to open %s. Are you in the right directory ? !\n", file_path);
            getchar();
            return 0;
        }

        return shaderCode;
}
void WidgetOGL::print_log(GLuint object)
{
  GLint log_length = 0;
  if (glIsShader(object))
    glGetShaderiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else if (glIsProgram(object))
    glGetProgramiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else {
    fprintf(stderr, "printlog: Não é um shader ou programa\n");
    return;
  }

  char* log = (char*)malloc(log_length);

  if (glIsShader(object))
    glGetShaderInfoLog(object, log_length, NULL, log);
  else if (glIsProgram(object))
    glGetProgramInfoLog(object, log_length, NULL, log);

  fprintf(stderr, "%s", log);
  free(log);
}
GLuint WidgetOGL::link(GLuint ProgramID, GLuint VertexShaderID, GLuint FragmentShaderID){


    ProgramID = glCreateProgram();
    printf("Linking program\n");
        glAttachShader(ProgramID, VertexShaderID);
        glAttachShader(ProgramID, FragmentShaderID);
        glLinkProgram(ProgramID);

        GLint link_ok = GL_FALSE;
        glGetProgramiv(ProgramID, GL_LINK_STATUS, &link_ok);
        if (!link_ok) {
          fprintf(stderr, "error in linkage:");
           print_log(ProgramID);
           return 0;
        }
    return ProgramID;
}
//-----------------------------------------------

// Matrix /-----------------------------------------
 glm::mat4 WidgetOGL::getMatrixModelViewProjection(){

return  Projection * View * Model;;

}
 glm::mat4 WidgetOGL::getMatrixModelView(){
   // Matrix
return  View * Model;

}
 glm::vec3  WidgetOGL::getCameraPosition(){

    glm::mat4 viewModel = glm::inverse(View*Model);
    return glm::vec3(viewModel[3]);

 }
//------------------------------------------------

// get data----------------------------------------
unsigned char* WidgetOGL::loadBMP_custom(const char * imagepath) {

    // Data read from the header of the BMP file
    unsigned char header[54]; // Each BMP file begins by a 54-bytes header
    unsigned int dataPos;     // Position in the file where the actual data begins
    unsigned int width, height;
    unsigned int imageSize;   // = width*height*3

    // Actual RGB data
   unsigned char * img_data;

    // Open the file
    FILE * file = fopen(imagepath,"rb");
    if (!file){printf("Image could not be opened\n"); return 0;}

    if ( fread(header, 1, 54, file)!=54 ){ // If not 54 bytes read : problem
        printf("Not a correct BMP file\n");
        return false;
    }
    if ( header[0]!='B' || header[1]!='M' ){
        printf("Not a correct BMP file\n");
        return 0;
    }

    // Read ints from the byte array
    dataPos    = *(int*)&(header[0x0A]);
    imageSize  = *(int*)&(header[0x22]);
    width      = *(int*)&(header[0x12]);
    height     = *(int*)&(header[0x16]);

    // Some BMP files are misformatted, guess missing information
    if (imageSize==0)    imageSize=width*height*3; // 3 : one byte for each Red, Green and Blue component
    if (dataPos==0)      dataPos=54; // The BMP header is done that way

    // Create a buffer
    img_data = new unsigned char [imageSize];

    fread(img_data,1,imageSize,file);

    //Everything is in memory now, the file can be closed
     fclose(file);

    return img_data;
}
std::vector<glm::vec3>  WidgetOGL::get_data_xyz(char* file_path){

   std::vector<float> vertex;
  // std::vector<glm::vec3> points;
   vertex.clear();points.clear(); m_points.clear();

    int aux=0;

    int size = 20000; int init = 0;

    std::string shaderCode;
    std::ifstream Stream(file_path, std::ios::in);
        if(Stream.is_open())
        {
            std::string Line = "";
            std::string word = "";

            while(getline(Stream, Line,'\n'))
            {
                std::stringstream line(Line);

                while(std::getline(line, word,' '))
                {
                    float value;
                    std::stringstream(word)>>value;
                    vertex.push_back(value);
                }

                if(aux!=0 /*&& (aux>init && aux<(init+size))*/)
                {
                    points.push_back(glm::vec3(vertex[0],vertex[1],vertex[2]));
                    m_points.push_back(Point_3(vertex[0],vertex[1],vertex[2]));
                    color.push_back(glm::vec3(vertex[4]/255,vertex[5]/255,vertex[6]/255)) ;
                }
               // if(aux>(init+size)){break;}

              aux++;
              vertex.clear();


            }
            Stream.close();
        }

        else{
            printf("Impossible to open %s. Are you in the right directory ? !\n", file_path);
            getchar();

        }
            std::cout<<points.size()<<" done!";

        /*    for(int i = 0;i<points.size() -1;i++){
          std::cout<<"point:"<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<
                    "color:"<<color[i].x<<" "<<color[i].y<<" "<<color[i].z<<"\n\n";
           }*/


 points_original = points;
 return points;

}
//--------------------------------------------------

// EVENTS-------------------------------------------
void WidgetOGL::mousePressEvent(QMouseEvent *e){
     lastPos = e->pos();

}
void WidgetOGL::mouseMoveEvent(QMouseEvent *e){
    GLfloat dx = ( e->pos().x() - lastPos.x())/( GLfloat)width ;
    GLfloat dy = (e->pos().y() - lastPos.y())/(GLfloat)height;


    if (e->buttons() == Qt::LeftButton ){
         angleX = angleX + 100*dy;
         angleY = angleY + 100*dx;

    }

    else if(e->buttons() == Qt::RightButton){
         angleZ = angleZ +  100*dy;
    }
    else if(e->buttons() == Qt::MiddleButton){
         transZ = transZ + 10*dy;
    }


    View = glm::rotate(View,angleX,glm::vec3(1,0,0));
    View = glm::rotate(View,angleY,glm::vec3(0,1,0));
    View = glm::rotate(View,angleZ,glm::vec3(0,0,1));
    // translate view

    View = glm::translate(View,glm::vec3(0,0,transZ));

 // update
 WidgetOGL::update();
    lastPos = e->pos();
     angleX =0;angleY =0;angleZ =0;
    /* transX =0;transY =0;*/transZ = 0;


}
void WidgetOGL::keyPressEvent(QKeyEvent *k){

      // navigation
      if (k->key() == Qt::Key_Z ){
          transY = transY - 0.1;
      }
      if (k->key() == Qt::Key_S ){
          transY = transY + 0.1;
      }
      if (k->key() == Qt::Key_D ){
          transX = transX + 0.1;
      }
      if (k->key() == Qt::Key_Q ){
          transX = transX - 0.1;
      }

      //changement
      View = glm::translate(View,glm::vec3(transX,transY,0));
      // update
      WidgetOGL::update();
      transY = 0;transX = 0;

      // to change the ray
      if (k->key() == Qt::Key_Up ){
        R = R + boost*100; qDebug()<<R;
      }
      if (k->key() == Qt::Key_Down ){
        R = R - boost*100; qDebug()<<R;
      }
      if (k->key() == Qt::Key_Right ){
        boost++; qDebug()<<boost;
      }
      if (k->key() == Qt::Key_Left ){
        boost--; qDebug()<<boost;
      }

      // coordenates camera
      if (k->key() == Qt::Key_C ){
         qDebug()<<cameraPos.x<<","<<cameraPos.y<<","<<cameraPos.z<<",";
      }

      //Active/Desactive HPR temp réel
      if (k->key() == Qt::Key_Space ){
            HPR_toogle = !HPR_toogle;
            if(!HPR_toogle){
                points = points_original;
            }
      }

      // call HPR
      if (k->key() == Qt::Key_H ){
       HPR_action = !HPR_action;

       if(HPR_action){
        timer.start();
           HPR(m_points,-getCameraPosition(), R);
        timer.stop();qDebug()<<"time:"<<timer.time()<<"s";timer.reset();
       }
       else
         points = points_original;
    }

      // call HPR ray optmal
      if (k->key() == Qt::Key_O){
       HPR_action = !HPR_action;
       if(HPR_action){
           timer.start();
             HPR_ray_optmal(m_points,-getCameraPosition());
           timer.stop();qDebug()<<"time:"<<timer.time()<<"s";timer.reset();
        }
       else
           points = points_original;
    }

      // reade pixels
      if (k->key() == Qt::Key_I){
       //   toogle_read_pixels = !toogle_read_pixels;
          read_pixels();
      }


}
// ---------------------------------------------


// HPR ----------------------------------------------------------------------
    // Simpe HPR
    void WidgetOGL::HPR(std::vector<Point_3> points,glm::vec3 view_p, double R){

       // view
       Vector_3 view_point = Vector_3(view_p.x ,view_p.y,view_p.z);

       // transformation
       std::vector<Point_3> points_trans ;
       points_trans = transformation(points, view_point, R);

       // hull convex
       Polyhedron poly;
       CGAL::convex_hull_3(points_trans.begin(), points_trans.end(), poly);

       // clear points
       points_trans.erase(points_trans.begin() , points_trans.end());

       // points visibles
        this->points = visibility(poly,view_point,R);

       // clear polyhedron
        poly.erase_all();


    }
        std::vector<Point_3> WidgetOGL::transformation(std::vector<Point_3> points,Vector_3 view_point,double R){

        // transformation related to origin (0,0,0)
            std::vector<Point_3> points_transform;
            Vector_3 v_pt;
            Point_3 pt;
            int aux=0;
            std::vector<Point_3>::const_iterator it;
            double   norm_pt;

            for(it = points.begin();
                it != points.end();
                it++)
            {

                pt = *it;

                pt = transform(pt, view_point, R);
                points_transform.push_back(pt);

                /* clean the m_points for to clean memory?*/

                aux++;
            }
            // add view_point to the set of points transformeds
               points_transform.push_back(CGAL::ORIGIN -view_point);

         return  points_transform;
        }
        Point_3 WidgetOGL::transform(Point_3 p , Vector_3 view_point , double R ){
            // camera to origin
            p = p + view_point;
            // transform point to vector (preparation for product)
            Vector_3 v_p = p-CGAL::ORIGIN;

            // norm -> distance of point to centre (view_point)
             double norm_p = sqrt(v_p*v_p);

           // transformation
            v_p = v_p*((2*R - norm_p)/(norm_p));

            p = CGAL::ORIGIN + v_p;

            // store the point in list of points
            p = p - view_point;

            return p;
        }
        std::vector<glm::vec3>  WidgetOGL::visibility(Polyhedron poly, Vector_3 view_point, double R){

        Point_3 pt_inv;
        Vector_3 v_pt_inv;
        std::vector<glm::vec3> points_visibles;

        Vertex_iterator it;
        for(it = poly.vertices_begin();
            it != poly.vertices_end();
            it++)
        {

            // inverse transform

             pt_inv =  it->point();
             pt_inv = pt_inv + view_point;
             v_pt_inv = pt_inv-CGAL::ORIGIN;

             double norm_pt_inv = sqrt(v_pt_inv*v_pt_inv);

             Vector_3 v = v_pt_inv*((2*R - norm_pt_inv)/(norm_pt_inv));

             pt_inv = CGAL::ORIGIN + v;

             pt_inv = pt_inv - view_point;
            points_visibles.push_back(glm::vec3(pt_inv.x(),pt_inv.y(),pt_inv.z() ) );
        }

    return points_visibles;
    }

    // HPR Ray optmal
     void WidgetOGL::HPR_ray_optmal(std::vector<Point_3> points ,glm::vec3 view_p ){

         // view
         Vector_3 view_point = Vector_3(view_p.x ,view_p.y,view_p.z);
         //get view opposite
         Vector_3 view_point_opp = view_point_opposite(points,view_point);

         // get R minimum
         double Rmin;
         Rmin = calcul_R_min(points, view_point);

         // get R optmal
         double R_opt;
         R_opt = calcul_R_optmal(points, view_point,view_point_opp, Rmin);


         HPR(points,-getCameraPosition(), R_opt);

     }
        Vector_3 WidgetOGL::view_point_opposite(std::vector<Point_3> points, Vector_3 view_point ){

         std::vector<Point_3>::const_iterator  it ;
         int size = points.size();
         Vector_3 v;
         Point_3 p;
         Vector_3 v_centre = Vector_3(0,0,0);

         // calculating center of mass
         for (it= points.begin();
              it!=points.end();
              it++){

             p = *it;
             v = p - CGAL::ORIGIN;
             v_centre = v_centre + v;
         }

         v_centre =  v_centre/size; //centre de masse


     Vector_3 view_point_to_CM =  ((v_centre)- ( view_point ) );
     Vector_3 view_point_to_view_point_opp = 2* view_point_to_CM;
     Vector_3 view_point_opp = view_point_to_view_point_opp + view_point;

     return view_point_opp;

     }
        double WidgetOGL::calcul_R_min(std::vector<Point_3> points,Vector_3 view_point){

            Vector_3 v_pt;
            Point_3 pt;

            double norm_pt;

        /* smallest and largest should point for the same element */
            std::vector<Point_3>::const_iterator it;
            it = points.begin();
            pt = *it;
            // camera to origin
            pt = pt + view_point;
            // transform point to vector (preparation for product)
            v_pt = pt-CGAL::ORIGIN;
            // norm -> distance of point to centre (view_point)
             norm_pt = sqrt(v_pt*v_pt);

            // they point for first element
            double largest = norm_pt;
            double smallest = norm_pt;

            for(it++;
                it != points.end();
                it++)
            {
               pt =  *it;

               // camera to origin
               pt = pt + view_point;
               // transform point to vector (preparation for product)
               v_pt = pt-CGAL::ORIGIN;
               // norm -> distance of point to centre (view_point)
                norm_pt = sqrt(v_pt*v_pt);

               if (largest < norm_pt)
                   largest=norm_pt;
               else if(smallest>norm_pt)
                   smallest=norm_pt;

            }

            double Rmin_view = largest;
            double Rmin_view_opposite = 2*sqrt(view_point*view_point) - smallest;


            return std::max(Rmin_view,Rmin_view_opposite);

        }

            double WidgetOGL::calcul_R_optmal(std::vector<Point_3> points, Vector_3 view_point,
                                      Vector_3 view_point_opp, double Rmin){
        int aux=0;
        double R_opt=0;


        // calculating of the first two values (R1 = 10*R0)
        double R_before = Rmin;
        double F_R_before = function_disjoint(points,view_point,view_point_opp , R_before);
            double R_current = 10*Rmin ;
            double F_R_current = function_disjoint(points,view_point,view_point_opp , R_current);

        //intitialize parameters
        double lamb =  100*Rmin;//Rmin;
            qDebug()<<"lambda:"<<lamb;
        double criteria_stop=0;
        double grad=0;
        double R_prox = 0;
        double F_R_prox = 0;
        // ----------------------------------------------------

       //Descendu du Gradient
        do
        {
            qDebug()<<"Before > R:"<<R_before<<"-- Function_disj:"<<F_R_before;
            qDebug()<<"Current> R"<<R_current<<"-- Function_disj:"<<F_R_current;

        //gradient
        grad = grad_function_disjoint(F_R_current,F_R_before, R_current, R_before);
        grad = atan(grad)*2/PI ;
            qDebug()<<"gradiante:"<<grad; //

        // R prox et F_prox
        R_prox = R_current + lamb*grad;
        F_R_prox = function_disjoint(points,view_point,view_point_opp , R_prox);
            qDebug()<<"    Prox > R:"<<R_prox<<"-- Function_disj:"<<F_R_prox<<"\n";

        // update
        R_before = R_current;
        F_R_before = F_R_current;
        R_current = R_prox;
        F_R_current = F_R_prox;

        // update criteria
        criteria_stop =  abs(R_current -R_before)/R_before;
            qDebug()<<"criteria_stop"<<criteria_stop;

        // number of iterations
        aux++;

        }while( criteria_stop>0.10 && aux<20);

            qDebug()<<"\nnumberof iterations:"<<aux;
            qDebug()<<"criteria_stop"<<criteria_stop;
            qDebug()<<R_prox ;

        R_prox = R_prox*sqrt(Rmin) ;

           qDebug()<<"Prox > R:"<<R_prox;

        return R_prox;

        }
            double WidgetOGL::function_disjoint(std::vector<Point_3> points,Vector_3 view_point,
                                     Vector_3 view_point_opp, double R){

            std::vector<Point_3> points_transf;
            Polyhedron poly;
            std::vector<Point_3> points_transf_opp;
            Polyhedron poly_opp;

            // view point
            points_transf = transformation(points,view_point,R);
               CGAL::convex_hull_3(points_transf.begin(),points_transf.end(),poly);
            //   qDebug()<<poly.size_of_vertices()<<" vertex in gull convex -- pv1";

           // view point opposite
            points_transf_opp = transformation(points,view_point_opp ,R);
              CGAL::convex_hull_3(points_transf_opp.begin(),points_transf_opp.end(),poly_opp);
            // qDebug()<<poly_opp.size_of_vertices()<<" vertex in gull convex -- pv2";

            // get intersection
             int intersec = intersection(poly,view_point,poly_opp,view_point_opp,R);

            // numbers of points disjoints
             double disjoint = (poly.size_of_vertices() + poly_opp.size_of_vertices() - 2*intersec);
             poly.erase_all();poly_opp.erase_all();

            return disjoint;
        }
                int WidgetOGL::intersection(Polyhedron poly1,Vector_3 view_point ,
                                Polyhedron poly2,Vector_3 view_point_opp, double R){

                   int aux =0;
                   Vertex_iterator it = poly1.vertices_begin();
                   Vertex_iterator it2  = poly2.vertices_begin();
                    Point_3 p1  = it->point();
                    Point_3 p2  = it2->point();

                   for(it;
                       it != poly1.vertices_end();
                       it++)
                   {
                       p1 = transform(it->point(),view_point,R);
                       while (it2!=poly2.vertices_end())
                       {
                           p2 = transform(it2->point(),view_point_opp,R);

                           if(compare_point(p1,p2))
                           {
                               aux++;
                               break; // a match
                           }

                        it2++;
                       }

                    it2 = poly2.vertices_begin();
                   }

                  return aux;
            }
                    bool  WidgetOGL::compare_point(Point_3 p1,Point_3 p2){
                      double  prec = 100000000;
                        return ( (int)(p1.x()*prec) == (int)(p2.x()*prec)&&
                                 (int)(p1.y()*prec) == (int)(p2.y()*prec)&&
                                 (int)(p1.z()*prec) == (int)(p2.z()*prec)
                                );

            }
            double WidgetOGL::grad_function_disjoint(double F_R_current, double F_R_before,
                                                 double R_current, double R_before ){
                  return ((F_R_current - F_R_before)/((R_current-R_before)));
              }

// --------------------------------------------------------------------------------------
    void WidgetOGL::read_pixels(){


        int aux=0;
        int size =1*width*height; // 1 chanel -> GL_R
        // read pixels
        GLfloat myData[size];
        glReadPixels(0,0,width,height,GL_RED,GL_FLOAT,myData);
        qDebug()<<myData[i+0]<<" "<<myData[i+1]<<" "<<myData[i+2]<<"\n";

        // file
          std::ofstream myfile;

           myfile.open("../test.txt");


           if (myfile.is_open())
             {
              //
                for(int j=height-1;j>=0;j--){
                    for(int i=0;i<width;i++){
                        myfile<<myData[(j)*(width) + i]<<" ";
                    }
                    myfile<<"\n";
                }
                /*/

               for(int j=0;j<height;j++){
                   for(int i=0;i<width;i++){
                       myfile<<myData[aux]<<" ";
                        aux++;
                   }
                   myfile<<"\n";
               }
                //*/

               myfile.close();
             }
             else cout << "Unable to open file";


    }


// cicle of life
void WidgetOGL::initializeGL()
{

   initializeOpenGLFunctions();

    // intitializations --------------------
         glClearColor(1.0,1.0,1.0,1.0);
         glShadeModel(GL_SMOOTH);
         glEnable(GL_DEPTH_TEST);
         glDisable(GL_CULL_FACE);
         glFrontFace(GL_CCW);
        // enable Alpha in shader
         glDisable(GL_TEXTURE_2D);
         glDisable(GL_BLEND);
    //--------------------------------------

    // VAO
    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Shader -------------------------------------------
        // create
        GLuint vs, fs;
        if ((vs = create_shader("../loadShaders/shaders/obj.vert", GL_VERTEX_SHADER))   == 0) exit(0);
        if ((fs = create_shader("../loadShaders/shaders/obj.frag", GL_FRAGMENT_SHADER)) == 0) exit(0);
        // link
        program =  link(program,vs,fs);
        // Free
        glDetachShader(program, vs); glDeleteShader(vs);
        glDetachShader(program, fs);  glDeleteShader(fs);
    //------------------------------------------------------

    // Get Data
       //  std::vector<glm::vec3> points = get_data_xyz("../data/dataXYZ/test2.xyz");
          std::vector<glm::vec3> points = get_data_xyz("../data/dataXYZ/frog.xyz");
        //  std::vector<glm::vec3> points = get_data_xyz("../data/dataXYZ/kinect_norm.xyz");
       // std::vector<glm::vec3> points = get_data_xyz("../data/dataXYZ/data-leo/Prueba_10.pts");
         numbers_points =  points.size();


    // location of variables
        MVPID = glGetUniformLocation(program, "MVP");
        MVID = glGetUniformLocation(program, "MV");
}

void WidgetOGL::paintGL()
{
        glViewport(0, 0, this->width, this->height);

        // HPR
        if (HPR_toogle)
        {      
            HPR(m_points,-getCameraPosition(), R);
         }

    //VBO // --------------------------------------------------
        // data vertex in VBO (copy CPU to buffer in GPU)
        glGenBuffers(1, &vbo_attributes_vertex);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_attributes_vertex);
        glBufferData(GL_ARRAY_BUFFER, points.size()*3*sizeof(float), points.data() ,GL_STATIC_DRAW);
        // disable this VBO (it's not delete)
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // data color in VBO (copy CPU to buffer in GPU)
        glGenBuffers(1, &vbo_attributes_color);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_attributes_color);
        glBufferData(GL_ARRAY_BUFFER, color.size()*sizeof(glm::vec3), &color[0] ,GL_STATIC_DRAW);
        // disable this VBO (it's not delete)
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    // --------------------------------------------------------------

    // parameters screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(1.0, 1.0, 1.0, 1.0);

    // Matrix
        //Définition des matrices
        glm::mat4 MVP = Projection * View * Model;
        glm::mat4 MV =  View * Model;
        //send matrix
        glUniformMatrix4fv(MVPID, 1, GL_FALSE, &MVP[0][0]);
        glUniformMatrix4fv(MVID, 1, GL_FALSE, &MV[0][0]);

    // program indication
    glUseProgram(program);

    // enable buffer VBO for to give the attributes vertex
    glBindBuffer(GL_ARRAY_BUFFER, vbo_attributes_vertex);
        glEnableVertexAttribArray(0);
        // ({location},{number of elements},{type}, {its false}, {stride}, {offse })
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

   //enable buffer VBO for give the attributes color
    glBindBuffer(GL_ARRAY_BUFFER, vbo_attributes_color);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE, 0, (GLvoid*)0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glPointSize(2);
    // draw in screen
    glDrawArrays(GL_POINTS, 0, numbers_points);



  //  if(toogle_read_pixels){read_pixels();}

    // clean
     glDeleteBuffers(1, &vbo_attributes_vertex);
     // clean
      glDeleteBuffers(1, &vbo_attributes_color);


 }

void WidgetOGL::resizeGL(int width, int height)
{
    this->width = width;
    this->height = height;
    glViewport(0, 0, width, height);
}
