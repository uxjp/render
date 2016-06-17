
#include <QApplication>
#include <QGLFormat>

#include <QSurfaceFormat>


#include "WidgetOGL.h"


int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Specify an OpenGL 3.3 format using the Core profile.
   QGLFormat glFormat;
        glFormat.setVersion( 3, 3 );
        glFormat.setProfile( QGLFormat::CoreProfile ); // Requires >=Qt-4.8.0
        glFormat.setSampleBuffers( true );

    QSurfaceFormat format;
       format.setDepthBufferSize(24);
       format.setStencilBufferSize(8);
       format.setVersion(3, 2);
       format.setProfile(QSurfaceFormat::CoreProfile);
       QSurfaceFormat::setDefaultFormat(format);

    WidgetOGL exemple(glFormat);
        exemple.setWindowTitle(QObject::tr("Exemple with shader"));
        exemple.resize(256,256);
        exemple.show();

    return app.exec();
}
