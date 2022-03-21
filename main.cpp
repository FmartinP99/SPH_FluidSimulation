#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "cstdlib"
#include "cmath"
#include "cstdio"
#include <iostream>
#include "Classes\Container.cpp"
#include <ctime>
#include <chrono>

#define FPS_NUMBER 60
#define PARTICLES_NUMBER 1600 //37 FPS WITH 800 PARTICLE  at 10-10
double width = 50, height = 30, depth = 1, gravity = -10, density = 1;
double radius= 0.5;

using namespace std::chrono;
auto START_TIME = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
auto TIME_DELTA = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
float COLOR_ARRAY[3][PARTICLES_NUMBER];


int xOrigin = -1;
int yOrigin = -1;
double templz = 0, templx = 0, temply = 0;
double angle = 0.0f;
double x,z, y;
double lx, ly, lz;
double deltaAngleX, deltaAngleY;
Container* container;

int initial_time = time(NULL), final_time, frame_count = 0;




void processSpecialKeys(int key, int mouse_pos_x, int mouse_pos_Y) {
}

void processNormalKeys(unsigned char key, int mouse_pos_x, int mouse_pos_Y) {
    float fraction = 1.0f;
    switch (key) {
        case 119:
        case 87:
            z -=  fraction;
            break;
        case 83:
        case 115:
            z +=  fraction;
            break;
        case 97:
        case 65:
            x -=  fraction;
            break;
        case 68:
        case 100:
            x += fraction;
            break;
        case 67:
        case 99:
            y -= fraction;
            break;
        case 32:
            y += fraction;
            break;
        default:
            break;

    }
}



void mouseButton(int button, int state, int mouse_pos_x, int mouse_pos_y) {
    // only start motion if the left button is pressed
    if (button == GLUT_LEFT_BUTTON) {

        // when the button is released
        if (state == GLUT_UP) {

        }
        else {// state = GLUT_DOWN
            xOrigin = mouse_pos_x; // lenyomtam
            yOrigin = mouse_pos_y;
        }
    }
}



void mouseMove(int mouse_pos_x, int mouse_pos_y) {


    if (xOrigin >= 0) {
        deltaAngleX = (mouse_pos_x - xOrigin) * 0.001f; // kiszámolt
        lx = sin(angle  - deltaAngleX);
        lz = -cos(angle - deltaAngleX);
    }

    if(yOrigin >= 0){
        deltaAngleY = (mouse_pos_y - yOrigin) * 0.001f;
        ly = -sin(angle - deltaAngleY);
    }


}



void changeSize(int w, int h) {

    if (h == 0)
        h = 1;

    float ratio =  w * 1.0 / h;

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    glViewport(0, 0, w, h);

    gluPerspective(45,ratio,1,1000);

    glMatrixMode(GL_MODELVIEW);
}


void drawParticles(int index, SPHParticle& particle){
    glColor3f(0.2, 0.6, 0.8);
   // glColor3f(std::min(std::abs(particle.get_vx()), 1.0), std::min(std::abs(particle.get_vy()), 1.0), std::min(std::abs(particle.get_vz()), 1.0));
    glTranslatef(0.0f ,0.0f, 0.0f);
    glutSolidSphere(radius,10,10);
}

void idle(int){
    glutPostRedisplay();
    glutTimerFunc(1000/FPS_NUMBER, idle, 0);
}

void renderScene() {


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(	x, y, z,  // eye position
                  x+lx+templx, y+ly+temply,  z+lz+templz, // hova nézek
                  0.0f, 1.0f,  0.0f); // up direction


    glRotatef(angle, 0.0f, 1.0f, 0.0f);

    glColor3f(0.0f, 0.1f, 0.1f);

    //EZ RAJZOLJA KI AZ ALAPOT
    glBegin(GL_QUADS);
    glVertex3f(-100.0f, -1.0f, -100.0f);
    glVertex3f(-100.0f, -1.0f,  100.0f);
    glVertex3f( 100.0f, -1.0f,  100.0f);
    glVertex3f( 100.0f, -1.0f, -100.0f);
    glEnd();

    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(0, height, 0);
    glVertex3d(width, 0, 0);
    glVertex3d(width, height, 0);
    glVertex3d(0, 0, depth);
    glVertex3d(0, height, depth);
    glVertex3d(width, 0, depth);
    glVertex3d(width, height, depth);

    glVertex3d(width, height, depth);
    glVertex3d(0, height, depth);
    glVertex3d(width, height, depth);
    glVertex3d(width, height, 0);
    glVertex3d(0, height, depth);
    glVertex3d(0, height, 0);
    glVertex3d(0, height, 0);
    glVertex3d(width, height, 0);

    glVertex3d(width, 0, depth);
    glVertex3d(0, 0, depth);
    glVertex3d(width, 0, depth);
    glVertex3d(width, 0, 0);
    glVertex3d(0, 0, depth);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 0, 0);
    glVertex3d(width, 0, 0);
    glEnd();

// ---------------------------------------------------------------------------------------------------------------------------------------------------------

    //TIME_DELTA = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - START_TIME;
   // std::cout<<TIME_DELTA<<std::endl;
    //if (TIME_DELTA <= 100) container->calculate_physics(TIME_DELTA);
    //START_TIME = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    container->calculate_physics(30.0/1000.0);
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

    for(int i=0; i < container->get_current_particles(); i++) {
        glPushMatrix();
        auto particle = container->get_particles()[i];
        glTranslatef(particle.get_px(), particle.get_py(), particle.get_pz());
   //     std::cout<<particle.get_vy()<<std::endl;
        drawParticles(i, particle);
        glPopMatrix();
    }

    glutSwapBuffers();

    frame_count++;
    final_time = time(NULL);

    if(final_time - initial_time >= 1){
        std::cout<<"FPS: "<<frame_count / (final_time - initial_time)<<std::endl;
        frame_count = 0;
        initial_time = final_time;
    }

}


int main(int argc, char **argv) {

    container = new Container(width, height, depth, gravity, density, PARTICLES_NUMBER, radius);
    x = width/2;
    y = height/2;
    z = 1 + width/2 + height/2 + depth + 5 ;



    lx=( sin(angle -(x - xOrigin) * 0.001f));
    lz=(-cos(angle -(x - xOrigin) * 0.001f));
    ly=(-sin(angle -(y - yOrigin) * 0.001f));

    for (int i = 0; i < PARTICLES_NUMBER * 3; i++){
        float random_number = ((float)rand()/(RAND_MAX));
        COLOR_ARRAY[0][i] = random_number;
    }

    // init GLUT and create window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(1280,720);
    glutCreateWindow("SPH_FluidSimulation");

    // register callbacks
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);

    glutTimerFunc(1000/FPS_NUMBER, idle, 0);

    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);

    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);


    // OpenGL init
    glEnable(GL_DEPTH_TEST);


    START_TIME = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    glutMainLoop();


    return 1;
}
