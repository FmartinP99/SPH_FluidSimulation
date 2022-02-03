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


float red=1.0f, blue=1.0f, green=1.0f;
float angle = 0.0f;
float lx=0.5f,lz=-3.0f, ly=0.0f;
float x=0.0f,z=2.0f, y=1.0f;
float deltaAngle = 0.0f;
float deltaMove = 0;
int xOrigin = -1;
Container* container;


void processNormalKeys(unsigned char key, int x, int y) {

    float fraction = 0.1f;

    switch (key) {
        case 87:
            x += lx * fraction;
            z += lz * fraction;
            break;
        case 119:
            x += lx * fraction;
            z += lz * fraction;
            break;
        case 83:
            x -= lx * fraction;
            z -= lz * fraction;
            break;
        case 115:
            x -= lx * fraction;
            z -= lz * fraction;
            break;
    }
}

void processSpecialKeys(int key, int x, int y) {

    float fraction = 0.1f;

    std::cout<<key;
    switch (key) {
        case GLUT_KEY_LEFT :
            angle -= 0.01f;
            lx = sin(angle);
            lz = -cos(angle);
            break;
        case GLUT_KEY_RIGHT :
            angle += 0.01f;
            lx = sin(angle);
            lz = -cos(angle);
            break;
            /*
            case GLUT_KEY_UP :
                v_angle += 0.01f;
                lx = cos(v_angle);
                ly = -sin(v_angle);
                break;
            case GLUT_KEY_DOWN :
                break;
                */
    }
}

void pressKey(int key, int xx, int yy) {

    switch (key) {
        case GLUT_KEY_LEFT : deltaAngle = -0.01f; break;
        case GLUT_KEY_RIGHT : deltaAngle = 0.01f; break;
        case 87 : deltaMove = 0.5f; break;
        case 119 : deltaMove = 0.5f; break;
        case 83 : deltaMove = -0.5f; break;
        case 115 : deltaMove = -0.5f; break;
    }
}

void releaseKey(int key, int x, int y) {

    switch (key) {
        case GLUT_KEY_LEFT :
        case GLUT_KEY_RIGHT : deltaAngle = 0.0f;break;
        case 87 :
        case 119 :
        case 83 : deltaMove = 0;break;
        case 115 : deltaMove = 0;break;
    }
}


void mouseButton(int button, int state, int x, int y) {

    // only start motion if the left button is pressed
    if (button == GLUT_LEFT_BUTTON) {

        // when the button is released
        if (state == GLUT_UP) {
            angle -= deltaAngle;
            xOrigin = -1;
        }
        else {// state = GLUT_DOWN
            xOrigin = x;
        }
    }
}

void mouseMove(int x, int y) {
    if (xOrigin >= 0) {
        float deltaAngleX = (x - xOrigin) * 0.001f;

        lx = sin(angle  - deltaAngleX);
        lz = -cos(angle - deltaAngleX);
    }
}



void changeSize(int w, int h) {

// Prevent a divide by zero, when window is too short
// (you cant make a window of zero width).
    if (h == 0)
        h = 1;

    float ratio =  w * 1.0 / h;

    // Use the Projection Matrix
    glMatrixMode(GL_PROJECTION);

    // Reset Matrix
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);

    // Set the correct perspective.
    gluPerspective(45,ratio,1,100);

    // Get Back to the Modelview
    glMatrixMode(GL_MODELVIEW);
}

void drawParticles(){
    glColor3f(1.0f, 1.0f, 1.0f);
    glTranslatef(0.0f ,0.0f, 0.0f);
    glutSolidSphere(0.01f,10,20);
}

void drawSnowMan() {

    glColor3f(1.0f, 1.0f, 1.0f);

    // Draw Body
    glTranslatef(0.0f ,0.75f, 0.0f);
    glutSolidSphere(0.75f,20,20);

    // Draw Head
    glTranslatef(0.0f, 1.0f, 0.0f);
    glutSolidSphere(0.25f,20,20);

    // Draw Eyes
    glPushMatrix();
    glColor3f(0.0f,0.0f,0.0f);
    glTranslatef(0.05f, 0.10f, 0.18f);
    glutSolidSphere(0.05f,10,10);
    glTranslatef(-0.1f, 0.0f, 0.0f);
    glutSolidSphere(0.05f,10,10);
    glPopMatrix();

    // Draw Nose
    glColor3f(1.0f, 0.5f , 0.5f);
    glutSolidCone(0.08f,0.5f,10,2);
}


void computePos(float deltaMove) {

    x += deltaMove * lx * 0.1f;
    z += deltaMove * lz * 0.1f;
}

void computeDir(float deltaAngle) {

    angle += deltaAngle;
    lx = sin(angle);
    lz = -cos(angle);
}


void renderScene() {

    if (deltaMove) computePos(deltaMove);
    if (deltaAngle) computeDir(deltaAngle);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(	x, 1.0f, z,
                  x+lx, 1.0f,  z+lz,
                  0.0f, 1.0f,  0.0f);

    glRotatef(angle, 0.0f, 1.0f, 0.0f);

    glColor3f(.9f, 0.9f, 0.9f);
    /*
    //EZ RAJZOLJA KI AZ ALAPOT
    glBegin(GL_QUADS);
    glVertex3f(-100.0f, 0.0f, -100.0f);
    glVertex3f(-100.0f, 0.0f,  100.0f);
    glVertex3f( 100.0f, 0.0f,  100.0f);
    glVertex3f( 100.0f, 0.0f, -100.0f);
    glEnd();
    */
        for(int i=0; i < container->get_particles_number(); i++) {
            glPushMatrix();
            //drawSnowMan();
            auto particle = container->get_particles()[i];
            glTranslatef(particle.get_pos_x(), particle.get_pos_y(), particle.get_pos_z());
            drawParticles();
            glPopMatrix();
        }

    glutSwapBuffers();
}

int main(int argc, char **argv) {

    container = new Container(1, 1, 1, 500);

    // init GLUT and create window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(1920,1080);
    glutCreateWindow("SPH_FluidSimulation");

    // register callbacks
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutIdleFunc(renderScene);

    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);
    glutIgnoreKeyRepeat(1);
    glutSpecialUpFunc(releaseKey);

    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);

    // OpenGL init
    glEnable(GL_DEPTH_TEST);

    // enter GLUT event processing loop
    glutMainLoop();

    return 1;
}