#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <fstream>
#include <windows.h>
#include <GL\glut.h>
#endif
#include "cstdlib"
#include "cstdio"
#include <iostream>
#include "Classes\OpenCLContainer.h"
#include "Classes\Timer.cpp"
#include <ctime>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <CL/cl.hpp>

//COMMAND LINE ARGUMENTS
int PARTICLES_NUMBER;
float WIDTH = 30, HEIGHT = 30, DEPTH = 5;
bool IS_NAIVE = false;
bool IS_GPU = false;
std::string KERNEL_PATH = "Kernel/kernel.cl";
int FPS_NUMBER = 60;
//


float gravity = -10, density = 1;
float radius= 0.40;

int xOrigin = -1;
int yOrigin = -1;
float templz = 0, templx = 0, temply = 0;
float angle = 0.0f;
float x,z, y;
float lx, ly, lz;
float deltaAngleX, deltaAngleY;

const double dtime = 30.0 / 1000.0;

double enviroment_array[21];

Container* container;
OpenCLContainer* openClContainer;

int initial_time = time(NULL), final_time, frame_count = 0;
std::vector<long long>* timer_vec_algo = new std::vector<long long>();
std::vector<long long>* timer_vec_drawing = new std::vector<long long>();
long long global_tick_count = 0;
long long START_TIME = duration_cast<seconds>(system_clock::now().time_since_epoch()).count();

void passToGPU(OpenCLContainer &helper, Container* container) {
    int err = CL_SUCCESS;
    auto asd = CL_DEVICE_EXTENSIONS;
    auto program = helper.getProgram();
    auto context = helper.getContext();
    auto devices = helper.getDevices();
    auto kernelUpdate = helper.updateKernel;
    auto dtimeBuffer = helper.dtimeBuffer;
    auto particleBuffer = helper.particleBuffer;
    auto enviromentarrayBuffer = helper.enviromentarrayBuffer;



    container->set_current_iteration(container->get_current_iteration() + 1);
    if(container->get_current_particles() < container->get_particles_number() && container->get_current_iteration() % container->get_new_particles_every_iteration() == 0)
        container->fill_container_gradually(radius);

    enviroment_array[3] = container->get_current_particles();

    SPHParticle* particles = container->get_particles();


    SPHParticle host_buffer[PARTICLES_NUMBER];
    for (int i=0; i<container->get_current_particles(); i++){
        host_buffer[i] = particles[i];
    }


    try {
        cl::Event event;
        cl::CommandQueue myQueue(context, devices[0], 0, &err);

        //err = myQueue.enqueueWriteBuffer(testBuffer, true, 0,  10 * sizeof(double) , doublearray);
        err = myQueue.enqueueWriteBuffer(particleBuffer, true, 0,  PARTICLES_NUMBER * sizeof(SPHParticle), particles);
        err = myQueue.enqueueWriteBuffer(dtimeBuffer, true, 0, sizeof(double), &dtime);
        //err = myQueue.enqueueWriteBuffer(containerBuffer, true, 0, sizeof(Container), container);
        err = myQueue.enqueueWriteBuffer(enviromentarrayBuffer, true, 0, sizeof(double)*21, enviroment_array);

        err = myQueue.enqueueNDRangeKernel(kernelUpdate,
                                           cl::NullRange,
                                           cl::NDRange(enviroment_array[3], 1),
                                           cl::NullRange,
                                           NULL,
                                           &event);
        event.wait();
        err = myQueue.enqueueReadBuffer(particleBuffer, true, 0, PARTICLES_NUMBER * sizeof(SPHParticle), particles);


    }
    catch (cl::BuildError &e) {
        std::cout << "build error" << std::endl;
        std::cout << e.getBuildLog().front().second << std::endl;
    }
    catch (cl::Error &e) {
        std::cout << e.what() << std::endl;
        std::cout << e.err() << std::endl;
    }
}

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

void write_to_csv(bool failed){
    std::ofstream file;
    file.open("results.csv", std::ios_base::out | std::ios_base::app);

    int failed_num = 0;
    if (failed) failed_num = 1;
    for(int i = 0; i < timer_vec_drawing->size(); i++){
        file << PARTICLES_NUMBER << "," << WIDTH << "," << HEIGHT << "," << DEPTH << "," << IS_NAIVE << "," << IS_GPU << "," << timer_vec_algo->at(i) << "," << timer_vec_drawing->at(i) << "," <<failed_num << std::endl;
    }
    file.close();
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
    glVertex3d(0, HEIGHT, 0);
    glVertex3d(WIDTH, 0, 0);
    glVertex3d(WIDTH, HEIGHT, 0);
    glVertex3d(0, 0, DEPTH);
    glVertex3d(0, HEIGHT, DEPTH);
    glVertex3d(WIDTH, 0, DEPTH);
    glVertex3d(WIDTH, HEIGHT, DEPTH);

    glVertex3d(WIDTH, HEIGHT, DEPTH);
    glVertex3d(0, HEIGHT, DEPTH);
    glVertex3d(WIDTH, HEIGHT, DEPTH);
    glVertex3d(WIDTH, HEIGHT, 0);
    glVertex3d(0, HEIGHT, DEPTH);
    glVertex3d(0, HEIGHT, 0);
    glVertex3d(0, HEIGHT, 0);
    glVertex3d(WIDTH, HEIGHT, 0);

    glVertex3d(WIDTH, 0, DEPTH);
    glVertex3d(0, 0, DEPTH);
    glVertex3d(WIDTH, 0, DEPTH);
    glVertex3d(WIDTH, 0, 0);
    glVertex3d(0, 0, DEPTH);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 0, 0);
    glVertex3d(WIDTH, 0, 0);
    glEnd();

// ---------------------------------------------------------------------------------------------------------------------------------------------------------

    //TIME_DELTA = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count() - START_TIME;
   // std::cout<<TIME_DELTA<<std::endl;
    //if (TIME_DELTA <= 100) container->calculate_physics(TIME_DELTA);
    //START_TIME = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    {
        Timer timer(timer_vec_algo);
        if (IS_GPU)
                passToGPU(*openClContainer, container);
        else if (IS_NAIVE)
                container->calculate_physicsv2(dtime);
        else
                container->calculate_physics(dtime);
    }
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
    {
        Timer timer(timer_vec_drawing);
        for (int i = 0; i < container->get_current_particles(); i++) {
            glPushMatrix();
            auto particle = container->get_particles()[i];
            glTranslatef(particle.get_px(), particle.get_py(), particle.get_pz());
            //     std::cout<<particle.get_vy()<<std::endl;
            drawParticles(i, particle);
            glPopMatrix();
        }
    }
    //---------------------------------------------------------------------------------------------------------------------------------------------------
    glutSwapBuffers();

    frame_count++;
    final_time = time(NULL);

    if(final_time - initial_time >= 1){
        std::cout<<"FPS: "<<frame_count <<std::endl;
        frame_count = 0;
        initial_time = final_time;
    }
    global_tick_count += 1;
    long long end_time = duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
    /*
    if (global_tick_count > 1500 || end_time - START_TIME > 1200) { //maximum tickig megy vagy maximum 20 percig
        if (end_time - START_TIME > 1200) write_to_csv(true);
        else write_to_csv(false);
        glutDestroyWindow(glutGetWindow());
    }
     */
}

void set_variables(int argc, char** argv){
    if (argc > 1)
        PARTICLES_NUMBER = atoi(argv[1]);
    if (argc > 2)
        WIDTH = atof(argv[2]);
    if (argc > 3)
        HEIGHT = atof(argv[3]);
    if (argc > 4)
        DEPTH = atof(argv[4]);
    if (argc > 5) {
        auto _var = atoi(argv[5]);
        if (_var > 0)
        FPS_NUMBER = _var;
    }
    if (argc > 6) {
        auto _var = atof(argv[6]);
        if (_var == 1)
            IS_NAIVE = true;
    }
    if (argc > 7) {
        auto _var = atof(argv[7]);
        if (_var == 1)
            IS_GPU = true;
    }
    if (argc > 8)
        KERNEL_PATH = argv[8];

}

int main(int argc, char **argv) {

    set_variables(argc, argv);
    std::cout<<argv[0]<<std::endl;
    std::cout<<argv[1]<<std::endl;
    if (IS_GPU)
        std::cout<<"GPU ENABLED!"<<std::endl;
    else {
        std::cout << "CPU ENABLED!" << std::endl;
        if (IS_NAIVE)
            std::cout<< "NAIVE ALGORITHM IS ENABLED!" <<std::endl;
        else
            std::cout<< "GRID ALGORITHM IS ENABLED!" <<std::endl;
    }



    container = new Container(WIDTH, HEIGHT, DEPTH, gravity, density, PARTICLES_NUMBER, radius);
    if (IS_GPU) {
        enviroment_array[0] = container->get_damping_coeff();
        enviroment_array[1] = container->get_gas_constant();
        enviroment_array[2] = container->get_u();
        enviroment_array[3] = container->get_current_particles();
        enviroment_array[4] = container->get_current_iteration();
        enviroment_array[5] = container->get_new_particles_every_iteration();
        enviroment_array[6] = container->get_kernel_smoother_length();
        enviroment_array[7] = container->get_normalization_density();
        enviroment_array[8] = container->get_norm_pressure();
        enviroment_array[9] = container->get_norm_visc();
        enviroment_array[10] = container->get_width();
        enviroment_array[11] = container->get_height();
        enviroment_array[12] = container->get_depth();
        enviroment_array[13] = container->get_gravity();
        enviroment_array[14] = container->get_density();
        enviroment_array[15] = container->get_particles_number();
        enviroment_array[16] = container->get_particle_radius();
        enviroment_array[17] = container->get_max_number_of_grids();
        enviroment_array[18] = container->get_max_column_per_row_grid();
        enviroment_array[19] = container->get_max_row_per_layer_grid();
        enviroment_array[20] = container->get_max_layer_per_cube_grid();


        openClContainer = new OpenCLContainer(KERNEL_PATH,
                                              *container, PARTICLES_NUMBER);
    }
    x = WIDTH / 2;
    y = HEIGHT / 2;
    z = 1 + WIDTH / 2 + HEIGHT / 2 + DEPTH + 5 ;

    lx=( sin(angle -(x - xOrigin) * 0.001f));
    lz=(-cos(angle -(x - xOrigin) * 0.001f));
    ly=(-sin(angle -(y - yOrigin) * 0.001f));

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
    glDepthFunc(GL_LESS);

    glutMainLoop();




}


