/*
typedef struct Kernel_Particle_Distance{
    int index;
    float distance;

}Kernel_Particle_Distance;

typedef struct Kernel_Neighboring_Particles{
    Kernel_Particle_Distance* kernelParticleDistance;
    int items;
}Kernel_Neighboring_Particles;

typedef struct Kernel_Neighboring_Grids{
    struct Kernel_Neighboring_Particles* kernelNeighboringParticles;
    int items;
}Kernel_Neighboring_Grids;



typedef struct Kernel_SPHParticle{
    float weight;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float ax;
    float ay;
    float az;
    float radius;
    float mass;
    float pressure;
    float density;
    float fx;
    float fy;
    float fz;
} Kernel_SPHParticle;

/*
typedef struct Container{

    float damping_coeff;
    unsigned int gas_constant;
    float kernel_smoother_length;
    float constant_pressure;
    float U;

    unsigned int current_particles;
    unsigned int current_iteration ;
    unsigned int new_particles_every_iteration;

    float normalization_density;
    float norm_pressure;
    float norm_visc;

    float WIDTH;
    float HEIGHT;
    float DEPTH;
    float gravity;
    float density;
    unsigned int particles_number;
    float particle_radius;

    Kernel_SPHParticle* particles;

    unsigned int max_number_of_grids;
    unsigned int max_column_per_row_grid;
    unsigned int max_row_per_layer_grid;
    unsigned int max_layer_per_cube_grid;

}Container;
*/
unsigned int calculate_grid(Kernel_SPHParticle* elem, float kernel_smoother_length, unsigned int max_column_per_row_grid, unsigned int max_row_per_layer_grid, unsigned int max_layer_per_cube_grid){

    unsigned int column = elem->px / kernel_smoother_length + 1;
    unsigned int row = elem->py / kernel_smoother_length + 1;
    unsigned int layer = elem->pz / kernel_smoother_length + 1;


    if(column > max_column_per_row_grid) column = max_column_per_row_grid;
    if(row > max_row_per_layer_grid) row = max_row_per_layer_grid;
    if(layer > max_layer_per_cube_grid) layer = max_layer_per_cube_grid;

    unsigned int return_value = (row - 1) * max_column_per_row_grid + column + (max_row_per_layer_grid * max_column_per_row_grid * (layer - 1));


    return return_value;
}

float kernel_pow(float ertek, unsigned int hatvany){
    float vegeredmeny = 1;
    for (int i = 0; i < hatvany; ++i) vegeredmeny *= ertek;
    return vegeredmeny;

}

static float calculate_absolute_distance_power2(__global Kernel_SPHParticle* elem, __global Kernel_SPHParticle* elem2){
return (elem->px - elem2->px)*(elem->px - elem2->px) +
(elem->py - elem2->py)*(elem->py - elem2->py) +
(elem->pz - elem2->pz)*(elem->pz - elem2->pz);

}

void check_boundaries(float WIDTH, float HEIGHT, float DEPTH, float damping_coeff,  __global Kernel_SPHParticle* elem){
if ((*elem).py + (*elem).radius  >= HEIGHT) {
(*elem).vy = (*elem).vy  * damping_coeff * -1;
(*elem).py =HEIGHT - (*elem).radius;
}else if ((*elem).py - (*elem).radius <= 0){
(*elem).vy =(*elem).vy  * damping_coeff * -1;
(*elem).py = 0 + (*elem).radius ;
}

if ((*elem).px + (*elem).radius  >= WIDTH) {
(*elem).vx = (*elem).vx  * damping_coeff * -1;
(*elem).px = WIDTH - (*elem).radius;
}else if ((*elem).px - (*elem).radius  <= 0){
(*elem).vx =(*elem).vx  * damping_coeff * -1;
(*elem).px = 0 + (*elem).radius ;
}

if ((*elem).pz + (*elem).radius  >= DEPTH) {
(*elem).vz = (*elem).vz  * damping_coeff * -1;
(*elem).pz = DEPTH - (*elem).radius;
}else if ((*elem).pz - (*elem).radius  <= 0){
(*elem).vz = (*elem).vz  * damping_coeff * -1;
(*elem).pz = 0 + (*elem).radius;
}

}


void __calculate_pressure_density(__global Kernel_SPHParticle* sphParticle, __global float* enviromentArray){

float kernel_smoother_length = enviromentArray[6];
float current_particles = enviromentArray[3];
float norm_density = enviromentArray[7];
float density = enviromentArray[14];
float gas_constant = enviromentArray[1];


    for (int i = 0; i < current_particles; i++){
        __global Kernel_SPHParticle *elem = &sphParticle[i];
        for (int j = 0; j < current_particles; j++) {
            if( i != j){
                __global Kernel_SPHParticle* elem2 = &sphParticle[j];
                float abs_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
                if (abs_distance_pow2 < kernel_smoother_length){
        elem->density += norm_density * kernel_pow((kernel_pow(kernel_smoother_length, 2) - kernel_pow( abs_distance_pow2,
        2)), 3) * elem2->mass;
        }
                }
            }
        if (elem->density < density) elem->density = density;
        elem->pressure = gas_constant * (elem->density - density);
    }
}

void __calculate_pressure_velocity(__global Kernel_SPHParticle* sphParticle, __global float* enviromentArray){
float kernel_smoother_length = enviromentArray[6];
float current_particles = enviromentArray[3];
float norm_pressure = enviromentArray[8];
float norm_visc = enviromentArray[9];
float U = enviromentArray[2];

for(int i = 0; i < current_particles; i++){

    __global Kernel_SPHParticle* elem = &sphParticle[i];

    for(int j = 0; j < current_particles; j++) {

        if (i != j) {
            __global Kernel_SPHParticle* elem2 = &sphParticle[j];
            float neighboring_particle_distance = calculate_absolute_distance_power2(elem, elem2);

            if (elem->pressure != 0 && elem2->pressure != 0 && neighboring_particle_distance < kernel_smoother_length) {

                float pressure_weight =
                        norm_pressure * kernel_pow((kernel_smoother_length - neighboring_particle_distance), 2);
                float pforce =
                        (elem->pressure + elem2->pressure) * elem2->mass / (2 * elem2->density) *
                        pressure_weight;
                float distance_x = elem2->px - elem->px;
                float distance_y = elem2->py - elem->py;
                float distance_z = elem2->pz - elem->pz;
                elem->fx += distance_x * pforce;
                elem->fy += distance_y * pforce;
                elem->fz += distance_z * pforce;

            }

            if (elem2->density != 0 && neighboring_particle_distance < kernel_smoother_length) {
                float viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
                float vel_difference_x = elem2->vx - elem->vx;
                float vel_difference_y = elem2->vy - elem->vy;
                float vel_difference_z = elem2->vz - elem->vz;
                float v_force_x =
                        vel_difference_x * (elem2->mass / elem2->density * viscosity_weight);
                float v_force_y =
                        vel_difference_y * (elem2->mass / elem2->density * viscosity_weight);
                float v_force_z =
                        vel_difference_z * (elem2->mass / elem2->density * viscosity_weight);
                elem->fx += U * (v_force_x);
                elem->fy += U * (v_force_y);
                elem->fz += U * (v_force_z);
                }
            }
        }
    }
}


void __move_particles(__global Kernel_SPHParticle* sphParticle, __global float* dtime, __global float* enviromentArray){
     float current_particles = enviromentArray[3];
     float gravity = enviromentArray[13];
     float WIDTH = enviromentArray[10];
     float HEIGHT = enviromentArray[11];
     float DEPTH = enviromentArray[12];
     float damping_coeff = enviromentArray[0];



     for(int i = 0; i < current_particles; i++){

      __global Kernel_SPHParticle* elem = &sphParticle[i];

     elem->fy +=  gravity;

     if (elem->density != 0) {

         elem->vx = elem->vx + *dtime * elem->fx / elem->density;
         elem->vy = elem->vy + *dtime * elem->fy / elem->density;
         elem->vz = elem->vz + *dtime * elem->fz / elem->density;
     }
     elem->px = elem->px + *dtime * elem->vx;
     elem->py = elem->py + *dtime * elem->vy;
     elem->pz = elem->pz + *dtime * elem->vz;

     check_boundaries(WIDTH, HEIGHT, DEPTH, damping_coeff, elem);

     }
}

__kernel void calculate_physics(__global Kernel_SPHParticle* sphParticle, __global float* dtime, __global float* enviromentArray){

__calculate_pressure_density(sphParticle, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
__calculate_pressure_velocity(sphParticle, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
__move_particles(sphParticle, dtime, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
printf("vegigment a kernel!\n");

}

*/