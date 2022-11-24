#define EndianSwap(n) (rotate(n & 0x00FF00FF, 24U)|(rotate(n, 8U) & 0x00FF00FF)

typedef struct SPHParticle{
    double weight;
    double px;
    double py;
    double pz;
    double vx;
    double vy;
    double vz;
    double radius;
    double mass;
    double pressure;
    double density;
    double fx;
    double fy;
    double fz;

} SPHParticle;

double kernel_pow(double ertek, unsigned int hatvany){
    double vegeredmeny = 1;
    for (int i = 0; i < hatvany; ++i) vegeredmeny *= ertek;
    return vegeredmeny;

}

static double calculate_absolute_distance_power2(__global SPHParticle* elem, __global SPHParticle* elem2){
return (elem->px - elem2->px)*(elem->px - elem2->px) +
(elem->py - elem2->py)*(elem->py - elem2->py) +
(elem->pz - elem2->pz)*(elem->pz - elem2->pz);

}

void check_boundaries(double WIDTH, double HEIGHT, double DEPTH, double damping_coeff,  __global SPHParticle* elem){

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


void __calculate_pressure_density(__global SPHParticle* sphParticle, __global double* enviromentArray){


double kernel_smoother_length = enviromentArray[6];
double current_particles = enviromentArray[3];
double norm_density = enviromentArray[7];
double density = enviromentArray[14];
double gas_constant = enviromentArray[1];

int i = get_global_id(0);
__global SPHParticle *elem = &sphParticle[i];
for (int j = 0; j < current_particles; j++) {
if (i != j){
__global SPHParticle* elem2 = &sphParticle[j];
double abs_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
double distance_rooted = sqrt(abs_distance_pow2);
if (distance_rooted <= kernel_smoother_length && distance_rooted >= 0){
elem->density += norm_density * kernel_pow((kernel_pow(kernel_smoother_length, 2) - kernel_pow( distance_rooted,
2)), 3) * elem2->mass;
}
}
}
if (elem->density < density) elem->density = density;
elem->pressure = gas_constant * (elem->density - density);
}

void __calculate_viscosity(__global SPHParticle* sphParticle, __global double* enviromentArray){

double kernel_smoother_length = enviromentArray[6];
double current_particles = enviromentArray[3];
double norm_pressure = enviromentArray[8];
double norm_visc = enviromentArray[9];
double U = enviromentArray[2];

int i = get_global_id(0);


__global SPHParticle* elem = &sphParticle[i];

for(int j = 0; j < current_particles; j++) {

if (i != j) {
__global SPHParticle* elem2 = &sphParticle[j];
double abs_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
double neighboring_particle_distance = sqrt(abs_distance_pow2);

if (neighboring_particle_distance <= kernel_smoother_length && neighboring_particle_distance >= 0 && elem->pressure != 0 && elem2->density != 0) {

double pressure_weight =
        norm_pressure *
         kernel_pow((kernel_smoother_length - neighboring_particle_distance), 2);
double pforce =
        (elem->pressure + elem2->pressure) *
        elem2->mass /
        (2 * elem2->density) *
        pressure_weight;
double distance_x = elem2->px - elem->px;
double distance_y = elem2->py - elem->py;
double distance_z = elem2->pz - elem->pz;
elem->fx += distance_x * pforce;
elem->fy += distance_y * pforce;
elem->fz += distance_z * pforce;

}

if (neighboring_particle_distance <= kernel_smoother_length && neighboring_particle_distance >= 0 && elem2->density != 0) {
double viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
double vel_difference_x = elem2->vx - elem->vx;
double vel_difference_y = elem2->vy - elem->vy;
double vel_difference_z = elem2->vz - elem->vz;
double v_force_x =
        vel_difference_x * (elem2->mass / elem2->density * viscosity_weight);
double v_force_y =
        vel_difference_y * (elem2->mass / elem2->density * viscosity_weight);
double v_force_z =
        vel_difference_z * (elem2->mass / elem2->density * viscosity_weight);
elem->fx += U * (v_force_x);
elem->fy += U * (v_force_y);
elem->fz += U * (v_force_z);
}
}
}

}


void __move_particles(__global SPHParticle* sphParticle, __global double* dtime, __global double* enviromentArray){
double current_particles = enviromentArray[3];
double gravity = enviromentArray[13];
double WIDTH = enviromentArray[10];
double HEIGHT = enviromentArray[11];
double DEPTH = enviromentArray[12];
double damping_coeff = enviromentArray[0];


int i = get_global_id(0);


__global SPHParticle* elem = &sphParticle[i];
elem->fy +=  gravity;

if (elem->density != 0){
elem->vx = elem->vx + *dtime * elem->fx / elem->density;
elem->vy = elem->vy + *dtime * elem->fy / elem->density;
elem->vz = elem->vz + *dtime * elem->fz / elem->density;
}
elem->px = elem->px + *dtime * elem->vx;
elem->py = elem->py + *dtime * elem->vy;
elem->pz = elem->pz + *dtime * elem->vz;

check_boundaries(WIDTH, HEIGHT, DEPTH, damping_coeff, elem);
elem->pressure = 0;
elem->fx = 1;
elem->fy = 1;
elem->fz = 1;
elem->density = 0;
}

__kernel void calculate_physics(__global SPHParticle* sphParticle, __global double* dtime, __global double* enviromentArray){
__calculate_pressure_density(sphParticle, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
__calculate_viscosity(sphParticle, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
__move_particles(sphParticle, dtime, enviromentArray);
barrier(CLK_GLOBAL_MEM_FENCE);
}
