#include "Container.h"

Container::Container(double width, double height, double depth, double gravity, double density, unsigned int particlesNumber, double particle_radius) :
        width(width),
        height(height),
        depth(depth),
        gravity(gravity),
        density(density),
        particles_number(
                particlesNumber),
                particle_radius(particle_radius),
        particles(new SPHParticle[particlesNumber]){
    if(this->width <= 0 || this->height <=0 || this->depth <=0 || particles_number <=0){
        throw std::invalid_argument("Every value of the container's attribute must be bigger than 0!");
    }
    max_number_of_grids = calculate_max_grid_number();
    std::cout<<normalization_density<<std::endl;
    std::cout<<norm_pressure<<std::endl;
    std::cout<<norm_visc<<std::endl;

    for(int i = 1; i <=max_number_of_grids; i++){
        neighbouring_grids.push_back(calculate_neighboring_grids(i));
        particles_on_grid.emplace_back();
    }

}

template<typename T> void Container::print(T asd){
    std::cout<<asd<<std::endl;
}

Container::~Container() {
    delete[] particles;
}

double Container::get_width() const {
    return width;
}

void Container::set_width(unsigned int width) {
    if(width <=0) throw std::invalid_argument("The WIDTH of the container must be bigger than 0!");
    width = width;
}

double Container::get_height() const {
    return height;
}

void Container::set_height(unsigned int height) {
    if(height <=0) throw std::invalid_argument("The HEIGHT of the container must be bigger than 0!");
    height = height;
}

double Container::get_depth() const {
    return depth;
}

void Container::set_depth(unsigned int depth) {
    if(depth <=0) throw std::invalid_argument("The DEPTH of the container must be bigger than 0!");
    depth = depth;
}

unsigned int Container::get_particles_number() const {
    return particles_number;
}

SPHParticle *Container::get_particles() const {
    return particles;
}

unsigned int Container::get_current_particles() const {
    return current_particles;
}

double Container::get_damping_coeff() const {
    return damping_coeff;
}

unsigned int Container::get_gas_constant() const {
    return gas_constant;
}

double Container::get_u() const {
    return U;
}

double Container::get_kernel_smoother_length() const {
    return kernel_smoother_length;
}

double Container::get_normalization_density() const {
    return normalization_density;
}

double Container::get_norm_pressure() const {
    return norm_pressure;
}

double Container::get_norm_visc() const {
    return norm_visc;
}

double Container::get_gravity() const {
    return gravity;
}

double Container::get_density() const {
    return density;
}

double Container::get_particle_radius() const {
    return particle_radius;
}

unsigned int Container::get_max_number_of_grids() const {
    return max_number_of_grids;
}

unsigned int Container::get_max_column_per_row_grid() const {
    return max_column_per_row_grid;
}

unsigned int Container::get_max_row_per_layer_grid() const {
    return max_row_per_layer_grid;
}

unsigned int Container::get_max_layer_per_cube_grid() const {
    return max_layer_per_cube_grid;
}

const std::vector<std::set<int>> &Container::get_neighbouring_grids() const {
    return neighbouring_grids;
}

const std::vector<std::set<int>> &Container::get_particles_on_grid() const {
    return particles_on_grid;
}

void Container::calculate_physics(double dtime){

    auto* particle_pressures = new double[particles_number]{ 0 }; // the other particles pressure
    auto* particle_forces = new double[particles_number][3] { { 1 } }; //the force on the particle x,y,z
    auto* particle_densities = new double[particles_number]{ 0 }; // store each density at the particles location

    std::vector<std::vector<std::pair<int, double>>> neighbor_particles;


    current_iteration++;
    if(current_particles < particles_number && current_iteration % new_particles_every_iteration == 0)
        fill_container_gradually(particle_radius);


    for(int i = 0; i < current_particles; i++) {
        std::vector<std::pair<int, double>> _vector_of_neighbour_pairs;
        SPHParticle &elem = particles[i];

        std::set<int> _nieghbour_grids = neighbouring_grids[calculate_grid(elem) - 1];

        for (auto n_grid: _nieghbour_grids) {


            for (auto j: particles_on_grid[n_grid - 1]) {

                if (i != j) {

                    SPHParticle &elem2 = particles[j];
                    double particle_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
                    double particle_distance = sqrt(particle_distance_pow2);
                    if (particle_distance <= kernel_smoother_length  &&
                        particle_distance >= 0) {
                        _vector_of_neighbour_pairs.emplace_back(j, sqrt(particle_distance));
                    }
                }
            }
        }
        neighbor_particles.push_back(_vector_of_neighbour_pairs);

        for(int j = 0; j < neighbor_particles[i].size(); j++){
            particle_densities[i] += normalization_density * kernel_pow(
                    (kernel_pow(kernel_smoother_length, 2) - kernel_pow(neighbor_particles[i][j].second, 2)), 3) *
                            particles[neighbor_particles[i][j].first].get_mass();
        }

        if(particle_densities[i] < density) particle_densities[i] = density;
        particle_pressures[i] = gas_constant * (particle_densities[i] - density); //pressure = k (p - p0)
    }




    for(int i = 0; i < current_particles; i ++) {
        SPHParticle elem = particles[i];
        std::vector<int> list_of_neighbours;
        for(auto elem_vec : neighbor_particles[i]){
            list_of_neighbours.push_back(elem_vec.first);
        }

        for(int j = 0; j < list_of_neighbours.size(); j++) {

            SPHParticle &elem2 = particles[list_of_neighbours[j]];
            int neighbor_index = list_of_neighbours[j];
            double neighboring_particle_distance = neighbor_particles[i][j].second;

            if (particle_pressures[i] != 0 && particle_densities[neighbor_index] != 0){

                double pressure_weight = norm_pressure *
                        kernel_pow((kernel_smoother_length - neighboring_particle_distance), 2);
                double pforce =  (particle_pressures[i] + particle_pressures[neighbor_index]) *
                        particles[neighbor_index].get_mass() /
                        (2 * particle_densities[neighbor_index]) * pressure_weight;
                double distance_x = elem2.get_px() - elem.get_px();
                double distance_y = elem2.get_py() - elem.get_py();
                double distance_z = elem2.get_pz() - elem.get_pz();
                particle_forces[i][0] += distance_x * pforce;
                particle_forces[i][1] += distance_y * pforce;
                particle_forces[i][2] += distance_z * pforce;
            }

            if (particle_densities[neighbor_index] != 0) {

                double viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
                double vel_difference_x = particles[neighbor_index].get_vx() - particles[i].get_vx();
                double vel_difference_y = particles[neighbor_index].get_vy() - particles[i].get_vy();
                double vel_difference_z = particles[neighbor_index].get_vz() - particles[i].get_vz();
                double v_force_x = vel_difference_x * (particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                double v_force_y = vel_difference_y * (particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                double v_force_z = vel_difference_z * (particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                particle_forces[i][0] += U * (v_force_x );
                particle_forces[i][1] += U * (v_force_y );
                particle_forces[i][2] += U * (v_force_z );
            }

        }
    }

    for(int i = 0; i < current_particles; i++){

        SPHParticle &elem = particles[i];
        particle_forces[i][1] += gravity;
        if (particle_densities[i] != 0) {
            elem.set_vx(elem.get_vx() + dtime * particle_forces[i][0] / particle_densities[i]);
            elem.set_vy(elem.get_vy() + dtime * particle_forces[i][1] / particle_densities[i]);
            elem.set_vz(elem.get_vz() + dtime * particle_forces[i][2] / particle_densities[i]);
        }

        particles_on_grid[calculate_grid(elem) - 1].erase(i);

        elem.set_px(elem.get_px() + dtime * elem.get_vx());
        elem.set_py(elem.get_py() + dtime * elem.get_vy());
        elem.set_pz(elem.get_pz() + dtime * elem.get_vz());


        check_boundaries(elem);
        particles_on_grid[calculate_grid(elem) - 1].insert(i);


    }
    delete[] particle_densities;
    delete[] particle_forces;
    delete[] particle_pressures;
}

void Container::fill_container_gradually2(const double radius){
    int row = 0;
    int column = 0;
    double x_cor;
    double z_cor;
    int max_particle_per_row = this->width / kernel_smoother_length;
    int max_particle_per_column = this->depth / kernel_smoother_length;

    for(int i = 0; i < 200; i++) {
        this->particles[current_particles] = SPHParticle(1,  column * kernel_smoother_length/2, height - 2 * radius, row * kernel_smoother_length/2,
                                                        10, -10, -13, radius, 1);



        this->current_particles++;
        if(this->current_particles == this->particles_number) break;

        this->particles[current_particles] = SPHParticle(1,  width - (column * kernel_smoother_length/2), height - 2 * radius, depth - (row * kernel_smoother_length/2),
                                                        -10, -10, 13, radius, 1);
        this->current_particles++;
        if(this->current_particles == this->particles_number) break;

        column++;
        if (i != 0 && i % max_particle_per_row == 0) {
            row++;
            column = 0;
        }
        if(row > max_particle_per_column) break;
    }
}

unsigned int Container::get_current_iteration() const {
    return current_iteration;
}

void Container::set_current_iteration(unsigned int iteration){
    current_iteration = iteration;
}

unsigned int Container::get_new_particles_every_iteration() const {
    return new_particles_every_iteration;
}



 unsigned int Container::calculate_grid(SPHParticle elem){

     unsigned int column = elem.get_px() / kernel_smoother_length + 1;
     unsigned int row = elem.get_py() / kernel_smoother_length + 1;
     unsigned int layer = elem.get_pz() / kernel_smoother_length + 1;


  if(column > max_column_per_row_grid) column = max_column_per_row_grid;
  if(row > max_row_per_layer_grid) row = max_row_per_layer_grid;
  if(layer > max_layer_per_cube_grid) layer = max_layer_per_cube_grid;

     unsigned int return_value = (row - 1) * max_column_per_row_grid + column + (max_row_per_layer_grid * max_column_per_row_grid * (layer - 1));


     return return_value;
}

unsigned int Container::calculate_max_grid_number(){

    unsigned int column = width / kernel_smoother_length;
    unsigned int row = height / kernel_smoother_length;
    unsigned int layer = depth / kernel_smoother_length;

    //IF IT HAS A REMAINDER IT ADDS +1. FOR EXAMPLE  5/2 = 3 BECAUSE int(5/2) == 2

    if(column != width / kernel_smoother_length || column == 0) column++;
    if(row != height / kernel_smoother_length || row == 0) row++;
    if(layer != depth / kernel_smoother_length || layer == 0) layer++;

    max_column_per_row_grid = column;
    max_row_per_layer_grid = row;
    max_layer_per_cube_grid = layer;

    return column * row * layer;
}

double Container::calculate_absolute_distance_power2(SPHParticle& elem, SPHParticle& elem2){
    return (elem.get_px() - elem2.get_px())*(elem.get_px() - elem2.get_px()) +
           (elem.get_py() - elem2.get_py())*(elem.get_py() - elem2.get_py()) +
           (elem.get_pz() - elem2.get_pz())*(elem.get_pz() - elem2.get_pz());
}

void Container::check_boundaries(SPHParticle& elem) const{

    if (elem.get_py() + elem.get_radius()  >= height) {
        elem.set_vy(elem.get_vy()  * damping_coeff * -1);
        elem.set_py(height - elem.get_radius());
    }else if (elem.get_py() - elem.get_radius()  <= 0){
        elem.set_vy(elem.get_vy()  * damping_coeff * -1);
        elem.set_py(0 + elem.get_radius() );
    }

    if (elem.get_px() + elem.get_radius()  >= width) {
        elem.set_vx(elem.get_vx()  * damping_coeff * -1);
        elem.set_px(width - elem.get_radius());
    }else if (elem.get_px() - elem.get_radius()  <= 0){
        elem.set_vx(elem.get_vx()  * damping_coeff * -1);
        elem.set_px(0 + elem.get_radius() );
    }

    if (elem.get_pz() + elem.get_radius()  >= depth) {
        elem.set_vz(elem.get_vz()  * damping_coeff * -1);
        elem.set_pz(depth - elem.get_radius());
    }else if (elem.get_pz() - elem.get_radius()  <= 0){
        elem.set_vz(elem.get_vz()  * damping_coeff * -1);
        elem.set_pz(0 + elem.get_radius() );
    }

}



void Container::fill_container_gradually(double radius){
    int row = 0;
    int column = 0;
    int max_particle_per_row = width / kernel_smoother_length;
    int max_particle_per_column = depth / kernel_smoother_length;

    for(int i = 0; i < 20; i++) {

        SPHParticle elem = SPHParticle(1,  column * kernel_smoother_length/2, height - 2 * radius, row * kernel_smoother_length/2,
                    10, -10, -5, radius, 1);

        particles[current_particles] = elem;
        unsigned int grid_number = calculate_grid(elem);

        particles_on_grid[grid_number - 1].insert(current_particles);

        current_particles++;
        if(current_particles == particles_number) break;
        elem = SPHParticle(1,  width - (column * kernel_smoother_length/2), height - 2 * radius, depth - (row * kernel_smoother_length/2),
                                                         -10, -10, 5, radius, 1);

        particles[current_particles] = elem;
        grid_number = calculate_grid(elem);

        particles_on_grid[grid_number - 1].insert(current_particles);

        current_particles++;
        if(current_particles == particles_number) break;

        column++;
        if (i != 0 && i % max_particle_per_row == 0) {
            row++;
            column = 0;
        }
        if(row > max_particle_per_column) break;
    }
}

void Container::fill_container(double radius){
    unsigned z_sor = 0;
    unsigned y_sor = 1;
    double offset = radius * 3;
    bool full = false;
    int modulus = width / offset;
    for (int i = 0; i<particles_number; i++){
        double shift = i % 2 == 1 ? radius : 0;

        if(i % modulus == modulus - 2) z_sor++;
        if(z_sor * offset + shift + radius >= get_depth() && i % modulus == modulus - 2) {
            y_sor++;
            z_sor = 0;
        }

        if (y_sor * radius * 3 >= height - radius){
            particles_number = i - 1;
            std::cout<<"Rendered particles (because the container is full):  "<<particles_number <<std::endl;
            full = true;
            break;
        }

        double random_x = double(rand())/double((RAND_MAX)) * 11 - 6 ;
        double random_z = double(rand())/double((RAND_MAX)) * 11 - 6 ;
        SPHParticle particle_to_add = SPHParticle(1, i % modulus * offset + shift + offset/2, y_sor * offset, z_sor * offset + shift,
                                                  2, -10, 0, radius, 1);
        particles[i] = particle_to_add;
    }

    if (!full)std::cout<<"Rendered particles:  "<<particles_number <<std::endl;
}

std::set<int> Container::calculate_neighboring_grids(int given_grid){
    //INDEXING STARTS FROM ONE

    std::set<int> return_set;
    std::vector<int> vector_to_compare;

    int right = 1;
    int left = -1;
    int up = max_column_per_row_grid;
    int down = -up;
    int front = max_column_per_row_grid * max_row_per_layer_grid;
    int back = -front;

    int up_right = up + right;
    int up_left = up + left;
    int up_front = up + front;
    int up_back = up + back;
    int down_right = down + right;
    int down_left = down + left;
    int down_front = down + front;
    int down_back = down + back;
    int right_front = right + front;
    int left_front = left + front;
    int right_back = right + back;
    int left_back = left + back;
    int up_right_front = up + right + front;
    int up_left_front = up + left + front;
    int up_right_back = up + right + back;
    int up_left_back = up + left + back;
    int down_right_front = down + right + front;
    int down_left_front = down + left + front;
    int down_right_back = down + right + back;
    int down_left_back = down + left + back;

    vector_to_compare.push_back(given_grid);

    if(given_grid == 1){

        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + up_right_front);

    }else if(given_grid == max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + up_left_front);

    }else if (given_grid == (front - max_column_per_row_grid + 1)){

        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + down_right_front);

    }else if(given_grid == front){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + down_left_front);

    }else if (given_grid == max_number_of_grids - back + 1){

        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + up_right_back);

    }else if (given_grid == max_number_of_grids - back + max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + up_left_back);

    }else if (given_grid == max_number_of_grids - max_column_per_row_grid + 1){

        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + down_right_back);

    }else if (given_grid == max_number_of_grids){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + down_left_back);

    }else if(given_grid % max_column_per_row_grid == 1 && given_grid < front){

        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + up_right_front);
        vector_to_compare.push_back(given_grid + down_right_front);

    }else if(given_grid % max_column_per_row_grid == 0 && given_grid < front){

        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + up_left_front);
        vector_to_compare.push_back(given_grid + down_left_front);

    }else if(given_grid < front &&
             given_grid > (front - max_column_per_row_grid)){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + down_left_front);
        vector_to_compare.push_back(given_grid + down_right_front);

    }else if(given_grid > 1 && given_grid < max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + up_left_front);
        vector_to_compare.push_back(given_grid + up_right_front);

    }else if(given_grid % max_column_per_row_grid == 1 && given_grid > max_number_of_grids - front){

        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + up_right_back);
        vector_to_compare.push_back(given_grid + down_right_back);

    }else if(given_grid % max_column_per_row_grid == 0 && given_grid > max_number_of_grids - front){

        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + up_left_back);
        vector_to_compare.push_back(given_grid + down_left_back);

    }else if(given_grid > max_number_of_grids - front && given_grid < max_number_of_grids - front + max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + down_left_back);
        vector_to_compare.push_back(given_grid + down_right_back);

    }else if(given_grid > max_number_of_grids - back + 1 && given_grid < max_number_of_grids - back + max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + up_left_back);
        vector_to_compare.push_back(given_grid + up_right_back);

    }else if(given_grid % front == 1){

        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + up_right_front);
        vector_to_compare.push_back(given_grid + up_right_back);

    }else if(given_grid % front == max_column_per_row_grid){

        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + up_left_front);
        vector_to_compare.push_back(given_grid + up_left_back);

    }else if(given_grid % front == front - max_column_per_row_grid + 1){

        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + down_right_front);
        vector_to_compare.push_back(given_grid + down_right_back);

    }else if(given_grid % front == 0){

        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + down_left_front);
        vector_to_compare.push_back(given_grid + down_left_back);

    }else{

        vector_to_compare.push_back(given_grid + right);
        vector_to_compare.push_back(given_grid + left);
        vector_to_compare.push_back(given_grid + up);
        vector_to_compare.push_back(given_grid + down);
        vector_to_compare.push_back(given_grid + front);
        vector_to_compare.push_back(given_grid + back);
        vector_to_compare.push_back(given_grid + up_right);
        vector_to_compare.push_back(given_grid + up_left);
        vector_to_compare.push_back(given_grid + up_front);
        vector_to_compare.push_back(given_grid + up_back);
        vector_to_compare.push_back(given_grid + down_right);
        vector_to_compare.push_back(given_grid + down_left);
        vector_to_compare.push_back(given_grid + down_front);
        vector_to_compare.push_back(given_grid + down_back);
        vector_to_compare.push_back(given_grid + right_front);
        vector_to_compare.push_back(given_grid + left_front);
        vector_to_compare.push_back(given_grid + right_back);
        vector_to_compare.push_back(given_grid + left_back);
        vector_to_compare.push_back(given_grid + up_right_front);
        vector_to_compare.push_back(given_grid + up_left_front);
        vector_to_compare.push_back(given_grid + up_right_back);
        vector_to_compare.push_back(given_grid + up_left_back);
        vector_to_compare.push_back(given_grid + down_right_front);
        vector_to_compare.push_back(given_grid + down_left_front);
        vector_to_compare.push_back(given_grid + down_right_back);
        vector_to_compare.push_back(given_grid + down_left_back);

    }

    for(auto elem : vector_to_compare){
        if(elem <= max_number_of_grids && elem > 0){
            return_set.insert(elem);
        }
    }

    return return_set;
}

double Container::kernel_pow(double ertek, unsigned int order){
    double vegeredmeny = 1;
    for (int i = 0; i < order; ++i) vegeredmeny *= ertek;
    return vegeredmeny;

}


void Container::calculate_physicsv2(const double dtime){

    auto* particle_pressures = new double[this->particles_number]{ 0 }; // the other particles pressure
    auto* particle_forces = new double[this->particles_number][3] { { 1 } }; //the force on the particle x,y,z
    auto* particle_densities = new double[this->particles_number]{ 0 }; // store each density at the particles location
    std::vector<std::vector<std::pair<int, double>>> neighbor_particles;

    current_iteration++;
    if(this->current_particles < this->particles_number && current_iteration % new_particles_every_iteration == 0)
        fill_container_gradually(this->particle_radius);


    for(int i = 0; i < this->current_particles; i++){
        std::vector<std::pair<int, double>> _vector_of_neighbour_pairs;
        SPHParticle& elem = this->particles[i];

        for(int j = 0; j < this->current_particles; j++ ){
            if (i != j) {
                SPHParticle &elem2 = this->particles[j];
                double particle_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
                double particle_distance = sqrt(particle_distance_pow2);
                if (particle_distance <= kernel_smoother_length  && particle_distance >= 0){
                    _vector_of_neighbour_pairs.emplace_back(j, sqrt(particle_distance));
                }
            }
        }

        neighbor_particles.push_back(_vector_of_neighbour_pairs);



        for(auto & _vector_of_neighbour_pair : _vector_of_neighbour_pairs){
          particle_densities[i] += normalization_density *
                   kernel_pow((kernel_pow(kernel_smoother_length, 2) -
                kernel_pow(_vector_of_neighbour_pair.second, 2)), 3) * this->particles[_vector_of_neighbour_pair.first].get_mass();
        }
        if(particle_densities[i] < this->density) particle_densities[i] = this->density;

        particle_pressures[i] = gas_constant * (particle_densities[i] - this->density); //pressure = k (p - p0)
    }



    // r - position, p - pressure /p/ - density

    for(int i = 0; i < this->current_particles; i++){

        SPHParticle elem = this->particles[i];
        std::vector<int> list_of_neighbours;
        for(auto elem_vec : neighbor_particles[i]){
            list_of_neighbours.push_back(elem_vec.first);
        }

        for(int j = 0; j < list_of_neighbours.size(); j++) {

            SPHParticle &elem2 = this->particles[list_of_neighbours[j]];
            int neighbor_index = list_of_neighbours[j];
            double neighboring_particle_distance = neighbor_particles[i][j].second;

            if (particle_pressures[i] != 0 && particle_densities[neighbor_index] != 0){

                double pressure_weight =
                        norm_pressure *
                                kernel_pow((kernel_smoother_length - neighboring_particle_distance), 2);
                double pforce =
                        (particle_pressures[i] + particle_pressures[neighbor_index]) *
                                elem2.get_mass() /
                        (2 * particle_densities[neighbor_index]) *
                        pressure_weight;
                double distance_x = elem2.get_px() - elem.get_px();
                double distance_y = elem2.get_py() - elem.get_py();
                double distance_z = elem2.get_pz() - elem.get_pz();
                particle_forces[i][0] += distance_x * pforce;
                particle_forces[i][1] += distance_y * pforce;
                particle_forces[i][2] += distance_z * pforce;
            }

            if (particle_densities[neighbor_index] != 0) {
                double viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
                double vel_difference_x = elem2.get_vx() - elem.get_vx();
                double vel_difference_y = elem2.get_vy() - elem.get_vy();
                double vel_difference_z = elem2.get_vz() - elem.get_vz();
                double v_force_x = vel_difference_x * (elem2.get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                double v_force_y = vel_difference_y * (elem2.get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                double v_force_z = vel_difference_z * (elem2.get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                particle_forces[i][0] += U * (v_force_x );
                particle_forces[i][1] += U * (v_force_y );
                particle_forces[i][2] += U * (v_force_z );
            }
        }
    }

    for(int i = 0; i < this->current_particles; i++){

        particle_forces[i][1] += this->gravity;
        SPHParticle &elem = this->particles[i];

        if (particle_densities[i] != 0) {

            elem.set_vx(elem.get_vx() + dtime * particle_forces[i][0] / particle_densities[i]);
            elem.set_vy(elem.get_vy() + dtime * particle_forces[i][1] / particle_densities[i]);
            elem.set_vz(elem.get_vz() + dtime * particle_forces[i][2] / particle_densities[i]);
        }
        elem.set_px(elem.get_px() + dtime * elem.get_vx());
        elem.set_py(elem.get_py() + dtime * elem.get_vy());
        elem.set_pz(elem.get_pz() + dtime * elem.get_vz());



        check_boundaries(elem);

    }

    delete[] particle_densities;
    delete[] particle_forces;
    delete[] particle_pressures;
    //   delete[] particle_distances;


}


void Container::calculate_physicsv3(double dtime, double* enviromentArray){

    if(this->current_particles < this->particles_number && current_iteration % new_particles_every_iteration == 0)
        fill_container_gradually(this->particle_radius);
    current_iteration++;
    for (int i = 0; i < current_particles;i++) {
        SPHParticle &elem = this->particles[i];
        for (int j = 0; j < current_particles; j++) {
            SPHParticle &elem2 = this->particles[j];
            double abs_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
            double neighboring_particle_distance = sqrt(abs_distance_pow2);
            if (neighboring_particle_distance <= kernel_smoother_length  && neighboring_particle_distance >= 0) {
                elem.density += normalization_density *
                                 kernel_pow((kernel_pow(kernel_smoother_length, 2) - kernel_pow(sqrt(abs_distance_pow2),
                                                                                                2)), 3) *
                                 elem2.mass;
            }
        }
        if (elem.density < density) elem.density = density;
        elem.pressure = gas_constant * (elem.density - density);
    }


    for (int i = 0; i < current_particles;i++) {
    SPHParticle& elem = this->particles[i];

    for(int j = 0; j < current_particles; j++) {

        if (i != j) {
            SPHParticle& elem2 = this->particles[j];
            double abs_distance_pow2 = calculate_absolute_distance_power2(elem, elem2);
            double neighboring_particle_distance = sqrt(abs_distance_pow2);

            if (neighboring_particle_distance <= kernel_smoother_length  && neighboring_particle_distance >= 0 && elem.pressure != 0 && elem2.density != 0) {

                double pressure_weight =
                        norm_pressure *
                        kernel_pow((kernel_smoother_length - neighboring_particle_distance), 2);
                double pforce =
                        (elem.pressure + elem2.pressure) *
                        elem2.mass /
                        (2 * elem2.density) *
                        pressure_weight;
                double distance_x = elem2.px - elem.px;
                double distance_y = elem2.py - elem.py;
                double distance_z = elem2.pz - elem.pz;
                elem.fx += distance_x * pforce;
                elem.fy += distance_y * pforce;
                elem.fz += distance_z * pforce;
            }

            if (neighboring_particle_distance <= kernel_smoother_length  && neighboring_particle_distance >= 0 && elem2.density != 0) {
                double viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
                double vel_difference_x = elem2.vx - elem.vx;
                double vel_difference_y = elem2.vy - elem.vy;
                double vel_difference_z = elem2.vz - elem.vz;
                double v_force_x =
                        vel_difference_x * (elem2.mass / elem2.density * viscosity_weight);
                double v_force_y =
                        vel_difference_y * (elem2.mass / elem2.density * viscosity_weight);
                double v_force_z =
                        vel_difference_z * (elem2.mass / elem2.density * viscosity_weight);
                elem.fx += U * (v_force_x);
                elem.fy += U * (v_force_y);
                elem.fz += U * (v_force_z);
            }
        }
    }
     }



    for (int i = 0; i < current_particles;i++) {


    SPHParticle& elem = this->particles[i];


    elem.fy +=  gravity;

    if (elem.density != 0) {
        elem.vx = elem.vx + dtime * elem.fx / elem.density;
        elem.vy = elem.vy + dtime * elem.fy / elem.density;
        elem.vz = elem.vz + dtime * elem.fz / elem.density;
    }
    elem.px = elem.px + dtime * elem.vx;
    elem.py = elem.py + dtime * elem.vy;
    elem.pz = elem.pz + dtime * elem.vz;

    check_boundaries(elem);

    elem.pressure = 0;
    elem.fx = 1;
    elem.fy = 1;
    elem.fz = 1;
    elem.density = 0;

    }
}