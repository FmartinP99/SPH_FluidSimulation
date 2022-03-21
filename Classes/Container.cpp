#include <iostream>
#include <string>
#include "SPHParticle.cpp"
#include <cmath>
#include <vector>
#include <set>


class Container{
private:

    double damping_coeff = 0.25;
    unsigned int gas_constant = 400;
    double kernel_smoother_length = 1;
    double U = 1;

    unsigned int current_particles = 0;
    unsigned int current_iteration = 0;
    unsigned int new_particles_every_iteration = 3;

    double normalization_density;
    double norm_pressure;
    double norm_visc;

    double width;
    double height;
    double depth;
    double gravity;
    double density;
    unsigned int particles_number;
    double particle_radius;
    SPHParticle* particles;
    unsigned int max_number_of_grids;
    unsigned int max_column_per_row_grid;
    unsigned int max_row_per_layer_grid;
    unsigned int max_layer_per_cube_grid;
    std::vector<std::set<int>> neighbouring_grids;  // megadja melyik  gridnek melyek a szomszédjai
    std::vector<std::set<int>> particles_on_grid; // megadja melyik griden melyik indexű particle-ek vannak

public:

    Container(double width, double height, double depth, double gravity, double density, unsigned int particlesNumber, double particle_radius) :
            width(width),
            height(height),
            depth(depth),
            gravity(gravity),
            density(density),
            particles_number(
                    particlesNumber),
                    particle_radius(particle_radius),
            particles(new SPHParticle[particlesNumber]){
        if(this->width <= 0 || this->height <=0 || this->depth <=0 || this->particles_number <=0){
            throw std::invalid_argument("Every value of the container's attribute must be bigger than 0!");
        }
       // fill_container(particle_radius);
        this->max_number_of_grids = calculate_max_grid_number();
        normalization_density = (315 / (64 * M_PI * pow(kernel_smoother_length, 9)));
        norm_pressure = (-45 / (M_PI * pow(kernel_smoother_length, 6)));
        norm_visc = (45 / ( M_PI * pow(kernel_smoother_length, 6)));
        std::cout<<normalization_density<<std::endl;
        std::cout<<norm_pressure<<std::endl;
        std::cout<<norm_visc<<std::endl;

        for(int i = 1; i <=max_number_of_grids; i++){
            neighbouring_grids.push_back(calculate_neighboring_grids(i));
            particles_on_grid.emplace_back();
        }

    }

    virtual ~Container() {
        delete[] particles;
    }

    double get_width() const {
        return width;
    }

    void set_width(unsigned int width) {
        if(width <=0) throw std::invalid_argument("The width of the container must be bigger than 0!");
        Container::width = width;
    }

    double get_height() const {
        return height;
    }

    void set_height(unsigned int height) {
        if(height <=0) throw std::invalid_argument("The height of the container must be bigger than 0!");
        Container::height = height;
    }

    double get_depth() const {
        return depth;
    }

    void set_depth(unsigned int depth) {
        if(depth <=0) throw std::invalid_argument("The depth of the container must be bigger than 0!");
        Container::depth = depth;
    }

    double get_particles_number() const {
        return particles_number;
    }

    SPHParticle *get_particles() const {
        return particles;
    }

    unsigned int get_current_particles() const {
        return current_particles;
    }


    void calculate_physics(double dtime){

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
            std::set<int> _nieghbour_grids = this->neighbouring_grids[calculate_grid(elem) - 1]; // 1 ->  {2, 6, 7, 26, 27, 31, 32}  5X5 kocka esetén ami a tömb 0. indexén van


            for(auto n_grid : _nieghbour_grids) {

                for (auto j : this->particles_on_grid[n_grid - 1]) {

                   if (i != j) {

                        SPHParticle &elem2 = this->particles[j];
                        double particle_distance = calculate_absolute_distance_power2(elem, elem2);
                        if (particle_distance < kernel_smoother_length * kernel_smoother_length &&
                            particle_distance > 0) {
                            _vector_of_neighbour_pairs.emplace_back(j, sqrt(particle_distance));
                        }
                   }
                }
            }
            neighbor_particles.push_back(_vector_of_neighbour_pairs);

            for(int j = 0; j < neighbor_particles[i].size(); j++){
                particle_densities[i] += normalization_density * pow((pow(kernel_smoother_length, 2) - pow(neighbor_particles[i][j].second, 2)), 3) * this->particles[neighbor_particles[i][j].first].get_mass();
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

                    double pressure_weight = norm_pressure * pow((kernel_smoother_length - neighboring_particle_distance), 2);
                    double pforce =  (particle_pressures[i] + particle_pressures[neighbor_index]) * this->particles[neighbor_index].get_mass() / (2 * particle_densities[neighbor_index]) * pressure_weight;
                    double distance_x = elem2.get_px() - elem.get_px();
                    double distance_y = elem2.get_py() - elem.get_py();
                    double distance_z = elem2.get_pz() - elem.get_pz();
                    particle_forces[i][0] += distance_x * pforce;
                    particle_forces[i][1] += distance_y * pforce;
                    particle_forces[i][2] += distance_z * pforce;

                }

                if (particle_densities[neighbor_index] != 0) {
                    double viscosity_weight = norm_visc * (kernel_smoother_length - neighboring_particle_distance);
                    double vel_difference_x = this->particles[neighbor_index].get_vx() - this->particles[i].get_vx();
                    double vel_difference_y = this->particles[neighbor_index].get_vy() - this->particles[i].get_vy();
                    double vel_difference_z = this->particles[neighbor_index].get_vz() - this->particles[i].get_vz();
                    double v_force_x = vel_difference_x * (this->particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                    double v_force_y = vel_difference_y * (this->particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
                    double v_force_z = vel_difference_z * (this->particles[neighbor_index].get_mass() / particle_densities[neighbor_index] * viscosity_weight);
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

private:

     unsigned int calculate_grid(SPHParticle elem){

         unsigned int column = elem.get_px() / kernel_smoother_length + 1;
         unsigned int row = elem.get_py() / kernel_smoother_length + 1;
         unsigned int layer = elem.get_pz() / kernel_smoother_length + 1;


      if(column > this->max_column_per_row_grid) column = this->max_column_per_row_grid;
      if(row > this->max_row_per_layer_grid) row = this->max_row_per_layer_grid;
      if(layer > this->max_layer_per_cube_grid) layer = this->max_layer_per_cube_grid;

         unsigned int return_value = (row - 1) * this->max_column_per_row_grid + column + (this->max_row_per_layer_grid * this->max_column_per_row_grid * (layer - 1));


         return return_value;
    }

    unsigned int calculate_max_grid_number(){

        unsigned int column = this->width / kernel_smoother_length;
        unsigned int row = this->height / kernel_smoother_length;
        unsigned int layer = this->depth / kernel_smoother_length;

        //IF IT HAS A REMAINDER IT ADDS +1. FOR EXAMPLE  5/2 = 3 BECAUSE int(5/2) == 2

        if(column != this->width / kernel_smoother_length || column == 0) column++;
        if(row != this->height / kernel_smoother_length || row == 0) row++;
        if(layer != this->depth / kernel_smoother_length || layer == 0) layer++;

        this->max_column_per_row_grid = column;
        this->max_row_per_layer_grid = row;
        this->max_layer_per_cube_grid = layer;

        return column * row * layer;
    }

    static double calculate_absolute_distance_power2(SPHParticle& elem, SPHParticle& elem2){
        return (elem.get_px() - elem2.get_px())*(elem.get_px() - elem2.get_px()) +
               (elem.get_py() - elem2.get_py())*(elem.get_py() - elem2.get_py()) +
               (elem.get_pz() - elem2.get_pz())*(elem.get_pz() - elem2.get_pz());
    }

    void check_boundaries(SPHParticle& elem) const{

        if (elem.get_py() + elem.get_radius()  >= this->height) {
            elem.set_vy(elem.get_vy()  * damping_coeff * -1);
            elem.set_py(this->height - elem.get_radius());
        }else if (elem.get_py() - elem.get_radius()  <= 0){
            elem.set_vy(elem.get_vy()  * damping_coeff * -1);
            elem.set_py(0 + elem.get_radius() );
        }

        if (elem.get_px() + elem.get_radius()  >= this->width) {
            elem.set_vx(elem.get_vx()  * damping_coeff * -1);
            elem.set_px(this->width - elem.get_radius());
        }else if (elem.get_px() - elem.get_radius()  <= 0){
            elem.set_vx(elem.get_vx()  * damping_coeff * -1);
            elem.set_px(0 + elem.get_radius() );
        }

        if (elem.get_pz() + elem.get_radius()  >= this->depth) {
            elem.set_vz(elem.get_vz()  * damping_coeff * -1);
            elem.set_pz(this->depth - elem.get_radius());
        }else if (elem.get_pz() - elem.get_radius()  <= 0){
            elem.set_vz(elem.get_vz()  * damping_coeff * -1);
            elem.set_pz(0 + elem.get_radius() );
        }

    }

    void fill_container_gradually(double radius){
        int row = 0;
        int column = 0;
        int max_particle_per_row = this->width / kernel_smoother_length; // so they are chill
        int max_particle_per_column = this->depth / kernel_smoother_length; // so they are chill

        for(int i = 0; i < 20; i++) {

            SPHParticle elem = SPHParticle(1,  column * kernel_smoother_length/2, this->height - 2 * radius, row * kernel_smoother_length/2,
                        10, -10, -13, radius, 1);

            this->particles[current_particles] = elem;
            unsigned int grid_number = calculate_grid(elem);
            particles_on_grid[grid_number - 1].insert(current_particles);
            this->current_particles++;
            if(this->current_particles == this->particles_number) break;

            elem = SPHParticle(1,  this->width - (column * kernel_smoother_length/2), this->height - 2 * radius, this->depth - (row * kernel_smoother_length/2),
                                                             -10, -10, 13, radius, 1);

            this->particles[current_particles] = elem;
            grid_number = calculate_grid(elem);
            particles_on_grid[grid_number - 1].insert(current_particles);
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

    void fill_container(double radius){
        unsigned z_sor = 0;
        unsigned y_sor = 1;
        double offset = radius * 3;
        bool full = false;
        int modulus = this->width / offset;
        for (int i = 0; i<this->particles_number; i++){
            double shift = i % 2 == 1 ? radius : 0;

            if(i % modulus == modulus - 2) z_sor++;
            if(z_sor * offset + shift + radius >= this->get_depth() && i % modulus == modulus - 2) {
                y_sor++;
                z_sor = 0;
            }

            if (y_sor * radius * 3 >= this->height - radius){
                this->particles_number = i - 1;
                std::cout<<"Rendered particles (because the container is full):  "<<this->particles_number <<std::endl;
                full = true;
                break;
            }

            float random_x = float(rand())/float((RAND_MAX)) * 11 - 6 ;
            float random_z = float(rand())/float((RAND_MAX)) * 11 - 6 ;
            SPHParticle particle_to_add = SPHParticle(1, i % modulus * offset + shift + offset/2, y_sor * offset, z_sor * offset + shift,
                                                      2, -10, 0, radius, 1);
            this->particles[i] = particle_to_add;
        }

        if (!full)std::cout<<"Rendered particles:  "<<this->particles_number <<std::endl;
    }

    std::set<int> calculate_neighboring_grids(int given_grid){
        //INDEXING STARTS FROM ONE
        //TODO IMPLEMENTING THIS

        std::set<int> return_set;
        std::vector<int> vector_to_compare;

        int right = 1;
        int left = -1;
        int up = this->max_column_per_row_grid;
        int down = -up;
        int front = this->max_column_per_row_grid * max_row_per_layer_grid;
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

        }else if(given_grid == this->max_column_per_row_grid){

            vector_to_compare.push_back(given_grid + left);
            vector_to_compare.push_back(given_grid + up);
            vector_to_compare.push_back(given_grid + up_left);
            vector_to_compare.push_back(given_grid + front);
            vector_to_compare.push_back(given_grid + up_front);
            vector_to_compare.push_back(given_grid + left_front);
            vector_to_compare.push_back(given_grid + up_left_front);

        }else if (given_grid == (front - this->max_column_per_row_grid + 1)){

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

        }else if (given_grid == this->max_number_of_grids - back + 1){ // bal felső hátsó

            vector_to_compare.push_back(given_grid + right);
            vector_to_compare.push_back(given_grid + up);
            vector_to_compare.push_back(given_grid + up_right);
            vector_to_compare.push_back(given_grid + back);
            vector_to_compare.push_back(given_grid + up_back);
            vector_to_compare.push_back(given_grid + right_back);
            vector_to_compare.push_back(given_grid + up_right_back);

        }else if (given_grid == this->max_number_of_grids - back + this->max_column_per_row_grid){

            vector_to_compare.push_back(given_grid + left);
            vector_to_compare.push_back(given_grid + up);
            vector_to_compare.push_back(given_grid + up_left);
            vector_to_compare.push_back(given_grid + back);
            vector_to_compare.push_back(given_grid + up_back);
            vector_to_compare.push_back(given_grid + left_back);
            vector_to_compare.push_back(given_grid + up_left_back);

        }else if (given_grid == this->max_number_of_grids - this->max_column_per_row_grid + 1){

            vector_to_compare.push_back(given_grid + right);
            vector_to_compare.push_back(given_grid + down);
            vector_to_compare.push_back(given_grid + down_right);
            vector_to_compare.push_back(given_grid + back);
            vector_to_compare.push_back(given_grid + down_back);
            vector_to_compare.push_back(given_grid + right_back);
            vector_to_compare.push_back(given_grid + down_right_back);

        }else if (given_grid == this->max_number_of_grids){

            vector_to_compare.push_back(given_grid + left);
            vector_to_compare.push_back(given_grid + down);
            vector_to_compare.push_back(given_grid + down_left);
            vector_to_compare.push_back(given_grid + back);
            vector_to_compare.push_back(given_grid + down_back);
            vector_to_compare.push_back(given_grid + left_back);
            vector_to_compare.push_back(given_grid + down_left_back);

        }else if(given_grid % this->max_column_per_row_grid == 1 && given_grid < front){

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

        }else if(given_grid % this->max_column_per_row_grid == 0 && given_grid < front){

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
                 given_grid > (front - this->max_column_per_row_grid)){

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

        }else if(given_grid > 1 && given_grid < this->max_column_per_row_grid){

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

        }else if(given_grid % this->max_column_per_row_grid == 1 && given_grid > this->max_number_of_grids - front){

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

        }else if(given_grid % this->max_column_per_row_grid == 0 && given_grid > this->max_number_of_grids - front){

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

        }else if(given_grid > this->max_number_of_grids - front && given_grid < this->max_number_of_grids - front + this->max_column_per_row_grid){ //ez volt hibás

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

        }else if(given_grid > this->max_number_of_grids - back + 1 && given_grid < this->max_number_of_grids - back + this->max_column_per_row_grid){

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

        }else if(given_grid % front == this->max_column_per_row_grid){

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

        }else if(given_grid % front == front - this->max_column_per_row_grid + 1){

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
            if(elem <= this->max_number_of_grids && elem > 0){
                return_set.insert(elem);
            }
        }

        return return_set;
    }
};