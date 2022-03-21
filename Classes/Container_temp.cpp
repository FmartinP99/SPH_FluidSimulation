#include <iostream>
#include <string>
#include "SPHParticle.cpp"
#include <cmath>
#include <vector>
#include <set>


class Container{
private:

    double damping_coeff = 0.3;
    unsigned int gas_constant = 400;
    double kernel_smoother_length = 1;
    double constant_pressure = 1;
    double U = 1;

    unsigned int current_particles = 0;
    unsigned int current_iteration = 0;
    unsigned int new_particles_every_iteration = 10;

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
            for(int j = 0; j < this->current_particles; j++ ){
                if (i != j) {
                    SPHParticle &elem2 = this->particles[j];
                    double particle_distance = calculate_absolute_distance_power2(elem, elem2);
                    if (particle_distance < kernel_smoother_length * kernel_smoother_length && particle_distance > 0){
                        _vector_of_neighbour_pairs.emplace_back(j, sqrt(particle_distance));
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
            elem.set_px(elem.get_px() + dtime * elem.get_vx());
            elem.set_py(elem.get_py() + dtime * elem.get_vy());
            elem.set_pz(elem.get_pz() + dtime * elem.get_vz());

            for(int j = 0; j < this->current_particles; j++){
                if(i != j){
                    SPHParticle& elem2 = this->particles[j];
                    //   calculate_collided_particles(elem, elem2);
                }
            }


            check_boundaries(elem);

        }

        delete[] particle_densities;
        delete[] particle_forces;
        delete[] particle_pressures;
        //   delete[] particle_distances;


    }

private:

    unsigned int calculate_grid(SPHParticle elem){

        int current_column = int(elem.get_px() / this->kernel_smoother_length) + 1;
        int current_row = int(elem.get_py() / this->kernel_smoother_length) + 1;
        int current_layer = int(elem.get_pz() / this->kernel_smoother_length) + 1;

        return (this->max_column_per_row_grid * this->max_row_per_layer_grid * (current_layer - 1)) + (this->max_column_per_row_grid * (current_row - 1) + current_column);
    }

    unsigned int calculate_max_grid_number(){

        unsigned int column = this->width / kernel_smoother_length;
        unsigned int row = this->height / kernel_smoother_length;
        unsigned int layer = this->depth / kernel_smoother_length;

        //IF IT HAS A REMAINDER IT ADDS +1. FOR EXAMPLE  5/2 = 3 BECAUSE int(5/2) == 2

        if(column != this->width / kernel_smoother_length) column++;
        if(row != this->height / kernel_smoother_length) row++;
        if(layer != this->depth / kernel_smoother_length) layer++;

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
        double x_cor;
        double z_cor;
        int max_particle_per_row = this->width / kernel_smoother_length; // so they are chill
        int max_particle_per_column = this->depth / kernel_smoother_length; // so they are chill

        for(int i = 0; i < 20; i++) {
            this->particles[current_particles] = SPHParticle(1,  column * radius, this->height - 2 * radius, row * radius,
                                                             5, -10, 5, radius, 1);

            //   std::cout<<current_particles<<std::endl;

            this->current_particles++;
            if(this->current_particles == this->particles_number) break;

            this->particles[current_particles] = SPHParticle(1,  this->width - (column * radius), this->height - 2 * radius, this->depth - (row * radius),
                                                             -5, -10, -5, radius, 1);
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
};