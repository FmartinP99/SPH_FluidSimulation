#ifndef SPH_FLUIDSIMULATION_CONTAINER_H
#define SPH_FLUIDSIMULATION_CONTAINER_H
#include <iostream>
#include <string>
#include "SPHParticle.h"
#include <vector>
#include <set>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
class Container {
private:

    double damping_coeff = 0.25;
    unsigned int gas_constant = 400; //400
    double U = 1;

    unsigned int current_particles = 0;
    unsigned int current_iteration = 0;
    unsigned int new_particles_every_iteration = 3;

    double kernel_smoother_length = 1;
    double normalization_density = (315 / (64 * M_PI * kernel_pow(kernel_smoother_length, 9)));
    double norm_pressure = (-45 / (M_PI * kernel_pow(kernel_smoother_length, 6)));
    double norm_visc = (45 / (M_PI * kernel_pow(kernel_smoother_length, 6)));

    double width;
    double height;
    double depth;
    double gravity;
    double density;
    unsigned int particles_number;
    double particle_radius;

    SPHParticle *particles;

    unsigned int max_number_of_grids;
    unsigned int max_column_per_row_grid;
    unsigned int max_row_per_layer_grid;
    unsigned int max_layer_per_cube_grid;
    std::vector<std::set<int>> neighbouring_grids;  // megadja melyik  gridnek melyek a szomszédjai
    std::vector<std::set<int>> particles_on_grid; // megadja melyik griden melyik indexű particle-ek vannak


public:
    Container(double width, double height, double depth, double gravity, double density, unsigned int particlesNumber, double particle_radius);
    template<typename T> void print(T asd);
    virtual ~Container();
    double get_width() const;
    void set_width(unsigned int width);
    double get_height() const;
    void set_height(unsigned int height);
    double get_depth() const;
    void set_depth(unsigned int depth);
    unsigned int get_particles_number() const;
    SPHParticle *get_particles() const;
    unsigned int get_current_particles() const;
    double get_damping_coeff() const;
    unsigned int get_gas_constant() const;
    double get_u() const;
    double get_kernel_smoother_length() const;
    double get_normalization_density() const;
    double get_norm_pressure() const;
    double get_norm_visc() const;
    double get_gravity() const;
    double get_density() const;
    double get_particle_radius() const;
    unsigned int get_max_number_of_grids() const;
    unsigned int get_max_column_per_row_grid() const;
    unsigned int get_max_row_per_layer_grid() const;
    unsigned int get_max_layer_per_cube_grid() const;
    const std::vector<std::set<int>> &get_neighbouring_grids() const;
    const std::vector<std::set<int>> &get_particles_on_grid() const;
    void calculate_physics(double dtime);
    void fill_container_gradually2(const double radius);
    unsigned int get_current_iteration() const;
    void set_current_iteration(unsigned int iteration);
    unsigned int get_new_particles_every_iteration() const;

private:
    unsigned int calculate_grid(SPHParticle elem);
    unsigned int calculate_max_grid_number();
    static double calculate_absolute_distance_power2(SPHParticle& elem, SPHParticle& elem2);
    void check_boundaries(SPHParticle& elem) const;
public:
    void fill_container_gradually(double radius);
    void fill_container(double radius);
    std::set<int> calculate_neighboring_grids(int given_grid);
    double kernel_pow(double ertek, unsigned int order);
    void calculate_physicsv2(const double dtime);
    void calculate_physicsv3(double dtime, double* enviromentArray);

};

#endif //SPH_FLUIDSIMULATION_CONTAINER_H

