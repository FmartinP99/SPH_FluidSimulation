#include <iostream>
#include <string>
#include "SPHParticle.cpp"


class Container{
private:
    unsigned width;
    unsigned height;
    unsigned depth;
    unsigned particles_number;
    float SCALE = 0.01;
    SPHParticle* particles;
public:

    Container(unsigned int width, unsigned int height, unsigned int depth, unsigned int particlesNumber) : width(width),
                                                                                                           height(height),
                                                                                                           depth(depth),
                                                                                                           particles_number(
                                                                                                                   particlesNumber),
                                                                                                                   particles(new SPHParticle[particlesNumber]){
        if(this->width <= 0 || this->height <=0 || this->depth <=0 || this->particles_number <=0){
            throw std::invalid_argument("Every value of the container's attribute must be bigger than 0!");
        }
        fill_container();
    }

    virtual ~Container() {
        delete[] particles;
    }

    unsigned int get_width() const {
        return width;
    }

    void set_width(unsigned int width) {
        if(width <=0) throw std::invalid_argument("The width of the container must be bigger than 0!");
        Container::width = width;
    }

    unsigned int get_height() const {
        return height;
    }

    void set_height(unsigned int height) {
        if(height <=0) throw std::invalid_argument("The height of the container must be bigger than 0!");
        Container::height = height;
    }

    unsigned int get_depth() const {
        return depth;
    }

    void set_depth(unsigned int depth) {
        if(depth <=0) throw std::invalid_argument("The depth of the container must be bigger than 0!");
        Container::depth = depth;
    }

    unsigned int get_particles_number() const {
        return particles_number;
    }

    SPHParticle *get_particles() const {
        return particles;
    }

private:
    void fill_container(){
        unsigned z_sor = 0;
        unsigned y_sor = 0;
        int modulus = this->width * (1/SCALE*0.1);
        for (int i = 0; i<this->particles_number; i++){
            SPHParticle particle_to_add = SPHParticle(1, i % modulus * 0.1, y_sor * 0.1, z_sor * 0.1, 0, 0, 0, 0.01);
            std::cout<<particle_to_add.get_pos_y()<<std::endl;
            this->particles[i] = particle_to_add;
            if(i % modulus == modulus - 1) z_sor++;
            if(z_sor % modulus == modulus - 1 && i % modulus == modulus - 1) {
                y_sor++;
                z_sor = 0;
            }
        }
    }
};