#include "SPHParticle.h"

    SPHParticle& SPHParticle::operator =(SPHParticle elem2){
        this->set_weight(elem2.get_weight());
        this->set_px(elem2.get_px());
        this->set_py(elem2.get_py());
        this->set_pz(elem2.get_pz());
        this->set_vx(elem2.get_vx());
        this->set_vy(elem2.get_vy());
        this->set_vz(elem2.get_vz());
        this->radius = elem2.get_radius();
        this->mass = elem2.get_mass();
        this->pressure = elem2.get_pressure();
        this->density = elem2.get_density();
        this->fx = elem2.get_fx();
        this->fy = elem2.get_fy();
        this->fz = elem2.get_fz();

        return *this;
    }

SPHParticle::SPHParticle() {
    }

SPHParticle::SPHParticle(double weight, double posX, double posY, double posZ, double speedX, double speedY, double speedZ,
                double radius, double mass ) : weight(weight), px(posX), py(posY), pz(posZ), vx(speedX), vy(speedY), vz(speedZ), radius(radius), mass(mass) {
        pressure = 0;
        density = 0;
        fx = 1.0;
        fy = 1.0;
        fz = 1.0;
    }

    double SPHParticle::get_radius() const {
        return radius;
    }

    void SPHParticle::set_radius(double radius) {
        SPHParticle::radius = radius;
    }

    double SPHParticle::get_mass() const {
        return mass;
    }

    void SPHParticle::set_mass(double mass) {
        SPHParticle::mass = mass;
    }

    double SPHParticle::get_pressure() const {
        return pressure;
    }

    double SPHParticle::get_density() const {
        return density;
    }

    double SPHParticle::get_fx() const {
        return fx;
    }

    double SPHParticle::get_fy() const {
        return fy;
    }

    double SPHParticle::get_fz() const {
        return fz;
    }

    double SPHParticle::get_weight() const {
        return weight;
    }

    void SPHParticle::set_weight(double weight) {
        this->weight = weight;
    }

    double SPHParticle::get_px() const {
        return px;
    }

    void SPHParticle::set_px(double posX) {
        px = posX;
    }

    double SPHParticle::get_py() const {
        return py;
    }

    void SPHParticle::set_py(double posY) {
        py = posY;
    }

    double SPHParticle::get_pz() const {
        return pz;
    }

    void SPHParticle::set_pz(double posZ) {
        pz = posZ;
    }

    double SPHParticle::get_vx() const {
        return vx;
    }

    void SPHParticle::set_vx(double speedX) {
        vx = speedX;
    }

    double SPHParticle::get_vy() const {
        return vy;
    }

    void SPHParticle::set_vy(double speedY) {
        vy = speedY;
    }

    double SPHParticle::get_vz() const {
        return vz;
    }

    void SPHParticle::set_vz(double speedZ) {
        vz = speedZ;
    }


