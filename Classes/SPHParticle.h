#ifndef SPH_FLUIDSIMULATION_SPHPARTICLE_H
#define SPH_FLUIDSIMULATION_SPHPARTICLE_H
#include <iostream>

class SPHParticle {
public:
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


public:
    SPHParticle& operator =(SPHParticle elem2);
    SPHParticle();
    SPHParticle(double weight, double posX, double posY, double posZ, double speedX, double speedY, double speedZ,
                double radius, double mass);
    double get_radius() const;
    void set_radius(double radius);
    double get_mass() const;
    void set_mass(double mass);
    double get_pressure() const;
    double get_density() const;
    double get_fx() const;
    double get_fy() const;
    double get_fz() const;
    double get_weight() const;
    void set_weight(double weight);
    double get_px() const;
    void set_px(double posX);
    double get_py() const;
    void set_py(double posY);
    double get_pz() const;
    void set_pz(double posZ);
    double get_vx() const;
    void set_vx(double speedX);
    double get_vy() const;
    void set_vy(double speedY);
    double get_vz() const;
    void set_vz(double speedZ);
    double get_ax() const;
    void set_ax(double _ax);
    double get_ay() const;
    void set_ay(double _ay);
    double get_az() const;
    void set_az(double _az);
};

#endif //SPH_FLUIDSIMULATION_SPHPARTICLE_H
