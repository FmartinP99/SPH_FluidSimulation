#include "SPHObject.cpp"

class SPHParticle : public SPHObject{
private:
    double radius;
    double mass;
public:
    SPHParticle() {
    }

    SPHParticle(double weight, double posX, double posY, double posZ, double speedX, double speedY, double speedZ,
                double radius, double mass) : SPHObject(weight, posX, posY, posZ, speedX, speedY, speedZ), radius(radius), mass(mass) {
        if (this->radius  <= 0){
            throw std::invalid_argument("The radius of the object must be bigger than 0!");
        }
    }

    virtual ~SPHParticle() {
    }

    double get_radius() const {
        return radius;
    }

    void set_radius(double radius) {
        if (radius <= 0){
            throw std::invalid_argument("The radius of the object must be bigger than 0!");
        }
        SPHParticle::radius = radius;
    }

    double get_mass() const {
        return mass;
    }

    void set_mass(double mass) {
        SPHParticle::mass = mass;
    }


};
