#include "SPHObject.cpp"

class SPHParticle : public SPHObject{
private:
    float radius;
public:
    SPHParticle() {
    }

    SPHParticle(float weight, float posX, float posY, float posZ, float speedX, float speedY, float speedZ,
                float radius) : SPHObject(weight, posX, posY, posZ, speedX, speedY, speedZ), radius(radius) {
        if (this->radius  <= 0){
            throw std::invalid_argument("The radius of the object must be bigger than 0!");
        }
    }

    virtual ~SPHParticle() {
    }

    float get_radius() const {
        return radius;
    }

    void set_radius(float radius) {
        if (radius <= 0){
            throw std::invalid_argument("The radius of the object must be bigger than 0!");
        }
        SPHParticle::radius = radius;
    }
};
