#include <iostream>

class SPHObject{
private:
    double weight;
    double px;
    double py;
    double pz;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    unsigned int grid;

public:

    SPHObject() {

    }

    SPHObject(double weight, double posX, double posY, double posZ, double speedX, double speedY, double speedZ) : weight(
            weight), px(posX), py(posY), pz(posZ), vx(speedX), vy(speedY), vz(speedZ) {
        if (this->weight  <= 0){
            throw std::invalid_argument("The weight of the object must be bigger than 0!");
        }
    }

    virtual ~SPHObject() {
    }


    double get_weight() const {
        return weight;
    }

    void set_weight(double weight) {
        if (weight <= 0){
            throw std::invalid_argument("The weight of the object must be bigger than 0!");
        }
        SPHObject::weight = weight;
    }

    double get_px() const {
        return px;
    }

    void set_px(double posX) {
        px = posX;
    }

    double get_py() const {
        return py;
    }

    void set_py(double posY) {
        py = posY;
    }

    double get_pz() const {
        return pz;
    }

    void set_pz(double posZ) {
        pz = posZ;
    }

    double get_vx() const {
        return vx;
    }

    void set_vx(double speedX) {
        vx = speedX;
    }

    double get_vy() const {
        return vy;
    }

    void set_vy(double speedY) {
        vy = speedY;
    }

    double get_vz() const {
        return vz;
    }

    void set_vz(double speedZ) {
        vz = speedZ;
    }

    double get_ax() const {
        return ax;
    }

    void set_ax(double ax) {
        SPHObject::ax = ax;
    }

    double get_ay() const {
        return ay;
    }

    void set_ay(double ay) {
        SPHObject::ay = ay;
    }

    double get_az() const {
        return az;
    }

    void set_az(double az) {
        SPHObject::az = az;
    }


};

