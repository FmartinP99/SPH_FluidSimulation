#include <iostream>

class SPHObject{
private:
    float weight;
    float pos_x;
    float pos_y;
    float pos_z;
    float speed_x;
    float speed_y;
    float speed_z;

public:

    SPHObject() {

    }

    SPHObject(float weight, float posX, float posY, float posZ, float speedX, float speedY, float speedZ) : weight(
            weight), pos_x(posX), pos_y(posY), pos_z(posZ), speed_x(speedX), speed_y(speedY), speed_z(speedZ) {
        if (this->weight  <= 0){
            throw std::invalid_argument("The weight of the object must be bigger than 0!");
        }
    }

    virtual ~SPHObject() {
    }


    float get_weight() const {
        return weight;
    }

    void set_weight(float weight) {
        if (weight <= 0){
            throw std::invalid_argument("The weight of the object must be bigger than 0!");
        }
        SPHObject::weight = weight;
    }

    float get_pos_x() const {
        return pos_x;
    }

    void set_pos_x(float posX) {
        pos_x = posX;
    }

    float get_pos_y() const {
        return pos_y;
    }

    void set_pos_y(float posY) {
        pos_y = posY;
    }

    float get_pos_z() const {
        return pos_z;
    }

    void set_pos_z(float posZ) {
        pos_z = posZ;
    }

    float get_speed_x() const {
        return speed_x;
    }

    void set_speed_x(float speedX) {
        speed_x = speedX;
    }

    float get_speed_y() const {
        return speed_y;
    }

    void set_speed_y(float speedY) {
        speed_y = speedY;
    }

    float get_speed_z() const {
        return speed_z;
    }

    void set_speed_z(float speedZ) {
        speed_z = speedZ;
    }
};

