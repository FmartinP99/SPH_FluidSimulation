#include <iostream>
#include <chrono>
#include <vector>

using namespace std::chrono;

class Timer{
private:
    long long START_TIME = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
    std::vector<long long>* t1;
    std::string name;

public:
    Timer(std::vector<long long>* t1){
        this->t1 = t1;
        START_TIME = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
        name = "";
    }

    Timer(std::string name){
        START_TIME = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
        this->name = name;
    }
    ~Timer(){
        long long END_TIME = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count() - START_TIME;
        //std::cout<<"( "<<this->name<<" )"<<"The time it took in microseconds: "<< END_TIME<<std::endl;
        t1->push_back(END_TIME);
    }
};
