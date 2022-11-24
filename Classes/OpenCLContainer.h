#ifndef SPH_FLUIDSIMULATION_OPENCLCONTAINER_H
#define SPH_FLUIDSIMULATION_OPENCLCONTAINER_H
#define __CL_ENABLE_EXCEPTIONS
#include "../CL/opencl.hpp"
#include <iostream>
#include "Container.h"

class OpenCLContainer {
private:
    cl::Program program;
    cl::Context context;
    std::vector<cl::Device> devices;
    cl::Device device;

    void programCreate(const std::string &file, const Container &hostBuffer, float dtime, unsigned int max_particles);
    void getError(const cl::Program&, int err);
    static auto getSource(std::string const& fileName);

public:
    unsigned int max_particles;
    OpenCLContainer(const std::string& file, Container &hostBuffer, float dtime, unsigned int max_particles);
    cl::Program& getProgram();
    cl::Context& getContext();
    std::vector<cl::Device>& getDevices();
    cl::Device& getDevice();

    cl::Kernel updateKernel;
    cl::Buffer inputBuffer;
    cl::Buffer dtimeBuffer;
    cl::Buffer particleBuffer;
    cl::Buffer enviromentarrayBuffer;
    cl::Buffer testBuffer;
    cl::CommandQueue queue;
};


#endif //SPH_FLUIDSIMULATION_OPENCLCONTAINER_H
