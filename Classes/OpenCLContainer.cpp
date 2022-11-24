//
// Created by treki on 2022. 09. 20..
//
#include <fstream>
#include <sstream>
#include "OpenCLContainer.h"

void OpenCLContainer::getError(const cl::Program &_program, int err) {
    if (err != 0)
        std::cout << _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
}

cl::Program &OpenCLContainer::getProgram() {
    return program;
}

cl::Context &OpenCLContainer::getContext() {
    return context;
}

std::vector<cl::Device> &OpenCLContainer::getDevices() {
    return devices;
}

cl::Device &OpenCLContainer::getDevice() {
    return device;
}


OpenCLContainer::OpenCLContainer(const std::string &file, Container& hostBuffer, float dtime, unsigned int max_particles) {
    programCreate(file, hostBuffer, dtime, max_particles);
}

auto OpenCLContainer::getSource(std::string const &fileName) {
    std::ifstream f{fileName};
    if (!f.is_open()) {
        std::cout << "Can't open file: " << fileName << std::endl;
        throw std::runtime_error{"Cannot open file"};
    }
    return std::string{std::istreambuf_iterator<char>{f}, std::istreambuf_iterator<char>{}};
}

void OpenCLContainer::programCreate(const std::string &file, const Container& hostBuffer, float dtime, unsigned int _max_particles){
    try {
        int err = CL_SUCCESS;
        int bufferSize = hostBuffer.get_particles_number();
        int particlesize = 52;
        int containersize = 92;

        std::vector<cl::Platform> platforms;
        err = cl::Platform::get(&platforms);
        std::cout << err << std::endl;

        if (platforms.empty()) {
            std::cout << "Unable to find suitable platform." << std::endl;
        } else {
            std::cout << "Platform count : " << platforms.size() << std::endl;
        }

        cl_context_properties properties[] =
                {CL_CONTEXT_PLATFORM, (cl_context_properties) (platforms[0])(), 0};
        cl::Context _context(CL_DEVICE_TYPE_GPU, properties);

        std::vector<cl::Device> _devices = _context.getInfo<CL_CONTEXT_DEVICES>();
        std::cout << "DEVICES : " << _devices.size() << std::endl;
        for (const auto &dev : _devices) {
            std::cout << "NAME : " << dev.getInfo<CL_DEVICE_NAME>() << std::endl;
        }

        auto programSource = getSource(file);
        cl::Program _program = cl::Program(_context, programSource, false, &err);

        getError(_program, err);

        // ERROR ITT ALATTA


        err = _program.build(_devices);

        getError(_program, err);
        std::cout << err << std::endl;


        this->program = _program;
        this->context = program.getInfo<CL_PROGRAM_CONTEXT>();
        this->devices = context.getInfo<CL_CONTEXT_DEVICES>();
        this->device = devices.front();

        cl::Kernel kernelUpdate{program, "calculate_physics", &err};


        cl::Buffer _particleBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, _max_particles * sizeof(SPHParticle), NULL,
                                                &err);

        cl::Buffer _dtimeBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double), NULL,
                                              &err);
        cl::Buffer _enviromentArrayBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double) * 21, NULL,
                                                       &err);

        cl::Buffer _containerBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(Container), NULL,
                                                 &err);

        cl::Buffer _testbuffer = cl::Buffer(context, CL_MEM_READ_WRITE, 10 * sizeof(double) , NULL,
                                                &err);


       //err = kernelUpdate.setArg(0, _testbuffer);
       err = kernelUpdate.setArg(0, _particleBuffer);
        err = kernelUpdate.setArg(1, _dtimeBuffer);
        //err = kernelUpdate.setArg(2, _containerBuffer);
        err = kernelUpdate.setArg(2, _enviromentArrayBuffer);



        this->updateKernel = kernelUpdate;
        this->inputBuffer = _containerBuffer;
        this->dtimeBuffer = _dtimeBuffer;
        this->particleBuffer = _particleBuffer;
        this->testBuffer = _testbuffer;
        this->enviromentarrayBuffer = _enviromentArrayBuffer;

        getError(program, err);

    }catch (cl::BuildError &e) {
        std::cout << "cl error was thrown" << std::endl;
        for (auto elem : e.getBuildLog()){
            std::cout << elem.second << std::endl;
        }
    }
    catch (cl::Error &e) {
        std::cout << "cl error was thrown" << std::endl;
        std::cout << e.what() << " : " << e.err();
    }
}