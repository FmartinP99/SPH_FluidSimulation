cmake_minimum_required(VERSION 3.21)
project(SPH_FluidSimulation)

set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_C_FLAGS -m64)
set(CMAKE_CXX_FLAGS -m64)

set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")

link_directories("C:\\Program\ Files\ (x86)\\AMD APP SDK\\3.0\\lib\\x86_64")
include_directories("C:\\Program\ Files\ (x86)\\AMD APP SDK\\3.0\\include")
include_directories("I:\\Szakdoga\\MINGW-Freeglut\\freeglut\\lib\\x64")

link_directories("I:/Szakdoga/MSVC-Freeglut/freeglut")
include_directories("I:/Szakdoga/MSVC-Freeglut/freeglut/include")



add_executable(SPH_FluidSimulation main.cpp Classes/Container.h Classes/Container.cpp Classes/SPHParticle.cpp Classes/Container_temp.cpp Classes/Timer.cpp CL/cl.hpp CL/opencl.hpp Classes/OpenCLContainer.cpp Classes/OpenCLContainer.h Classes/SPHParticle.h)

target_link_libraries(SPH_FluidSimulation -lOpenGL32 -lfreeglut -lglu32 I:\\Szakdoga\\GLUT\\glut-amd\\lib\\x86_64\\OpenCL.lib)
