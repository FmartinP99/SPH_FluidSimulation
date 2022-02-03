cmake_minimum_required(VERSION 3.21)
project(SPH_FluidSimulation)

set(CMAKE_CXX_STANDARD 20)

add_executable(SPH_FluidSimulation main.cpp Classes/Container.cpp Classes/SPHObject.cpp Classes/SPHParticle.cpp)
target_link_libraries(SPH_FluidSimulation -lOpenGL32 -lfreeGLUT -lglu32)
