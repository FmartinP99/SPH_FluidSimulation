
# REFERENCE ARTICLE

https://matthias-research.github.io/pages/publications/sca03.pdf <br>
**D. C. a. M. G. Matthias Müller, „Particle-Based Fluid Simulation for Interactive Applications,” ParticleBased Fluid Simulation for Interactive Applications, 2003.**
# PC 

I've ran this program on a **Ryzen 5 3600X** CPU and on a **RX 580 8GB** GPU.
Results are in the *results.csv* file.

#Program

This SPH Fluid Simulation is a bit different, than the others.<br>
The particles have no edges; they are merely points in space.<br>  So this means, since they can get really close to eachother, <br>
the numbers in the calcualtions can be really large or really small. To combat this, I had to use *double* types instead of *float*.<br><br>
The problem with *double* types is, that the modern consumable GPU's are not optimized for *double* operations, so it is way slower than it'd be with float operations.<br>

The non-naive CPU algorithm divides the 3D space into cubes.<br>
When the program calculates the force exerted on a particular particle by other particles,<br> it only considers the particles within its own cube and those in the neighboring cubes.<br>
Otherwise if the naive CPU algorithm (or GPU algorithm) is running, then it considers every other particle.
<br><br>
To run:

SPH_FluidSimulation.exe <particle_numbers> <container_width> <container_height>
<container_depth> <maximum_fps>
<naiv_algorithm> <gpu_algorithm> <opencl_kernel_path> <br>
<br>
Where:
- the <naiv_algorithm> and the <gpu_algorithm> are 1 or 0.
- the <opencl_kernel_path> is a string.
- the others are integers.

If the <naiv_algorithm> is a 1 it uses a primitive "bruteforce" algorithm, if 0 it uses a more advanced one.<br>
The more advanced one is only available on the CPU. <br>
If the <gpu_algorithm> is set to 1, then the physics calculation runs on the GPU and the <naiv_algorithm> is ignored.<br>
If the <opencl_kernel_path> is empty, then the default path will be "../Kernel/kernel.cl" (only necessary if the program<br>
runs on the GPU).