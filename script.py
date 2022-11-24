import os
import subprocess
import sys

particles = [2000, 2800, 4000, 5200, 6400, 8000, 12000]

for idx, particle in enumerate(particles): #400-asával növelni
    for is_gpu in range(2):
        if is_gpu == 0:
            for is_naive in range(2):
                command = f"start /wait I:\Szakdoga\SPH_FluidSimulation\SPH_FluidSimulation\cmake-build-release\SPH_FluidSimulation.exe {particle} {30 + idx * 10 } 20 5 120 {is_naive} {is_gpu}"
                subprocess.Popen(command, shell=True).wait()
        else:
            command = f"start /wait I:\Szakdoga\SPH_FluidSimulation\SPH_FluidSimulation\cmake-build-release\SPH_FluidSimulation.exe {particle} {30 + idx * 10 } 20 5 120 0 {is_gpu}"
            subprocess.Popen(command, shell=True).wait()

