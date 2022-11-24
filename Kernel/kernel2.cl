

__kernel void calculate_physics(__global double* data)
{
double da = data[8] + data[9];
printf("%f\n", da);
barrier(CLK_GLOBAL_MEM_FENCE);
}