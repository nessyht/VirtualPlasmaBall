// CUDA Runtime, Interop, and includes
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>

// CUDA utilities
#include <helper_cuda.h>

// Helper functions
#include <helper_functions.h>
#include <helper_timer.h>

// CUDA Poisson functions
// Plasma ball
#include "Plasma_Ball.h"

Plasma_Ball Ball;

int *pArgc;
char **pArgv;

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	
	pArgc = &argc;
	pArgv = argv;
	
	double real_size = 0.15;
	double Vapp = -10000;

	int number_of_filaments = 10;
	double branching = 0.05;
	int stepsize = 3;
	int res = 128;
	int iterations;
	std::string name;
	int n;

	if (checkCmdLineFlag(argc, (const char **)argv, "fils"))
	{
		n = getCmdLineArgumentInt(argc, (const char **)argv, "fils");
		number_of_filaments = n;
	}

	if (checkCmdLineFlag(argc, (const char **)argv, "branching"))
	{
		n = getCmdLineArgumentInt(argc, (const char **)argv, "branching");
		branching = n;
	}

	if (checkCmdLineFlag(argc, (const char **)argv, "stepsize"))
	{
		n = getCmdLineArgumentInt(argc, (const char **)argv, "stepsize");
		stepsize = n;
	}

	if (checkCmdLineFlag(argc, (const char **)argv, "res"))
	{
		n = getCmdLineArgumentInt(argc, (const char **)argv, "res");
		res = n;
	}

	if (checkCmdLineFlag(argc, (const char **)argv, "name"))
	{
		name = getCmdLineArgumentInt(argc, (const char **)argv, "name");
	}

	if (checkCmdLineFlag(argc, (const char **)argv, "its"))
	{
		n = getCmdLineArgumentInt(argc, (const char **)argv, "its");
		iterations = n;
	}

	Ball = Plasma_Ball(CIRCULAR, real_size, res, Vapp);
	Ball.set_grad_descent(true); 
	Ball.name_dir(name);
	Ball.init_DBM(number_of_filaments, stepsize, branching, 3, 10000);

	for (int i = 0; i < iterations; ++i)
		Ball.step(false);

return 0;
}
