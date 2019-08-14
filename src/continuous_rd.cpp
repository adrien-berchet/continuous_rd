#include "continuous_rd_device.h"

int main( int argc , char* argv[] )
{
	Parameters params;
	params.read(argc, argv);
	cout<<params<<endl;

	// std::filesystem::create_directories(param.result_folder);

	simulate_rd(params);

}

/*int main(int argc, const char *argv[])
{
	// Read parameters

	// Define variables
	int nx, ny;

	// Create initial state

	// Compute evolution

	// Export results
}*/
