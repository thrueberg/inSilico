/* FILE MANAGEMENT DEFINITIONS */

#include "filem.hpp"

// Write VTK File
void filem::writeVTK ( Pool* pool, int step )
{
	//Create file with proper numeration
	char filename[64];
	std::ofstream file_step;
	if ( step < 10 )
	{
		std::snprintf ( filename, sizeof(char) * 64, 
						"./results/cells/cells_000%i.vtk",	step );
	}
	if (( step < 100 ) && ( step >= 10 ))
	{
		std::snprintf ( filename, sizeof(char) * 64, 
						"./results/cells/cells_00%i.vtk", step );
	}
	if (( step < 1000 ) && ( step >= 100 ))
	{
		std::snprintf ( filename, sizeof(char) * 64, 
						"./results/cells/cells_0%i.vtk", step );
	}
	if ( step >= 1000 )
	{
		std::snprintf ( filename, sizeof(char) * 64, 
						"./results/cells/cells_%i.vtk", step );
	}
	file_step.open ( filename );

	//Create header and set points positions
	int numCells;
	numCells = pool -> cells.size();
	file_step << "# vtk DataFile Version 2.0" << "\n"
			  << "Beta Cells Discrete model - Step:" << step << "\n"
			  << "ASCII" << "\n"
			  << "DATASET UNSTRUCTURED_GRID" << "\n"
         	  << "POINTS " << numCells << " float" << "\n";
	
	for ( int i=0; i < numCells; ++i )
	{
		double* position = pool -> cells[i].getPosition();
		for ( int j=0; j < 3; ++j )
		{
			file_step << position[j] << " ";
		}
		file_step << "\n";
	}
	file_step << "\n";
	
	//Scalars and Vectors
	file_step << "POINT_DATA " << numCells << "\n";
	filem::writeVTKScalar ( pool, file_step, numCells ); 
	filem::writeVTKVector ( pool, file_step, numCells ); 
	

	file_step.close();
    return;
};

void filem::writeVTKScalar ( Pool* pool, std::ofstream &file_step, 
							 int numCells )
{
	//ID
	file_step << "SCALARS id int 1" << "\n" 
			  << "LOOKUP_TABLE default" << "\n";
	for ( int i=0; i < numCells; ++i )
	{
		file_step << pool -> cells[i].getID() << "\n";
	}; 
	
	//Radius
	file_step << "SCALARS radius float 1" << "\n" 
			  << "LOOKUP_TABLE default" << "\n";

	for ( int i=0; i < numCells; ++i )
	{
		file_step << pool -> cells[i].getRadius() << "\n"; 
	};
	
	//Cluster
	file_step << "SCALARS cluster int 1" << "\n"
			  << "LOOKUP_TABLE default" << "\n";
	for ( int i=0; i < numCells; ++i )
	{
		file_step << pool -> cells[i].getCluster() << "\n"; 
	};
	
	//Status
	file_step << "SCALARS status int 1" << "\n"
			  << "LOOKUP_TABLE default" << "\n";
	for ( int i=0; i < numCells; ++i )
	{
		file_step << pool -> cells[i].getStatus() << "\n";
	};

	//Health
	file_step << "SCALARS health float 1" << "\n" 
			  << "LOOKUP_TABLE default" << "\n";

	for ( int i=0; i < numCells; ++i )
	{
		file_step << pool -> cells[i].getHealth() << "\n"; 
	};
	
};
void filem::writeVTKVector ( Pool* pool, std::ofstream &file_step, 
							 int numCells )
{
	//Forces
	file_step << "VECTORS forces double" << "\n";
	for ( int i=0; i < numCells; ++i )
	{
		double* force = pool -> cells[i].getForce();
		for ( int j=0; j < 3; ++j )
		{
			file_step << force[j] << " ";
		}
	file_step << "\n";	
	} 
	
	//Velocity
	file_step << "VECTORS velocity double" << "\n";
	for ( int i=0; i < numCells; ++i )
	{
		double* velocity = pool -> cells[i].getVelocity();
		for ( int j=0; j < 3; ++j )
		{
			file_step << velocity[j] << " ";
		}
	file_step << "\n";
	} 
	
};

