/* SCAFFOLD CLASS DEFINITION */ 

// Header
#include "scaffold.hpp"



/* CONSTRUCTOR & DESTRUCTOR */

// Constructor
Scaffold::Scaffold()
{
	size[0] = 5000;
	size[1] = 5000;
	size[2] = 1000;
	porosity = 0.5;
	Dmin = 1.0;
	Dmax = 10.0;
};

// Destructor
Scaffold::~Scaffold(){};



/* DATA MANIPULATION METHODS */ 

// Print Scaffold Info
void Scaffold::printInfo()
{
	std::cout << "\n·Scaffold:" << std::endl;		
	// Print Size 
		std::cout << "\tSize: [  ";
		for(int i=0; i < 3 ; i++) {
			std::cout << size[i] << "  ";
		};
		std::cout << "] µm" << std::endl;
	
	std::cout << "\tPorosity: "<< porosity << std::endl;		
	std::cout << "\tMin. Diff. Constant:" << Dmin << std::endl;
	std::cout << "\tMax. Diff. Constant:" << Dmax << std::endl;
}

//Getters
//Get Diffusion Constants
double Scaffold::getDiffConst ( int choice )
{
	switch ( choice )
	{
		case 1:
			return Dmax;
		case 2:
			return Dmin;
	};
};

// Parameter Modification Methods
// Set Size
void Scaffold::setSize( int x, int y, int z )
{
	size[0]= x;
	size[1]= y;
	size[2]= z;
}; 

// Set Porosity
void Scaffold::setPorosity ( float n_porosity )
{
	porosity = n_porosity;
};
// Set Diffusion constants
void Scaffold::setDiffConst ( double nDmin, double nDmax )
{
	Dmin = nDmin;
	Dmax = nDmax;
}; 
