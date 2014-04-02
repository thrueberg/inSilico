/* MEDIUM CLASS DEFINITION */ 

// Header
#include "medium.hpp"


// Constructor
Medium::Medium()
{
	volume = 5.0;
	c_cs = 1.0;
	c_o2 = 1.0; 
};
// Destructor
Medium::~Medium(){};

void Medium::printInfo()
{
	std::cout << "·Medium:" << std::endl;
	std::cout << "\tVolume: " << volume << " mL" << std::endl;
	std::cout << "\t[Carbon Source]: " << c_cs << " mg/mL" << std::endl;
	std::cout << "\t[Oxygen]: " << c_o2 << " mg/mL" << std::endl;
};

// Getters
// Get Concentrations
double Medium::getConc ( int comp )
{
	switch ( comp )
	{
		case 1:
			return c_cs; // Carbon source
		case 2:
			return c_o2; // Oxygen
	};
};

// Parameters modification methods
// Set Concentration
void Medium::setConc ( double conc_cs, double conc_o2 )
{
	c_cs = conc_cs;
	c_o2 = conc_o2;
};
