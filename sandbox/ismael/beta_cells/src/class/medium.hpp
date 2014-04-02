/* MEDIUM CLASS HEADER */

// Guard
#ifndef _H_MEDIUM_
#define _H_MEDIUM_

// Standard Libraries
#include <string>
#include <iostream>


class Medium {

private:
// Private Variables
	double volume;	// Medium volume (mL)
	double c_cs; 	// Carbon Source Concentration (mg/mL)
	double c_o2;	// Oxygen Concentration (mg/mL)
	
public:
// Public Methods
	// Constructor
	Medium();
	// Destructor
	virtual ~Medium();

	// Print Medium Info
	void printInfo(); 

	// Getters
		// Get concentration
		double getConc ( int comp ); 	// Carbon Source, Oxygen

	// Parameters modification methods
		// Set Concentrations
		void setConc ( double conc_cs, double conc_o2 );
};


#endif
