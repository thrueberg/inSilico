/* SCAFFOLD CLASS HEADER */

// Guard
#ifndef _H_SCAFFOLD_
#define _H_SCAFFOLD_

// Standard Libraries
#include <string>
#include <iostream>


class Scaffold {

private:
// Private Variables
	int size[3];			// Scaffold Size µm
	float porosity;			// Scaffold Porosity
	double Dmin, Dmax;		// Max/Min Diff. Constant
public:
// Public Methods
	// Constructor
	Scaffold();
	// Destructor
	virtual ~Scaffold();

	// Print Scaffold Info
	void printInfo(); 

	// Getters
		// Get Size
		int* getSize () { return size; };
		// Get Diffusion constants
		double getDiffConst ( int choice ); // 1 = Dmin, 2 = Dmax

	// Parameter Modification Methods
		// Set size
		void setSize ( int x, int y, int z );
		// Set porosity
		void setPorosity ( float n_porosity );
		// Set diffusion constants
		void setDiffConst ( double nDmin, double nDmax );

};


#endif
