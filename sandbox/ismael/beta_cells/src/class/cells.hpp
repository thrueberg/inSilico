/* CELLS CLASS HEADER */

// Guard
#ifndef _H_CELLS_
#define _H_CELLS_

// Standard Libraries
#include <string>
#include <iostream>
#include <cmath>

class Cell {
private:
	// Cell Parameters 
	int id; 				// Cell id
	int radius;				// Cell radius
	int age;				// Cell age
	int cluster;			// Cluster	
	std::string status;		// Cell status(alive, dead, duplicating)
	double position[3];		// Cell position vector [x y z]
	double force[3];		// Force applied to this cell
	double velocity[3];		// Velocity
	double health;			// Health

	// Force Model Constants
	static double A;			//dyne	
	static double B;			//dyne
	static double Eps_1;		//µm;
	static double Eps_2;		//µm
	static double fd;			//dyne*s/µm²

	// Health Paremeters	
	static double concT;		// Concentration Threshold
	static double hRate;		// Healing Rate (hp/s)
	static double dRate;		// Damage Rate (hp/s)

	//Comsuption rate
	static double cRate;
	// Cell counter
	static int counter;
	// Dead cells counter
	static int dcounter;

public:
	// Constructor
	Cell();
	// Destructor
	virtual ~Cell();

	// Print Cell Info
	void printInfo();
	
	// Getters
		int getID() 			{ return id; };
		int getRadius() 		{ return radius; };
		int getAge() 			{ return age; };
		int getCluster() 		{ return cluster; };
		int getStatus();
		int getCounter() 		{ return counter; };
		int getDCounter() 		{ return dcounter; };
		double getHealth() 		{ return health; };
		double getRate() 		{ return cRate; };
		double* getPosition() 	{ return position; };
		double* getForce() 		{ return force; };
		double* getVelocity() 	{ return velocity; };

	// Data manipulation Methods
		// Set position
		void setPosition( double x, double y, double z );
		// Set size
		void setRadius ( int n_radius ) { radius = n_radius; };
		// Kill Cell
		void killCell() { status = "dead"; };
		// Set Cluster
		void setCluster( int n_cluster ) { cluster = n_cluster; };
		// Set Comsuption rate
		static void setRate ( double tRate ) { cRate = tRate; };
		// Set Nutrient Threshold
		static void setThreshold ( double nutrT ) { concT = nutrT; };
		// Set Health Rates
		static void setHPR ( double hR, double dR ) { hRate = hR; dRate = dR; };
		// Duplicate Cell
		void duplicateCell();
		// Size growth
		void growCell();

	// Calculation Methods
		// Reset Force value to zero
		void resetForce();
		// Forces Calculation 
		void calculateForce( Cell* otherCell );
		// Velocity Calculation
		void calculateVelocity();
		// Calculate new position 
		void newPosition( double time_step, int* size );

	// Health System
		// Update Health	
		void updateHealth ( double time, double concentration );
		// Check viability
		void checkHealth ();
};


#endif
