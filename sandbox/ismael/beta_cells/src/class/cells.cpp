/* CELLS CLASS DEFINITION */

// Header
#include "cells.hpp"

/* VARIABLES */

double Cell::cRate = 0;
// Force Model Constants
double Cell::A = -9.3184e-2;		//dyne	
double Cell::B = 6.10542e-3;		//dyne
double Cell::Eps_1 = 2.94;		//µm;
double Cell::Eps_2 = 5.88;		//µm
double Cell::fd = 10e-2;			//dyne*s/µm²

// Health Paremeters	
double Cell::concT = 0.5;		// Concentration Threshold
double Cell::hRate = 3;		// Healing Rate (hp/h)
double Cell::dRate = 5;		// Damage Rate (hp/h)

// Counters
int Cell::counter = 1;
int Cell::dcounter = 0;

/* CONSTRUCTOR & DESTRUCTOR */

// Constructor
Cell::Cell()
{
	// Set Default Parameters
	id = counter;
	radius = 10;
	cluster = 0;
	status = "alive";
	// Position
	position[0] = 0;	//x
	position[1] = 0;	//y	
	position[2] = 0;	//z
	// Forces
	force[0] = 0;		//x
	force[1] = 0;		//y
	force[2] = 0;		//z
	// Velocity
	velocity[0] = 0;
	velocity[1] = 0;
	velocity[2] = 0;
	// Cell Health
	health = 100.0;

	// Update Counter
	counter++;		
};

// Destructor
Cell::~Cell(){};


/* DATA MANIPULATION METHODS */

// Print Cell Info
void Cell::printInfo()
{
	std::cout << "\n·Cell: " << id << std::endl;
	std::cout << "\tRadius: " << radius << " µm" << std::endl;
	std::cout << "\tCluster: " << cluster << std::endl;
	std::cout << "\tStatus: " << status <<  std::endl;
	std::cout << "\tHealth: " << health << std::endl;
	// Print Position
		std::cout << "\tPosition: [  ";
		for(int i=0; i < 3 ; i++) {
			std::cout << position[i] << "  ";
		};
		std::cout << "] (µm)" << std::endl;
	// Print Force
		std::cout << "\tForce:    [  ";
		for(int i=0; i < 3 ; i++) {
			std::cout << force[i] << "  ";
		};
		std::cout << "] (dyne)" << std::endl;
	// Print Velocity
		std::cout << "\tVelocity: [  ";
		for(int i=0; i < 3 ; i++) {
			std::cout << velocity[i] << "  ";
		};
		std::cout << "] (µm/s)" << std::endl;
};

// Get Status
int Cell::getStatus()
{
	if ( this -> status == "dead" ) { return 0; };
	if ( this -> status == "alive" ) { return 1; }; 
	if ( this -> status == "replicating" ) { return 2; }; 
};


// Set Position
void Cell::setPosition ( double x, double y, double z )
{
	position[0]= x;
	position[1]= y;
	position[2]= z;
};  


// Duplicate Cell
// Size growth




/* CALCULATION METHODS  */

// Reset Force value to zero
void Cell::resetForce()
{
	for ( int i=0; i < 3; ++i )
	{
		force[i] = 0;
	};
};

// Force Calculation
void Cell::calculateForce ( Cell* otherCell )
{
	double 	force_cc;			// Cell to Cell interaction Force
	//double 	force_p				// Propulsion Force
	double 	force_module;		// Force module
	double	rPos_module;		// Relative Position module
	double 	rPos[3];			// Relative Position vector
	double 	n_vector[3];		// Normalized vector
	
	// Relative Position between cells
	for ( int i=0; i < 3; ++i )
	{
			rPos[i] = ( otherCell -> position[i] ) - ( this -> position[i] );
	};
	
	// Normalize relative position vector
	rPos_module = sqrt ( ( rPos[0] * rPos[0] ) + ( rPos[1] * rPos[1] ) 
						+ ( rPos[2] * rPos[2] ) );

	if ( rPos_module == 0 ) 
	{
		for ( int i=0; i < 3; ++i ){
			n_vector[i] = 0;
		};
	}
	else 
	{
		for ( int i=0; i < 3; ++i )
		{
			n_vector[i] = rPos[i] / rPos_module;
		};
		
		// Check if cell is in a cluster
		if ( rPos_module <= ( 1.1 * ( this -> radius ) + 
			 ( otherCell -> radius ) ) )
		{
			this -> cluster += 1;
		}

	}

	// Force Calculation
	force_cc = ( A * exp ( - ( rPos_module  / Eps_1 ) ) + 
				 B * exp ( - ( rPos_module  / Eps_2 ) ) );  
	force_module = force_cc;
	
	// Store force & velocity
	for ( int i=0; i < 3; ++i )
	{
		this -> force[i] =  force[i] + ( force_module * n_vector[i] );
	};
};

// Velocity calculation
void Cell::calculateVelocity()
{
	for ( int i=0; i < 3; ++i )
	{
		velocity[i] = force[i] / ( fd * radius );
	};
};

// New position
void Cell::newPosition( double time_step, int* size)
{
	double posTemp[3];
	for ( int i=0; i < 3; ++i )
	{
		posTemp[i] =  velocity[i] * time_step + position[i];
		if ( ( (posTemp[i]) - double (size[i]) ) > 1e-8 )
		{
			position[i] = size[i] - radius;
		}
		else
		{
			if ( posTemp[i] < 1e-8 )
			{
				position[i] = 0 + radius;
			}
			else
			{
				position[i] = posTemp[i];
			}
		}
	};
};


/* HEALTH SYSTEM */
// Update Health
void Cell::updateHealth ( double time, double concentration )
{
	if ( concentration >= concT )
	{
		health = health + ( time / 3600 ) * hRate; 
		if ( health >= ( 100.0 + 1E-7 ) ) { health = 100.0; }
	}
	else
	{
		health = health - ( time / 3600 ) * dRate;
	}

}

// Check viability
void Cell::checkHealth()
{
	double cHealth = this -> health;
	if ( cHealth <= 1E-7 )
	{
		this -> killCell();
		dcounter++;
	}
}


