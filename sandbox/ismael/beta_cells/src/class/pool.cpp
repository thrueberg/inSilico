/* POOL CLASS DEFINITION */

// Header
#include "pool.hpp"


/* CONSTRUCTOR & DESTRUCTOR */

// Constructor
Pool::Pool(){};
// Destructor
Pool::~Pool(){};


/* DATA MANIPULATION METHODS */

// Print Pool Info
void Pool::printInfo()
{
	// Cells Info
	int vectSize = cells.size();
	for ( int i=0; i < vectSize; ++i)
	{
		cells[i].printInfo();
	};
	// Scaffold Info
	scaffold.printInfo();
	// Medium Info
	medium.printInfo();
	// Pool Info
	std::cout << "\n·Pool:"<< std::endl;
	std::cout << "\tSize: "<< cells.size() << " cells" << std::endl;
};


// Add cell to the container 
void Pool::addCell ( int quantity )
{	
	for (int i=0; i < quantity; ++i){
		cells.push_back( Cell() );
	}
};


/* Randomize position of selected cells */

	// All cells in the container
	void Pool::randomPosition ()
	{		
		
		int vectSize = cells.size();
		
		for ( int i=0; i < vectSize; ++i )
		{
		// Increase thread duration for better randomization
			std::this_thread::sleep_for ( std::chrono::nanoseconds ( 100 ) );
		// Call randomizer for this cell
			this -> randomPosition ( i );
		};

	};

	// Cell selected
	void Pool::randomPosition ( int index )
	{	
		
		double x,y,z;
		int* size = scaffold.getSize();

		// Define generator and randomize seed with epoch time
		std::default_random_engine generator;
		
		typedef std::chrono::high_resolution_clock myclock;
		generator.seed ( myclock::now().time_since_epoch().count() );


		// Get numbers
		std::uniform_real_distribution<double> dist_x ( 0, size[0] );
		std::uniform_real_distribution<double> dist_y ( 0, size[1] );
		std::uniform_real_distribution<double> dist_z ( 0, size[2] );

		x = dist_x ( generator );
		y = dist_y ( generator ); 
		z = dist_z ( generator );
		
		// Set values in cell
		cells[index].setPosition ( x, y, z );
	};



/* Randomize radius of cells */

	// All cells in the container
	void Pool::randomSize ()
	{		
		int vectSize = cells.size();
		
		for ( int i=0; i < vectSize; ++i )
		{
			// Increase thread duration for better randomization
			std::this_thread::sleep_for ( std::chrono::nanoseconds ( 100 ) );
			// Call randomizer for this cell
			this -> randomSize ( i );
		};

	};

	// Cell selected
	void Pool::randomSize ( int index )
	{
		int n_radius;

		// Define generator and randomize seed with epoch time
		std::default_random_engine generator;
		
		typedef std::chrono::high_resolution_clock myclock;
		generator.seed ( myclock::now().time_since_epoch().count() );


		// Get numbers
		std::uniform_int_distribution<int> dist_r ( 6, 12 );

		n_radius = dist_r ( generator );
		
		// Set values in cell
		cells[index].setRadius ( n_radius );
	};



/* CALCULATION METHODS */

// Calculate force
void Pool::calculateForce ( )
{
	int vectSize = cells.size();
	
	#pragma omp parallel for num_threads(10) 
	for ( int i=0; i < vectSize; ++i )
	{
		this -> calculateForce ( i );
	};
};
void Pool::calculateForce ( int index )
{
	int vectSize = cells.size();

	this -> cells[index].resetForce();
	this -> cells[index].setCluster ( 0 );

	if ( cells[index].getStatus() != 0 )
	{
	for ( int i=0; i < vectSize; ++i )
	{
		if ( ( index != i ) && ( this -> cells[i].getStatus() != 0 ) )
		{
		this -> cells[index].calculateForce( &this -> cells[i] );	
		}
	};
	// Calculate velocity
	this -> cells[index].calculateVelocity();
	};
};

// Calculate new position
void Pool::newPosition()
{
	int vectSize = cells.size();
	
	#pragma omp parallel for num_threads(10)
	for ( int i=0; i < vectSize; ++i )
	{
		if ( this -> cells[i].getStatus() != 0 ) 
		{ 
			this -> newPosition ( i ); 
		}
	};
};

void Pool::newPosition( int index )
{
	this -> cells[index].newPosition ( time_step, scaffold.getSize() );
};




/* HEALTH SYSTEM */
	void Pool::procHealth ( int index, double concentration )
	{
		this -> cells[index].updateHealth ( time_step, concentration );
		this -> cells[index].checkHealth ();
	}


