#include <string>
#include <vector>

// Cells Class
namespace cells {

class Cells {

private:
// Private Variables
	const long int id; 		//Cell id
	double radius;			//Cell radius
	long int age;			//Cell age
	//long int lineage;		//Cell lineage(string??)
	int cluster;			//Cluster	
	std::string status;		//Cell status(alive, dead, duplicating)
	
	long int position[3];		//Cell position vector
	//Relative Position
	//Force


public:
// Public Methods
	// Constructor
	Cells();
	// Destructor
	virtual ~Cells();

};

};//End of Namespace Cells
