#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include "collapseMesh.hpp"

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>

int main( int argc, char* argv[] )
{
    const unsigned dim     = 2;
    const unsigned geomDeg = 1;
    const base::Shape shape = base::SimplexShape<dim-1>::value;
    const double tol  = 1.e-8;

    const std::string inFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string outFile = boost::lexical_cast<std::string>( argv[2] );

    std::ifstream smfIn(   inFile.c_str() );
    std::ifstream smfOut( outFile.c_str() );

    typedef base::Unstructured<shape,geomDeg,dim> Mesh;
    Mesh meshIn, meshOut;

    base::io::smf::readMesh(  smfIn,  meshIn );
    tfm::collapseMesh( meshIn, meshOut, tol );
    base::io::smf::writeMesh( meshOut, smfOut );
    
    return 0;
}
