#ifndef fluidsolid_generatemesh_h
#define fluidsolid_generatemesh_h

//------------------------------------------------------------------------------
#include <fstream>
#include <base/verify.hpp>
#include <base/io/smf/Reader.hpp>

#include "../generateMesh.hpp"


//------------------------------------------------------------------------------
//! \ingroup thomas
template<typename MESH, typename UIN>
void generateMesh( MESH& mesh, const UIN& userInput ) 
{
    if ( userInput.readMeshFromFile ) {
        std::ifstream smf( userInput.meshFile.c_str() );
        VERIFY_MSG( smf.is_open(),
                    x2s("Cannot find input file ") + userInput.meshFile );
        base::io::smf::readMesh( smf, mesh );
    }
    else {
        ::generateMesh<UIN::dim>( mesh, userInput.N, userInput.bbmin, userInput.bbmax );
    }
}


#endif
