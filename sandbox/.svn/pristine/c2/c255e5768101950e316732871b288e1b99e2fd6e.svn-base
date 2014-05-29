#ifndef generatemesh_h
#define generatemesh_h

//------------------------------------------------------------------------------
#include <sstream>
#include <base/verify.hpp>
#include <base/io/smf/Reader.hpp>

#include <tools/meshGeneration/unitCube/unitCube.hpp>
#include <tools/converter/smfAffine/smfAffine.hpp>

//------------------------------------------------------------------------------
// Generate a mesh with N elements on the interval [a,b]
template<typename MESH>
void generateMesh( MESH& mesh, const unsigned N, const double a, const double b )
{
    VERIFY_MSG( b > a, "Right end must be right of left end" );
    
    static const unsigned dim = MESH::Node::dim;

    // number of elements per direction
    const unsigned e1 = N;
    const unsigned e2 = (dim > 1? N : 1);
    const unsigned e3 = (dim > 2? N : 1);

    typedef tools::meshGeneration::unitCube::SMF<dim,false,MESH::Element::GeomFun::degree> SMF;
    
    std::stringstream buffer;
    SMF::apply( e1, e2, e3, buffer );

    
    // apply an affine transformation to convert [0,1]^d to [a,b]^d
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            tools::converter::smfAffine::A(d1,d2) = (d1 == d2 ?  b-a : 0.0);
        }
        tools::converter::smfAffine::c[d1] = a;
    }

    std::stringstream buffer2;
    tools::converter::smfAffine::Converter<SMF::shape,
                                           MESH::Element::GeomFun::degree>::apply( buffer,
                                                                                   buffer2 );
    
    base::io::smf::readMesh( buffer2, mesh );

    return;
}


#endif
