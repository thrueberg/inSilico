#ifndef thomas_generatemesh_h
#define thomas_generatemesh_h

//------------------------------------------------------------------------------
#include <sstream>
#include <base/verify.hpp>
#include <base/io/smf/Reader.hpp>

#include <tools/meshGeneration/unitCube/unitCube.hpp>
#include <tools/converter/smfMap/smfAffine.hpp>

//------------------------------------------------------------------------------
/** \ingroup thomas
 *  Generate a structured hyper-cube mesh in the domain
 *  \f$ [a_i, b_i]^d \f$ with
 *  \f$ N_i \f$ elements per direction \f$ x_i \f$. 
 */
template<unsigned DIM, typename MESH>
void generateMesh(
    MESH& mesh,
    const typename base::Vector<DIM,unsigned>::Type&  N,
    const typename base::Vector<DIM,double  >::Type&  a,
    const typename base::Vector<DIM,double  >::Type&  b )
{
    static const unsigned dim = DIM;

    for ( unsigned d = 0; d < dim; d++ )
        VERIFY_MSG( b[d] > a[d], "Right end must be right of left end" );
    
    VERIFY_MSG( MESH::Element::GeomFun::degree == 1,
                "Only linear elements currently supported by unitCube" );

    // number of elements per direction
    const unsigned e1 = N[0];
    const unsigned e2 = (dim > 1? N[1] : 1);
    const unsigned e3 = (dim > 2? N[2] : 1);

    // Type of hypercube mesh generator for [0,1]^DIM
    static const bool simplex = false;
    typedef tools::meshGeneration::unitCube::SMF<dim,simplex> SMF;

    // use buffer
    std::stringstream buffer;
    SMF::apply( e1, e2, e3, buffer );

    // apply an affine transformation to convert [0,1]^DIM to [a_i,b_i]^DIM
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            tools::converter::smfMap::A(d1,d2) = (d1 == d2 ?  b[d1]-a[d1] : 0.0);
        }
        tools::converter::smfMap::c[d1] = a[d1];
    }

    // another buffer
    std::stringstream buffer2;
    tools::converter::smfMap::Affine<SMF::shape,1>::apply( buffer, buffer2 );

    // read mesh from second buffer
    base::io::smf::readMesh( buffer2, mesh );

    return;
}
#endif
