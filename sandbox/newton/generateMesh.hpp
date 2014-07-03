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
    
    VERIFY_MSG( MESH::Element::GeomFun::degree == 1,
                "Only linear elements currently supported by unitCube" );

    static const unsigned dim = MESH::Node::dim;

    // number of elements per direction
    const unsigned e1 = N;
    const unsigned e2 = (dim > 1? N : 1);
    const unsigned e3 = (dim > 2? N : 1);

    // number of nodes per direction
    const unsigned n1 = N + 1;
    const unsigned n2 = (dim > 1? e2 + 1 : 1);
    const unsigned n3 = (dim > 2? e3 + 1 : 1);

    typedef tools::meshGeneration::unitCube::SMF<dim,false> SMF;
    
    std::stringstream buffer;
    SMF::apply( n1, n2, n3, e1, e2, e3, buffer );

    
    // apply an affine transformation to convert [0,1]^d to [a,b]^d
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            tools::converter::smfAffine::A(d1,d2) = (d1 == d2 ?  b-a : 0.0);
        }
        tools::converter::smfAffine::c[d1] = a;
    }

    std::stringstream buffer2;
    tools::converter::smfAffine::Converter<SMF::shape,1>::apply( buffer, buffer2 );
    
    base::io::smf::readMesh( buffer2, mesh );

    return;
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const FIELD&   displacement,
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet,
                   const std::vector<double>& supports )
{
    const std::string vtkFile = baseName + "." + 
        base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );

    std::vector<double> distances;
    std::transform( levelSet.begin(), levelSet.end(),
                    std::back_inserter( distances ),
                    boost::bind( &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );

    vtkWriter.writeUnstructuredGrid( mesh );
    vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
    base::io::vtk::writePointData( vtkWriter, mesh, displacement, "disp" );

    typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDof;
    std::vector<VecDof> nodalValues2;
    base::post::evaluateHistoryAtNodes<1>( mesh, displacement, nodalValues2 );
    vtkWriter.writePointData( nodalValues2.begin(), nodalValues2.end(), "prevDisp" );

    std::vector<double> doFStatus;
    typename FIELD::DoFPtrConstIter dIter = displacement.doFsBegin();
    typename FIELD::DoFPtrConstIter dEnd  = displacement.doFsEnd();
    for ( ; dIter != dEnd; ++dIter ) {

        bool isConstrained = false;
        for ( unsigned d = 0; d < FIELD::DegreeOfFreedom::size; d++ ) {
            if ( (*dIter) -> isConstrained(d) ) isConstrained = true;
        }

        if      (  isConstrained          ) doFStatus.push_back(  0. );
        else if ( (*dIter) -> isActive(0) ) doFStatus.push_back(  1. );
        else                                doFStatus.push_back( -1. );

    }

    if ( FIELD::Element::FEFun::degree == 1 ) {
        vtkWriter.writePointData( doFStatus.begin(), doFStatus.end(), "status" );
        vtkWriter.writePointData( supports.begin(), supports.end(), "support" );
    }

    
    vtk.close();
}

#endif
