#ifndef generatemesh_h
#define generatemesh_h

//------------------------------------------------------------------------------
#include <sstream>
#include <base/verify.hpp>
#include <base/io/smf/Reader.hpp>


//------------------------------------------------------------------------------
// Generate a mesh with N elements on the interval [a,b]
template<typename MESH>
void generateMesh( MESH& mesh, const unsigned N, const double a, const double b )
{
    VERIFY_MSG( b > a, "Right end must be right of left end" );
    
    static const unsigned degree = MESH::Element::GeomFun::degree;
    
    std::stringstream buffer;
    // write header
    buffer << "! elementShape line \n"
           << "! elementNumPoints  " << degree+1 << "\n";

    // total number of elements
    const unsigned numElements = N;
    // total number of nodes
    const unsigned numNodes = 1 + degree * numElements;

    buffer << numNodes << "  " << numElements << "\n";

    // element lengths
    const double h = (b-a) / static_cast<double>( N );

    // generate coordinates
    double x = a;
    for ( unsigned n = 0; n < N; n++ ) {
        for ( unsigned d = 0; d < degree; d++ ) {
            buffer << x << "  0.  0. \n";
            x += h;
        }
    }

    // last node
    buffer << x << "  0.  0. \n";
    
    // generate connectivity
    unsigned v = 0;
    for ( unsigned n = 0; n < numElements; n++ ) {
        for ( unsigned d = 0; d < degree+1; d++ ) {
            buffer << v+d << " ";
        }
        buffer << "\n";
        v += degree;
    }

    base::io::smf::readMesh( buffer, mesh );

    return;
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const FIELD&   displacement,
                   //const FIELD&   velocity,
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
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
    //base::io::vtk::writePointData( vtkWriter, mesh, velocity,     "veloc" );
    vtk.close();
}

#endif
