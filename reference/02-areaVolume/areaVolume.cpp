#include <iostream>
#include <fstream>
#include <string>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/SurfaceElement.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/io/smf/Reader.hpp>

//
#include <base/Quadrature.hpp>
#include <base/geometry.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/Format.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/Measure.hpp>

//------------------------------------------------------------------------------
namespace ref02{

    template<unsigned QUADDEG, typename MESH> 
    std::pair<double,double> computeAreaAndVolume(
        const MESH&  mesh,
        const typename base::mesh::BoundaryMeshBinder<typename MESH::Element>::Type&
        boundaryMesh );

    
    int areaVolume();
}

//------------------------------------------------------------------------------
/** Compute the mesh volume and surface area for a given degree of quadrature.
 *  Given a volume mesh and the corresponding surface mesh, this function
 *  computes the measures of the volume and the surface area.
 *  The function demonstrates how define the right quadrature for a volume or
 *  a surface integral and how to perform a simple (geometric) integration by
 *  using the base::kernel::Measure measure-kernel.
 *
 *  \tparam QUADDEG  Degree of quadrature
 *  \tparam MESH     Type of volume mesh
 *  \param[in] mesh         Reference to a volume mesh
 *  \param[in] boundaryMesh Reference to a boundary mesh
 *  \return    A pair of the volume and the surface measure
 */
template<unsigned QUADDEG, typename MESH> 
std::pair<double,double> ref02::computeAreaAndVolume(
    const MESH&  mesh,
    const typename base::mesh::BoundaryMeshBinder<typename MESH::Element>::Type&
    boundaryMesh )
{
    const unsigned degEstimate = QUADDEG;
    // Introspect element shape
    const base::Shape shape    = MESH::Element::shape;
    
    // compute volume of the mesh
    double volume = 0.;
    {
        // define a volume quadrature
        typedef base::Quadrature<degEstimate,shape> Quadrature;
        Quadrature quadrature;

        // Bind the mesh
        typedef base::asmb::FieldBinder<const MESH> Binder;
        Binder binder( mesh );
        typedef typename Binder::template TupleBinder<>::Type TB;
        
        // call simple integration
        base::asmb::simplyIntegrate<TB>( quadrature, volume, binder,
                                         base::kernel::Measure<typename TB::Tuple>() );
    }

    // compute area of the surface mesh
    double area = 0.;
    {
        typedef base::SurfaceQuadrature<degEstimate,shape> SurfaceQuadrature;
        SurfaceQuadrature surfaceQuadrature;

        typedef 
            typename base::mesh::BoundaryMeshBinder<typename MESH::Element>::Type
            BoundaryMesh;
        typedef base::asmb::SurfaceFieldBinder<const BoundaryMesh>
            SurfaceBinder;
        SurfaceBinder surfaceBinder( boundaryMesh );
        typedef typename SurfaceBinder::template TupleBinder<>::Type STB;
        base::asmb::simplyIntegrate<STB>( surfaceQuadrature, area, surfaceBinder,
                                          base::kernel::Measure<typename STB::Tuple>() );
        
    }

    return std::make_pair( volume, area );
}

//------------------------------------------------------------------------------
/** Numerical integration.
 *  Other than in ref01::unstructured, no command line arguments are given here.
 *  It is assumed that a file called 'input.dat' is present in the working
 *  directory. Then the base::io::PropertiesParser is used to read and validate
 *  this file. Based on the read data, a mesh is created (base::Unstructured)
 *  and the function ref02::computeAreaAndVolume is called for various degrees
 *  of quadrature. The results are written to a base::io::Table.
 *
 *  Error in volume and area computation for various orders of polynomial
 *  exactness of the quadrature rules:
 *  \include measure.ref.dat
 *
 *  It is assumed that a mesh of a sphere with quadratic tetrahedra is used.
 *  Why does the error of the computed values not go to zero?
 */
int ref02::areaVolume()
{
    const unsigned    geomDeg  = 2; 
    const base::Shape shape    = base::TET;

    std::string meshFile, surfFile;
    double radius;
    {    
        // Feed properties parser with the variables to read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile", meshFile );
        prop.registerPropertiesVar( "radius",   radius   );

        // Read variables from the input.dat file
        const std::string inputFile = "./input.dat";
        std::ifstream inp( inputFile.c_str()  );
        VERIFY( inp.is_open() );
        VERIFY( prop.readValuesAndCheck( inp ) );
        inp.close( );
    }
    
    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    
    // Input volume mesh
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // extract the boundary of the volume mesh
    typedef
        base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create list of <Element,faceNo> pairs
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }


    std::cout << "# Computation of volume and area with a quadrature rule \n";

    // compute for various quadrature rules
    const std::pair<double,double> result1 = computeAreaAndVolume<1>( mesh, boundaryMesh );
    const std::pair<double,double> result2 = computeAreaAndVolume<2>( mesh, boundaryMesh );
    const std::pair<double,double> result3 = computeAreaAndVolume<3>( mesh, boundaryMesh );
    const std::pair<double,double> result4 = computeAreaAndVolume<4>( mesh, boundaryMesh );
    const std::pair<double,double> result5 = computeAreaAndVolume<5>( mesh, boundaryMesh );

    // write error to a table
    const unsigned L = 15;
    base::io::Table<6>::WidthArray widths = {{ 0, L, L, L, L, L }};
    base::io::Table<6> table( widths );
    table % "Degree" % 1 % 2 % 3 % 4 % 5;
    std::cout << table;

    const double vol = (4. /3. * M_PI * radius * radius * radius);

    table % "Volume"
        % std::abs(vol-result1.first) % std::abs(vol-result2.first)
        % std::abs(vol-result3.first) % std::abs(vol-result4.first)
        % std::abs(vol-result5.first);
    std::cout << table;

    const double area = (4. * M_PI * radius * radius);
    
    table % "Area"
        % std::abs(area-result1.second) % std::abs(area-result2.second)
        % std::abs(area-result3.second) % std::abs(area-result4.second)
        % std::abs(area-result5.second);
    std::cout << table;
    
    return 0;
}

//------------------------------------------------------------------------------
// Delegate the call 
int main(  )
{
    return ref02::areaVolume();
}
