//[tut01Includes]{
#include <iostream>
#include <fstream>
#include <string>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/SurfaceElement.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/io/smf/Reader.hpp>
//[tut01Includes]}

//[tut2Includes]{
#include <base/Quadrature.hpp>
#include <base/geometry.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/Format.hpp>
//[tut2Includes]}

//------------------------------------------------------------------------------
// Integrate Area/Volume of element and add result to sum storage
template<typename ELEMENT, typename QUAD>
class Metric
{
public:
    
    Metric( const QUAD&   quad,
            double& result )
    : quad_( quad ),
      result_( result )
    { }

    void operator()( const ELEMENT* ep )
    {
        typename QUAD::Iter qIter = quad_.begin();
        typename QUAD::Iter qEnd  = quad_.end();

        for ( ; qIter != qEnd; ++qIter ) {

            // quadrature weight and point
            const double            weight = qIter -> first;
            const typename QUAD::VecDim xi = qIter -> second;
            
            // Get Jacobian of element
            const double detJ =
                base::Jacobian<ELEMENT>()( ep, xi );

            result_ += detJ * weight;
        }
    }


private:
    const QUAD& quad_;
    double&     result_; 
};

//------------------------------------------------------------------------------
template<unsigned QUADDEG, typename MESH> 
std::pair<double,double> computeAreaAndVolume( const MESH&  mesh )
{
    const unsigned degEstimate = QUADDEG;
    const base::Shape shape    = MESH::Element::shape;
    
    // compute volume of the mesh
    double volume = 0.;
    {
        typedef base::Quadrature<degEstimate,shape> Quadrature;
        Quadrature quadrature;

        Metric<typename MESH::Element,Quadrature> volumeInt( quadrature, volume );

        std::for_each( mesh.elementsBegin(), mesh.elementsEnd(), volumeInt );
    }


    typedef typename
        base::mesh::BoundaryMeshBinder<typename MESH::Element>::Type BoundaryMesh;
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


    // compute area of the surface mesh
    double area = 0.;
    {
        typedef base::SurfaceQuadrature<degEstimate,shape> SurfaceQuadrature;
        SurfaceQuadrature surfaceQuadrature;
    
        Metric<typename BoundaryMesh::Element,
               SurfaceQuadrature> areaInt( surfaceQuadrature, area );
    
        std::for_each( boundaryMesh.elementsBegin(), boundaryMesh.elementsEnd(), areaInt );
    
    }

    return std::make_pair( volume, area );
}

//------------------------------------------------------------------------------
int main()
{
    const unsigned    geomDeg  = 2; 
    const base::Shape shape    = base::TET;

    //PP
    std::string meshFile, surfFile;
    double radius;
    {    
        //! Feed properties parser with the variables to read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile", meshFile );
        prop.registerPropertiesVar( "radius",   radius   );

        //! Read variables from the input.dat file
        const std::string inputFile = "./input.dat";
        std::ifstream inp( inputFile.c_str()  );
        VERIFY( inp.is_open() );
        prop.readValues( inp );
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

    std::cout << "# Computation of volume and area with a quadrature rule \n";
    
    const std::pair<double,double> result1 = computeAreaAndVolume<1>( mesh );
    const std::pair<double,double> result2 = computeAreaAndVolume<2>( mesh );
    const std::pair<double,double> result3 = computeAreaAndVolume<3>( mesh );
    const std::pair<double,double> result4 = computeAreaAndVolume<4>( mesh );
    const std::pair<double,double> result5 = computeAreaAndVolume<5>( mesh );

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
