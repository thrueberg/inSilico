//[tut01Includes]{
#include <iostream>
#include <fstream>
#include <string>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/SurfaceElement.hpp>
#include <base/io/smf/Reader.hpp>
//[tut01Includes]}

//[tut2Includes]{
#include <base/Quadrature.hpp>
#include <base/geometry.hpp>
#include <base/io/PropertiesParser.hpp>
//[tut2Includes]}

//------------------------------------------------------------------------------
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
int main( int argc, char * argv[] )
{
    const unsigned    geomDeg  = 2; 
    const base::Shape shape    = base::TET;
    const unsigned degEstimate = 4;

    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    
    //PP
    std::string meshFile, surfFile;
    double radius;
    {    
        //! Feed properties parser with the variables to read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile", meshFile );
        prop.registerPropertiesVar( "surfFile", surfFile );
        prop.registerPropertiesVar( "radius",   radius   );

        //! Read variables from the input.dat file
        const std::string inputFile = "./input.dat";
        std::ifstream inp( inputFile.c_str()  );
        VERIFY( inp.is_open() );
        prop.readValues( inp );
        inp.close( );
    }
    
    // Input volume mesh
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // compute volume of the mesh
    {
        typedef base::Quadrature<degEstimate,shape> Quadrature;
        Quadrature quadrature;

        double volume = 0.;
        Metric<Mesh::Element,Quadrature> volumeInt( quadrature, volume );

        std::for_each( mesh.elementsBegin(), mesh.elementsEnd(), volumeInt );
        
        std::cout << "Computed volume = " << volume
                  << "  ---  " 
                  << "Exact volume    = " << 4. /3. * M_PI * radius * radius * radius
                  << '\n';
    }



    // Input surface mesh
    typedef base::mesh::SurfaceElement<Mesh::Element>   SurfaceElement;
    typedef base::mesh::Unstructured<SurfaceElement>    SurfaceMesh;
    SurfaceMesh surfaceMesh;
    {
        std::ifstream smf( surfFile.c_str() );
        base::io::smf::readMesh( smf, surfaceMesh );
        smf.close();
    }

    // compute area of the surface mesh
    {
        typedef base::SurfaceQuadrature<degEstimate,shape> SurfaceQuadrature;
        SurfaceQuadrature surfaceQuadrature;
    
        double area = 0.;
        Metric<SurfaceElement,SurfaceQuadrature> areaInt( surfaceQuadrature, area );
    
        std::for_each( surfaceMesh.elementsBegin(), surfaceMesh.elementsEnd(), areaInt );
    
        std::cout << "Computed area   = " << area
                  << "  ---  " 
                  << "Exact area      = " << 4. * M_PI * radius * radius
                  << '\n';
    }

    return 0;
}
