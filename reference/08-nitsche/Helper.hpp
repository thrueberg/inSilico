//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Helper.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef apps_nitsche_helper_h_hpp
#define apps_nitsche_helper_h_hpp

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <string>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/shape.hpp>
#include <base/Structured.hpp>
#include <base/Unstructured.hpp>
// base/mesh includes
#include <base/mesh/MeshBoundary.hpp>
// base/io includes
#include <base/io/sgf/Reader.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

//------------------------------------------------------------------------------
namespace apps{
    namespace nitsche{
    
        template<unsigned DIM,unsigned GEOMDEG,unsigned FIELDDEG>
        struct StructuredHelper;

        template<unsigned DIM,unsigned GEOMDEG,unsigned FIELDDEG>
        struct UnstructuredHelper;
    }
}


//------------------------------------------------------------------------------
// Helper for structured grid based on B-splines
template<unsigned DIM,unsigned GEOMDEG,unsigned FIELDDEG>
struct apps::nitsche::StructuredHelper
{
    // Shape has to be a hypercube
    static const base::Shape shape = base::HyperCubeShape<DIM>::value;

    // A structured mesh
    typedef base::Structured<DIM,GEOMDEG> Mesh;

    // Input file suffix
    static std::string suffix() { return ".sgf"; }

    // Read input file and fill grid
    static void readFromFile( std::istream& in, Mesh& mesh )
    {
        base::io::sgf::readGrid( in, mesh );
    }

    // Basis based on B-Splines
    typedef base::fe::Basis<shape,FIELDDEG,base::BSPLINE> FEBasis;

    // Box-boundary of the grid
    static void meshBoundary( base::mesh::MeshBoundary &meshBoundary,
                              const Mesh& mesh )
    {
        meshBoundary.create<DIM>( mesh.gridSizes() );
    }

    // Write a VTK file
    template<typename FIELD>
    static void writeVTK( std::ostream & vtk, Mesh& mesh, FIELD& field )
    {
        static const unsigned doFSize = FIELD::DegreeOfFreedom::size;
        
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeStructuredGrid( mesh );
        {
            // Evaluate the solution field at every geometry node
            std::vector<typename base::Vector<doFSize>::Type> nodalValues;

            base::post::sampleField( mesh, field, 
                                     std::back_inserter( nodalValues ) );
            vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "heat" );
        }

        {
            // Evaluate the solution field at every geometry node
            std::vector<typename base::Matrix<DIM,doFSize>::Type> cellValues;

            base::post::sampleFieldGradient( mesh, field, 
                                             std::back_inserter( cellValues ), 0 );
            vtkWriter.writeCellData( cellValues.begin(), cellValues.end(), "flux" );
        }
    }
};

//------------------------------------------------------------------------------
// Helper for structured grid based on Lagrangian Simplex elements
template<unsigned DIM,unsigned GEOMDEG,unsigned FIELDDEG>
struct apps::nitsche::UnstructuredHelper
{
    // Use simplex elements
    static const base::Shape shape = base::SimplexShape<DIM>::value;

    // Unstructured mesh
    typedef base::Unstructured<shape,GEOMDEG> Mesh;

    // SMF is the suffix of the input file
    static std::string suffix() { return ".smf"; }

    // Read mesh file and fill mesh
    static void readFromFile( std::istream& in, Mesh& mesh )
    {
        base::io::smf::readMesh( in, mesh );
    }

    // Basis with Lagrangian shape functions
    typedef base::fe::Basis<shape,FIELDDEG,base::LAGRANGE> FEBasis;

    // Extract faces which are not duplicates
    static void meshBoundary( base::mesh::MeshBoundary &meshBoundary,
                              const Mesh& mesh )
    {
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
    }

    // Write a VTK file
    template<typename FIELD>
    static void writeVTK( std::ostream & vtk, Mesh& mesh, FIELD& field )
    {
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );

        base::io::vtk::writePointData( vtkWriter, mesh, field, "heat" );

        const typename base::Vector<DIM>::Type xi =
            base::ShapeCentroid<shape>::apply();
        
        base::io::vtk::writeCellData( vtkWriter, mesh, field,
                                      boost::bind(
                                          base::post::template
                                          evaluateFieldGradient<typename Mesh::Element,
                                          typename FIELD::Element>, _1, _2, xi ), "flux" );
    }

};


#endif
