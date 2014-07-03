#ifndef writevtk_h
#define writevtk_h

#include <string>
#include <vector>
#include <fstream>

#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>
#include <base/shape.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/post/evaluateAtNodes.hpp>

//------------------------------------------------------------------------------
//! Combine basename, a step number and the vtk-suffix to a string
std::string makeVTKFileName( const std::string& baseName,
                             const unsigned step )
{
    const std::string vtkFile = baseName + "." + 
        base::io::leadingZeros( step ) + ".vtk";
    return vtkFile;
}

//------------------------------------------------------------------------------
//! Maintain stream in order to allow for initialiser construction of VTK-Writer
struct HoldStream
{
    HoldStream( const std::string& name ) : bla( name.c_str() ) { }

    std::ofstream& access() { return bla; }
    
    std::ofstream bla;
};

//------------------------------------------------------------------------------
//! Convenience structure for VTK file writing
template<typename MESH>
class VTKWriter
{
public:
    //--------------------------------------------------------------------------
    VTKWriter( const MESH& mesh,
               const std::string& baseName,
               const unsigned     step )
        : mesh_( mesh ),
          holdStream_( makeVTKFileName( baseName, step ) ),
          vtkWriter_(  holdStream_.access() ) 
    {
        vtkWriter_.writeUnstructuredGrid( mesh_ );
    }

    //--------------------------------------------------------------------------
    template<typename FIELD>
    void writeField( const FIELD& field, const std::string& name )
    {
        base::io::vtk::writePointData( vtkWriter_, mesh_, field, name );
    }

    //--------------------------------------------------------------------------
    template<typename ITER>
    void writePointData( ITER first, ITER last, const std::string& name )
    {
        vtkWriter_.writePointData( first, last, name );
    }

    //--------------------------------------------------------------------------
    template<typename ITER>
    void writeCellData( ITER first, ITER last, const std::string& name )
    {
        vtkWriter_.writeCellData( first, last, name );
    }

    //--------------------------------------------------------------------------
    template<typename FIELD>
    void writePreviousField( const FIELD& field, const std::string& name )
    {
        typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDof;
        std::vector<VecDof> nodalValues;
        base::post::evaluateFieldHistoryAtNodes<1>( mesh_, field, nodalValues );
        vtkWriter_.writePointData( nodalValues.begin(), nodalValues.end(), name );
    }


    //--------------------------------------------------------------------------
    void writeDistances(
        const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet,
        const std::string& name = "distance" )
    {
        
        std::vector<double> distances;
        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( distances ),
                        boost::bind(
                            &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );
        vtkWriter_.writePointData( distances.begin(), distances.end(), name );
    }

    //--------------------------------------------------------------------------
    template<typename FIELD>
    void writeDoFStatus( const FIELD& field, const std::string& name )
    {
        typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDof;
        std::vector<VecDof> nodalValues;
        
        base::post::evaluateAtNodes(
            mesh_, field, 
            boost::bind( base::dof::Status<typename FIELD::DegreeOfFreedom>::apply, _1, _2 ),
            nodalValues );
        
        vtkWriter_.writePointData( nodalValues.begin(), nodalValues.end(), name );
    }

    
private:
    const MESH&                 mesh_;
    HoldStream                  holdStream_;
    base::io::vtk::LegacyWriter vtkWriter_;
};



#endif
