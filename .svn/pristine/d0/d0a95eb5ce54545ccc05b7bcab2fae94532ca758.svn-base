//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfMap.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_converter_smfmap_h
#define tools_converter_smfmap_h

//------------------------------------------------------------------------------
// std includes
#include <iostream>
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
// base/io includes
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smfMap{

            //------------------------------------------------------------------
            // 3D vector
            typedef base::Vector<3,double>::Type   Vec;

            // Identity map
            struct Identity
            {
                Vec operator()( const Vec& x ) const
                {
                    return x;
                }
            };

            static boost::function<Vec( const Vec& )> coordinateMap = Identity();

            //------------------------------------------------------------------
            /** Read smf file, apply coordinate map and write new smf file
             *  \tparam SHAPE  Type of element shape
             *  \tparam DEGREE Polynomial degree of the element
             */
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Converter
            {
                static void apply( std::istream& smfIn,
                                   std::ostream& smfOut )
                {
                    // Attributes of the mesh
                    static const base::Shape shape   = SHAPE;
                    static const unsigned degree     = DEGREE;
                    static const unsigned    dim     = 3;
        
                    // Mesh type and object
                    typedef base::Unstructured<shape,degree,dim>  Mesh;
                    Mesh mesh;

                    // SMF input
                    base::io::smf::readMesh( smfIn, mesh );
                    
                    // go through nodes and modify coordinates
                    typename Mesh::NodePtrIter nIter = mesh.nodesBegin();
                    typename Mesh::NodePtrIter nEnd  = mesh.nodesEnd();
                    for (; nIter != nEnd; ++nIter ) {
                        Vec xOld;
                        (*nIter) -> getX( &(xOld[0]) );
                        const Vec xNew = coordinateMap( xOld );
                        (*nIter) -> setX( &(xNew[0]) );
                    }

                    // SMF output
                    base::io::smf::writeMesh( mesh, smfOut );
                }

            };


        } // namespace smfAffine
    } // namespace converter
} // namespace tools

#endif
