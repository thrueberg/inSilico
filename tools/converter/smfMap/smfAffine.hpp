//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfAffine.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_converter_smfaffine_h
#define tools_converter_smfaffine_h

//------------------------------------------------------------------------------
// std includes
#include <iostream>
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
            // 3D transformation matrix and translation vector
            typedef base::Matrix<3,3,double>::Type Mat;
            typedef base::Vector<3,double>::Type   Vec;

            static Mat A = Mat::Identity();
            static Vec c = Vec::Zero();
    
            //------------------------------------------------------------------
            /** Read smf file, apply affine transformation and write new smf file
             *  \tparam SHAPE  Type of element shape
             *  \tparam DEGREE Polynomial degree of the element
             */
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Affine
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
                        const Vec xNew = A * xOld + c;
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
