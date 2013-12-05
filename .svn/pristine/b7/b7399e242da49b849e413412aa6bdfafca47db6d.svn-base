//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LagrangeElement.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fe_lagrangeelement_hpp
#define base_fe_lagrangeelement_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>
#include <base/meta.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace fe{

        //----------------------------------------------------------------------
        //! Lagrangian DIM-Simplex of polynomial DEGREE.
        template<unsigned DEGREE, unsigned DIM>
        struct LagrangeSimplex
        {
            static const unsigned numDoFsPerEdge =
                base::Binomial<DEGREE-1,base::EDGE>::value;
            static const unsigned numDoFsPerFace =
                (DIM > 1? base::Binomial<DEGREE-1,base::FACE>::value : 0);
            static const unsigned numDoFsPerCell =
                (DIM > 2? base::Binomial<DEGREE-1,base::CELL>::value : 0);

        };

        //----------------------------------------------------------------------
        //! Lagrangian DIM-HyperCube of polynomial DEGREE.
        template<unsigned DEGREE, unsigned DIM>
        struct LagrangeHyperCube
        {
            static const unsigned numDoFsPerEdge =
                base::MToTheN<DEGREE-1,base::EDGE>::value;
            static const unsigned numDoFsPerFace =
                ( DIM > 1 ? base::MToTheN<DEGREE-1,base::FACE>::value : 0 );
            static const unsigned numDoFsPerCell =
                ( DIM > 2 ? base::MToTheN<DEGREE-1,base::CELL>::value : 0 );
            
        };

        //----------------------------------------------------------------------
        /** Representation of attributes of a Lagrangian element with given
         *  polynomial degree and underlying geometric shape.
         *  \tparam SHAPE  Geometric shape of element
         *  \tparam DEGREE Polynomial degree of the shape functions
         */
        template<base::Shape SHAPE, unsigned DEGREE>
        struct LagrangeElement
            : public base::IfElse< SHAPE == base::HyperCubeShape<
                                       base::ShapeDim<SHAPE>::value
                                       >::value, //!< Condition 'is HyperCube'
                                   LagrangeHyperCube<  //!< IF-result
                                           DEGREE, base::ShapeDim<SHAPE>::value>, 
                                   LagrangeSimplex<    //!< ELSE-result
                                           DEGREE, base::ShapeDim<SHAPE>::value>
                                   >::Type
        {
            //! Number of dofs per vertex = 1
            static const unsigned numDoFsPerVertex = 1;

            //------------------------------------------------------------------
            // For legibilty redefine the base class type
            static const base::Shape shape = SHAPE;
            static const unsigned dim   = base::ShapeDim<shape>::value;
            static const bool condition = (shape == base::HyperCubeShape<dim>::value);
            typedef LagrangeHyperCube<DEGREE,dim> IfType;
            typedef LagrangeSimplex<  DEGREE,dim> ElseType;

            typedef typename
            base::IfElse< condition, IfType, ElseType >::Type  LagrangeShape;

            //------------------------------------------------------------------
            //! @name Number of basic geometric entities
            //@{
            static const unsigned numVertices =
                base::NumNFaces<shape,base::VERTEX>::value;
            static const unsigned numEdges    =
                base::NumNFaces<shape,base::EDGE>::value;
            static const unsigned numFaces    =
                base::NumNFaces<shape,base::FACE>::value;
            static const unsigned numCells    = (dim == 3 ? 1 : 0);
            //@}

            //! @name Total number of DoFs on every sub-face category
            //@{
            static const unsigned numVertexDoFs =
                numDoFsPerVertex * numVertices;
            
            static const unsigned numEdgeDoFs =
                LagrangeShape::numDoFsPerEdge * numEdges;

            static const unsigned numFaceDoFs =
                LagrangeShape::numDoFsPerFace * numFaces;

            static const unsigned numCellDoFs =
                LagrangeShape::numDoFsPerCell * numCells;
            //@}

            //! Total number of DoFs
            static const unsigned numTotalDoFs =
                numVertexDoFs + numEdgeDoFs + numFaceDoFs + numCellDoFs;

            //! Shape function associated with this element
            typedef base::LagrangeShapeFun<DEGREE,shape>    ShapeFun;
        };

        //----------------------------------------------------------------------
        //! Specialisation for a constant shape function: just one DoF
        template<base::Shape SHAPE>
        struct LagrangeElement<SHAPE,0>
        {
            static const base::Shape shape = SHAPE;
            
            static const unsigned dim   = base::ShapeDim<shape>::value;

            static const unsigned numDoFsPerVertex = (dim==0 ? 1 : 0);
            static const unsigned numDoFsPerEdge   = (dim==1 ? 1 : 0);
            static const unsigned numDoFsPerFace   = (dim==2 ? 1 : 0);
            static const unsigned numDoFsPerCell   = (dim==3 ? 1 : 0);

            //------------------------------------------------------------------
            //! @name Number of basic geometric entities
            //@{
            static const unsigned numVertices =
                base::NumNFaces<shape,base::VERTEX>::value;
            static const unsigned numEdges    =
                base::NumNFaces<shape,base::EDGE>::value;
            static const unsigned numFaces    =
                base::NumNFaces<shape,base::FACE>::value;
            static const unsigned numCells    = (dim == 3 ? 1 : 0);
            //@}

            //! @name Total number of DoFs on every sub-face category
            //@{
            static const unsigned numVertexDoFs =
                numDoFsPerVertex * numVertices;
            
            static const unsigned numEdgeDoFs =
                numDoFsPerEdge * numEdges;

            static const unsigned numFaceDoFs =
                numDoFsPerFace * numFaces;

            static const unsigned numCellDoFs =
                numDoFsPerCell * numCells;
            //@}

            //! Total number of DoFs
            static const unsigned numTotalDoFs =
                numVertexDoFs + numEdgeDoFs + numFaceDoFs + numCellDoFs;

            //! Shape function associated with this element
            typedef base::LagrangeShapeFun<0,shape>    ShapeFun;
                                                      
        };

    }
}

#endif
