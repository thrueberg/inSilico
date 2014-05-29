//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   generateCutCells.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_generatecutcells_hpp
#define base_cut_generatecutcells_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base  includes
#include <base/linearAlgebra.hpp>
#include <base/mesh/HierarchicOrder.hpp>
// base/cut includes
#include <base/cut/LevelSet.hpp>
#include <base/cut/Cell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename MESH, typename CELL>
        void generateCutCells(
            const MESH& mesh,
            const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet,
            std::vector<CELL> & cutCells,
            const bool   update   = false, 
            const double epsilon  = 1.e-12 );

        //----------------------------------------------------------------------
        namespace detail_{

            //! Proxy for access of signed distance function
            template<bool ISSURFACE> struct GetSignedDistance;

            //! Get from nodal data (volume element)
            template<>
            struct GetSignedDistance<false>
            {
                template<typename ELEMENT,typename LSVEC>
                static double apply( const ELEMENT* ep,
                                     const unsigned v,
                                     const LSVEC& levelSet )
                {
                    const std::size_t nodeID =
                        (ep -> nodePtr(v)) -> getID();
                    
                    return levelSet[ nodeID ].getSignedDistance();
                }
            };

            //! Evaluate at parameter values of the nodes (surface element)
            template<>
            struct GetSignedDistance<true>
            {
                template<typename ELEMENT,typename LSVEC>
                static double apply( const ELEMENT* ep,
                                     const unsigned v,
                                     const LSVEC& levelSet )
                {
                    typename ELEMENT::ParamConstIter pIter = ep -> parametricBegin();
                    std::advance( pIter, v );

                    return
                        base::cut::signedDistance( ep -> getDomainElementPointer(),
                                                   *pIter, 
                                                   levelSet );
                }                

            };


        } // end namespace detail_
    }
}

//------------------------------------------------------------------------------
/** Function to generate cut-cell structures based on a given level set.
 *  Get for every element of the domain mesh the (nodal) values of the level set
 *  function. These are then used in order to construct the volume and surface
 *  triangulations which recover the implicit domain in the cut cells.
 *
 *  \warning The level set data is retrieved via the index of the nodes of the
 *           domain mesh. This will not work in case of a structured mesh with
 *           higher-order B-splines and a multi-index approach would be
 *           necessary.
 *
 *  \tparam MESH  Type of domain mesh
 *  \tparam CELL  Type of cell for the cut-cell structure
 *  \param[in]  mesh     Access to the domain mesh
 *  \param[in]  levelSet All the level set data
 *  \param[out] cutCells The cells representing elements (inside, outside, cut)
 *  \param[in]  update   If true, the level set data is added to existing cell
 *  \param[in]  epsilon  Small interval around zero to avoid degenerate cases
 */
template<typename MESH, typename CELL>
void base::cut::generateCutCells(
    const MESH& mesh,
    const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet,
    std::vector<CELL> & cutCells,
    const bool   update, 
    const double epsilon  )
{
    // number of vertices of a cell
    static const unsigned numNodes = CELL::numElementNodes;

    // needs re-ordering if these conditions are fulfilled
    static const bool reorder =
        ( (CELL::shape ==
           base::HyperCubeShape<base::ShapeDim<CELL::shape>::value>::value) and
          (MESH::Element::GeomFun::ordering == base::sfun::LEXICOGRAPHIC) );

    // treatment of a surface mesh (actually the level-set dim should be used)
    static const bool isSurface = (MESH::Node::dim != CELL::dim);

    typename MESH::ElementPtrConstIter eIter = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter eEnd  = mesh.elementsEnd();

    // make space for the cells
    if ( update )
        VERIFY_MSG( (std::distance(eIter,eEnd) == cutCells.size()),
                    "Expect existing cut-cell structures" );
    else cutCells.resize( std::distance( eIter, eEnd ) );

    // go through all elements of the mesh
    for ( std::size_t elemNum = 0; eIter != eEnd; ++eIter, elemNum++ ) {

        // this element's array of signed distances
        boost::array<double,numNodes> signedDistances;

        // get signed distances via node IDs
        for ( unsigned v = 0; v < numNodes; v++ ) {
            
            signedDistances[ v ] =
                detail_::GetSignedDistance<isSurface>::apply( *eIter, v,
                                                              levelSet );
        }

        // a distance of exactly (!) zero is not good
        for ( std::size_t s = 0; s < signedDistances.size(); s++ ) {
            double sd = signedDistances[s];
            if ( std::abs( sd ) < epsilon ) {
                if ( sd < 0. ) sd -= epsilon;
                else           sd += epsilon;

                signedDistances[ s ] = sd;
            }
        }

        // if asked to update and cell was cut 
        if ( update ) {
            if ( cutCells[ elemNum ].isCut() )
                cutCells[ elemNum ].update( signedDistances );
        }
        else cutCells[ elemNum ].destroy();


        // if necessary, reorder from lexicographic to hierarchic 
        if ( reorder ) {

            static const base::Shape hcShape =
                base::HyperCubeShape<base::ShapeDim<CELL::shape>::value>::value;
            typedef base::mesh::HierarchicOrder<hcShape,CELL::geomDegree> Hierarchic;
            
            boost::array<double,numNodes> tmp;

            for ( unsigned v = 0; v < numNodes; v++ ) {
                const unsigned hier = Hierarchic::apply( v );
                //tmp[ v ] = signedDistances[ hier ];
                tmp[ hier ] = signedDistances[ v ];
            }

            signedDistances = tmp;
        }

        // create the cut cell structure
        if ( cutCells[ elemNum ].isInside() )
            cutCells[ elemNum ].create( signedDistances );

    }

    return;
}

#endif
