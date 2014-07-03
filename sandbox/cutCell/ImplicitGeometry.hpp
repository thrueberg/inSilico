#ifndef cutcell_implicitgeometry_h
#define cutcell_implicitgeometry_h

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/bruteForce.hpp> 

//------------------------------------------------------------------------------
/** Implicit geometry representation.
 *  
 */
template<typename MESH>
class ImplicitGeometry
{
public:
    typedef MESH Mesh;
    typedef typename base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;

    static const unsigned    dim       = Mesh::Node::dim;
    static const base::Shape shape     = Mesh::Element::shape;
    static const base::Shape surfShape = BoundaryMesh::Element::shape;

    typedef base::cut::LevelSet<dim>       LevelSet;
    typedef base::cut::Cell<shape>         Cell;
    typedef base::cut::Cell<surfShape>     SurfCell;
    typedef typename base::Vector<dim,double>::Type VecDim;

    //--------------------------------------------------------------------------
    //! Store mesh references and create a default level set state (all in)
    ImplicitGeometry( const Mesh& mesh, const BoundaryMesh& boundaryMesh )
        : mesh_( mesh ), boundaryMesh_( boundaryMesh )
    {
        // generate a default level set state
        const std::size_t numNodes = std::distance( mesh.nodesBegin(), mesh.nodesEnd() );
        levelSet_.resize( numNodes );
        for ( std::size_t n = 0; n < numNodes; n++ ) {
            typename Mesh::Node* nodePtr = mesh_.nodePtr( n );
            VecDim x;
            nodePtr -> getX( &(x[0]) );
            levelSet_[n].setDefault( x, true );
        }
        
        this -> cutCells_( levelSet_, base::cut::CREATE );
    }

    //--------------------------------------------------------------------------
    //! Compute new level set, update the cell structures and unit the level sets
    template<typename SURF>
    void intersectAnalytical( const SURF& surface, const double cutThreshold )
    {
        // compute new level set
        std::vector<LevelSet> newLevelSet;
        base::cut::analyticLevelSet( mesh_, surface, true, newLevelSet );
        // update the cut-cell containers
        this -> cutCells_( newLevelSet, base::cut::INTERSECT );

        // compute intersection of level set functions
        std::transform( levelSet_.begin(), levelSet_.end(),
                        newLevelSet.begin(),
                        levelSet_.begin(),
                        boost::bind( base::cut::setIntersection<dim>, _1, _2 ) );
        
        // compress = remove degenerate cells
        std::for_each( cells_.begin(), cells_.end(),
                       boost::bind( &Cell::compress, _1, cutThreshold ) );

    }

    //--------------------------------------------------------------------------
    //! Compute new level set, update the cell structures and unit the level sets
    template<typename SMESH>
    void intersectSurfaceMesh( const SMESH& surfMesh,
                               const double cutThreshold )
    {
        std::vector<LevelSet> newLevelSet;
        base::cut::bruteForce( mesh_, surfMesh, true, newLevelSet );
        this -> cutCells_( newLevelSet, base::cut::INTERSECT );
        
        // compute intersection of level set functions
        std::transform( levelSet_.begin(), levelSet_.end(),
                        newLevelSet.begin(),
                        levelSet_.begin(),
                        boost::bind( base::cut::setIntersection<dim>, _1, _2 ) );
        
        // compress = remove degenerate cells
        std::for_each( cells_.begin(), cells_.end(),
                       boost::bind( &Cell::compress, _1, cutThreshold ) );
    }

    //--------------------------------------------------------------------------
    // Accessors
    const std::vector<LevelSet>& getLevelSet()  const { return levelSet_; }
    const std::vector<Cell>&     getCells()     const { return cells_; }
    const std::vector<SurfCell>& getSurfCells() const { return surfCells_; }

private:
    void cutCells_( const std::vector<LevelSet>& newLevelSet, 
                    const base::cut::SetOperation& so )
    {
        base::cut::generateCutCells( mesh_,         newLevelSet, cells_,     so );
        base::cut::generateCutCells( boundaryMesh_, newLevelSet, surfCells_, so );
    }
    
private:
    const Mesh&         mesh_;
    const BoundaryMesh& boundaryMesh_;

    std::vector<LevelSet> levelSet_;
    std::vector<Cell>     cells_;
    std::vector<SurfCell> surfCells_;
};

#endif 
