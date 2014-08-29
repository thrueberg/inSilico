//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   InternalTraction.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_internaltraction_hpp
#define tfm_internaltraction_hpp

//------------------------------------------------------------------------------
#include <iterator>
#include <fstream>
#include <string>

#include <boost/random.hpp>

#include <base/shape.hpp>
#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/linearAlgebra.hpp>
#include <base/dof/generate.hpp>
#include <base/post/evaluateField.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

//------------------------------------------------------------------------------
namespace tfm{

    template<typename FIELD, typename SURFACEMESH>
    class InternalTraction;

    namespace detail_{

        // random number generator, seeded with time 
        static boost::random::mt19937 rng( static_cast<unsigned>( std::time(0)) );

        // return a random number from [-1,1]
        static double randomNum()
        {
            // use boost to get a random number
            boost::random::uniform_int_distribution<> randInt;
            const int x = randInt( rng );

            // random double from [0,1]
            const double y =
                static_cast<double>( x ) / static_cast<double>( randInt.max() );

            // map from [0,1] to [-1,1]
            return 2. * y - 1.;
        }
            
    }
}

//------------------------------------------------------------------------------
template<typename FIELD, typename SURFACEMESH>
class tfm::InternalTraction
{
public:
    typedef FIELD       Field;
    typedef SURFACEMESH SurfaceMesh;
    
    typedef typename base::Vector<Field::DegreeOfFreedom::size>::Type VecDoF;
    typedef typename base::Vector<SurfaceMesh::Node::dim>::Type       VecDim;

    static const base::Shape surfaceShape = SurfaceMesh::Element::shape;
    static const unsigned    degree       = SurfaceMesh::Element::GeomFun::degree;
    typedef base::fe::Basis<surfaceShape,degree,base::LAGRANGE,-1> FEBasis;
    typedef base::Field<FEBasis,Field::DegreeOfFreedom::size>  SurfField;
        

    InternalTraction( Field& adjoint, SurfaceMesh& surfaceMesh,
                      const double delta )
        : adjoint_(     adjoint ),
          surfaceMesh_( surfaceMesh ),
          delta_(       delta )
    {
        base::dof::generate<FEBasis>( surfaceMesh_, traction_ );
    }

    //----------------------------------------------------------------------
    // set new traction field
    void update( const double stepSize )
    {
        std::size_t numSurfElem = std::distance( surfaceMesh_.elementsBegin(),
                                                 surfaceMesh_.elementsEnd() );
            
        for ( std::size_t s = 0; s < numSurfElem; s++ ) {

            typename SurfaceMesh::Element* surfEp = surfaceMesh_.elementPtr( s ); 

            const typename SurfaceMesh::Element::DomainElement* domainEp =
                surfEp ->  getDomainElementPointer();
            const std::size_t elemID = domainEp -> getID();
            typename Field::Element* adjPtr = adjoint_.elementPtr( elemID );

            typename SurfField::Element* trElemPtr = traction_.elementPtr( s );

            typename SurfaceMesh::Element::ParamIter pBegin = surfEp -> parametricBegin();
            typename SurfaceMesh::Element::ParamIter pEnd   = surfEp -> parametricEnd();
            typename SurfField::Element::DoFPtrIter  dIter  = trElemPtr -> doFsBegin();
            for ( ; pBegin != pEnd; ++pBegin, ++dIter ) {

                const VecDoF adjointValue =
                    base::post::evaluateField( domainEp, adjPtr, *pBegin );
                    
                for ( unsigned d = 0; d < Field::DegreeOfFreedom::size; d++ ) {
                    const double oldT = (*dIter) -> getValue( d );

                    const double newT =
                        oldT - stepSize * ( adjointValue[d] + delta_ * oldT);

                    (*dIter) -> setValue( d, newT );

                }

            }


        }
    }
        
    //----------------------------------------------------------------------
    // set to constant
    void setToConstant( const VecDoF t ) 
    {
        std::size_t numSurfElem = std::distance( surfaceMesh_.elementsBegin(),
                                                 surfaceMesh_.elementsEnd() );
            
        for ( std::size_t s = 0; s < numSurfElem; s++ ) {

            typename SurfField::Element* trElemPtr = traction_.elementPtr( s );

            typename SurfField::Element::DoFPtrIter  dIter  = trElemPtr -> doFsBegin();
            typename SurfField::Element::DoFPtrIter  dEnd   = trElemPtr -> doFsEnd();  
            for ( ; dIter!= dEnd; ++dIter ) {

                for ( unsigned d = 0; d < Field::DegreeOfFreedom::size; d++ ) {

                    (*dIter) -> setValue( d, t[d] );
                        
                }
            }
        }

        return;
    }

    //----------------------------------------------------------------------
    // set to bla
    void setToBla( ) 
    {
        std::size_t numSurfElem = std::distance( surfaceMesh_.elementsBegin(),
                                                 surfaceMesh_.elementsEnd() );
            
        for ( std::size_t s = 0; s < numSurfElem; s++ ) {

            typename SurfField::Element* trElemPtr = traction_.elementPtr( s );

            typename SurfField::Element::DoFPtrIter  dIter  = trElemPtr -> doFsBegin();
            typename SurfField::Element::DoFPtrIter  dEnd   = trElemPtr -> doFsEnd();  
            for ( ; dIter!= dEnd; ++dIter ) {

                for ( unsigned d = 0; d < Field::DegreeOfFreedom::size; d++ ) {
                    const double value =
                        ( d== 0?
                          std::cos( 0.01 * s ) :
                          std::sin( 0.01 * s ) );
                    (*dIter) -> setValue( d, value );
                        
                }
            }
        }

        return;
    }

    //----------------------------------------------------------------------
    void randomize()
    {
        std::size_t numSurfElem = std::distance( surfaceMesh_.elementsBegin(),
                                                 surfaceMesh_.elementsEnd() );

        // 2D only!!
        VERIFY_MSG( SurfaceMesh::Node::dim == 2,
                    "This function only works in 2D ");
        VERIFY_MSG( degree == 1, "Only for linear elements" );
        
        std::vector<VecDoF> tractions;
            
        for ( std::size_t s = 0; s < numSurfElem; s++ ) {
            VecDoF t;
            for ( unsigned d = 0; d < 2; d++ ) t[d] = detail_::randomNum();
            tractions.push_back( t );
        }
        tractions.push_back( tractions[0] ); // make ring


        for ( std::size_t s = 0; s < numSurfElem; s++ ) {
            typename SurfField::Element* trElemPtr = traction_.elementPtr( s );

            typename SurfField::Element::DoFPtrIter  dIter  = trElemPtr -> doFsBegin();
            typename SurfField::Element::DoFPtrIter  dEnd   = trElemPtr -> doFsEnd();  
            for ( unsigned i=0; dIter!= dEnd; ++dIter, i++ ) {

                VecDoF t = tractions[s+i];
                                         

                for ( unsigned d = 0; d < Field::DegreeOfFreedom::size; d++ ) {
                    (*dIter) -> setValue( d, t[d] );
                }
            }
        }

        return;

    }

    //----------------------------------------------------------------------
    // return traction field
    VecDoF apply( const typename SurfaceMesh::Element* ep,
                  const typename SurfaceMesh::Element::GeomFun::VecDim& eta )
    {
        const std::size_t surfID = ep -> getID();
        typename SurfField::Element* tracPtr = traction_.elementPtr( surfID );

        return base::post::evaluateField( ep, tracPtr, eta );
    }

    //----------------------------------------------------------------------
    void writeVTKFile( const std::string& baseName,
                       const unsigned stepNum )
    {
        const std::string name = baseName + "." + 
            base::io::leadingZeros( stepNum ) + ".vtk";

        //
        std::ofstream vtk( name.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( surfaceMesh_ );
        base::io::vtk::writePointData( vtkWriter, surfaceMesh_, traction_, "traction" );
    }

    //--------------------------------------------------------------------------
    VecDoF monitor()
    {
        const std::size_t numElements = std::distance( surfaceMesh_.elementsBegin(),
                                                       surfaceMesh_.elementsEnd() );
        return
            this -> apply( surfaceMesh_.elementPtr( numElements/2 ),
                           base::ShapeCentroid<SurfaceMesh::Element::shape>::apply() );
    }

    //--------------------------------------------------------------------------
    SurfaceMesh& getMesh() { return surfaceMesh_; }

private:
    Field&       adjoint_;
    SurfaceMesh& surfaceMesh_;

    SurfField    traction_;
        
    const double delta_;
};

#endif
