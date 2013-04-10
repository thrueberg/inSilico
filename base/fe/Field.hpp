//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   fe/Field.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_fe_field_hpp
#define base_fe_field_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// boost includes
#include <boost/utility.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace fe{
        
        template<typename ELEMENT, typename COEFF,
                 bool EXTERNALCOEFF = false>
        class Field;
        
    }
}

//------------------------------------------------------------------------------
/** Container for representing a field with finite elements.
 *  Independent of the type of field (geometry, solution, auxiliary variable,
 *  post-processing), in a classical FE approach it will be represented by
 *  a piecewise approximation based on the decomposition of the domain (or its
 *  approximation) into elements. On every element, a local interpolation is
 *  used for the field representation. This class provides a container for such
 *  a field. The interface is to add elements and coefficients (e.g., nodes or
 *  degrees of freedom), give acceess via iterators and random access, and
 *  deallocate the container. In some circumstances, the field is subset of an
 *  already existing field but with a possibly different element topology.
 *  This is for instance the case for a surface mesh derived from a volume mesh.
 *  In such case, the coefficients of the field representation are already
 *  stored in a different object and this class will hold only copies thereof.
 *  In that case the externalCoefficients flag is set to true and duplicate
 *  allocation and deallocation is prevented.
 *  \tparam ELEMENT    Type of element for domain decomposition
 *  \tparam COEFF      Type of field coefficient
 *  \tparam EXTERNALCOEFF  Flag for sharing the coefficients 
 */
template<typename ELEMENT, typename COEFF, bool EXTERNALCOEFF>
class base::fe::Field
    : public boost::noncopyable
{
public:
    //! @name Template parameter: type of element and node, node-sharing
    //@{
    typedef ELEMENT Element;
    typedef COEFF   Coefficient;
    static const bool externalCoefficients = EXTERNALCOEFF;
    //@}

    //! @name Container typedefs
    //@{
    typedef std::vector<Coefficient*>    CoeffPtrContainer;
    typedef std::vector<Element*>        ElementPtrContainer;
    //@}

    //! Destructor destroys the containers
    virtual ~Field()
    {
        if ( not externalCoefficients ) 
            this -> deallocateCoefficients_();
        this -> deallocateElements_();
    }
    

    //! @name Iterator access typedefs
    //@{
    typedef typename CoeffPtrContainer::iterator         CoeffPtrIter;
    typedef typename CoeffPtrContainer::const_iterator   CoeffPtrConstIter;
    typedef typename ElementPtrContainer::iterator       ElementPtrIter;
    typedef typename ElementPtrContainer::const_iterator ElementPtrConstIter;
    //@}

    //! @name Iterator access functions
    //@{
    CoeffPtrIter      coefficientsBegin()        { return coefficientPtrs_.begin(); }
    CoeffPtrIter      coefficientsEnd()          { return coefficientPtrs_.end();   }
    CoeffPtrConstIter coefficientsBegin()  const { return coefficientPtrs_.begin(); }
    CoeffPtrConstIter coefficientsEnd()    const { return coefficientPtrs_.end();   }

    ElementPtrIter      elementsBegin()       { return elementPtrs_.begin(); }
    ElementPtrIter      elementsEnd()         { return elementPtrs_.end();   }
    ElementPtrConstIter elementsBegin() const { return elementPtrs_.begin(); }
    ElementPtrConstIter elementsEnd()   const { return elementPtrs_.end();   }
    //@}

    //! @name Random access
    //@{
    Coefficient* coefficientPtr( const std::size_t c) const { return coefficientPtrs_[c]; }
    Element*     elementPtr(     const std::size_t e) const { return elementPtrs_[e]; }
    //@}

protected:
    
    //! @name Allocate functions for coefficients and elements
    //@{
    void addCoefficients_( const std::size_t toBeAdded )
    {
        // Current size of node container
        const std::size_t currentSize = coefficientPtrs_.size();

        // Increase size to new size
        coefficientPtrs_.reserve( currentSize + toBeAdded );

        // Allocate new nodes or just space
        for ( unsigned n = 0; n < toBeAdded; n++ ) {
            if ( externalCoefficients )
                coefficientPtrs_.push_back( static_cast<Coefficient*>( NULL ) );
            else
                coefficientPtrs_.push_back( new Coefficient );
        }
    }
    
    void addElements_( const std::size_t toBeAdded )
    {
        // Current size of element container
        const std::size_t currentSize = elementPtrs_.size();

        // Increase size to new size
        elementPtrs_.reserve( currentSize + toBeAdded );

        // Allocate new elements
        for ( unsigned e = 0; e < toBeAdded; e++ )
            elementPtrs_.push_back( new Element );
    }
    //@}

private:
    //! @name Deallocate functions for nodes and elements
    //@{
    void deallocateCoefficients_( )
    {
        for ( unsigned n = 0; n < coefficientPtrs_.size(); n++ )
            delete coefficientPtrs_[n];
        coefficientPtrs_.resize( 0 );
    }
    
    void deallocateElements_( )
    {
        for ( unsigned e = 0; e < elementPtrs_.size(); e++ )
            delete elementPtrs_[e];
        elementPtrs_.resize( 0 );
    }
    //@}

private:
    //! @name Data which defines the field representation
    //@{
    CoeffPtrContainer   coefficientPtrs_; //!< Pointers to coefficients
    ElementPtrContainer elementPtrs_;     //!< Pointers to elements
    //@}
};

#endif
