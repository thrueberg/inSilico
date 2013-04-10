#include <iostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/mesh/ExtractMeshFaces.hpp>
#include <base/mesh/GenerateMeshFromFaces.hpp>
#include <base/mesh/FaceIterator.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>
#include <base/LagrangeShapeFun.hpp>

#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/evaluateAtNodes.hpp>
#include <base/dof/scaleConstraints.hpp>
#include <base/dof/generate.hpp>

#include <base/StiffnessMatrix.hpp>
#include <base/ForceIntegrator.hpp>
#include <base/NeumannForce.hpp>


#include <base/solver/Eigen3.hpp>

#include <base/fe/Basis.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/CompNeoHookean.hpp>

#include <mat/Lame.hpp>

#include <solid/Elasticity.hpp>
#include <solid/Stress.hpp>

//------------------------------------------------------------------------------
template<unsigned DIM>
class PulledSheetProblem
{
public:
    typedef typename base::VectorType<DIM>::Type VecDim;

    template<typename DOF>
    static void dirichletBC( const VecDim& x, DOF* doFPtr, const double value ) 
    {
        const double tol = 1.e-5;
    
        const bool onLeftBdr =
            ( std::abs( x[0] -  0. ) < tol );

        const bool onRightBdr =
            ( std::abs( x[0] -  1. ) < tol );

        
        if ( onLeftBdr ) {
            for ( unsigned d = 0; d < DOF::size; d++ ) {
                if ( doFPtr -> isActive(d) )
                    doFPtr -> constrainValue( d, 0.0 );
            }
        }

        if (  onRightBdr ) {
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
            
        }
        return;
    }

    static VecDim neumannBC( const VecDim& x,
                             const VecDim& normal )
    {
        VecDim result = VecDim::Constant( 0. );

        const double tol = 1.e-5;
        
        const bool onTractionBdr = false
            ( std::abs( x[0] -  1. ) < tol );

        if ( onTractionBdr ) result[1] = - 10.0;
        
        return result;
    }
    
};


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );
    
    std::string meshFile;
    double E, nu, pull;
    unsigned maxIter, loadSteps;
    {    
        //! Feed properties parser with the variables to read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",  meshFile );
        prop.registerPropertiesVar( "E",         E );
        prop.registerPropertiesVar( "nu",        nu );
        prop.registerPropertiesVar( "pull",      pull );
        prop.registerPropertiesVar( "maxIter",   maxIter );
        prop.registerPropertiesVar( "loadSteps", loadSteps );

        //! Read variables from the input.dat file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY( inp.is_open() );
        prop.readValues( inp );
        inp.close( );
    }
    
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const base::Shape shape    = base::QUAD;

    //--------------------------------------------------------------------------
    const unsigned    dim     = base::ShapeDim<shape>::value;
    typedef base::mesh::Node<dim>                 Node;
    typedef base::LagrangeShapeFun<geomDeg,shape> SFun;
    typedef base::mesh::Element<Node,SFun>        Element;
    typedef base::mesh::Unstructured<Element>     Mesh;

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, smf ); 
        smf.close();
    }

    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // DOF handling
    const unsigned    doFSize = dim;
    typedef base::dof::DegreeOfFreedom<doFSize>    DoF;
    typedef base::LagrangeShapeFun<fieldDeg,shape> FieldFun;
    typedef base::dof::Element<DoF,FieldFun>       FieldElement;
    typedef base::dof::Field<FieldElement>         Field;
    Field field;

    // generate DoFs from mesh
    typedef base::fe::Basis<shape,fieldDeg> FEBasis;
    base::dof::generate<FEBasis>( mesh, field );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary
    const double firstPull = pull / static_cast<double>( loadSteps );
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, field, 
                                           boost::bind( &PulledSheetProblem<dim>::dirichletBC<DoF>,
                                                        _1, _2, firstPull ) );

    typedef mat::hypel::CompNeoHookean Material;
    Material material( mat::lambda( E, nu), mat::mu( E, nu ) );

    typedef solid::Elasticity<Material,Element,FieldElement> Elasticity;
    Elasticity elasticity( material );
            
    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << " Number of dofs " << numDofs << std::endl;

    for ( unsigned step = 0; step < loadSteps; step++ ) {
        std::cout << "* Load step: " << step << std::endl;
        std::cout << "   Non-linear iterations ";
        unsigned iter = 0;
        
        while ( iter < maxIter ) {

            std::cout <<  iter << "  " << std::flush;
    
            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( numDofs );


            //typedef mat::hypel::StVenant Material;
            //--------------------------------------------------------------------------
            typedef base::InternalForceIntegrator<Quadrature,Solver,Element,
                                                  FieldElement,FieldElement>
                InternalForceIntegrator;
            
            InternalForceIntegrator::ForceKernel internalForce =
                boost::bind( &Elasticity::internalForce, &elasticity, _1, _2, _3, _4, _5, _6 );
        
            InternalForceIntegrator forceInt( internalForce, quadrature, solver );
            base::aux::forEach3( mesh.elementsBegin(), mesh.elementsEnd(),
                                 field.elementsBegin(), field.elementsBegin(),
                                 forceInt );
        
            // Compute element stiffness matrices and assemble them
            base::stiffnessMatrixComputation( quadrature, solver,
                                              mesh, field, field,
                                              boost::bind( &Elasticity::tangentStiffness,
                                                           &elasticity,
                                                           _1, _2, _3, _4, _5, _6),
                                              iter > 0 );

            // Finalise assembly
            solver.finishAssembly();

            std::cout << "|F| = " << solver.norm() << " " << std::flush;

            // Solve
            solver.choleskySolve();

            
            // distribute results back to dofs
            base::dof::Distribute<DoF,Solver,base::dof::ADD>
                distributeDoF( solver, iter > 0 );
            std::for_each( field.doFsBegin(), field.doFsEnd(), distributeDoF );
    
            std::cout << "|x| = " << solver.norm() << " " << std::flush;
    
            iter++;
        }

        // OUTPUT
        {
            // Evaluate the solution field at every geometry node
            std::vector<base::VectorType<doFSize>::Type> nodalValues;
            base::dof::evaluateAtNodes( mesh, field, nodalValues );
            const std::string vtkFile =
                "test." + base::io::leadingZeros( step ) + ".vtk";
            std::ofstream vtk( vtkFile.c_str() );
            base::io::vtk::LegacyWriter vtkWriter( vtk );
            vtkWriter.writeUnstructuredGrid( mesh );
            vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "disp" );

            // compute Cauch stress at element mid points
            {
                std::vector<mat::Tensor> cauchyStress;
                std::transform( mesh.elementsBegin(), mesh.elementsEnd(),
                                field.elementsBegin(), std::back_inserter( cauchyStress ),
                                boost::bind( solid::cauchy<Element,FieldElement,Material>,
                                             _1, _2, material ) );
                vtkWriter.writeCellData( cauchyStress.begin(), cauchyStress.end(), "sigma" );

                
            }

            
            vtk.close();

        }


        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    return 0;
}
