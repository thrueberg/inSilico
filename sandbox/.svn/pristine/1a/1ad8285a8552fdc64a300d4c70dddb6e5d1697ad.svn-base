#include <iostream>
#include <fstream>
#include <string>

#include <base/io/Format.hpp>
#include <base/post/ErrorNorm.hpp>


#include "Helper.hpp"
#include "Solid.hpp"
#include "ReferenceSolution.hpp"

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //typedef mat::hypel::StVenant Material;
    typedef mat::hypel::NeoHookeanCompressible Material;

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile, baseName;
    double E, nu, pull, tolerance;
    unsigned maxIter, loadSteps;
    bool updated;

    if ( not userInput( argc, argv,
                        meshFile, baseName,
                        E, nu, pull, tolerance,
                        maxIter, loadSteps, updated ) )
        return 0;


    
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const unsigned    dim = 2;

    typedef ::Solid<Material,dim,fieldDeg,geomDeg> Solid;
    Solid solid( meshFile, E, nu );

    const double firstPull = pull / static_cast<double>( loadSteps );

    typedef ReferenceSolution<dim> RefSol;
    typedef BoundaryValueProblem<RefSol> BVP;
    
    RefSol refSol( firstPull );
    BVP     bvp( mat::Lame::lambda( E, nu ), mat::Lame::mu(     E, nu ),
                 refSol, true ); //(not updated) );
    
    solid.constrainBoundary( boost::bind( &BVP::dirichletBC<Solid::DoF>,
                                          &bvp, _1, _2 ) );
    
    solid.numberDoFs();

    // create table for writing the convergence behaviour of the nonlinear solves
    base::io::Table<4>::WidthArray widths = {{ 2, 10, 10, 10 }};
    base::io::Table<4> table( widths );
    //table % "Load step" % "iteration" % "|F|"  % "|x|";
    //std::cout << "#" << table;

    // write a vtk file
    solid.writeVTKFile( baseName, 0 );


    
    
    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < loadSteps; step++ ) {

        refSol.setFactor( (step+1.) * firstPull );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        //----------------------------------------------------------------------
        
        unsigned iter = 0;
        while ( iter < maxIter ) {

            std::cout << step << " " << iter << "\n"
                      << "-------------------------\n";

            //table % step % iter;

            const std::pair<double,double> conv =
                solid.iterate( iter, updated, boost::bind( &BVP::bodyForce, &bvp, _1 ) );

            
            //table % conv.first % conv.second;
            //std::cout << table;

            //if ( updated ) solid.updateGeometry();

            if ( (conv.first < tolerance) or (conv.second < tolerance ) )
                break;

            iter++;
            
        }

        //if ( updated ) solid.updateGeometry();
        
        // Finished non-linear iterations
        //----------------------------------------------------------------------

        // warning
        if ( iter == maxIter ) {
            std::cout << "# (WW) Step " << step << " has not converged within "
                      << maxIter << " iterations \n";
        }

        // write a vtk file
        solid.writeVTKFile( baseName, step+1 );

        //----------------------------------------------------------------------
        {
            // compute L2-error
            std::cout //<< "L2-error = "
                      << base::post::errorComputation<0>(
                          solid.accessQuadrature(),
                          solid.accessMesh(),
                          solid.accessDisplacement(),
                          boost::bind( &BVP::solution, &bvp, _1 ) ) << "  ";
            //<< '\n';
            
            // compute H1-error
            std::cout //<< "H1-error = "
                      << base::post::errorComputation<1>(
                          solid.accessQuadrature(),
                          solid.accessMesh(),
                          solid.accessDisplacement(),
                          boost::bind( &BVP::solutionGradient, &bvp, _1 ) )
                << '\n';
        }



    }
    // Finished load steps
    //--------------------------------------------------------------------------
    
    return 0;
}
