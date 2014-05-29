#include <string>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/linearAlgebra.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/auxi/FundamentalSolution.hpp>

#include <mat/Lame.hpp>
#include <Eigen/Cholesky>

//------------------------------------------------------------------------------
// Give the minimal distance between mesh a and mesh b.
// Address check is done in order to change policy for identical meshes.
template<typename MESH>
double minimalDistance( const MESH& a, const MESH& b )
{
    const std::size_t numA = std::distance( a.nodesBegin(), a.nodesEnd() );
    const std::size_t numB = std::distance( b.nodesBegin(), b.nodesEnd() );

    const bool identical = (&a == &b);

    double result = std::numeric_limits<double>::max();
    for ( std::size_t nA = 0; nA < numA; nA++ ) {

        const typename MESH::Node* aPtr = a.nodePtr( nA );
        const typename MESH::Node::VecDim xA = aPtr -> getX( );

        const std::size_t lowerB = (identical ? nA+1 : 0);

        for ( std::size_t nB = lowerB; nB < numB; nB++ ) {

            const typename MESH::Node* bPtr = b.nodePtr( nB );
            const typename MESH::Node::VecDim xB = bPtr -> getX();

            const double dist = (xA - xB).norm();
            if ( dist < result ) result = dist;
        }

    }

    return result;
}

//------------------------------------------------------------------------------
template<typename MESH>
bool geomFilter( const MESH& cell, const typename MESH::Node* nodePtr,
                 const double thresh )
{
    const typename MESH::Node::VecDim x = nodePtr -> getX();
    
    const std::size_t numNodes =
        std::distance( cell.nodesBegin(), cell.nodesEnd() );

    double result = std::numeric_limits<double>::max();
    for ( std::size_t n = 0; n < numNodes; n++ ) {

        const typename MESH::Node* aPtr = cell.nodePtr( n );
        const typename MESH::Node::VecDim xA = aPtr -> getX( );

        const double dist = (xA - x).norm();
        if ( dist < result ) result = dist;
    }

    return ( result > thresh );
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // Read input from command line
    if ( argc != 4 ) {
        std::cerr << "Usage:  " << argv[0] << "forceLocation.smf "
                  << " dispLocation.smf  dispValues \n\n";
        return -1;
    }
    const std::string forcePointFile = boost::lexical_cast<std::string>( argv[1] );
    const std::string  dispPointFile = boost::lexical_cast<std::string>( argv[2] );
    const std::string       dispFile = boost::lexical_cast<std::string>( argv[3] );

    // material
    const double E = 1.;
    const double nu = 0.3;
    const double thresh = 10.;

    // attributes of the point clouds
    static const unsigned    dim    = 3;
    static const base::Shape shape  = base::POINT;
    static const unsigned    degree = 1;
    
    typedef base::Unstructured<shape,degree,dim> PointCloud;
    typedef PointCloud::Node::VecDim VecDim;
    typedef base::auxi::FundSolElastoStatic<dim> FundSol;

    // read points of force location
    PointCloud forcePoints;
    {
        std::ifstream smf( forcePointFile.c_str() );
        base::io::smf::readMesh( smf, forcePoints );
    }

    // read points of displacement location
    PointCloud dispPoints;
    {
        std::ifstream smf( dispPointFile.c_str() );
        base::io::smf::readMesh( smf, dispPoints );
    }

    // Check the distances of the point clouds
    {
        std::cout << "Minimal distance among force points: "
                  << minimalDistance( forcePoints, forcePoints )
                  << std::endl;
        std::cout << "Minimal distance among displacement points: "
                  << minimalDistance( dispPoints, dispPoints )
                  << std::endl;
        std::cout << "Minimal distance between point clouds: "
                  << minimalDistance( dispPoints, forcePoints )
                  << std::endl;
    }

    // numbers involved
    const std::size_t numDisp = std::distance( dispPoints.nodesBegin(),
                                               dispPoints.nodesEnd() );
    const std::size_t numForc = std::distance( forcePoints.nodesBegin(),
                                               forcePoints.nodesEnd() );

    VERIFY_MSG( numDisp >= numForc, "System is not solvable (not enough data)" );

    // Read displacements from file
    std::vector<VecDim> displacements, activeDisplacements;
    {
        // read double values first
        std::ifstream dis( dispFile.c_str() );
        std::istream_iterator<double> eos;       
        std::istream_iterator<double> iit (dis);
        std::vector<double> dummy;
        while ( iit != eos ) {
            dummy.push_back( *iit );
            ++iit;
        }

        VERIFY_MSG( dummy.size() == numDisp * dim,
                    "Not enough displacements provided" );

        // convert to VecDims
        for ( std::size_t n = 0; n < numDisp; n++ ) {
            VecDim disp = base::constantVector<dim>( 0. );
            //for ( unsigned d = 0; d < dim; d++ ) {
            for ( unsigned d = 0; d < 2; d++ ) {
                disp[d] = dummy[ n * dim + d];
            }
            displacements.push_back( disp );
        }
    }

    // Write displacements
    {
        const std::string vtkFile = "displacement.vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( dispPoints );
        vtkWriter.writePointData( displacements.begin(), displacements.end(), "u" );
    }

    // filter geometry
    std::size_t numActive = 0;
    std::vector<bool> active;
    PointCloud::NodePtrConstIter nIter = dispPoints.nodesBegin();
    PointCloud::NodePtrConstIter nEnd  = dispPoints.nodesEnd();
    for ( std::size_t ctr = 0; nIter != nEnd; ++nIter, ctr++ ) {
        const bool isActive = geomFilter( forcePoints,  *nIter, thresh );

        if ( isActive ) numActive++;
        active.push_back( isActive );
        if ( isActive ) activeDisplacements.push_back( displacements[ctr] );
    }

    PointCloud dispPointsFiltered;
    dispPointsFiltered.allocate( numActive, 0 );
    nIter  = dispPoints.nodesBegin();
    nEnd   = dispPoints.nodesEnd();
    PointCloud::NodePtrIter      nIter2 = dispPointsFiltered.nodesBegin();
    for ( std::size_t ctr = 0; nIter != nEnd; ++nIter, ctr++ ) {
        const bool isActive = active[ ctr ];
        if ( isActive ) {
            (*nIter2) -> deepCopy( *nIter );
            nIter2++;
        }
    }
    
    // Write active displacements
    {
        const std::string vtkFile = "activeDisplacement.vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( dispPointsFiltered );
        vtkWriter.writePointData( activeDisplacements.begin(),
                                  activeDisplacements.end(), "u" );
    }


    VERIFY_MSG( numActive >= numForc, "System is not solvable (filter to strict)" );

    std::cout << "Using " << numActive << " out of " << numDisp << " Points \n";

#if 1
    // Fundamental solution of linear elastostatics
    FundSol fundSol( mat::Lame::lambda( E, nu ), mat::Lame::mu( E, nu ) );
        
    // System matrices
    base::MatrixD A = base::MatrixD::Zero( numActive * dim, numForc * dim );
    base::VectorD b = base::VectorD::Zero( numActive * dim );

    // Assemble system
    for ( std::size_t nD = 0; nD < numActive; nD++ ) {

        const PointCloud::Node* xPtr = dispPointsFiltered.nodePtr( nD );
        const VecDim x = xPtr -> getX( );

        // LHS
        for ( std::size_t nF = 0; nF < numForc; nF++ ) {
            const PointCloud::Node* yPtr = forcePoints.nodePtr( nF );
            VecDim y;
            yPtr -> getX( &(y[0]) );
            
            FundSol::MatDimDim U = fundSol.U( x, y );
            
            A.block( nD*dim, nF*dim, dim, dim ) = U;
        }

        // RHS
        b.segment( nD*dim, dim ) = displacements[nD];

    }
        

    // Solve least squares system
    std::vector<VecDim> forces, computedDisplacements;
    {
        base::MatrixD AtA; AtA.noalias() = A.transpose() * A;
        base::VectorD Atb; Atb.noalias() = A.transpose() * b;
        Eigen::LDLT<base::MatrixD> cholesky( AtA );
        base::VectorD x = cholesky.solve( Atb );

        // VERIFY_MSG( cholesky.isPositive(),
        //             "Matrix is bad! " );
        
        for ( std::size_t nF = 0; nF < numForc; nF++ ) {
            const VecDim F = x.segment( nF*dim, dim );
            forces.push_back( F );
        }

        // compute displacement
        const base::VectorD u = A * x;
        for ( std::size_t nD = 0; nD < numActive; nD++ ) {
            const VecDim uD = u.segment( nD*dim, dim );
            computedDisplacements.push_back( uD );
        }

        // Compute error
        {
            base::VectorD e1 = b - u;
            base::VectorD e2 = Atb - AtA * x;
            std::cout << "Errors " << e1.norm()
                      << " and " << e2.norm()
                      << std::endl;
        }
    }

    // Write forces
    {
        const std::string vtkFile = "forces.vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( forcePoints );
        vtkWriter.writePointData( forces.begin(), forces.end(), "F" );
    }

    // Write computed displacements
    {
        const std::string vtkFile = "computedDisplacement.vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( dispPointsFiltered );
        vtkWriter.writePointData( computedDisplacements.begin(),
                                  computedDisplacements.end(), "u" );
    }

#endif
    
    return 0;
}
