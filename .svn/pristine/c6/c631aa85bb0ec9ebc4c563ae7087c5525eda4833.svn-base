//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   vtkComp.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// System includes
#include <string>
#include <limits>
// Boost includes
#include <boost/lexical_cast.hpp>
// VTK includes
#include <vtkSmartPointer.h>
#include <vtkXMLDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
// base includes
#include <base/verify.hpp>

// approximate relative comparison due to Ericsson
bool relativelyEqual( const double& a, const double& b, const double& rtol )
{
    return
        std::abs( a - b ) / rtol * ( std::abs(a) + std::abs(b) + 1.0 );
}

//------------------------------------------------------------------------------
/** Tool to compare two vtu/vts files and quit upon difference
 */
int main( int argc, char *argv[] )
{
    // Check proper call
    if ( (argc != 3) and (argc != 4)) {
        std::cerr << "Usage: " << argv[0]
                  << " file1.vt{u|s}  file2.vt{u|s}  [rtol] " << std::endl;
        return 0;
    }

    // files to compare
    const std::string filename1 = boost::lexical_cast<std::string>( argv[1] );
    const std::string filename2 = boost::lexical_cast<std::string>( argv[2] );

    // tolerance
    const double rtol =
        ( argc == 4 ? boost::lexical_cast<double>( argv[3] ) :
          std::sqrt( std::numeric_limits<double>::min() ) );
    
    // check if the type of files is VTS
    bool isVTS;
    if (      filename1.find( ".vts" ) != std::string::npos ) isVTS = true;
    else if ( filename1.find( ".vtu" ) != std::string::npos ) isVTS = false;
    else 
        VERIFY_MSG( false, "Neither \'vtu\' nor \'vts\' file suffix found" );

    // Sanity check concerning other file
    if (     isVTS ) VERIFY_MSG( filename2.find(".vts") != std::string::npos,
                                 "Second file is not a \'vts\' file" );
    if ( not isVTS ) VERIFY_MSG( filename2.find(".vtu") != std::string::npos,
                                 "Second file is not a \'vtu\' file" );

    // output file name
    const std::string suffix = ( isVTS ? ".vts" : ".vtu" );
    const std::string outFilename = "diff" + suffix;

    // Set up readers and writer according to VTU/VTS choice
    vtkSmartPointer<vtkXMLDataReader> reader1, reader2;
    if ( isVTS ) {
        reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        reader2 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
    }
    else {
        reader1 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader2 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    }

    // Read
    {
        reader1 -> SetFileName( filename1.c_str() );
        reader1 -> Update();

        reader2 -> SetFileName( filename2.c_str() );
        reader2 -> Update();
    }
 
    // Get the grids from the reader
    vtkSmartPointer<vtkPointSet> grid1, grid2;
    if ( isVTS ) {
        vtkXMLStructuredGridReader * tmp1 
            = static_cast<vtkXMLStructuredGridReader*>( reader1.GetPointer() );
        grid1 = tmp1 -> GetOutput();
        vtkXMLStructuredGridReader * tmp2 
            = static_cast<vtkXMLStructuredGridReader*>( reader2.GetPointer() );
        grid2 = tmp2 -> GetOutput();
    }
    else {
        vtkXMLUnstructuredGridReader * tmp1 
            = static_cast<vtkXMLUnstructuredGridReader*>( reader1.GetPointer() );
        grid1 = tmp1 -> GetOutput();
        vtkXMLUnstructuredGridReader * tmp2 
            = static_cast<vtkXMLUnstructuredGridReader*>( reader2.GetPointer() );
        grid2 = tmp2 -> GetOutput();
    }

    //--------------------------------------------------------------------------
    // GRID COMPARISON

    // number of points
    const unsigned nPoints1 = grid1 -> GetNumberOfPoints();
    const unsigned nPoints2 = grid2 -> GetNumberOfPoints();
    VERIFY_MSG( nPoints1 == nPoints2, "Number of points is not equal" );

    // point coordinates
    for ( unsigned n = 0; n < nPoints1; n ++ ) {
        double * p1 = grid1 -> GetPoint( n );
        double * p2 = grid2 -> GetPoint( n );
        for ( unsigned d = 0; d < 3; d ++ ) {
            VERIFY_MSG( relativelyEqual( p1[d], p2[d], rtol ),
                        "Coordinates differ" );
        }
    }

    // number of cells
    const unsigned nCells1 = grid1 -> GetNumberOfCells();
    const unsigned nCells2 = grid2 -> GetNumberOfCells();
    VERIFY_MSG( nCells1 == nCells2, "Number of cells is not equal" );

    // cell structure
    for ( unsigned c = 0; c < nCells1; c ++ ) {
        // cell types
        VERIFY_MSG( (grid1 -> GetCellType( c )) == (grid2 -> GetCellType( c ) ),
                    "Cell types are not equal" ); 
        // individual cells
        vtkCell * cell1 = grid1 -> GetCell( c );
        vtkCell * cell2 = grid2 -> GetCell( c );
        // number of points
        const unsigned nPointsPerCell1 = cell1 -> GetNumberOfPoints();
        const unsigned nPointsPerCell2 = cell2 -> GetNumberOfPoints();
        VERIFY_MSG( nPointsPerCell1 == nPointsPerCell2,
                "Number of points per cell are not equal" );
        // individual points
        for ( unsigned n = 0; n < nPointsPerCell1; n ++ ) {
            VERIFY_MSG( (cell1 -> GetPointId(n)) == (cell2 -> GetPointId(n) ),
                    "Connectivity is not equal" );
        }
    }
        
    //--------------------------------------------------------------------------
    // DIFF POINTDATA
    {
        // get point data
        vtkSmartPointer<vtkPointData> pointData1 = grid1 -> GetPointData();
        vtkSmartPointer<vtkPointData> pointData2 = grid2 -> GetPointData();

        // number of arrays
        const unsigned nArrays1 = pointData1 -> GetNumberOfArrays();
        const unsigned nArrays2 = pointData2 -> GetNumberOfArrays();
        VERIFY_MSG( nArrays1 == nArrays2,
                    "Number of fields per point data is not equal");

        std::cout << "Comparing " << nArrays1 << " arrays of point data: "
                  << std::endl;
    
        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = pointData1 -> GetArray( n );
            vtkDataArray * array2 = pointData2 -> GetArray( n );
            // array names
            const std::string name1( pointData1 -> GetArrayName( n ) );
            const std::string name2( pointData2 -> GetArrayName( n ) );
            std::cout << n << ": " << name1 << "  vs. " << name2 << std::endl;
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            VERIFY_MSG(  nComponents1 == nComponents2,
                         "Number of components are not equal" );
            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            VERIFY_MSG( nTuples1 == nTuples2,
                        "Number of tuples are not equal" );
            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // diff each component
                for ( unsigned c = 0; c < nComponents1; c ++ ) {
                    if ( not relativelyEqual( tuple1[c], tuple2[c], rtol ) ) {
                        std::cerr << "Point datum does not match for "
                                  << name1 << "/" << name2 << "\n";
                        return 1;
                    }
                    
                }
            }
        }
    } 
    //--------------------------------------------------------------------------
    // CHECK CELLDATA
    {
        vtkSmartPointer<vtkCellData> cellData1 = grid1 -> GetCellData();
        vtkSmartPointer<vtkCellData> cellData2 = grid2 -> GetCellData();

        // number of arrays
        const unsigned nArrays1 = cellData1 -> GetNumberOfArrays();
        const unsigned nArrays2 = cellData2 -> GetNumberOfArrays();
        VERIFY_MSG( nArrays1 == nArrays2,
                    "Number of fields per cell data is not equal");
    
        std::cout << "Comparing " << nArrays1 << " arrays of cell data: "
                  << std::endl;

        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = cellData1 -> GetArray( n );
            vtkDataArray * array2 = cellData2 -> GetArray( n );
            // array names
            const std::string name1( cellData1 -> GetArrayName( n ) );
            const std::string name2( cellData2 -> GetArrayName( n ) );
            std::cout << n << ": " << name1 << "  vs. " << name2 << std::endl;
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            VERIFY_MSG(  nComponents1 == nComponents2,
                         "Number of components are not equal" );

            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            VERIFY_MSG( nTuples1 == nTuples2,
                        "Number of tuples are not equal" );

            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // check each component
                if ( not relativelyEqual( tuple1[t], tuple2[t], rtol ) ) {
                    std::cerr << "Cell datum does not match for "
                              << name1 << "/" << name2 << "\n";
                    return 1;
                }
            }
        }
    }    

    
    return 0;
}
