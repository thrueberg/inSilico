// System includes
#include <string>
#include <iostream>
// Boost includes
#include <boost/lexical_cast.hpp>
// VTK includes
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
// Corlib includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
// Converts ASCII-XML vtu/s-files to binary 
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // print usage message
    if( argc != 2 ) {
        std::cout << "Usage: " << argv[0] << " filename.vtk" << std::endl
                  << std::endl
                  << "Output: filename.bin.vtu" << std::endl;
            
        return 0;
    }

    // get input file name
    const std::string filename = boost::lexical_cast<std::string>( argv[1] );
 
    // get basename
    const std::string basename = filename.substr( 0, filename.find( ".vtk" ) );
    const std::string outputFilename = basename + ".bin.vtk";

    // create reader and writer
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
        vtkSmartPointer<vtkUnstructuredGridReader>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    
    // read file
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // write XML binary file
    writer -> SetFileName( outputFilename.c_str() );
    writer -> SetInputConnection( reader->GetOutputPort() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();
 
    return 0;
}

 
