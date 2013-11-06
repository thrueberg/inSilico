//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   vtk2vts.cpp
//! @author Thomas Rueberg
//! @date   2013

// System includes
#include <string>
#include <iostream>
// Boost includes
#include <boost/lexical_cast.hpp>
// VTK includes
#include <vtkSmartPointer.h>
#include <vtkStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
// Corlib includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
// Converts ascii VTK file (with structured grid) to an XML-binary VTS file
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    // print usage message
    if( argc != 2 ) {
        std::cout << "Usage: " << argv[0] << " filename.vtk" << std::endl
                  << std::endl
                  << "Output: filename.bin.vts" << std::endl;
            
        return 0;
    }

    // get input file name
    const std::string filename = boost::lexical_cast<std::string>( argv[1] );
 
    // get basename
    const std::string basename = filename.substr( 0, filename.find( ".vtk" ) );
    const std::string outputFilename = basename + ".bin.vts";

    // create reader and read file
    vtkSmartPointer<vtkStructuredGridReader> reader =
        vtkSmartPointer<vtkStructuredGridReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // create writer and write file
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer -> SetFileName( outputFilename.c_str() );
    writer -> SetInputConnection( reader->GetOutputPort() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();
 
    return 0;
}

 
