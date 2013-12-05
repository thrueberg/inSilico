//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   vtk2vtu.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
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
// Converts ascii VTK file (with unstructured grid) to an XML-binary VTU file
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
    const std::string outputFilename = basename + ".bin.vtu";

    // create reader and read file
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
        vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // create writer and write file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer -> SetFileName( outputFilename.c_str() );
    writer -> SetInputConnection( reader->GetOutputPort() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();
 
    return 0;
}

 
