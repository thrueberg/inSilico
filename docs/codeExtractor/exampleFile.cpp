//! \cond SKIPDOX

//[std]{  std includes
#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
//[std]}

//[boost] boost includes
#include <boost/lexical_cast.hpp> // lexical cast between objects


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    //[inputCheck]{
    // Check the number of input arguments
    if ( argc != 2 ) { 
        std::cout << "Usage:  " << argv[0] << " pimmel \n";
        return 0;
    }
    //[inputCheck]}

    // Some pointless remarks:
    // absulotely pointless

    //[return] 
    return 0;
}
//! \endcond
