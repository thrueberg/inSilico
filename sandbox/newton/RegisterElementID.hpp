#ifndef register_hpp
#define register_hpp

#include <set>
#include <iterator>

struct RegisterElementID
{
    static std::set<std::size_t> idRegister;

    static void registerID( const std::size_t id )
    {
        idRegister.insert( id );
    }

    static void printIDs()
    {
        if ( idRegister.size() ) {
            std::cout << "Registered: ";
            std::copy( idRegister.begin(), idRegister.end(),
                       std::ostream_iterator<std::size_t>( std::cout, "  " ) );
            std::cout << std::endl;
        }
    }
};


std::set<std::size_t> RegisterElementID::idRegister = std::set<std::size_t>();

#endif
