//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Memory.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_auxi_memory_hpp
#define base_auxi_memory_hpp

//------------------------------------------------------------------------------
// std   includes
#include <sys/sysinfo.h>

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        //----------------------------------------------------------------------
        /** Return the number of bytes currently in use
         *  Note that this function only returns the difference between the
         *  total RAM available and the free RAM. So, it will return also the
         *  memory occupied by other processes at the time of call. It can only
         *  be used in order to find out the memory difference before and after
         *  some program parts (e.g., before and after assembly). Even then the
         *  results will be imprecise and can serve only as an indication.
         *  \return Number of bytes as difference betweeen total and free RAM
         */
        unsigned long memoryUsageInBytes()
        {
            struct sysinfo memInfo;
            sysinfo (&memInfo);
            
            unsigned long physMemUsed = memInfo.totalram - memInfo.freeram;
            //convert to bytes
            physMemUsed *= memInfo.mem_unit;
            
            return physMemUsed;
        }

        //! Return the number of MegaBytes currently in use
        double memoryUsageInMegaBytes( )
        {
            // get number of bytes
            const unsigned long bytes = memoryUsageInBytes();

            // convert to mega-bytes
            const double MB =
                ( static_cast<double>( bytes ) / 1000000. );
            
            return MB;
        }

    }
}

#endif
