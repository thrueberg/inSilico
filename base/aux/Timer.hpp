//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Timer.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_aux_timer_hpp
#define base_aux_timer_hpp

//------------------------------------------------------------------------------
// std   includes
#include <sys/time.h>
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace aux{
        typedef unsigned long long TimeUnit;
        class Timer;

        namespace detail_{

            //! Helper function to get time from the system (in micro-s)
            TimeUnit getTimeInMicroSeconds()
            {
                timeval tv;
                gettimeofday(&tv, NULL);
                TimeUnit ret = tv.tv_usec + (tv.tv_sec * 1000000);
                return ret;
            }
        }
        
    }
}

//------------------------------------------------------------------------------
/** Simple object to store the current time and return the elapsed time.
 *  Upon creation or call of the \p reset() function, the current system time
 *  is stored in micro-seconds. Interfaces give back the elapsed time in
 *  micro-seconds, milli-seconds, seconds or integral minutes. A \p print()
 *  function generates a string of type XhYYmZZ.ZZZs for time display.
 */
class base::aux::Timer
{
public:
    //! Constructor sets seconds field with and stores the start value
    Timer() : width_( 5 )
    {
        reset();
    }

    //! Reset the timer
    void reset()
    {
        start_ = detail_::getTimeInMicroSeconds();
    }

    //! Return number of elapsed micro-seconds
    TimeUnit microSeconds() const
    {
        return detail_::getTimeInMicroSeconds() - start_;
    }

    //! Return elapsed milli-seconds
    double milliSeconds() const
    {
        return static_cast<double>( microSeconds() ) / 1000.;
    }

    //! Return elapsed seconds
    double seconds() const
    {
        return milliSeconds() / 1000.;
    }

    //! Return elapsed minutes (integral number)
    int minutes() const
    {
        return static_cast<int>( std::floor( seconds()/60. ) );
    }

    //! Generate and return a string for time display
    std::string print() const
    {
        double sec = seconds();
        const int hours   = static_cast<int>( sec ) / 3600; sec -= hours * 3600;
        const int minutes = static_cast<int>( sec ) / 60;   sec -= minutes * 60;

        const std::string h = boost::lexical_cast<std::string>( hours )   + "h";
        const std::string m = boost::lexical_cast<std::string>( minutes ) + "m";
        const std::string s =
            boost::lexical_cast<std::string>( sec ).substr(0,width_) + "s";

        return h + m + s;
    }

private:
    TimeUnit start_; //!< start time in micro-seconds
    int      width_; //!< print width
};

#endif
