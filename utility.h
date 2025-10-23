#pragma once

#include <utility>
#ifndef __cpp_lib_unreachable

namespace std
{

inline void unreachable()
{
#    ifdef __GNUC__
    __builtin_unreachable();
#    elif defined( _MSC_VER )
    __assume( false );
#    endif
}

}

#endif
