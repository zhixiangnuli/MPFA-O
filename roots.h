#pragma once

#include "utility.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <numeric>
#include <type_traits>

namespace resim
{
namespace roots
{

enum Status
{
    SUCCESS                       = 0x00,
    BISECT_NO_CHANGE_OF_SIGN      = 0x10,
    BISECT_CANNOT_REACH_TOLERANCE = 0x11,
    NEWTON_GET_INVALID_DERIVATIVE = 0x20,
    NEWTON_CANNOT_REACH_TOLERANCE = 0x21
};

template <typename T, typename F>
int bisect( F f, T* x, T min, T max, T* res = nullptr, T tol = std::numeric_limits<T>::epsilon() )
{
    T r, *y = res != nullptr ? res : &r;
    assert( x != nullptr && max > min );

    T ymin = std::invoke( f, min );
    if ( std::abs( ymin ) <= tol )
    {
        *x = min, *y = ymin;
        return Status::SUCCESS;
    }

    T ymax = std::invoke( f, max );
    if ( std::abs( ymax ) <= tol )
    {
        *x = max, *y = ymax;
        return Status::SUCCESS;
    }

    if ( std::signbit( ymin ) ^ !std::signbit( ymax ) )
    {
        if ( std::abs( ymin ) < std::abs( ymax ) )
        {
            *x = min, *y = ymin;
            return Status::BISECT_NO_CHANGE_OF_SIGN;
        }
        else
        {
            *x = max, *y = ymax;
            return Status::BISECT_NO_CHANGE_OF_SIGN;
        }
    }

    const int its = 1 - std::floor( std::log2( tol / ( max - min ) ) );
    assert( its > 0 );

    *y = std::numeric_limits<T>::infinity();
    for ( int i = 0; i < its; ++i )
    {
        *x = ( min + max ) / 2, *y = std::invoke( f, *x );
        if ( std::abs( *y ) <= tol )
            return Status::SUCCESS;
        else if ( std::signbit( *y ) ^ std::signbit( ymin ) )
            max = *x;
        else if ( std::signbit( *y ) ^ std::signbit( ymax ) )
            min = *x;
        else
            return Status::SUCCESS;
    }

    return Status::BISECT_CANNOT_REACH_TOLERANCE;
}

template <typename T, typename F>
int inverse_quadratic_interp( F f, T* x, T min, T max, T* res = nullptr,
                              T tol = std::numeric_limits<T>::epsilon() )
{
    T r, *y = res != nullptr ? res : &r;
    assert( x != nullptr && max > min );

    T ymin = std::invoke( f, min );
    if ( std::abs( ymin ) <= tol )
    {
        *x = min, *y = ymin;
        return Status::SUCCESS;
    }

    T ymax = std::invoke( f, max );
    if ( std::abs( ymax ) <= tol )
    {
        *x = max, *y = ymax;
        return Status::SUCCESS;
    }

    if ( std::signbit( ymin ) ^ !std::signbit( ymax ) )
    {
        if ( std::abs( ymin ) < std::abs( ymax ) )
        {
            *x = min, *y = ymin;
            return Status::BISECT_NO_CHANGE_OF_SIGN;
        }
        else
        {
            *x = max, *y = ymax;
            return Status::BISECT_NO_CHANGE_OF_SIGN;
        }
    }

    T xhst[3] = { min, max, min + ( max - min ) / ( ymin - ymax ) * ymin };
    T yhst[3] = { ymin, ymax, std::invoke( f, *std::rbegin( xhst ) ) };

    const int its = 1 - std::floor( std::log2( tol / ( max - min ) ) );
    assert( its > 0 );

    *y = std::numeric_limits<T>::infinity();
    for ( int i = 0;; ++i )
    {
        *x = ( yhst[1] * yhst[2] * xhst[0] ) / ( ( yhst[0] - yhst[1] ) * ( yhst[0] - yhst[2] ) ) +
             ( yhst[2] * yhst[0] * xhst[1] ) / ( ( yhst[1] - yhst[2] ) * ( yhst[1] - yhst[0] ) ) +
             ( yhst[0] * yhst[1] * xhst[2] ) / ( ( yhst[2] - yhst[0] ) * ( yhst[2] - yhst[1] ) );

        if ( *x < min || *x > max || std::isnan( *x ) || i >= its )
        {
            return bisect<T>( f, x, min, max, y, tol );
        }
        else if ( *y = std::invoke( f, *x ); std::abs( *y ) <= tol )
        {
            return Status::SUCCESS;
        }
        else
        {
            int q = std::max_element( std::begin( yhst ), std::end( yhst ),
                                      []( T l, T r ) { return std::abs( l ) < std::abs( r ); } ) -
                    std::begin( yhst );
            xhst[q] = *x;
            yhst[q] = *y;
        }
    }

    std::unreachable();
}

template <typename T, typename F, typename N = std::size_t>
int newton_iterate( F f, T* x, T min = std::numeric_limits<T>::lowest(),
                    T max = std::numeric_limits<T>::max(), T* res = nullptr,
                    T tol = std::numeric_limits<T>::epsilon(), N its = 1000 )
{
    T r, *y = res != nullptr ? res : &r;

    std::deque<std::pair<T, T>> itv;
    for ( int i = 0;; ++i )
    {
        T g = 0.0;
        if constexpr ( std::is_invocable_r_v<T, F, T, T*> )
        {
            *y = std::invoke( f, *x, &g );
        }
        else if constexpr ( std::is_invocable_r_v<T, F, T> )
        {
            *y = std::invoke( f, *x );

            T e = std::copysign( std::min( std::abs( *x ), 1.0 ), *x );
            T t = e *= std::sqrt( std::numeric_limits<T>::epsilon() );
            T u = std::invoke( f, *x + e );
            T v = std::invoke( f, *x - t );

            if ( std::isnan( u ) || std::isinf( u ) ) u = *y, e = 0.0;
            if ( std::isnan( v ) || std::isinf( v ) ) v = *y, t = 0.0;
            g = ( u - v ) / ( e + t );
        }
        else
        {
            static_assert( !sizeof( F* ), "INVOKABLE_FUNCTION_EXPECTED" );
        }

        if ( itv.size() > 1 )
        {
            T u = itv.front().first, v = itv.back().first;
            if ( i >= 1 - std::floor( std::log2( tol / ( v - u ) ) ) )
            {
                if constexpr ( std::is_invocable_r_v<T, F, T, T*> )
                {
                    using std::placeholders::_1;
                    return inverse_quadratic_interp( std::bind( f, _1, &g ), x, u, v, res, tol );
                }
                else
                {
                    return inverse_quadratic_interp( f, x, u, v, res, tol );
                }
            }
        }

        if ( std::abs( g ) < std::numeric_limits<T>::epsilon() )
        {
            return Status::NEWTON_CANNOT_REACH_TOLERANCE;
        }
        else if ( std::isnan( g ) || std::isinf( g ) || std::isnan( *y ) || std::isinf( *y ) )
        {
            if ( itv.size() > 1 )
            {
                T u = itv.front().first, v = itv.back().first;
                if constexpr ( std::is_invocable_r_v<T, F, T, T*> )
                {
                    using std::placeholders::_1;
                    return inverse_quadratic_interp( std::bind( f, _1, &g ), x, u, v, res, tol );
                }
                else
                {
                    return inverse_quadratic_interp( f, x, u, v, res, tol );
                }
            }
            else
            {
                return Status::NEWTON_GET_INVALID_DERIVATIVE;
            }
        }
        else
        {
            if ( std::abs( *y ) <= tol )
            {
                return Status::SUCCESS;
            }
            else if ( i > its )
            {
                return Status::NEWTON_CANNOT_REACH_TOLERANCE;
            }
            else if ( itv.empty() )
            {
                itv.emplace_back( *x, *y );
            }
            else if ( itv.size() == 1 )
            {
                if ( std::signbit( *y ) ^ std::signbit( itv.back().second ) )
                {
                    if ( *x < itv.back().first )
                    {
                        itv.emplace_front( *x, *y );
                    }
                    else
                    {
                        itv.emplace_back( *x, *y );
                    }
                }
                else
                {
                    if ( std::abs( *y ) < std::abs( itv.back().second ) )
                    {
                        itv.pop_back();
                        itv.emplace_back( *x, *y );
                    }
                }
            }
            else
            {
                if ( std::signbit( *y ) ^ std::signbit( itv.back().second ) )
                {
                    assert( std::signbit( *y ) ^ !std::signbit( itv.front().second ) );
                    if ( *x > itv.front().first )
                    {
                        itv.pop_front();
                        itv.emplace_front( *x, *y );
                    }
                }
                else
                {
                    assert( std::signbit( *y ) ^ !std::signbit( itv.back().second ) );
                    if ( *x < itv.back().first )
                    {
                        itv.pop_back();
                        itv.emplace_back( *x, *y );
                    }
                }
            }

            T t = *x - *y / g;
            if ( t < min )
            {
                *x = ( *x + min ) / 2.0;
            }
            else if ( t > max )
            {
                *x = ( *x + max ) / 2.0;
            }
            else
            {
                *x = t;
            }
        }
    }

    std::unreachable();
}

}
}
