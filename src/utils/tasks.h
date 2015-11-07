#include "ale/utils/zip.h"

using ale::utils::zip;



template<class Func, class Iterator, std::size_t... I>
void simple_task_impl( Func & func, Iterator &tup, std::index_sequence<I...> )
{
  for ( auto it : tup ) {
    func( it.template get<I>()... );
  }
}

template<class Func, class... Args>
void simple_task( Func & func, Args&&... args )
{
  constexpr size_t N = sizeof...(Args);
  auto tup = zip( args... );
  simple_task_impl( func, tup, std::make_index_sequence<N>() );
}



// sample usage:
// 
//   auto internal_energy = [&](const real_t & a, const real_t & b, auto & c) 
//     { c = eos.compute_internal_energy(a, b); };
//   simple_task( internal_energy, d, p, e );
