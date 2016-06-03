namespace detail {
 
template<typename T, typename... U>
using Array = std::array<T, 1+sizeof...(U)>;

template<typename T, typename... U, std::size_t... I>
inline Array<T, U...>
to_array(const std::tuple<T, U...>& t, std::index_sequence<I...>)
{
  return Array<T, U...>{ get<I>(t)... };
}

} // detail

template<typename T, typename... U>
inline detail::Array<T, U...>
to_array(const tuple<T, U...>& t)
{
  using IndexTuple = std::make_index_sequence<1+sizeof...(U)>;
  return detail::to_array(t, IndexTuple());
} 
 
