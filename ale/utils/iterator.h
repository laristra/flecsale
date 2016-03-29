// http://www.sj-vs.net/c-implementing-const_iterator-and-non-const-iterator-without-code-duplication/

/**
 * Our data structure (e.g., a linked list)
 */
template <class ValueType>;
class MyDataStructure {
  /**
   * Inner class that describes a const_iterator and 'regular' iterator at the same time, depending
   * on the bool template parameter (default: true - a const_iterator)
   */
  template<bool is_const_iterator = true>
  class const_noconst_iterator : public std::iterator<std::bidirectional_iterator_tag, ValueType> {
    /**
     * For const_iterator:   define DataStructurePointerType to be a   const MyDataStructure*
     * For regular iterator: define DataStructurePointerType to be a   MyDataStructure*
     */
    typedef typename std::conditional<is_const_iterator, const MyDataStructure*, MyDataStructure*>::type DataStructurePointerType;
 
    /**
     * For const_iterator:   define ValueReferenceType to be a   const ValueType&
     * For regular iterator: define ValueReferenceType to be a   ValueType&
     */
    typedef typename std::conditional<is_const_iterator, const ValueType&, ValueType&>::type ValueReferenceType;
 
    /**
     * Regular constructor: set up your iterator.
     */
    const_noconst_iterator(DataStructurePointerType pointer_to_datastructure) : _datastructure(pointer_to_datastructure) {
      // You might want to do something here, but not necessarily
    }
 
    /**
     * Copy constructor. Allows for implicit conversion from a regular iterator to a const_iterator
     */
    const_noconst_iterator(const const_noconst_iterator<false>& other) : _datastructure(other._datastructure) {
      // Copy constructor. Depending on your iterator, you might want to add something here.
    }
 
    /**
     * Equals comparison operator
     */
    bool operator== (const const_noconst_iterator& other) const {
      // Up to you to define
    }
 
    /**
     * Not-equals comparison operator
     * @see operator==(const const_noconst_iterator&) const
     */
    bool operator!= (const const_noconst_iterator& other) const {
      return !(*this == other);
    }
 
    /**
     * Dereference operator
     * @return the value of the element this iterator is currently pointing at
     */
    ValueReferenceType operator*() {
      // Up to you to define: get a reference to an element in your data structure
    }
 
    /**
     * Prefix decrement operator (e.g., --it)
     */
    const_noconst_iterator &operator--(){
      // Up to you to define: move iterator backwards
    }
 
    /**
     * Postfix decrement operator (e.g., it--)
     */
    const_noconst_iterator operator--(int){
      // Use operator--()
      const const_noconst_iterator old(*this);
      --(*this);
      return old;
    }
 
    /**
     * Prefix increment operator (e.g., ++it)
     */
    const_noconst_iterator &operator++(){
      // Up to you to define: move iterator forwards
    }
 
    /**
     * Postfix increment operator (e.g., it++)
     */
    const_noconst_iterator operator++(int){
      // Use operator++()
      const const_noconst_iterator old(*this);
      ++(*this);
      return old;
    }
 
    /**
     * Make const_noconst_iterator<true> a friend class of const_noconst_iterator<false> so
     * the copy constructor can access the private member variables.
     */
    friend class const_noconst_iterator<true>;
 
  private:
    DataStructurePointerType _list; // store a reference to MyDataStructure
 
  } // end of nested class const_noconst_iterator
 
  /**
   * Shorthand for a regular iterator (non-const) type for MyDataStructure.
   */
    typedef const_noconst_iterator<false> iterator;
 
  /**
   * Shorthand for a constant iterator (const_iterator) type for a MyDataStructure.
   */
  typedef const_noconst_iterator<true> const_iterator;
 
  // (...)
} // end of class MyDataStructure
