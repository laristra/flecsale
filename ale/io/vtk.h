////////////////////////////////////////////////////////////////////////////////
///
/// Functions to write binary files in vtk format
///
/// \date Friday, May 20 2011
////////////////////////////////////////////////////////////////////////////////
#pragma once

// system includes
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

// user includes
#include "write_binary.h"

namespace ale {
namespace io {


////////////////////////////////////////////////////////////////////////////////
//! \brief a vtk writer class for legacy files
////////////////////////////////////////////////////////////////////////////////
class vtk_writer {

public :

  /*! *************************************************************************
   * type map
   ****************************************************************************/
  using type_map_t = std::unordered_map<std::type_index, std::string>;
  static const type_map_t type_map;

  /*! *************************************************************************
   * element map
   ****************************************************************************/
  enum class cell_type_t
  {
    triangle = 2,
    polygon = 7,
    quad = 9,
    tetra = 10,
    hexahedron = 12,
    wedge = 13,
    pyramid = 14
  };

  /*! *************************************************************************
   * Open a tecplot file for writing
   ****************************************************************************/
  auto open( const char* filename, bool binary = true ) 
  {

    // open file
    if ( binary ) file_.open(filename, std::ofstream::binary);
    else          file_.open(filename);

    binary_ = binary;
    
    // check for errors
    return !file_.good();
  }


  /*! *************************************************************************
   * close the tecplot file once comleted
   ****************************************************************************/
  auto close( void ) 
  {   

    file_.close();
    return !file_.good();
  }

  /*! *************************************************************************
   * write the header
   ****************************************************************************/
  auto init( const char* title ) 
  {
    
    // write out version number in ascii
    file_ << "# vtk DataFile Version 3.0" << std::endl;
    file_ << title << std::endl;
    if ( binary_ ) file_ << "BINARY" << std::endl;
    else           file_ << "ASCII" << std::endl;
    file_ << "DATASET UNSTRUCTURED_GRID" << std::endl;
    
    // check for write errors
    auto ierr = !file_.good();
    
    return ierr;

  }


  /*! *************************************************************************
   * write nodes
   ****************************************************************************/
  template< 
    template<typename,typename...> typename C, 
    typename T, typename... Args 
  >
  auto write_points( const C<T,Args...> & data, std::size_t npoints, std::size_t ndims )
  {

    assert( data.size() == npoints*ndims && "dimension mismatch" );

    // points header
    file_ << "POINTS " << npoints;
    file_ << " " << type_map.at( typeid(T) ) << std::endl;


    //--------------------------------------------------------------------------
    // write the data
    if ( binary_ ) {

      if ( isBigEndian() )
        for (auto val : data ) 
          WriteBinary<T>( file_, val );
      else
        for (auto val : data ) 
          WriteBinarySwap<T>( file_, val );
      
    }
    //--------------------------------------------------------------------------
    // ascii
    else {

      std::size_t cnt = 0;
      for (std::size_t p=0; p<npoints; p++ ) {
        for (std::size_t d=0; d<ndims; d++ ) 
          file_ << data[cnt++] << " ";
        file_ << std::endl;
      }

    }
    //--------------------------------------------------------------------------


    // check for write errors
    return !file_.good();
               
  }  



  /*! *****************************************************************
   * write connectivity
   ********************************************************************/
  template< 
    template<typename...> typename C, 
    typename T,
    typename... Args 
  >
  auto write_elements( const C<Args...> & data, const T * cell_type ) 
  {

    using size_t = std::size_t;
    using value_type = std::decay_t< decltype( data[0][0] )>;

    // figure out the size
    // per element: points plus # of points
    size_t size = 0;
    for ( const auto & elem : data )
      size += elem.size() + 1;

    auto nelem = data.size();

    // check endienness, vtk needs big endian
    auto swap = !isBigEndian();

    // write the number of cells ( per element: 4 points plus # of points )
    file_ << "CELLS " << nelem << " " << size << std::endl;

    //--------------------------------------------------------------------------
    // write the data
    if ( binary_ ) {

      if (swap)
        for ( const auto & elem : data ) {
          WriteBinarySwap<value_type>(file_, static_cast<value_type>(elem.size()) );
          for ( auto val : elem )
            WriteBinarySwap<value_type>(file_, val);                
        }
      else
        for ( const auto & elem : data ) {
          WriteBinary<value_type>(file_, elem.size());
          for ( auto val : elem )
            WriteBinary<value_type>(file_, val);                
        }

    }
    //--------------------------------------------------------------------------
    // ascii
    else {

      for ( const auto & elem : data ) {
        file_ << elem.size() << " ";
        for ( auto val : elem )
          file_ << val << " ";
        file_ << std::endl;
      }

    } // binary
    //--------------------------------------------------------------------------
    

    // write the cell types
    file_ << "CELL_TYPES " << nelem << std::endl;
    std::size_t cell = 0;

    //--------------------------------------------------------------------------
    // write the data
    if ( binary_ ) {

      if (swap) 
        for ( const auto & elem : data )
          WriteBinarySwap<value_type>(
            file_, static_cast<value_type>(cell_type[cell++]) 
          );
      else
        for ( const auto & elem : data )
          WriteBinary<value_type>(
            file_, static_cast<value_type>(cell_type[cell++])
          );    

    } 
    //--------------------------------------------------------------------------
    // ascii
    else {

      for ( const auto & elem : data )
        file_ << static_cast<value_type>(cell_type[cell++]) << " ";
      file_ << std::endl;

    } // binary 
    //--------------------------------------------------------------------------
    
    // check for write errors
    return !file_.good();
  
  }  


  /*! *****************************************************************
   * mark the start of cell data
   ********************************************************************/
  auto start_cell_data( std::size_t ncells ) 
  {      
    file_ << "CELL_DATA " << ncells << std::endl;    
    return !file_.good();
  }
  
  /*! *****************************************************************
   * mark the start of cell data
   ********************************************************************/
  auto start_point_data( std::size_t npoints ) 
  {   
    file_ << "POINT_DATA " << npoints << std::endl;
    return !file_.good();
  }


  /*! *****************************************************************
   * write nodes
   ********************************************************************/
  template< 
    template<typename,typename...> typename C, 
    typename T, typename... Args 
  >
  auto write_field( const char* name, const C<T,Args...> & data, std::size_t ndims = 1 )
  {

    // header
    file_ << "SCALARS " << name;
    file_ << " " << type_map.at( typeid(T) );
    file_ << " " << ndims << std::endl;
    file_ << "LOOKUP_TABLE default" << std::endl;
  

    //--------------------------------------------------------------------------
    // write the data
    if ( binary_ ) {

      if ( isBigEndian() )
        for (auto val : data ) 
          WriteBinary<T>( file_, val );
      else
        for (auto val : data ) 
          WriteBinarySwap<T>( file_, val );    

    }
    //--------------------------------------------------------------------------
    // ascii
    else {

      auto n = data.size() / ndims;
      std::size_t cnt = 0;
      
      for (std::size_t p=0; p<n; p++ ) {
        for (std::size_t d=0; d<ndims; d++ ) 
          file_ << data[cnt++] << " ";
        file_ << std::endl;
      }
      
    } // binary   
    //--------------------------------------------------------------------------
  
    
    // check for write errors
    return !file_.good();
               
  }  


private :

  //! \brief file pointer
  std::ofstream file_;

  //! \brief a boolean defining whether this is binary
  bool binary_;

};


} // namespace
} // namespace
