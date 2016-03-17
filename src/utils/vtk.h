////////////////////////////////////////////////////////////////////////////////
///
/// Functions to write binary files in vtk format
///
/// \date Friday, May 20 2011
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "write_binary.h"

namespace ale {
namespace utils {


////////////////////////////////////////////////////////////////////////////////
//! \brief a vtk writer class
////////////////////////////////////////////////////////////////////////////////
class vtk_writer {

public :

  /*! *************************************************************************
   * element map
   ****************************************************************************/
  enum class cell_type_t
  {
    vtk_triangle = 2,
    vtk_polygon = 7,
    vtk_quad = 9,
    vtk_tetra = 10,
    vtk_hexahedron = 12,
    vtk_wedge = 13,
    vtk_pyramid = 14
  };

  /*! *************************************************************************
   * Open a tecplot file for writing
   ****************************************************************************/
  auto open( const char* filename ) 
  {
    
    // open file
    file_.open(filename, std::ofstream::binary);
    auto ierr = !file_.good();
    
    // check for errors
    if (ierr) return ierr;
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
    file_ << "BINARY" << std::endl;
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
    typename T = float, typename... Args 
  >
  auto write_points( const C<T,Args...> & data, std::size_t npoints, std::size_t ndims )
  {

    assert( data.size() == npoints*ndims && "dimension mismatch" );

    // check endienness, vtk needs big endian
    bool swap = !isBigEndian();

    // points header
    file_ << "POINTS " << npoints;
    file_ << " float" << std::endl;

    // write the data
    if (swap)
      for (auto val : data ) 
        WriteBinaryFloatSwap( file_, val );
    else
      for (auto val : data ) 
        WriteBinaryFloat( file_, val );

    // check for write errors
    return !file_.good();
               
  }  

  /*! *************************************************************************
   * write nodes
   ****************************************************************************/
  template< 
    template<typename,typename...> typename C, 
    typename... Args 
  >
  auto write_points( const C<double,Args...> & data, std::size_t npoints, std::size_t ndims )
  {

    assert( data.size() == npoints*ndims && "dimension mismatch" );

    // check endienness, vtk needs big endian
    bool swap = !isBigEndian();

    // points header
    file_ << "POINTS " << npoints;
    file_ << " double" << std::endl;

    // write the data
    if (swap)
      for (auto val : data ) 
        WriteBinaryDoubleSwap( file_, val );
    else
      for (auto val : data ) 
        WriteBinaryDouble( file_, val );

    // check for write errors
    return !file_.good();
               
  }  



  /*! *****************************************************************
   * write connectivity
   ********************************************************************/
  template< 
    template<typename...> typename C, 
    typename... Args 
    >
  auto write_elements( const C<Args...> & data, cell_type_t cell_type ) 
  {

    using size_t = std::size_t;

    // figure out the size
    // per element: points plus # of points
    size_t size = 0;
    for ( const auto & elem : data )
      size += elem.size() + 1;

    auto nelem = data.size();

    // check endienness, vtk needs big endian
    bool swap = !isBigEndian();

    // write the number of cells ( per element: 4 points plus # of points )
    file_ << "CELLS " << nelem << " " << size << std::endl;

    // write the data
    if (swap) {

      for ( const auto & elem : data ) {
        WriteBinaryIntSwap(file_, elem.size());
        for ( auto val : elem )
          WriteBinaryIntSwap(file_, val);                
      }

    } 
    else {

      for ( const auto & elem : data ) {
        WriteBinaryInt(file_, elem.size());
        for ( auto val : elem )
          WriteBinaryInt(file_, val);                
      }
      
    }

    // write the cell types
    file_ << "CELL_TYPES " << nelem << std::endl;

    auto type_id = static_cast<int>(cell_type);        

    // now write it
    if (swap) 
      for ( const auto & elem : data )
        WriteBinaryIntSwap(file_, type_id);

    else

      for ( const auto & elem : data )
        WriteBinaryInt(file_, type_id);    


    // check for write errors
    return !file_.good();
  
  }  

private :

  //! \brief file pointer
  std::ofstream file_;


};


#if 0


/*! *****************************************************************
 * mark the start of cell data
 ********************************************************************/
void vtk_start_cell_data( int n, int ier ) 
{   

  file_ << "CELL_DATA " << n << std::endl;

  ier = !file_.good();
}

/*! *****************************************************************
 * mark the start of cell data
 ********************************************************************/
void vtk_start_point_data( int n, int ier ) 
{   

  file_ << "POINT_DATA " << n << std::endl;

  ier = !file_.good();
}

/*! *****************************************************************
 * write nodes
 ********************************************************************/
void vtk_write_scalar( char* name, int n, void *data, int isDouble, int &ier ) 
{

  int i;
  

  // check endienness, vtk needs big endian
  bool swap = !isBigEndian();




  // header
  file_ << "SCALARS " << name;

  if ( isDouble == 0 )
    file_ << " float" << std::endl;
  else
    file_ << " double" << std::endl;
  
  file_ << "LOOKUP_TABLE default" << std::endl;
  

  // write the data
  if (isDouble == 0) {
  
    if (swap)
      for (i=0; i<n; i++) 
        WriteBinaryFloatSwap( file_, ((float*)data)[i] );
    else
      for (i=0; i<n; i++) 
        WriteBinaryFloat( file_, ((float*)data)[i] );
  
  } else {
    exit(-1);
    if (swap)
      for (i=0; i<n; i++)
        WriteBinaryDoubleSwap( file_, (double)i );
    else
      for (i=0; i<n; i++)
        WriteBinaryDouble( file_, ((double*)data)[i] );        
  
  }


  // check for write errors
  ier = !file_.good();
               
}  

#endif

} // namespace
} // namespace
