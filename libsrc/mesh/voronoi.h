/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some functionality for generating vornoi meshes.
/// \see Kenamond, Shashkov, Mikhail, Berndt, "ShaPo Remeshing Requirements",
///      LANL Technical Report LA-UR-13-27946.
////////////////////////////////////////////////////////////////////////////////

#pragma once


// user includes
#include "detail/voronoi_impl.h"

// uncomment this to dump shapo data
#define DUMP_DIAGNOSTICS


namespace ale {
namespace mesh {

//! a unique pointer type for shapo
//! \tparam T the shapo type
template<class T>
using shapo_unique_ptr = 
  std::unique_ptr<T, std::function<void(T*)> >;

//! \brief Constructs a unique pointer for shapo of non-array type T
//! \tparam T the shapo type
//! \param [in]  args  the arguments to the constructor
//! \return a unique_ptr
template<class T, class... Args>
shapo_unique_ptr<T>
shapo_make_unique(Args&&... args) {
  return shapo_unique_ptr<T>( new T(std::forward<Args>(args)...), 
                              [](T* t){ t->Delete(); } );
}



////////////////////////////////////////////////////////////////////////////////
//! \brief Create a voronoi mesh
//!
//! \param [in] src_mesh  the source mesh to use for generators
//! \param [in] merge_tol the merge tolerance
//! \param [in] num_steps the number of relaxation steps
//! \param [in] use_clipping  
//! \return a new mesh object
////////////////////////////////////////////////////////////////////////////////
template< typename T >
T voronoi( 
  const T & src_mesh,
  typename T::real_t merge_tol,
  typename T::size_t num_steps,
  bool use_clipping ) 
{

#ifdef HAVE_SHAPO


  // the tesselator
  using Tessellator = shapo::Tessellator2D;

  //----------------------------------------------------------------------------
  // shapo setup
  //----------------------------------------------------------------------------

  // Initialize the tessellator
  decltype( shapo_make_unique<Tessellator>( 0, SHAPO_TESSELLATOR_CLIPPED ) ) tess;
  if ( use_clipping ) {
    cout << "Using clipping" << endl;
    tess = shapo_make_unique<Tessellator>( 0, SHAPO_TESSELLATOR_CLIPPED );
  }
  else {
    cout << "Using constrained triangulation" << endl;
    tess = shapo_make_unique<Tessellator>( 0, SHAPO_TESSELLATOR_CONSTRAINED );
  }

  // Ask for all connectivity information to be computed
  tess->SetConnectivityToAdvanced(); 
  
  // Ask for mergiang closed points
  tess->SetMerging(true);
  tess->SetMergingTolerance( merge_tol );
  cout << "Merge tolerance = " << merge_tol << std::endl;

  // Create generators
  detail::create_generators( src_mesh, *tess );

  // Create boundary points
  if ( use_clipping ) 
    detail::create_boundaries_clipped( src_mesh, *tess );
  else
    detail::create_boundaries_constrained( src_mesh, *tess );
 

  //----------------------------------------------------------------------------
  // Solve
  //----------------------------------------------------------------------------
  
  // Perform a check of the provided input
  auto ret = tess->PerformInputDiagnostic();
  if (ret != SHAPO_NO_ERROR) {
    cerr << endl << tess->GetInputDiagnosticString() << endl;
    raise_runtime_error("ShaPo Input Error");
  }

  // Run parameters: merging, nb_of_Lloyd_steps 
  cout << "Relaxing " << num_steps << " steps..." << std::flush;
  tess->RunCentroidal( num_steps );
  cout << "done." << endl;

  // Let's check that no error happen during the tessellation process
  ret = tess->PerformOutputDiagnostic();
  if (ret != SHAPO_NO_ERROR) {
    cerr << endl << tess->GetOutputDiagnosticString() << endl;
    //raise_runtime_error("ShaPo Error Tessellating");
  }
 

  // get the mesh
  auto shapo_mesh = tess->GetOutputMesh();
                 
  // Perform an output check of the obtained tessellation
  shapo_mesh->PrintLastError();
  
#ifdef DUMP_DIAGNOSTICS

  // Save tessellator output to file
  shapo_mesh->Save( "voromesh.vtk", true );

#endif

  //----------------------------------------------------------------------------
  // Transfer to FLAG
  //----------------------------------------------------------------------------
  
  // extract the tesselation into a new delfi mesh object
  auto new_mesh = detail::transfer_mesh<T>( *tess );

  return new_mesh;

#endif // HAVE_SHAPO  

} // voronoi

} // namespace mesh
} // namespace ale
