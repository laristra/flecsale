/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some utilities for using python.
////////////////////////////////////////////////////////////////////////////////

// use python
#include <Python.h>

namespace ale {
namespace utils {

//! \brief The python char type.
#if PY_MAJOR_VERSION < 3
using python_char_t = char;
#else
using python_char_t = wchar_t;
#endif

///////////////////////////////////////////////////////////////////////////////
//! \brief Given a character array, return a python string.
///////////////////////////////////////////////////////////////////////////////
auto make_python_string( const char * input ) 
{
#if PY_MAJOR_VERSION < 3
  return std::string(input);
#else
  return utils::to_wstring(input);
#endif
}

///////////////////////////////////////////////////////////////////////////////
//! \brief Given a character array, return a python string.
///////////////////////////////////////////////////////////////////////////////
auto make_python_string( std::string input ) 
{
	return make_python_string(input.c_str());
}

///////////////////////////////////////////////////////////////////////////////
//! \brief Given a python string, return a character array pointer.
///////////////////////////////////////////////////////////////////////////////
template< typename T >
auto as_python_char( T && input ) 
{
	return const_cast<python_char_t*>( std::forward<T>(input).c_str() ) ;
}



void python_set_program_name( const char * name ) 
{
	auto program = Py_DecodeLocale(name, NULL);
  if (program == nullptr) 
  	raise_runtime_error("Cannot decode program \"" << name << "\"");
  Py_SetProgramName(program);
  PyMem_RawFree(program);
}

void python_initialize( void ) 
{
	Py_Initialize();
}

void python_finalize( void ) 
{
	Py_Finalize();
}

auto python_import( const char * name ) 
{
	auto pname = PyUnicode_DecodeFSDefault(name);
  if (pname == nullptr) 
  	raise_runtime_error("Cannot decode module name \"" << name << "\"");
  auto pmodule = PyImport_Import(pname);
  if (pmodule == nullptr) 
  	raise_runtime_error("Cannot import module \"" << name << "\"");
  Py_DECREF(pname);
  return pmodule;
}


} // namespace utils
} // namespace ale