/////////////////////////////////////////////////////////////////////
/// \file write_binary.h
///
/// Utilities for writing binary files
///
/// \date Wednesday, March 14 2012
/////////////////////////////////////////////////////////////////////

#pragma once

namespace ale {
namespace utils {

//! data types
using int32 = int;
using float32 = float;
using float64 = double;

/*! *****************************************************************
 * Determine endienness
 ********************************************************************/
inline bool isBigEndian(void) {
  static unsigned long x(1);
  static bool result(reinterpret_cast<unsigned char*>(&x)[0] == 0);
  return result;
}

/*! *****************************************************************
 * The main binary writing functions
 ********************************************************************/
template <class T> inline
void WriteBinary(std::ostream &file, T buffer) {
  file.write(reinterpret_cast<char*>(&buffer), sizeof(T));
}

inline void WriteBinaryInt(std::ofstream &file, int32 buffer) 
{ WriteBinary(file, buffer); }

inline void WriteBinaryIntSwap(std::ofstream &file, int32 buffer) 
{ 
  union temp {
    int32  value;
    char   c[4];
  } in, out;

  in.value = buffer;

  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  
  file.write(out.c, sizeof(int32));
}

inline void WriteBinaryFloat(std::ofstream &file, float32 buffer) 
{ WriteBinary(file, buffer); }

inline void WriteBinaryFloatSwap(std::ofstream &file, float32 buffer) 
{ 
  union temp {
    float32  value;
    char     c[4];
  } in, out;

  in.value = buffer;

  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  
  file.write(out.c, sizeof(float32));
}


inline void WriteBinaryDouble(std::ofstream &file, float64 buffer) 
{ WriteBinary(file, buffer); }

inline void WriteBinaryDoubleSwap(std::ofstream &file, float64 buffer) 
{ 
  union temp {
    float64  value;
    char     c[8];
  } in, out;

  in.value = buffer;

  out.c[0] = in.c[7];
  out.c[1] = in.c[6];
  out.c[2] = in.c[5];
  out.c[3] = in.c[4];
  out.c[4] = in.c[3];
  out.c[5] = in.c[2];
  out.c[6] = in.c[1];
  out.c[7] = in.c[0];
  
  file.write(out.c, sizeof(float64));
}

inline void WriteString(std::ofstream &file, const char *S) {
  int L = 0;
  while (S[L] != '\0')
    WriteBinary<int32>(file,int(S[L++]));
  WriteBinary<int32>(file,0);
}

} // namspace
} // namspace
