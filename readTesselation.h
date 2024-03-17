// $Id$
//==============================================================================
//!
//! \file readTesselation.h
//!
//! \date Mar 17 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief %IFEM wrapper for geometry file reader.
//!
//==============================================================================

#ifndef UTL_READ_TESSELATION_H
#define UTL_READ_TESSELATION_H

#include "Vec3.h"
#include <vector>

using IntVec = std::vector<int>; //!< Convenience type


namespace utl
{
  int readTesselation(const std::string& fileName,
                      double scale,
                      std::vector<Vec3>& vertices,
                      std::vector<IntVec>& faces);
}

#endif
