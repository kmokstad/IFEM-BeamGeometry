// $Id$
//==============================================================================
//!
//! \file SIMBeamGeometry.h
//!
//! \date Mar 17 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for beams with separate geometry for visualization.
//!
//==============================================================================

#ifndef _SIM_BEAM_GEOMETRY_H
#define _SIM_BEAM_GEOMETRY_H

#include "SIMElasticBar.h"
#include "ElementBlock.h"


/*!
  \brief Driver class for beams with separate geometry for visualization.
*/

class SIMBeamGeometry : public SIMElasticBar
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMBeamGeometry(const char* hd) : SIMElasticBar(hd) {}
  //! \brief Empty destructor.
  virtual ~SIMBeamGeometry() {}

  //! \brief Writes current model geometry to the VTF-file.
  //! \details This method is overridden to also write out
  //! the tesselated geometry of the beam.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear);

protected:
  using SIMElasticBar::parse;
  //! \brief Parses a data section from an XML element.
  //! \details This method is overridden to also read in
  //! the tesselated geometry used for visualizing the beam.
  virtual bool parse(const tinyxml2::XMLElement* elem);

private:
  ElementBlock* myGeometry = nullptr; //!< Tesselated geometry
};

#endif
