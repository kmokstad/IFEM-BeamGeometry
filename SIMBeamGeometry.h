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

class ASMs1D;


/*!
  \brief Driver class for beams with separate geometry for visualization.
*/

class SIMBeamGeometry : public SIMElasticBar
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMBeamGeometry(const char* hd, unsigned char n = 1)
    : SIMElasticBar(hd,n) {}
  //! \brief Empty destructor.
  virtual ~SIMBeamGeometry() {}

  //! \brief Writes current model geometry to the VTF-file.
  //! \details This method is overridden to also write out
  //! the tesselated geometry of the beam.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear);

  //! \brief Writes primary solution for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolCmps Optional number of primary solution components
  //!
  //! \details This method is overridden to also write out the deformation
  //! on the tesselated geometry based on interpolation of the beam solution.
  virtual int writeGlvS1(const Vector& psol, int iStep, int& nBlock, double,
                         const char* pvecName, int idBlock, int psolCmps, bool);

protected:
  using SIMElasticBar::parse;
  //! \brief Parses a data section from an XML element.
  //! \details This method is overridden to also read in
  //! the tesselated geometry used for visualizing the beam.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Evaluates the interpolated beam solution at a specified point.
  //! \param[in] patch The patch containing the point to evaluate at
  //! \param[in] X Spatial coordinates of the point to evaluate at
  //! \param[out] U Global displacements at the point
  //! \param[out] T Global rotation matrix at the point
  bool getBeamSolution(const ASMs1D& patch, const Vec3& X,
                       Vec3& U, Tensor& T) const;

private:
  ElementBlock* myGeometry = nullptr; //!< Tesselated geometry
};

#endif
