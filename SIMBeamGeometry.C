// $Id$
//==============================================================================
//!
//! \file SIMBeamGeometry.C
//!
//! \date Mar 17 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for beams with separate geometry for visualization.
//!
//==============================================================================

#include "SIMBeamGeometry.h"
#include "readTesselation.h"
#include "Utilities.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <numeric>


/*!
  \brief Sub-class of ElementBlock with specialized constructor.
*/

class GeoBlock : public ElementBlock
{
public:
  //! \brief The constructor initializes the element block.
  //! \param vertices Array of vertex coordinates, swapped with \a coord member
  //! \param[in] elms Array of Element-node connectivities
  //! \param[in] nenod Number of nodes per element
  GeoBlock(std::vector<Vec3>& vertices,
           const std::vector<IntVec>& elms, size_t nenod) : ElementBlock(nenod)
  {
    IntVec nodes(nen*elms.size());
    for (size_t i = 0; i < elms.size(); i++)
      for (size_t j = 0; j < nen; j++)
        nodes[nen*i+j] = elms[i][j];

    coord.swap(vertices);
    MMNPC.swap(nodes);
    MINEX.resize(MMNPC.size()/nen,0);
    std::iota(MINEX.begin(),MINEX.end(),1);
  }
};


bool SIMBeamGeometry::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"postprocessing"))
  {
    const tinyxml2::XMLElement* geo = elem->FirstChildElement("geometry");
    if (geo && geo->FirstChild())
    {
      double scale = 1.0;
      utl::getAttribute(geo,"scale",scale);
      std::vector<Vec3> vertices;
      std::vector<IntVec> faces;
      int nen = utl::readTesselation(geo->FirstChild()->Value(),
                                     scale,vertices,faces);
      if (nen > 2 && !vertices.empty() && !faces.empty())
        myGeometry = new GeoBlock(vertices,faces,nen);
    }
  }

  return this->SIMElasticBar::parse(elem);
}


bool SIMBeamGeometry::writeGlvG (int& nBlock, const char* inpFile, bool doClear)
{
  if (!this->SIMElasticBar::writeGlvG(nBlock,inpFile,doClear))
    return false;
  else if (!myGeometry)
    return true;

  if (msgLevel > 1)
    IFEM::cout <<"Writing tesselated geometry to VTF ("
               << myGeometry->getNoNodes() <<")"<< std::endl;

  return this->getVTF()->writeGrid(myGeometry,"Tesselated geometry",++nBlock);
}
