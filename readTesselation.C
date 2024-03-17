// $Id$
//==============================================================================
//!
//! \file readTesselation.C
//!
//! \date Mar 17 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief %IFEM wrapper for geometry file reader.
//!
//==============================================================================

#include "readTesselation.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#ifdef HAS_FFALIB
#include "FFaLib/FFaAlgebra/FFaBody.H"
#endif
#include <fstream>


/*!
  \brief Wrapper for geometry file reader.
  \param[in] fileName Name of geometry file
  \param[in] scale Vertex coordinate scaling factor
  \param[out] vertices List of vertices
  \param[out] faces List of face topology

  \details This function uses the FEDEM foundation class FFaBody to parse a
  tesselated geometry from specified file. Various file formats are supported.
*/

int utl::readTesselation (const std::string& fileName,
                          double scale,
                          std::vector<Vec3>& vertices,
                          std::vector<IntVec>& faces)
{
  std::ifstream is(fileName,std::ios::in);
  if (!is)
  {
    std::cerr <<" *** Failed to open file "<< fileName << std::endl;
    return 0;
  }

  size_t nen = 0;
#ifdef HAS_FFALIB
  size_t slashPos = fileName.find_last_of("/\\");
  size_t prevSPos = fileName.size();
  while (slashPos == prevSPos-1 && slashPos > 0)
  {
    prevSPos = slashPos;
    slashPos = fileName.find_last_of("/\\",slashPos-1);
  }
  FFaBody::prefix = fileName.substr(0,slashPos+1);
  FFaBody* body = FFaBody::readFromCAD(is);
  if (!body)
  {
    std::cerr <<" *** Invalid file or unsupported file format"<< std::endl;
    return 0;
  }

  size_t nVert = body->getNoVertices();
  size_t nFace = body->getNoFaces();
  IFEM::cout <<"Read tesselated geometry from "<< fileName
	     <<"\n\t# Vertices: "<< nVert
	     <<"\n\t# Faces   : "<< nFace << std::endl;

  vertices.resize(nVert);
  for (size_t v = 0; v < nVert; v++)
    vertices[v] = Vec3(body->getVertex(v).getPt())*scale;

  faces.resize(nFace);
  for (size_t f = 0; f < nFace; f++)
  {
    int inod;
    IntVec& mnpc = faces[f];
    mnpc.reserve(3);
    for (int i = 0; (inod = body->getFaceVtx(f,i)) >= 0; i++)
      mnpc.push_back(inod);
    if (nen == 0)
      nen = mnpc.size();
    else if (nen > 1 && nen != mnpc.size())
      nen = 1;
  }

  delete body;
#endif

  return nen;
}
