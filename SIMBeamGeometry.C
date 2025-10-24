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
#include "ASMs1D.h"
#include "IntegrandBase.h"
#include "LocalIntegral.h"
#include "Utilities.h"
#include "Vec3Oper.h"
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
      double scale = 1.0, tol = -1.0;
      utl::getAttribute(geo,"scale",scale);
      utl::getAttribute(geo,"tolerance",tol);
      std::vector<Vec3> vertices;
      std::vector<IntVec> faces;
      int nen = utl::readTesselation(geo->FirstChild()->Value(),
                                     scale,tol,vertices,faces);
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


bool SIMBeamGeometry::getBeamSolution (const ASMs1D& patch, const Vec3& X,
                                       Vec3& U, Tensor& T) const
{
  class ElementVecs : public LocalIntegral
  {
  public:
    size_t size() const { return vec.front().size(); }
    double operator[](int i) { return vec.front()[i]; }
  };

  // Find element containing this point
  std::pair<int,double> elm = patch.findElement(X);
  if (elm.first < 1) return false; // Point is not in this patch

  // Extract the element displacement vector
  ElementVecs eV;
  if (!myProblem->initElement(*(patch.begin_elm()+elm.first-1),eV))
    return false;
  else if (eV.size() != 12)
  {
    // Consider support for higher-order elements later...
    std::cerr <<" *** SIMBeamGeometry::getBeamSolution: For linear element only"
              << std::endl;
    return false;
  }

  // Interpolate using first-order Lagrange polynomials.
  // TODO: Use Hermitian interpolation.
  auto&& N1 = [](double xi) { return 0.5*(1.0 - xi); };
  auto&& N2 = [](double xi) { return 0.5*(1.0 + xi); };

  Vec3 R;
  for (int i = 0; i < 3; i++)
  {
    U[i] = N1(elm.second)*eV[i]   + N2(elm.second)*eV[6+i];
    R[i] = N1(elm.second)*eV[3+i] + N2(elm.second)*eV[9+i];
  }
  T = Tensor(R.x,R.y,R.z);

  // TODO: Transform from local to global (assumed identity transform for now)
  return true;
}


int SIMBeamGeometry::writeGlvS1 (const Vector& psol, int iStep, int& nBlock,
                                 double, const char* pvecName, int idBlock,
                                 int psolCmps, bool)
{
  VTF* vtf = this->getVTF();
  if (psol.empty() || !vtf)
    return idBlock; // no primary solution

  IntVec vID;
  std::vector<IntVec> sID;
  sID.reserve(this->getNoFields());

  Matrix field;
  Vector lovec;

  int geomID = this->getStartGeo();
  for (const ASMbase* pch : myModel)
    if (!pch->empty())
    {
      if (msgLevel > 1)
        IFEM::cout <<"Writing primary solution for patch "
                   << pch->idx+1 << std::endl;

      // Evaluate primary solution variables
      pch->extractNodeVec(psol,lovec);
      if (!pch->evalSolution(field,lovec,opt.nViz))
        return -1;

      // Output as vector field
      if (!vtf->writeVres(field,++nBlock,++geomID,nsd))
        return -2;
      else
        vID.push_back(nBlock);

      // Output as scalar fields
      for (size_t j = 1; j <= field.rows(); j++)
        if (!vtf->writeNres(field.getRow(j),++nBlock,geomID))
          return -3;
        else if (j > sID.size())
          sID.push_back({nBlock});
        else
          sID[j-1].push_back(nBlock);
    }

  if (myGeometry && myProblem)
  {
    if (msgLevel > 1)
      IFEM::cout <<"Writing deformation for tesselated geometry"<< std::endl;

    Vec3 Xb, Ub, Upt;
    Tensor Tb(nsd);
    field.resize(nsd,myGeometry->getNoNodes());
    const ASMs1D* myPatch = nullptr;
    for (const ASMbase* pch : myModel)
      if (!pch->empty() && (myPatch = dynamic_cast<const ASMs1D*>(pch)))
      {
        myPatch->extractNodeVec(psol,myProblem->getSolution());
        for (size_t i = 0; i < field.cols(); i++)
        {
          const Vec3& X0 = myGeometry->getCoord(i);

          // Assume for now the beam is aligned with the global X-axis.
          // TODO: Project X0 onto the beam model, account for multi-patch, etc.
          Xb.x = X0.x;

          // Calculate deformation at this point from the beam solution
          if (this->getBeamSolution(*myPatch,Xb,Ub,Tb))
          {
            Upt = (Xb+Ub) + Tb*(X0-Xb) - X0;
            field.fillColumn(1+i,Upt.ptr());
          }
        }
      }

    // Output as vector field
    if (!vtf->writeVres(field,++nBlock,++geomID,nsd))
      return -2;
    else
      vID.push_back(nBlock);

    // Output as scalar fields
    for (size_t j = 1; j <= field.rows(); j++)
      if (!vtf->writeNres(field.getRow(j),++nBlock,geomID))
        return -3;
      else
        sID[j-1].push_back(nBlock);
  }

  // Write result block identifications

  bool ok = true;
  if (!vID.empty())
  {
    if (pvecName)
      ok = vtf->writeVblk(vID,pvecName,idBlock,iStep);
    else
      ok = vtf->writeDblk(vID,"Solution",idBlock,iStep);
  }

  if (idBlock <= static_cast<int>(nsd))
    idBlock = nsd+1; // since we might have written BCs

  std::string pname(pvecName ? pvecName : "d"); pname += "_w";
  for (size_t i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
  {
    if (myProblem && (!pvecName || psolCmps <= 0))
      pname = myProblem->getField1Name(i);
    else if (i > 0 && i%nsd == 0)
    {
      pname.back() = 'r';
      pname += 'x';
    }
    else
      ++pname.back();
    ok = vtf->writeSblk(sID[i],pname.c_str(),idBlock++,iStep);
  }

  return ok ? idBlock : -4;
}
