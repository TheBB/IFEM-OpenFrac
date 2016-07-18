// $Id$
//==============================================================================
//!
//! \file PoroFracture.C
//!
//! \date Apr 15 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for elasticity problems with fracture.
//!
//==============================================================================

#include "PoroFracture.h"
#include "PoroMaterial.h"
#include "FractureElasticityVoigt.h"
#include "FiniteElement.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


class PoroFractureElasticityVoigt : public FractureElasticityVoigt
{
public:
  PoroFractureElasticityVoigt(IntegrandBase *parent, unsigned short int n)
    : FractureElasticityVoigt(parent, n)
  {}

  void setMode(SIM::SolutionMode mode)
  {
    this->FractureElasticityVoigt::setMode(mode);
    iS = PoroElasticity::Fu + 1;
  }
};


PoroFracture::PoroFracture (unsigned short int n) : PoroElasticity(n)
{
  fracEl = new PoroFractureElasticityVoigt(this,n);

  L_per = 0.01;
  d_min = 0.1;
  Kc    = 83.0;
  eps   = 50.0;
}


PoroFracture::~PoroFracture()
{
  delete fracEl;
}


bool PoroFracture::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"crack"))
    return this->PoroElasticity::parse(elem) & fracEl->parse(elem);

  IFEM::cout <<"\tCrack parameters:";
  if (utl::getAttribute(elem,"Kc",Kc))
    IFEM::cout <<" Kc = "<< Kc;
  if (utl::getAttribute(elem,"eps",eps))
    IFEM::cout <<" eps = "<< eps;
  if (utl::getAttribute(elem,"L_per",L_per))
    IFEM::cout <<" L_per = "<< L_per;
  if (utl::getAttribute(elem,"d_min",d_min))
    IFEM::cout <<" d_min = "<< d_min;
  IFEM::cout << std::endl;

  return true;
}


Material* PoroFracture::parseMatProp (const TiXmlElement* elem, bool)
{
  this->PoroElasticity::parseMatProp(elem,true);
  fracEl->setMaterial(material);
  return material;
}


void PoroFracture::setMaterial (Material* mat)
{
  this->PoroElasticity::setMaterial(mat);
  fracEl->setMaterial(mat);
}


void PoroFracture::setMode (SIM::SolutionMode mode)
{
  this->PoroElasticity::setMode(mode);
  fracEl->setMode(mode);
  primsol.resize(fracEl->getNoSolutions());
}


void PoroFracture::initIntegration (size_t nGp, size_t nBp)
{
  fracEl->initIntegration(nGp,nBp);
}


LocalIntegral* PoroFracture::getLocalIntegral (size_t nen,
                                               size_t, bool neumann) const
{
  LocalIntegral* elmInt = this->PoroElasticity::getLocalIntegral(nen,0,neumann);
  fracEl->setVar(nsd+1);
  return elmInt;
}


LocalIntegral* PoroFracture::getLocalIntegral (const std::vector<size_t>& nen,
                                               size_t, bool neumann) const
{
  LocalIntegral* elmInt = this->PoroElasticity::getLocalIntegral(nen,0,neumann);
  fracEl->setVar(nsd);
  return elmInt;
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                LocalIntegral& elmInt)
{
  // Allocating three more vectors on the element level compared to global level
  // (1 = pressure, 2 = pressure rate, 3 = phase field)
  elmInt.vec.resize(primsol.size()+3);

  // Extract element displacement, velocity, acceleration and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elmInt))
    return false;

  // Extract element phase-field vector
  return fracEl->initElement(MNPC,elmInt);
}


bool PoroFracture::initElement (const std::vector<int>& MNPC,
                                const std::vector<size_t>& elem_sizes,
                                const std::vector<size_t>& basis_sizes,
                                LocalIntegral& elmInt)
{
  // Allocating three more vectors on the element level compared to global level
  // (1 = pressure, 2 = pressure rate, 3 = phase field)
  elmInt.vec.resize(primsol.size()+3);

  // Extract element displacement, velocity, acceleration and pressure vectors
  if (!this->PoroElasticity::initElement(MNPC,elem_sizes,basis_sizes,elmInt))
    return false;

  // Extract element phase-field vector
  std::vector<int>::const_iterator uEnd = MNPC.begin() + elem_sizes.front();
  return fracEl->initElement(std::vector<int>(MNPC.begin(),uEnd),elmInt);
}


const RealArray* PoroFracture::getTensileEnergy () const
{
  return fracEl->getTensileEnergy();
}

#define EXTRA 24

size_t PoroFracture::getNoFields(int fld) const
{
  size_t super = this->PoroElasticity::getNoFields(fld);
  if (fld > 1)
    super += nsd * (nsd + 1) / 2 + EXTRA;
  return super;
}


std::string PoroFracture::getField2Name(size_t i, const char *prefix) const
{
  if (i < this->PoroElasticity::getNoFields(2))
    return this->PoroElasticity::getField2Name(i, prefix);

  i -= this->PoroElasticity::getNoFields(2);
  std::string name;

  if (i < EXTRA) {
    static const char* t[] = {".ux", ".uy", ".phase", ".dphase_x", ".dphase_y",
                              ".invphase", ".|dphase|^2",
                              ".Fxx", ".Fxy", ".Fyx", ".Fyy",
                              ".Cxx", ".Cxy", ".Cyy",
                              ".Cinvxx", ".Cinvxy", ".Cinvyy",
                              ".Cinv dphase x", ".Cinv dphase y",
                              ".bracket_yy", ".lambda", ".detF", ".w",
                              "w"};
    name = t[i];
  }
  else {
    i -= EXTRA;
    static const char* s[][6] = {{"x", "y", "xy"},
                                 {"x", "y", "z", "yz", "xz", "xy"}};
    size_t ncomps = nsd * (nsd + 1) / 2;
    name = "perm" + std::string("_") + s[nsd-2][i % ncomps];
  }

  if (!prefix)
    return name;
  return prefix + std::string(" ") + name;
}


bool PoroFracture::evalSol(Vector& s, const FiniteElement& fe,
                           const Vec3& X, const std::vector<int>& MNPC) const
{
  if (!this->PoroElasticity::evalSol(s, fe, X, MNPC))
    return false;

  Vector eV;
  utl::gather(MNPC, 1, fracEl->getCVec(), eV);

  Vector eDtemp, eD(nsd * fe.N.size());
  utl::gather(MNPC, nsd+1, primsol.front(), eDtemp);
  for (size_t i = 0; i < nsd; i++)
    for (size_t bfun = 0; bfun < fe.N.size(); bfun++)
      eD[nsd*bfun+i] = eDtemp[(nsd+1)*bfun+i];

  SymmTensor K(nsd);
  Vector extra;
  double w = this->formCrackedPermeabilityTensor(K, eV, eD, fe, X, extra);

  s.insert(s.end(), extra.begin(), extra.end());
  s.insert(s.end(), w);

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) return false;

  Vec3 permeability = pmat->getPermeability(X);
  for (size_t i = 1; i <= K.dim(); i++)
    K(i,i) += permeability(i);

  if (nsd == 2) {
    s.insert(s.end(), K(1,1));
    s.insert(s.end(), K(2,2));
    s.insert(s.end(), K(1,2));
  }
  else if (nsd == 3) {
    s.insert(s.end(), K(1,1));
    s.insert(s.end(), K(2,2));
    s.insert(s.end(), K(3,3));
    s.insert(s.end(), K(2,3));
    s.insert(s.end(), K(1,3));
    s.insert(s.end(), K(1,2));
  }

  return true;
}


bool PoroFracture::evalElasticityMatrices (ElmMats& elMat, const Matrix&,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  // Evaluate tangent stiffness matrix and internal forces
  return fracEl->evalInt(elMat,fe,X);
}


/*!
  This method calculates the anisotropic permeability for the broken solid
  based on a Poiseuille flow approximation of the fluid flow in the crack.
  See Section 5.5.2 (eqs. (104)-(109)) of Miehe and Maute (2015)
  "Phase field modeling of fracture in multi-physics problems. Part III."
*/

double PoroFracture::formCrackedPermeabilityTensor (SymmTensor& Kcrack,
                                                    const Vector& eV,
                                                    const Vector& eD,
                                                    const FiniteElement& fe,
                                                    const Vec3& X,
                                                    Vector& extra) const
{
  for (size_t i = 0; i < EXTRA - 1; i++)
    extra.push_back(0.0);

  if (!eV.empty())
  {
    double ux = 0.0, uy = 0.0, p = 0.0;
    for (size_t i = 0; i < fe.N.size(); i++) {
      ux += eD[2*i] * fe.N[i];
      uy += eD[2*i+1] * fe.N[i];
      p += eV[i] * fe.N[i];
    }
    extra[0] = ux;
    extra[1] = uy;
    extra[2] = p;
  }

  Vec3 gradD; // Evaluate the phase field value and gradient
  double d = fracEl->evalPhaseField(gradD,eV,fe.N,fe.dNdX);

  gradD[0] = 0.0;               // FIXME

  extra[3] = gradD[0];
  extra[4] = gradD[1];

  if (d < 0.0)
  {
    std::cerr <<" *** PoroFracture::formCrackedPermeabilityTensor(X = "<< X
              <<")\n     Invalid phase field value: "<< d << std::endl;
    return d;
  }
  else if (d < d_min)
  {
    // The crack does not affect the permeability tensor at this point
    Kcrack.zero();
    return 0.0;
  }

  extra[5] = d;

  double d2 = gradD.length2();
  if (d2 <= 0.0) {              // FIXME
    gradD[1] = 1.0;
    d2 = gradD.length2();
  }
  // if (d2 <= 0.0)
  // {
  //   std::cerr <<" *** PoroFracture::formCrackedPermeabilityTensor(X = "<< X
  //             <<")\n     Zero phase field gradient: "<< gradD << std::endl;
  //   return -1.0;
  // }

  extra[6] = d2;

  Tensor F(nsd); // Calculate the deformation gradient
  if (!this->formDefGradient(eD,fe.N,fe.dNdX,X.x,F))
    return -2.0;

  F(1,1) = 1.0;                 // FIXME
  F(1,2) = 0.0;                 // FIXME
  F(2,1) = 0.0;                 // FIXME
  // F(2,2) = 1.0;                 // FIXME

  extra[7] = F(1,1);
  extra[8] = F(1,2);
  extra[9] = F(2,1);
  extra[10] = F(2,2);

  // Compute the inverse right Cauchy-Green tensor (C^-1)
  // if (Kcrack.rightCauchyGreen(F).inverse() == 0.0)
  //   return -3.0;

  // SymmTensor K(nsd);
  Kcrack.rightCauchyGreen(F);

  extra[11] = Kcrack(1,1);
  extra[12] = Kcrack(1,2);
  extra[13] = Kcrack(2,2);

  Kcrack.inverse();

  extra[14] = Kcrack(1,1);
  extra[15] = Kcrack(1,2);
  extra[16] = Kcrack(2,2);

  // std::cout << "gradD = " << gradD << " F = " << F << " Cinv = " << Kcrack;

  // Compute the symmetric tensor C^-1 - (C^-1*n0)otimes(C^-1*n0)
  // (the term in the bracket [] of eq. (108) in Miehe2015pfm3)
  Vec3 CigD = Kcrack*gradD; // C^-1*gradD

  extra[17] = CigD[0] / sqrt(d2);
  extra[18] = CigD[1] / sqrt(d2);

  // std::cout << "CigD = " << CigD << std::endl << std::endl;
  for (unsigned short int j = 1; j <= nsd; j++)
    for (unsigned short int i = 1; i <= j; i++)
      Kcrack(i,j) -= sqrt(CigD(i)*CigD(j)/d2);
  // std::cout << "Bracket term = " << Kcrack << std::endl;

  extra[19] = Kcrack(2,2);

  // Compute the perpendicular crack stretch
  // lambda = gradD*gradD / gradD*C^-1*gradD (see eq. (106))
  double lambda = sqrt(d2 / (gradD*CigD));

  extra[20] = lambda;

#if INT_DEBUG > 3
  std::cout <<"PoroFracture::formCrackedPermeabilityTensor(X = "<< X
            <<") lambda = "<< lambda << std::endl;
#endif
  if (lambda <= 1.0)
  {
    Kcrack.zero();
    return 0.0;
  }

  // Compute the permeability tensor, scale by d^eps*Kc*w^2*J (see eq. (108))
  double w = lambda*L_per - L_per; // Crack opening (see eq. (107))
  extra[21] = F.det();
  extra[22] = w;
  // Kcrack *= pow(d,eps) * (w*w/(1.002e-6)*9801.0 - 2.07e-9) * F.det();
  Kcrack *= pow(d,eps)*Kc*w*w*F.det();
  // std::cout << "Final = " << Kcrack << std::endl;
  return w < 0.0 ? 0.0 : w;
}


bool PoroFracture::formPermeabilityTensor (SymmTensor& K,
                                           const Vectors& eV,
                                           const FiniteElement& fe,
                                           const Vec3& X) const
{
  unsigned short int eC = fracEl->geteC();

  if (eV.size() <= eC)
  {
    std::cerr <<" *** PoroFracture::formPermeabilityTensor:"
              <<" Missing phase field solution vector."<< std::endl;
    return false;
  }
  else if (!eV[eC].empty() && eV[eC].size() != fe.N.size())
  {
    std::cerr <<" *** PoroFracture::formPermeabilityTensor:"
              <<" Invalid phase field vector.\n     size(eC) = "
              << eV[eC].size() <<"   size(N) = "<< fe.N.size() << std::endl;
    return false;
  }

  Vector temp;
  if (this->formCrackedPermeabilityTensor(K,eV[fracEl->geteC()],eV.front(),fe,X,temp) < 0.0)
    return false;

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) return false;

  Vec3 permeability = pmat->getPermeability(X);
  for (size_t i = 1; i <= K.dim(); i++)
    K(i,i) += permeability(i);

  return true;
}
