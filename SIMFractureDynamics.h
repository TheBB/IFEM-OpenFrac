// $Id$
//==============================================================================
//!
//! \file SIMFactureDynamics.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_DYNAMICS_H_
#define _SIM_FRACTURE_DYNAMICS_H_

#include "SIMCoupled.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSplineSurface.h"
#endif
#include <fstream>


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling=SIMCoupled>
class SIMFracture : public Coupling<SolidSolver,PhaseSolver>
{
public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile)
    : Coupling<SolidSolver,PhaseSolver>(s1,s2), infile(inputfile), aMin(0.0) {}
  //! \brief Empty destructor.
  virtual ~SIMFracture() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S1.registerDependency(&this->S2,"phasefield",1);
    // The tensile energy is defined on integration points and not nodal points.
    // It is a global buffer array across all patches in the model.
    // Use an explicit call instead of normal couplings for this.
    this->S2.setTensileEnergy(this->S1.getTensileEnergy());
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (!energFile.empty() && this->S1.getProcessAdm().getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      if (tp.step == 1)
        os <<"#t eps_e external_energy eps+ eps- eps_b |c|"
           <<" eps_d-eps_d(0) eps_d"<< std::endl;

      const Vector& n1 = this->S1.getGlobalNorms();
      const Vector& n2 = this->S2.getGlobalNorms();

      size_t i;
      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (i = 0; i < n1.size(); i++) os <<" "<< n1[i];
      for (i = 0; i < n2.size(); i++) os <<" "<< n2[i];
      os << std::endl;
    }

    return this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock);
  }

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Stores current solution state in an internal buffer.
  void saveState()
  {
    dsol = this->S1.getSolutions();
    psol = this->S2.getSolution();
    hsol = this->S2.getHistoryField();
  }

  //! \brief Refines the mesh on the initial configuration.
  bool initialRefine(double beta, double min_frac, size_t nrefinements)
  {
    TimeStep step0;
    int newElements = 1;
    while (newElements > 0)
      if (!this->S2.solveStep(step0))
        return false;
      else
        newElements = this->adaptMesh(1.0,min_frac,nrefinements);

    return newElements == 0;
  }

  //! \brief Refines the mesh with transfer of solution onto the new mesh.
  int adaptMesh(double beta, double min_frac, size_t nrefinements)
  {
#ifdef HAS_LRSPLINE
    ASMu2D* pch = dynamic_cast<ASMu2D*>(this->S1.getPatch(1));
    if (!pch)
      return -1;

    if (aMin <= 0.0) // maximum refinements per element
      aMin = pch->getBasis()->getElement(0)->area()/(nrefinements*nrefinements);

    // Fetch element norms to use as refinement criteria
    Vector eNorm;
    double gNorm = this->S2.getNorm(eNorm);

    // Sort element indices based on comparing values in eNorm
    IntVec idx(eNorm.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),
              [&eNorm](size_t i1, size_t i2) { return eNorm[i1] < eNorm[i2]; });

    double eMin = min_frac*gNorm/sqrt(idx.size()); // lower cap on norm value
    size_t eMax = beta < 0.0 ? idx.size() : idx.size()*beta/100.0;
    IFEM::cout <<"\n  Lowest element: "<< std::setw(8) << idx.front()
               <<"    |c| = "<< eNorm[idx.front()]
               <<"\n  Highest element:"<< std::setw(8) << idx.back()
               <<"    |c| = "<< eNorm[idx.back()]
               <<"\n  Minimum |c|-value for refinement: "<< eMin
               <<"\n  Minimum element area: "<< aMin << std::endl;

    IntVec elements; // Find the elements to refine
    for (size_t i = 0; i < idx.size() && elements.size() < eMax; i++)
      if (eNorm[idx[i]] > eMin)
        break;
      else if (pch->getBasis()->getElement(idx[i])->area() > aMin)
        elements.push_back(idx[i]);

    if (elements.empty())
      return 0;

    IFEM::cout <<"  Elements to refine: "<< elements.size()
               <<" (|c| = ["<< eNorm[elements.front()]
               <<","<< eNorm[elements.back()] <<"])\n"<< std::endl;

    std::shared_ptr<LR::LRSplineSurface> oldBasis;
    if (!hsol.empty())
      oldBasis.reset(pch->getBasis()->copy());

    // Do the mesh refinement
    LR::RefineData prm;
    prm.options = { 10, 1, 2 };
    prm.elements = pch->getFunctionsForElements(elements);
    if (!this->S1.refine(prm,dsol) || !this->S2.refine(prm,psol))
      return -2;

    // Re-initialize the simulators for the new mesh
    this->S1.clearProperties();
    this->S2.clearProperties();
    if (!this->S1.read(infile.c_str()) || !this->S2.read(infile.c_str()))
      return -3;

    if (!this->preprocess())
      return -4;

    if (!this->init(TimeStep()))
      return -5;

    if (!this->S1.initSystem(this->S1.opt.solver) ||
        !this->S2.initSystem(this->S2.opt.solver,1,1,false))
      return -6;

    // Transfer solution variables onto the new mesh
    if (!dsol.empty())
    {
      IFEM::cout <<"\nTransferring "<< dsol.size() <<"x"<< dsol.front().size()
                 <<" solution variables to new mesh for "<< this->S1.getName();
      this->S1.setSolutions(dsol);
    }
    if (!psol.empty())
    {
      IFEM::cout <<"\nTransferring "<< psol.size()
                 <<" solution variables to new mesh for "<< this->S2.getName();
      this->S2.setSolution(psol);
    }
    if (!hsol.empty())
    {
      IFEM::cout <<"\nTransferring "<< hsol.size()
                 <<" history variables to new mesh for "<< this->S2.getName()
                 << std::endl;
      this->S2.transferHistory2D(hsol,oldBasis.get());
    }

    return elements.size();
#endif

    return -1;
  }

private:
  std::string energFile; //!< File name for global energy output
  std::string infile;    //!< Input file parsed

  double    aMin; //!< Minimum element area
  Vectors   dsol; //!< Displacement state to transfer onto refined mesh
  Vector    psol; //!< Phase-field state to transfer onto refined mesh
  RealArray hsol; //!< History field to transfer onto refined mesh
};

#endif
