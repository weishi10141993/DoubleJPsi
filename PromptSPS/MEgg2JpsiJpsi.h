// -*- C++ -*-
//
// MEgg2JpsiJpsi.h is a plugin to Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef PLUGIN_MEgg2JpsiJpsi_H
#define PLUGIN_MEgg2JpsiJpsi_H
//
// This is the declaration of the MEgg2JpsiJpsi class.  It is
// compatible with Herwig++-2.5.0
//

#include "Herwig++/MatrixElement/HwMEBase.h"

namespace PLUGIN {
using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEgg2JpsiJpsi class implements the matrix element squared of
 * the process g g -> Jpsi Jpsi using hard coded formulae.  As such it
 * does not include spin correlations.  The non-perturbative Jpsi
 * wavefunction (_lambda) is included, and can be adjusted using the
 * interface provided
 *
 * @see \ref MEgg2JpsiJpsiInterfaces "The interfaces"
 * defined for MEgg2JpsiJpsi.
 */
class MEgg2JpsiJpsi: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEgg2JpsiJpsi();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element squared for \f$gg\to Jpsi \to l\bar{l}\f$
   */
  double gg2JpsiJpsiME() const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEgg2JpsiJpsi> initMEgg2JpsiJpsi;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEgg2JpsiJpsi & operator=(const MEgg2JpsiJpsi &);

private:

  /**
   *  pointer to the resonant meson
   */
  tcPDPtr _meson;

  /**
   *  overall gg->JpsiJpsi amplitude normalisation
   */
  double _lambda;

  /**
   *  Processes to include
   */
  unsigned int _process;

  /**
   *  Colour flow
   */
  mutable unsigned int _flow;

  /**
   *  Diagram
   */
  mutable unsigned int _diagram;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEgg2JpsiJpsi. */
template <>
struct BaseClassTrait<PLUGIN::MEgg2JpsiJpsi,1> {
  /** Typedef of the first base class of MEgg2JpsiJpsi. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEgg2JpsiJpsi class and the shared object where it is defined. */
template <>
struct ClassTraits<PLUGIN::MEgg2JpsiJpsi>
  : public ClassTraitsBase<PLUGIN::MEgg2JpsiJpsi> {
  /** Return a platform-independent class name */
  static string className() { return "PLUGIN::MEgg2JpsiJpsi"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEgg2JpsiJpsi is implemented. It may also include several, space-separated,
   * libraries if the class MEgg2JpsiJpsi depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "MEgg2JpsiJpsi.so"; }
};

/** @endcond */

}

#endif /* PLUGIN_MEgg2JpsiJpsi_H */
