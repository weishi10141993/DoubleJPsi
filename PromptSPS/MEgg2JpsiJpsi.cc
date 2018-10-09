// -*- C++ -*-
//
// MEgg2JpsiJpsi.cc is a plugin to Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEgg2JpsiJpsi class.  It is compatible with
// Herwig++-2.5.0
//

#include "MEgg2JpsiJpsi.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace PLUGIN;
using namespace Herwig;

MEgg2JpsiJpsi::MEgg2JpsiJpsi():_lambda(0.92),_process(0) {
  massOption(vector<unsigned int>(2,1));
}

void MEgg2JpsiJpsi::doinit() {
  // the final state meson can be changed here
  _meson=getParticleData(ParticleID::Jpsi);
}

IBPtr MEgg2JpsiJpsi::clone() const {
  return new_ptr(*this);
}

IBPtr MEgg2JpsiJpsi::fullclone() const {
  return new_ptr(*this);
}

double MEgg2JpsiJpsi::gg2JpsiJpsiME() const {
  //Energy2 u(uHat());
  Energy2 t(tHat()),s(sHat());
  Energy wave(1.*GeV);
  Energy2 mc2(sqr(0.5*_meson->mass()));
  double ta(acos( (-8.*mc2+s+2.*t)/(sqrt(s*s-16.*s*mc2)) ));
  double sb(s/mc2);
  double s1(sin(1.*ta));
  double c1(cos(1.*ta));
  double c2(cos(2.*ta));
  double c4(cos(4.*ta));
  double c6(cos(6.*ta));
  double c8(cos(8.*ta));
  double c10(cos(10.*ta));
  double c12(cos(12.*ta));

  // output is s^2*dsig/dt
  double output=
    (1332.*c2+243.*c4+1217.)*pow(sb,6)*pow(s1,8)
    -2.*(1365.*c2+2854.*c4+1899.*c6+729.*c8-2239.)*pow(sb,5)*pow(s1,4)
    +(18860.*c2-9913.*c4-14990.*c6+12805.*c8-990.*c10+3645.*c12+6967.)*pow(sb,4)
    -64.*(10208.*c2+8809.*c4+8629.*c6+7555.*c8+1035.*c10+1215.*c12+3509.)*pow(sb,3)
    +256.*(175240.*c2+127259.*c4+78808.*c6+27740.*c8+7200.*c10+3645.*c12+103372)*pow(sb,2)
    -32768.*pow(c1,2)*(21756.*c2+14212.*c4+3275.*c6+801.*c8+729.*c10+10939.)*pow(sb,1)
    +3145728.*pow(c1,4)*(1418.*c2-35.*c4+18.*c6+81.*c8+1406.)
    ;
  output *= 1./pow(sb,4)/pow((16.*c1*c1+sb*s1*s1),4);
  output *= wave*wave*wave*wave*wave*wave/mc2/mc2/mc2;

  double fact1=Constants::pi/81.;
  double fact2=16.*Constants::pi;  // used to convert |M|^2=16pi*s^2*dsig/dt, the desired output
  output *= fact1*fact2;

  _flow = 1;
  _diagram=1;
  return output;
}

Energy2 MEgg2JpsiJpsi::scale() const {
  // Energy2 s(sHat());
  // return s;
  return meMomenta()[2].m2()+meMomenta()[2].perp2();
}

void MEgg2JpsiJpsi::persistentOutput(PersistentOStream & os) const {
  os << _meson << _lambda << _process;
}

void MEgg2JpsiJpsi::persistentInput(PersistentIStream & is, int) {
  is >> _meson >> _lambda >> _process;
}

unsigned int MEgg2JpsiJpsi::orderInAlphaS() const {
  return 4;
}

unsigned int MEgg2JpsiJpsi::orderInAlphaEW() const {
  return 0;
}

ClassDescription<MEgg2JpsiJpsi> MEgg2JpsiJpsi::initMEgg2JpsiJpsi;
// Definition of the static class description member.

void MEgg2JpsiJpsi::Init() {

  static ClassDocumentation<MEgg2JpsiJpsi> documentation
    ("The MEgg2JpsiJpsi class implements the non-pertubative gg->JpsiJpsi process in hadron-hadron"
     " collisions");
  
  static Parameter<MEgg2JpsiJpsi,double> interfaceLambda
    ("Lambda",
     "The wavefunction^2 of the ccbar system at the origin",
     &MEgg2JpsiJpsi::_lambda, 0.92, 0.0, 9.2, false, false, false);
  
  // Implemented as a place-holder only
  static Switch<MEgg2JpsiJpsi,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEgg2JpsiJpsi::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all gg->Jpsi Jpsi modes",
     0);
}

Selector<MEBase::DiagramIndex>
MEgg2JpsiJpsi::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is trivial as we are directly computing |ME|^2
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if(diags[i]->id()==-int(_diagram)) sel.insert(1.0, i);
    else sel.insert(0., i);
  }
  return sel;
}

void MEgg2JpsiJpsi::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon=getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(2),gluon,gluon,1,gluon,
  	       3,_meson,3,_meson,-1)));
}

Selector<const ColourLines *>
MEgg2JpsiJpsi::colourGeometries(tcDiagPtr diag) const {
  // colour lines for g g -> Jpsi Jpsi
  static const ColourLines cggs("1 -2, 2 -1");
  // select the colour flow (only one option since the Jpsi's are
  // colour singlets)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case 1:
    sel.insert(1.0, &cggs);
    break;
  }
  return sel;
}

double MEgg2JpsiJpsi::me2() const {
  // total matrix element
  double me(0.);
  if (mePartonData()[0]->id()==ParticleID::g &&
      mePartonData()[1]->id()==ParticleID::g) {
    me = gg2JpsiJpsiME();
  }
  else {
    throw Exception() << "Unknown process in MEgg2JpsiJpsi::me2()" 
		      << Exception::abortnow;
  }

  // multiply by alpha_S^2 and return answer
  double alphas(SM().alphaS(scale()));
  // multiply for _lambda, the overall (non-pertubative) 
  // normalisation factor
  return pow(alphas,4)*sqr(_lambda)*me;
}
