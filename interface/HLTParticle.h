#ifndef HLTReco_HLTParticle_h
#define HLTReco_HLTParticle_h

/** \class HLTParticle
 *
 *
 *  The basic information (4-momentum etc.) on a reconstrcuted physics
 *  object (e/gamma/mu/jet/Met...) for those physics objects used in
 *  an HLT filter decision.
 *
 *
 *  $Date: 2006/05/18 17:00:05 $
 *  $Revision: 1.2 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include <cassert>

namespace reco
{

  class HLTParticle : public Particle {
    // Currently, we re-use the Particle class from the Candidate Model
    // Later, we may use our own (more compact) particle (base) class

  public:
    HLTParticle(): Particle(), id_() { }
    HLTParticle(const Particle& p, int i=0) : Particle(p), id_(i) {
      assert ( (-128<=i) && (i<=127) ) ; }

    int id() const { return (int) (id_); }

  protected:
    char id_; // id of particle (PDG scheme - fundamental particles only!)
  };

}

#endif
