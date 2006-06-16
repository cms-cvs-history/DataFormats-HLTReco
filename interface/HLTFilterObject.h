#ifndef HLTReco_HLTFilterObject_h
#define HLTReco_HLTFilterObject_h

/** \class HLTFilterObject
 *
 *
 *  If HLT cuts of intermediate or final HLT filters are satisfied,
 *  instances of this class hold the combination of reconstructed
 *  physics objects (e/gamma/mu/jet/MMet...) satisfying the cuts.
 *
 *  This implementation is not completely space-efficient as some
 *  physics object containers may stay empty. However, the big
 *  advantage is that the solution is generic, i.e., works for all
 *  possible HLT filters. Hence we accept the reasonably small
 *  overhead of empty containers.
 *
 *  $Date: 2006/05/25 16:52:20 $
 *  $Revision: 1.7 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Common/interface/HLTenums.h"
#include "DataFormats/HLTReco/interface/HLTParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <cassert>
#include <vector>
#include <map>

namespace reco
{
  using namespace std;

  class HLTFilterObjectBase {

  typedef edm::hlt::HLTScalar HLTScalar;

  private:
    bool accept_;                     // filter decision
    unsigned char module_;            // mdoule index of filter on path
    unsigned short int path_;         // path index of path in trigger tabge (cfg file)

  public:

    HLTFilterObjectBase(): accept_(), module_(), path_() { }

    bool getAccept() const { return accept_;}
    void setAccept(const bool accept) {accept_=accept;}

    unsigned int getModule() const {return (unsigned int)(module_);}
    void setModule(const unsigned int i) {assert(i<  256); module_=i;}

    unsigned int getPath()   const {return (unsigned int)(path_  );}
    void setPath  (const unsigned int i) {assert(i<65536); path_  =i;}
  };

  class HLTFilterObject : public HLTFilterObjectBase {

  typedef edm::hlt::HLTScalar HLTScalar;

  private:
    map<HLTScalar,float>  scalars_;   // scalar quantities used in filter (HT, ...)
    vector<HLTParticle>   particles_; // particles/MET (4-momentum vectors) used by filter

  public:

    HLTFilterObject(): HLTFilterObjectBase(), scalars_(), particles_() { }

    unsigned int numberScalars  () const { return   scalars_.size();}
    unsigned int numberParticles() const { return particles_.size();}

    bool getScalar(const HLTScalar scalar, float& value) const {
      if (scalars_.find(scalar)==scalars_.end()) {
        return false;
      } else {
        value = scalars_.find(scalar)->second;
        return true;
      }
    }
    void putScalar(const HLTScalar scalar, const float value) {
      scalars_[scalar] = value;
    }

    bool getParticle(const unsigned int i, HLTParticle& particle) const {
      if (i<particles_.size()) {
        particle = particles_[i];
        return true;
      } else {
        return false;
      }
    }

    void putParticle(const edm::RefToBase<Candidate>& ref) {
      particles_.push_back(HLTParticle(*ref));
    }

  };


  class HLTFilterObjectWithRefs : public HLTFilterObject {

  private:
    std::vector<edm::RefToBase<Candidate> > refs_;

  public:

    HLTFilterObjectWithRefs(): HLTFilterObject(), refs_() { }

    void putParticle(const edm::RefToBase<Candidate>& ref) {
      this->HLTFilterObject::putParticle(ref);
      refs_.push_back(ref);
    }

    bool getParticleRef(const unsigned int i, const Candidate* & candidate) const {
      if (i<refs_.size()) {
        candidate = (refs_[i]).get();
	return true;
      } else {
	return false;
      }
    }
 
  };
}

#endif
