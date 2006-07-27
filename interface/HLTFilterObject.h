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
 *  $Date: 2006/06/24 21:04:46 $
 *  $Revision: 1.13 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/Common/interface/RefToBase.h"
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

  private:
    unsigned short int path_;    // index of path in trigger table
    unsigned short int module_;  // index of module on trigger path

  public:

    HLTFilterObjectBase(int path=0, int module=0) {
      assert( (0<=path  ) && (path  < 65536) );
      assert( (0<=module) && (module< 65536) );
      path_   = path;
      module_ = module;
    }

    unsigned short int path()   const { return path_  ;}
    unsigned short int module() const { return module_;}

  };

  class HLTFilterObject : public HLTFilterObjectBase {

  private:
    vector<HLTParticle>   particles_; // particles/MET (4-momentum vectors) used by filter

  public:

    HLTFilterObject(int path=0, int module=0)
      : HLTFilterObjectBase(path,module), particles_() { }

    unsigned int size() const { return particles_.size();}

    const HLTParticle& getParticle(const unsigned int i) const {
      return particles_.at(i);
    }

    void putParticle(const edm::RefToBase<Candidate>& ref) {
      particles_.push_back(HLTParticle(*ref));
    }

  };


  class HLTFilterObjectWithRefs : public HLTFilterObject {

  private:
    std::vector<edm::RefToBase<Candidate> > refs_;

  public:

    HLTFilterObjectWithRefs(int path=0, int module=0)
      : HLTFilterObject(path,module), refs_() { }

    void putParticle(const edm::RefToBase<Candidate>& ref) {
      this->HLTFilterObject::putParticle(ref);
      refs_.push_back(ref);
    }

    const edm::RefToBase<Candidate>& getParticleRef(const unsigned int i) const {
      return refs_.at(i);
    }
 
  };
}

#endif
