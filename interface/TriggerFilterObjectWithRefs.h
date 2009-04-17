#ifndef HLTReco_TriggerFilterObjectWithRefs_h
#define HLTReco_TriggerFilterObjectWithRefs_h

/** \class trigger::TriggerFilterObjectWithRefs
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
 *  $Date: 2009/03/17 16:21:34 $
 *  $Revision: 1.13.2.1 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTagHash.h"
#include "FWCore/Utilities/interface/InputTagOps.h"

#include <string>
#include <vector>
#include <algorithm>
#include <ext/hash_set>

namespace trigger
{

  /// Transient book-keeping EDProduct filled by HLTFilter module to
  /// record physics objects firing the filter (never persistet in
  /// production; same functionality but different implementation
  /// compared to the old HLT data model's HLTFilterObjectWithRefs
  /// class)
  class TriggerFilterObjectWithRefs : public TriggerRefsCollections {
  public:
    typedef __gnu_cxx::hash_set<edm::InputTag, edm::InputTagHash, edm::InputTagEquality> InputTagSet;
    //typedef std::set<edm::InputTag, edm::InputTagCompare> InputTagSet;

  /// data members
  private:
    int  path_;
    int  module_;
    bool mask_;

    static InputTagSet s_collectionTags;

  /// methods
  public:
    /// constructors
    TriggerFilterObjectWithRefs():
      TriggerRefsCollections(),
      path_(-9),
      module_(-9),
      mask_(false) {}
    
    TriggerFilterObjectWithRefs(int path, int module):
      TriggerRefsCollections(),
      path_(path),
      module_(module),
      mask_(false) {}
    
    /// accessors
    int  path()   const { return path_; }
    int  module() const { return module_; }
    bool mask()   const { return mask_; }
    
    /// collectionTags
    void addCollectionTag(const edm::InputTag& collectionTag){
      s_collectionTags.insert(collectionTag);
      mask_ = true;
    }

    static    
    const InputTagSet & getCollectionTags() {
      return s_collectionTags;
    }

    /// utility
    void swap(TriggerFilterObjectWithRefs & other) {
      TriggerRefsCollections::swap(other);                  // swap base instance
      std::swap(path_,   other.path_);
      std::swap(module_, other.module_);
    }

  };

  // picked up via argument dependent lookup, e-g- by boost::swap()
  inline void swap(TriggerFilterObjectWithRefs & first, TriggerFilterObjectWithRefs & second) {
    first.swap(second);
  }

}

#endif
