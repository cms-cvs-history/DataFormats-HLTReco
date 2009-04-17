#ifndef FWCore_Utilities_InputTagHash_h
#define FWCore_Utilities_InputTagHash_h

#include <boost/functional/hash.hpp>

#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {

  inline
  std::size_t hash_value(const InputTag & tag) {
    std::size_t seed = 0;
    boost::hash_combine(seed, tag.label());

    // pattern copied from InputTag::encode()
    if (not tag.instance().empty() || not tag.process().empty()) {
      boost::hash_combine(seed, ':');
      boost::hash_combine(seed, tag.instance());
    }
    if (not tag.process().empty()) {
      boost::hash_combine(seed, ':');
      boost::hash_combine(seed, tag.process());
    }

    return seed;
  }

  typedef boost::hash<InputTag> InputTagHash;

}

#endif // FWCore_Utilities_InputTagHash_h
