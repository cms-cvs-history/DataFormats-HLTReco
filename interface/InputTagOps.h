#ifndef FWCore_Utilities_InputTagOps_h
#define FWCore_Utilities_InputTagOps_h

#include <string>

#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {

  class InputTagEquality {
  public:
    bool operator()(const InputTag & left, const InputTag & right) const {
      return left.process()  == right.process() && 
             left.label()    == right.label()   &&
             left.instance() == right.instance();
    }
  };

  class InputTagCompare {
  public:
    bool operator()(const InputTag & left, const InputTag & right) const {
      return left.process()  < right.process() && 
             left.label()    < right.label()   &&
             left.instance() < right.instance();
    }
  };

}

#endif // FWCore_Utilities_InputTagOps_h
