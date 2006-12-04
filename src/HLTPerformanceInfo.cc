// $Id: HLTPerformanceInfo.cc,v 1.1 2006/11/29 16:12:09 wittich Exp $
#include <functional>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <boost/lambda/lambda.hpp> 
#include <boost/lambda/bind.hpp> 
#include <boost/bind.hpp> 
using namespace boost;

#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"


// THIS IS MY BABY!
// iterator adaptor to allow you to do some operation to 
// the return value of a member function.  intended for
// things like std::accumulate
template<class V, class Res, class binary_op>
class BinaryOpMemFun {
 private:
  // pointer to a member funtion of V which has no arguments and returns
  // Res.
  Res (V::* _s)() const;
 public:
  BinaryOpMemFun(Res (V::*s)() const ) :
    _s(s) {}
  // this is for for_each
  Res operator()( V & v) const {
    // this is the function call
    return (v.*_s)();
  }
  // this is for accumulate
  Res operator()(Res b,const  V & v) const {
    // it applies binary_op(b, v->_s()) and returns the result
    // e.g. if you use plus it's b+ v->_s()
    return binary_op()((v.*_s)(), b);
  }
};


HLTPerformanceInfo::HLTPerformanceInfo()
{
  paths_.clear(); modules_.clear();
}

double HLTPerformanceInfo::Path::time() const
{
  double t = 0;
  for ( HLTPerformanceInfo::Path::const_iterator j = this->begin();
	j != this->end(); ++j ) {
    t += j->time();
  }
  return t;
}


void HLTPerformanceInfo::addModuleToPath(const char *mod, Path *p ) 
{
  Modules::const_iterator m = this->findModule(mod);
  if ( m != modules_.end() ) {
    size_t a = m - modules_.begin();
    p->addModuleRef(a);
  }
  else {
    // if we can't find the module, it probably just wasn't run. 
    // so no worries.
    Module newMod(mod, 0); // time = 0 since it wasn't run
    modules_.push_back(newMod);
    p->addModuleRef(modules_.size()-1); // last guy on the stack
  }
}

void HLTPerformanceInfo::addPath(Path & p )
{
  // need this to get back at the modules that we don't own
  p.setModules_(&modules_);
  paths_.push_back(p);
}

// class ugh {
//   HLTPerformanceInfo::* s_;
//   double operator()(double t, HLTPerformanceInfo::Module & a) {
//     t += (p.*s)();
//   }
// }


double HLTPerformanceInfo::totalTime() const
{
  double t = 0;
  t = std::accumulate(beginModules(), endModules(), 0.,
		      BinaryOpMemFun<HLTPerformanceInfo::Module, double,
		      std::plus<double> >(&HLTPerformanceInfo::Module::time));
//   double testt = t;
//   for ( Modules::const_iterator i = beginModules();
//         i != endModules(); ++i ) {
//     t += i->time();
//   }
//   assert(std::abs(t-testt)<0.001);
  return t;
}

HLTPerformanceInfo::Modules::const_iterator 
HLTPerformanceInfo::findModule(const char* moduleInstanceName) 
{
  return std::find(modules_.begin(), modules_.end(),
		     moduleInstanceName);
}

HLTPerformanceInfo::PathList::const_iterator 
HLTPerformanceInfo::findPath(const char* pathName) 
{
  PathList::const_iterator l = std::find(paths_.begin(), paths_.end(),
					 pathName);
  if ( l != endPaths() ) {
    return l;
  }
  else {
    return endPaths();
  }
}


double HLTPerformanceInfo::Path::lastModuleTime() const
{
  double prev_time = -1;
  for ( HLTPerformanceInfo::Path::const_iterator j = this->begin();
	j != this->end(); ++j ) {
    if ( j->status().wasrun() && !(j->status().accept()) )
      return prev_time;
    prev_time = j->time();
  }
  return -2; // no modules on the path
}


double HLTPerformanceInfo::longestModuleTime() const
{
  double t = -1;
//   t = std::accumulate(beginModules(), endModules(), -99, 
// 		      BinaryOpMemFun<HLTPerformanceInfo::Module, double,
// 		      std::max >(&HLTPerformanceInfo::Module::time));
  for ( Modules::const_iterator i = beginModules();
        i != endModules(); ++i ) {
    t = std::max(i->time(),t);
  }
  return t;
}

const char* HLTPerformanceInfo::longestModuleTimeName() const
{
  double t = -1;
  std::string slowpoke("unknown");
  for ( Modules::const_iterator i = beginModules();
        i != endModules(); ++i ) {
    if ( i->time() > t ) {
      slowpoke = i->name();
      t = i->time();
    }
  }
  return slowpoke.c_str();
}
