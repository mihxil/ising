
#include "spins.h"

namespace spins
{
#ifndef COMPILE_TEMPLATES 
  /*
  ostream& operator<<(ostream& os, const updown s)
    {	  
      os << ((s.value == -1) ? "-" : ((s.value == +1) ? "+" : "*") ); 
      return os;
    }
  */
  rnd::uniform* spinbase::draw = NULL;
#endif 

#ifdef COMPILE_TEMPLATES
  template <class floatt = float> 
    floatt  transverse<floatt>::multfac = 2; // for drawing purposes, the ratio y/x of a character
#endif
 
}
