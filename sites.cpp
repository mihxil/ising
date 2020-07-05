#include "sites.h"


namespace sites
{

  rnd::uniform * sitet::draw = NULL;

  unsigned sitet::maxindex = 0;

  template <int n, class st> unsigned cnD<n,st>::siz[n];
  

  /*
  ostream& operator<<(ostream& os, const simple3D s)
    {
      os << "(" << s.x[0] << "," << s.x[1] << "," << s.x[2] << ")"; 
      return os;
    }
  */
  /*  
  template <int n, class st>
  ostream& operator<<(ostream& os, const cnD<n, st>& s)
    {
      os << "(";
      for(int i = 0; i<n - 1; i++) os << s.x[i] << ",";
      os << s.x[n -1] <<  ")"; 
      return os;
    }
  */
}
