#include "rnd.h"
#include "spins.h"
#include "lattice.h"
#include "sites.h"
#include <iostream>


template <int n, class st> int sites::cnD<n,st>::siz[n];
main()
{
  rnd::exponential draw(- abs(1) -1 );
  // copy a pointer to this random number generator to the classes
  spins::spin::draw = &draw; 
  sites::sitet::draw = &draw;

  sites::simple3D::setsizei(0,4);
  sites::simple3D::setsizei(1,4);
  sites::simple3D::setsizei(2,4);
  
  lattice<sites::simple3D> l;
  l.init();

  typename sites::simple3D::pointer nb;

  for(int i = 0; i < sites::simple3D::max(); i++)
    {
      cout << l.r[i] <<  "      " << l.r[i].spin << endl;
      for(int j = 0; j < sites::simple3D::maxnb(); j++)
	{
	  nb = l.r[i].nb[j];
	  cout << "    " << *nb <<  "   " << nb->spin;
	  if(nb->spin == l.r[i].spin) cout << " zelfde" << endl;
	  else    cout << " anders" << endl;
	}
    }
  
    
  /*
  sites::transverse1D::setsizei(0,10);
  lattice<sites::transverse1D> l;
  l.init();
  cout << l;  
  
  const int test = 0;

  cout << "-------" << endl << l.r[test] << endl << "-----" << endl;
  for(int i = 0; i< sites::transverse1D::maxnb(); i++)
    cout << *(l.r[test].nb[i]) << endl;

  */
}
