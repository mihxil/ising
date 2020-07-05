#ifndef __LATTICE_H_
#define __LATTICE_H_

#include <cstddef>
#include <iostream>
#include "rnd.h"
#include <string>
#include "sites.h"


template <class site>
class lattice
{
  // data members:
 public: 
  site * r;

 public:
  //constructor:
  lattice();
  void  init();
  void  init(int i);
  void  print(ostream& os = cout) const;
  template <class st> friend ostream& operator<<(ostream&, const lattice<st>);
  typename site::spintype::spinsum  nup()  const;
  typename site::spintype::spinsum  size() const { return site::max() * site::size().getvalue(); }
  ostream& sizes(ostream&) const;

};

/*
  Implementation of this template-class probably should be in this header file,
  because the compiler cannot know for which types it should generate code.

  These _are not_ real implementations, of course. Because they are only compiled if 
   necessary.


 */


template <class site> lattice<site>::lattice()
{
  int i;
  r = new site[site::max()];
  // Make locations of the sites:
  for(i = 0; i< site::max(); i++) r[i].setindex(i);
  // Make neighbors:
  for(i = 0; i< site::max(); i++) 
    {
      for(int j = 0; j < site::maxnb(); j++)
	{
	  r[i].nb[j] = &(r[r[i].findnb(j)]);
	}
    }

}

template <class site> void lattice<site>::init()
// initialiseert het rooster r met willekeurige waarden.
{
  for(int i=0; i<site::max(); i++) r[i].spin.makerandom();
  return;  
}

template <class site> void lattice<site>::init(int v)
// initialiseert het rooster r met v
{
  for(int i=0; i<site::max(); i++) r[i].spin.makevalue(1);
  return;
}

template <class site> 
void lattice<site>::print(ostream& os) const
{
  for(int i=0; i<site::max(); i++) os << r[i].spin << (const char *) r[i].lprint();     
  return;
}

template <class site> 
ostream& operator<<(ostream& os, const lattice<site> l)
{
 l.print(os); 
 return os;
}

template <class site> 
typename site::spintype::spinsum  lattice<site>::nup() const
{ 
  typename site::spintype::spinsum hulp = 0;
  for(int i=0; i<site::max(); i++) hulp+= r[i].spin.upvalue();
  return hulp;
}


template <class site> 
ostream&  lattice<site>::sizes(ostream& os) const
{ 
  os << "(";
  for(int i=0; i<site::dimension() - 1; i++) os << site::siz[i]  <<  ",";
  os << site::siz[site::dimension()-1]  << ") "; 
  site::extrainform(os);
  return os;
}

#endif //__LATTICE_H_
