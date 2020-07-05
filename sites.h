// -*-C++-*-
#ifndef __SITES_H_
#define __SITES_H_

#include <string>
#include "rnd.h"
#include "spins.h"
#include "ising_error.h"

namespace sites // possible sitetypes
{

class sitet 
{
protected:
  static unsigned maxindex;
public:
  virtual unsigned index() const { return 0;}
  virtual unsigned operator()() const { return  index();}
  virtual const char * lprint() const { return  "";}    
  static unsigned max() { return maxindex; }
  static int maxnb() { return 0;}
  static unsigned random() { return draw->mdraw(maxindex); } // returns random site index.
  static rnd::uniform* draw;  
  static ostream&  extrainform(ostream& os){ return os; }
  static void  setmultfac(float f) { return;}
};
  
  // base class for sites:
template <class spinT>
class site : public sitet
{
 public:
  spinT spin;  
  typedef spinT spintype;
  static typename spintype::spinsum size() { return spintype::size(); }

};

// This one is undoubtable less efficient, but it's handy to have it of course.
template <int n, class st> // = typename spins::updown>
class cnD : public site <st>
{
  unsigned m_index;

public:
  typedef cnD<n, st> * pointer;
  
  //int x[n];         // if we have the neighbours, then it might be even not necessary to store the location
  
  pointer nb[2*n]; // pointers to the neighbours (should be filled from external)

  static unsigned siz[n];

   cnD(int i = 0)     
    {
     setindex(i);
    } 

 
  void setindex(unsigned i)
    {
     
      m_index = i;
    }
 
  unsigned getx(int i) const
    {
      unsigned x1[n];         // x's
      unsigned in = m_index;
      for(int j = 0; j< n; j++) 
	{
	  //	  x1[j] = x[j];
	  x1[j] = in % siz[j];
	  in   /= siz[j];
	}
      return x1[i];
    }

  unsigned findnb(int i) const // finds the index of a neighbor
    /* not a very efficient function, but is used only once 
       The answers should be stored in nb[2*n]
     */
    {
      // convert the index to a vector
      // this code is repeated a little, that's a pitty.
      int x1[n];         // x's
      unsigned in = m_index;
      for(int j = 0; j< n; j++) 
	{
	  //	  x1[j] = x[j];
	  x1[j] = in % siz[j];
	  in   /= siz[j];
	}

      // find the neigbor
      int dimension = i / 2; // in which dimension it is
      int plusormin = (i % 2) * 2 - 1; // +/- 1
      
      x1[dimension] += plusormin + siz[dimension];
      x1[dimension] %= siz[dimension]; // periodic bounday conditions      

      // convert it to an index:
      in = x1[n-1]; 
      for(int j = n - 2; j >= 0; j--)
	{ 
	  in = siz[j] * in + x1[j];
	} 
      return in;
      
    }

  unsigned index() const
    {
      return m_index;
      /*
      unsigned in = x[n-1]; 
      for(int i = n - 2; i >= 0; i--)
	{ 
	  in = siz[i] * in + x[i];
	} 
      return in;
      */
    }

  static setsizei(int i, int x) 
    {
      if(i >= n) throw ising_error("Specified index greater than dimension allows");
      siz[i] = x;
      // determine new maxsize
      int h = 1;
      for(int i = 0; i < n; i++) h*= siz[i];
      maxindex = h;
    }

  static int maxnb()     { return 2 * n; }
  static int dimension() { return n; }


  template <int tn, class tst>
  friend ostream& operator<<(ostream&, const cnD<tn,tst> &);  
    
};

template <int n, class st>
ostream& operator<<(ostream& os, const cnD<n, st>& s)
{
  os << "(";
  for(int i = 0; i<n - 1; i++) os << s.x[i] << ",";
  os << s.x[n -1] <<  ")"; 
  return os;
}

//--------------------------------------------------------


    
/* ---------------------------------------------------------------------------
   
   1D 

 */

class simple1D : public cnD<1, spins::updown>
{
public:
  simple1D(int n = 0) : cnD<1, spins::updown>(n) {}  
  
  //  friend ostream& operator<<(ostream&, const simple1D &); 
  
  const char * lprint() const { return getx(0) == siz[0]-1 ? "\n" : "";}
};
/*
ostream& operator<<(ostream& os, const simple1D& s)
{
 os << "(" << s.x[0] << ")"; 
 return os;
} 
*/



/* --------------------------------------------------------------------
class simple2D, is an example of which the lattice-class type can use for it sizes.

The idea is that if you want to create another type of lattice, you only need to redefine 
a new type for its sites, in which you should define at least the same memberfunctions.

 */
class simple2D : public cnD<2, spins::updown>
{    
 public:
  simple2D(int n = 0) : cnD<2, spins::updown>(n) {}  
  
  const char * lprint() const { return getx(0) == siz[0] - 1 ? "\n" : "";}
 
};


class  simple3D : public cnD<3, spins::updown>
{
public:
  simple3D(int n = 0) : cnD<3, spins::updown>(n) {}  
  const char * lprint()  const { return getx(0) == siz[0] - 1 ? getx(1) == siz[1] -1 ? "\n\n" : "\n" : "";}

};



/* -----------------------------------------------------------------------------------------
   Transverse sites.
   One of the coordinates is a floating point number

 */

template <class floatt>
class c_location
  {
  public:
    unsigned i;
    floatt   z;    
  public:
    c_location(rnd::uniform* draw, unsigned max, floatt f)
      {
	i = draw->mdraw(max);
	z = draw->fdraw() * f;
      }
    c_location(unsigned ii, floatt zz)
      {
	i = ii;
	z = zz;	
      }
   
    template <class f>
    friend ostream& operator<<(ostream&, const c_location<f> &);  

  };
  
template <class f>
ostream& operator<<(ostream& os, const c_location<f>& s)
{
  os << "(" << s.i << "," << s.z << ")";
  return os;
}



template <int n, class floatt>
class ctransverse : public cnD<n, typename spins::transverse<floatt> > 
{
  typedef spins::transverse<floatt>::interface interface;
  typedef typename spins::transverse<floatt> spintype;
  //  typedef c_location<n, floatt> location; 
  
public:
  ctransverse(int i = 0) : cnD<n, typename spins::transverse<floatt> >(i)
    //    :  site<typename spins::transverse<floatt> >::spin(length)
    {
    }  
  static void setsize(floatt x)
    {
      spins::transverse<floatt>::setlength(x);
    }
  static floatt getsize()
    {
      return spins::transverse<floatt>::getlength();
    }

  static ostream&  extrainform(ostream& os){ os <<  spins::transverse<floatt>::getlength(); return os; }
  static void  setmultfac(float f) { spintype::multfac = f; }
};


class transverse2D : public ctransverse<2, double>
{
public:
  transverse2D(int i = 0) : ctransverse<2, double>(i) {}  
  const char * lprint()  const { return getx(0) == siz[0] - 1 ? getx(1) == siz[1] -1 ? "\n\n" : "\n" : "";}
};

template <class floatt>
class transverse1D : public ctransverse<1, floatt>
{
public:
  transverse1D(int i = 0) : ctransverse<1, float> (i) {}  
};

}// namespace sites
#endif // __SITES_H_
