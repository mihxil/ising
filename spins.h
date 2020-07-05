// -*-C++-*-
#ifndef __SPINS_H_
#define __SPINS_H_

// march 1999
// MM 15 october 1998

#include <iostream.h>
#include <rnd.h>
#include <atree.h>
#include "ising_error.h"

/* In this file the classes for an actual system under study are described.

   A system is described by two classes, one which describes the 'spin' on a
   lattice site, and one which describes the properties of such a site.

 */

namespace spins // possible spintypes
{
  
class spinbase
{
public:  
  static rnd::uniform*  draw;
};

class spin : public spinbase
{
public:
  class spinsum
    {
      int value;
  public:
      spinsum(int i) { value = i;}
      int getvalue() const { return value;}
    };  // this is the type of the number of summed spins.
  friend  spinsum  operator+=(spinsum s1, spin s) { throw ising_error("+= operator not defined for this type");}
  spin*  flip() { throw ising_error("flip-operator not defined for this type");}
  friend  spinsum  operator*(spinsum s1, spin s) { throw ising_error("* operator not defined for this type");}
  bool operator==(spin s)  { throw ising_error("== operator not defined for this type");}
   
  
};
 

class updown : public spin
{
 private:
  signed value;
 public:
  // constructor
  updown(signed i = 1) { value = i; }

  // how it counts
  spinsum upvalue() const { return value == 1 ? 1 : 0; }

  // how to flip it  
  updown*  flip() { value = -value; return this; }

  // how to convert this spin to a float:
  //operator float() const{ return (float)value; }

  // how to make a random value:
  bool makerandom(    ) { value = (draw->mdraw(2) == 0) ? -1 : 1; }
  bool makevalue(int i) { value = i<0 ? -1 : 1; }
  // how to set it a value:
  
  // - operator
  updown   operator-() { return -value; } 
  friend  spinsum  operator+=(spinsum s1, updown s) { s1 +=s.value;}
  bool operator==(updown s) { return value == s.value; }

  // how to print this spin:
  friend ostream& operator<<(ostream& os, const updown s)
    {	  
      os << ((s.value == -1) ? "-" : ((s.value == +1) ? "+" : "*") ); 
      return os;
    }
  signed  getvalue() { return value;}
};



/* --------------------------------------------------------------------------------


 */
template <class floatt>
class transverse;


// interfacetype: a help type, it can not be regarded as a spin.
template <class floatt>
class interfacetype
{
public:
  floatt  location;   // where it is located
  
  typedef updown rightspin; // may it could be something else..
  rightspin  right;      // what spin is right from this interface
  bool       two;     // if this interface is a double interface (sign doesn't change on it)
                         // used for tempory interfaces, and for the interface on z=0
public:
  //  interfacetype() { location = 0; right = 0; two = false; } 
  interfacetype(floatt& loc, rightspin r, bool d = false) 
    {
      location = loc; 
      right = r; 
      two = d;
    }
  interfacetype(floatt loc = 0.0) // hopefully this is fast
    {
      location = loc; 
    }

  //  interfacetype(const interfacetype<floatt> &s) { location = s.location; right = s.right; }
  // a few operators have to be overloaded for this type (it must be comparable to fit in a avl::Tree)
  
  template <class f> 
  friend ostream& operator<<(ostream&, const interfacetype<f>&);

  bool operator<(const interfacetype<floatt> & i) const 
    {
      return location < i.location; 
    } 

  // it's better not to compare floatt on equalness.
  // but we need it, because it is used in the avl-tree algorithm
  bool operator==(const interfacetype<floatt> & i) const 
    { 
      return location == i.location; 
    }
    
  template <class f>
  friend ostream& operator<<(ostream&, const transverse<f>& );
  
  floatt gety() { return location; }
  rightspin getspin() { return right; } 

  friend class transverse<floatt>;
};

template<class floatt>
ostream& operator<<(ostream& os , const interfacetype<floatt>& a)
{
  os << a.location << "(" << a.right << ")"; 
  return os;
}
//------------------------------------------------------


#define LISTTYPE  atree::with_array<interfacetype<floatt>, 100 >
//#define LISTTYPE  atree::with_avl<interfacetype<floatt> >



template <class floatt> 
//class transverse : public spin, public avl::Tree<interfacetype<floatt> > 
class transverse : public spin, public  LISTTYPE
  {   
    template <class floattt>
    class  spinsum_ : public spinsum
    {
      floattt value; // the number of spins up now of course can be a floatt.
      unsigned ninterfaces;
    public:
      spinsum_(floattt x = 0, unsigned i = 0) { value = x; ninterfaces = i;}
      floattt getvalue() const { return value; }
      floattt getninterfaces() const { return ninterfaces; }
      friend  void  operator+=(spinsum_<floattt>& s1, floattt s) { s1.value += s; }
      friend  void  operator+=(spinsum_<floattt>& s1, spinsum_<floattt>& s2) { s1.value += s2.value; }
      void operator++() { ninterfaces++; }
      
      friend ostream& operator<<(ostream& os, const spinsum_<floattt>& s) { os << s.value; return os;} 
    };
    
    typedef spinsum_<floatt> spinsumt;
    typedef interfacetype<floatt> interface;
    
    typedef typename interface::rightspin  attype;

    static floatt  length;   // 'a3', all spins have the same length, therefore I made it static

  public:
    //interfacetype<floatt> left, right; // interfaces on the left and right hand of the area, so we can easily find where they are.

    // overload constructors
    transverse () : LISTTYPE()
    {
      floatt zero = 0;
      interfacetype<floatt> left(zero,1,true), right(length,1,true); // interfaces on the left and right hand of the a
      Insert(left);
      Insert(right);      
    };
    /*
    transverse (interface * i) : LISTTYPE() 
    {
      floatt zero = 0;
      interfacetype<floatt> left(zero,1,true), right(length,1,true); // interfaces on the left and right hand of the a
      Insert(left);
      Insert(right);      
    };
    */
    static void setlength(floatt x) 
    {
      length = x;     
    }
    static floatt getlength() { return length; }
    static floatt size() { return length;}

    attype at(floatt& y) 
      { 
	interface x(y);
	SSearch(x);
	return GetSmaller().right;
     } 
    /*
    void  GetNearestInterfaces(floatt y, interface ** l, interface **h)
      {
	SSearch(interface(y));
	*l = GetSmaller();
        *h = GetBigger();
	return;
      }
    */
    static floatt multfac; // for drawing purposes, the ratio y/x of a character

    spinsumt upvalue()
      {
	spinsumt hulp = 0;
	interface * f;
	interface * of;
        typename interface::rightspin cur;
	
	f = Lowest();
	while(f)
	  {
	    of = f;
	    cur = of->right;
	    f = Next();
	    if(f && cur.upvalue()) hulp += (f->location - of->location);
	    hulp++; // increase number of interfaces
	  }
	if(cur.upvalue()) hulp+= length - of->location;
        
	return hulp; 
      }// returns part which is up.
    floatt getvalue() { return 1;}

    bool  makerandom() 
    { 
      Smallest().right  = (draw->mdraw(2) == 0 ? -1 : 1);
      Biggest().right = Smallest().right;
      return true; 
    }
    bool makevalue(int i)
    {
      Smallest().right  = (i < 0 ? -1 : 1);
      Biggest().right = Smallest().right;
      return true;       
    }

    template <class f>
    friend ostream& operator<<(ostream&, const transverse<f>&);
  };



template<class floatt>
floatt transverse<floatt>::length = 1; // don't make it 0, because in that case the right interface won't be inserted


template<class floatt>
ostream& operator<<(ostream& os , const transverse<floatt>& a)
{
  

  float column = 0;

  typename transverse<floatt>::interface * f;
  typename transverse<floatt>::interface * of;

  f = a.Lowest();    
  while(f)
    {      
      of = f;
      //cur = of->right;
      f = a.Next();
      if(f) 
         {
	   os << '|'; column+=1; // write explicitily the interface, in this way we more easily remark bugs
	   // in principle interface always divide areas with different spin.
	   for(float i = column; i < (f->location) * a.multfac; i+=1) { os << of->right; column+=1;}
	 }
    }
  os << '|'; column +=1;
  for(float i = column; i<  (transverse<floatt>::length) * a.multfac; i+=1) os << of->right; 
  os << a.Count() - 2  << endl; // number of interfaces (but not the two on the ends)

  return os;
};

 


} // namespace spins

#endif // __SPINS_H_
