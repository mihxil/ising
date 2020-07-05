/* ------------------------------------------------------
   Ising Model 
   
   DESCRIPTION

// 	$Id: ising.cpp,v 1.1 1999/05/07 15:02:56 michiel Exp michiel $	

   
   ------------------------------------------------------
*/
    
#include <iostream>
#include "ising.h"
#include "spins.h"
#include "sites.h"
#include <astack.h>

template <class floatt, class sitetype>
ostream& operator<<(ostream& os , const area<floatt, sitetype>& a)
{	  
  os << "(";
  if(a.smaller) os << a.smaller;
  else    os << "NULL";
  os <<  ";";
  if(a.bigger) os << a.bigger;
  else   os << "NULL";
  os << ";" << a.i << ")";
  return os;
  // strange construction is left from the time that a.smaller and a.bigger were pointers.
  // for the case if I change it back...
}

// 
template <class floatt,  class sitetype> 
void  isingmodel<floatt, sitetype>::summarize(ostream& os)
{
  sizes(os) << " Ising model  at K = " << m_K << endl 
	    << "number of spins up: " << nup() << "\n";

}
/* ---------------------------
   Metropolis algoritme.
   'Ising Model'
   
   zie ook: Reif, F. - Fundamentals of Statistical and Thermal Physics 
                     - McGraw-Hill 1985
*/


template <class floatt, class sitetype> 
void isingmodel<floatt, sitetype>::MetropolisIteration()
{ 
  /*
  unsigned i;

  // chose random site (metropolis):
  i = m_draw->mdraw(sitetype::max());  

  // calculate energy difference if we would flip it:
  static floatt DeltaE;    //static to prevent realocation

  typename sitetype::spintype::spinsum sumnb = 0;
  //
  for(int j=0; j< sitetype::maxnb(); j++)
    {
      sumnb += r[i].nb[j]->spin.getvalue();
    }
    
  DeltaE = -2 * r[i].spin * ( m_Jx * sumnb.getvalue() - m_B );
  
  // assign new value for this spin:

  if(DeltaE < 0) // if DeltaE x> 0, we do nothing
    { 
      if(m_draw->draw()>=1.0/(1+exp(+2*DeltaE))) // misschien nog niet goed.
	{ 
	  r[i].spin.flip();
	}
    }
 
  m_nit++; // number of iterations
  */
}

/* ------------------------------------------------------------------------
   Swendsen / Wang 
   bond algoritme
   
  zie ook Phys. Rev. Letters {\bf 58} (1987) 86

*/

template <class floatt, class sitetype> 
void isingmodel<floatt, sitetype>::SwendsenWangIteration()
{
// not yet implemented
  /*
    There used to be an implementation which was recursive, but it doesn't fit in the schema now.
    
  */

}
/* ------------------------------------------------------------------------

   Wolff algorithm

   see also Phys. Rev. Letters {\bf 62} (1989) 361
   */

template <class floatt, class sitetype> 
void isingmodel<floatt,  sitetype>::WolffIteration()
{
  typename sitetype::pointer currentsite, neighbor;
  unsigned j;

  // chose random lattice site:  
  unsigned random = sitetype::random();
  currentsite = &(r[random]);
  
  //flip it:
  (currentsite->spin).flip();

  //We need a stack to put positions in.
  // Astack is a fast stack type, see astack.h  
  // by the way, it seems hardy lo matter which stacktype is used
  // even with 'faststack' the program doesn't run faster at all.
  static Astack<typename sitetype::pointer> thestack;

  thestack.push(currentsite); // we want to start with our start position on the stack, because it makes our loop simpler.   

   // !! Not yet fully optimized for efficiency.
  while(!thestack.empty())
    {
      currentsite = thestack.top(); thestack.pop();
      // looking at the neighbors
      for(j = 0; j < sitetype::maxnb(); j++)      
	{
	  neighbor = currentsite->nb[j];
	  if(!((neighbor->spin) == (currentsite->spin)))
	    { 
	      // decide to add it to the stack or not:	  	      
	      if(m_draw->draw() < m_one_minus_exp) // with probability 1 - exp(-K)
		{ 
		  neighbor->spin.flip(); //  flip this neighbor		  
		  thestack.push(neighbor);   //  and push it on the stack
		}
	    }
	} // for neighbors    
    } // while stack not empty
  m_nit++; // count iterations (not really necessary)
}// Wolff

/* ============================================================================
   Anisotropic/Transverse Wolff algorithm

   It seems to me impossible to do it with the same code, so that's why I create
   it again  here..
   ============================================================================
 */

/*
  JumpTransverseDirection
  This help function deals with the transverse part
*/

template <class floatt, class sitetype> 
area<floatt, sitetype> isingmodel<floatt, sitetype>::JumpsTransverseDirection
( sites::c_location<floatt>& i,  
  typename sitetype::spintype::attype& spinofcluster
  )
{
  //*debug */ static unsigned count = 0;
  //*debug */ cout << "going to handle (" << count++ << ") " << *i << endl;
   // get nearest interfaces:
  typename sitetype::interface iz(i.z);
  r[i.i].spin.SSearch(iz);

  area<floatt, sitetype> ca(r[i.i].spin.GetSmaller().gety(), 
			    r[i.i].spin.GetBigger().gety(), i.i);  // abreviation of 'current area'
  // we always find a smaller one. (namely on zero)
  // we also always find a bigger one. (namely on the other end)

  // We have to check on the boundaries.
  // boundary conditions:
  if(ca.smaller == r[i.i].spin.Smallest().gety()) // the smaller one is the left pseudo-interface
    {
      // look for the biggest interface
      r[i.i].spin.SSearchBiggest();     
      ca.smaller = r[i.i].spin.GetSmaller().gety(); // the really biggest is the other pseudo-interace
      // if there were no interfaces, it is still the left-inteface, otherwise it's the most right real one
    }
  // analogously for the other side:
  if(ca.bigger == r[i.i].spin.Biggest().gety())
    {
      r[i.i].spin.SSearchSmallest();
      ca.bigger = r[i.i].spin.GetBigger().gety();
    }
  //--------------------------------------------------------------------------
  // We need the distance to these two interfaces:
  floatt distl = i.z - ca.smaller; 
  floatt distr = ca.bigger - i.z;  
  // in case the boundary conditions applied then these distances are negative
  // we correct it here:
  if(distl < 0)  distl += sitetype::getsize(); 
  else  // they can't be both smaller than 0
    if(distr < 0) distr += sitetype::getsize();
  
  //---------------------------------------------------------------------------
  // a jump; an expontial distributed number:
  // static, because we only draw one, when we used one, and of course the first time, so it is initialised
  /*static*/ floatt jump = m_draw->edraw();  
  //it's not static in ani
 
  // first the left interface:
  if(jump < distl || ca.smaller == r[i.i].spin.Smallest().gety()) // add an interface
    {
      iz.location = i.z - jump;   // calculate new z:
      if(iz.location < 0)
	{
	  iz.location += sitetype::getsize();       
	  if(iz.location < i.z) // Jump was very big, this can only be if there where no interfaces
	    {
	      r[i.i].spin.Smallest().right = -r[i.i].spin.Smallest().right;
	      //* debug */ cout << " (left jump to big) handled and pushed on stack: " << ca << "(with "<< *i << ")" << endl;
	      //* debug */ cout << *this;
	      return ca; // jumped over whole iterface, no need to do the right interface or insert anything
	    }
	}
      iz.right     = -spinofcluster;
      ca.smaller =  r[i.i].spin.Insert_ReturnKey(iz).gety();
      // draw new random number:
      jump = m_draw->edraw();
    }
  else // remove the interface
    {
      iz.location = ca.smaller;
      r[i.i].spin.Delete(iz);	
      // no need to draw a new random number, the resting part can be reused:
      jump -= distl;
    }

  // and right interface analogously
  if(jump < distr || ca.bigger == r[i.i].spin.Biggest().gety()) // add an interface
    {
      iz.location = i.z + jump;   // calculate new z:
      if(iz.location > sitetype::getsize())
	{ 
	  iz.location -= sitetype::getsize();       	  
	  //r[i.i].spin.left.right = -r[i.i].spin.left.right;
	  if(iz.location > i.z) 
	    {
	      r[i.i].spin.Smallest().right = -r[i.i].spin.Smallest().right;
	      //* debug */ cout << " (right jump to big) handled and pushed on stack: " << ca << "(with "<< *i << ")" << endl;
	      //* debug */ cout << *this;
	      return ca; // no need to insert interface.
	    }
	}     
      iz.right     = spinofcluster;
      ca.bigger = r[i.i].spin.Insert_ReturnKey(iz).gety();
      jump = m_draw->edraw();
    }
  else // remove the interface
    {          
      iz.location = ca.bigger;
      r[i.i].spin.Delete(iz);
      // no need to draw a new random number, the resting part can be reused:
      jump -= distr;
    }  
  // if the new area goes over the boundary, then the spin on 0 must be flipped:
  if(ca.bigger  < ca.smaller)  r[i.i].spin.Smallest().right = -r[i.i].spin.Smallest().right;

  //* debug */ cout << "handled and pushed on stack: " << ca << "(with "<< *i << ")" << endl;
  //* debug */ cout << *this;

  return ca;   
}


template <class floatt, class sitetype> 
void isingmodel<floatt, sitetype>::TransverseWolffIteration()
{
  // We need a stack to put areas in.
  // Astack is a fast stack type, see astack.h  
  static Astack<area<floatt, sitetype>, 500 > thestack; // static to prevent reallocation of pretty large array.

  // draw random position in the lattice:
  // this can may be done more efficiently (with e.g. only one draw)
  sites::c_location<floatt> i(m_draw,  sitetype::max(), sitetype::getsize());   // constructor is random 

  // remember the spin of the starting point.
  typename sitetype::spintype::attype spinofcluster = r[i.i].spin.at(i.z);

  // treat the first region:
  //  - change it size
  //  - flip the spin
  //  - push it to the stack
 
  //*debug*/ cout << "random point : " << *i << endl;
  
  thestack.push(JumpsTransverseDirection(i, spinofcluster));

  //*debug*/cout << "found spin: " << spinofcluster << endl;
  
  // when reading from the stack, we get areas, they will be stored in 'ca' (current area)
  area<floatt, sitetype> ca;
  // jumps in the transverse direction are summed in 'step'
  floatt step; // it's static because that's little quicker
  floatt width; 

  //*debug*/unsigned sidejumps = 0;
   
  //step += m_t2p * m_draw->edraw();	  	  
  
  while(!thestack.empty())
    { 
      // get an area from the stack, and look at it's neighbors:
      ca = thestack.top(); thestack.pop();      
      //*debug*/cout << "size of stack: " << thestack.size() << endl;
      //*debug*/jump =  m_t2p * m_draw->edraw(); 
      //*debug*/step = jump;
      step = m_t2p * m_draw->edraw(); 
      
      //*debug*/ cout << "found ca: " << ca << endl;
      
      width = ca.bigger - ca.smaller;
      if(width < 0) width += sitetype::getsize();

      //*debug*/ cout << "breedte: " << width << endl;

      for(int j = 0 ; j < sitetype::maxnb(); j++)
	{	       	  	 	 
	  i.i = r[ca.i].nb[j]->index();
	  i.z = ca.smaller + step;

	  if(i.z > sitetype::getsize()) i.z -= sitetype::getsize(); 

	  //*debug*/cout << "nb: " << j <<  "  new step : " << step <<  "  jump : " << jump << "   width: " << width << "  i.z: " << i.z << endl; 

	  while(step < width)
	    {
	      if(r[i.i].spin.at(i.z) == spinofcluster) 
		{	
		  //*debug*/ cout << "ok" << endl;
		  thestack.push(JumpsTransverseDirection(i, spinofcluster));
		  //*debug*/ sidejumps++;
		}
	      	      
	      step += m_t2p * m_draw->edraw();
	      //*debug*/jump  = m_t2p * m_draw->edraw(); 
	      //*debug*/step += jump;

	      i.z = ca.smaller + step;
	      if(i.z > sitetype::getsize()) i.z -= sitetype::getsize();

	      //*debug*/cout << "nb: " << j <<  "  new step : " << step <<  "  jump : " << jump << "   width: " << width << "  i.z: " << i.z << endl; 
	    }
	  step -= width;
	}// loop over neighbors      
    }// while(!thestack.empty()

  //Wolff iteration completed  
  //*debug*/ cout << "Wolfiteration " << m_nit << "  sidejumps: " << sidejumps << endl;
   m_nit++; // count iterations..
  return;
}//Anisotropic Wolff
