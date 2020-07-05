#ifndef _ISING_H__
#define _ISING_H__

#include "lattice.h"
#include "rnd.h"
#include "sites.h"
#include "blitz/array.h"

/*

       
*/

template <class floatt, class sitetype>
class area
{
public:
  //  area(typename sitetype::interface * s = NULL, typename sitetype::interface * b = NULL, const unsigned i_ = 0)
  area(floatt s = 0, floatt b = 0, const unsigned i_ = 0)
    {
      smaller = s;
      bigger = b;
      i = i_;
    }
  area(const area<floatt, sitetype>& a) { smaller = a.smaller; bigger = a.bigger; i = a.i; }

  set(floatt s = 0, floatt b = 0, const unsigned i_ = 0)
    {
      smaller = s;
      bigger = b;
      i = i_;
    }
  //typename sitetype::interface *smaller,*  bigger;
  floatt smaller, bigger;
  unsigned i;

  template <class f, class s>
  friend ostream& operator<<(ostream&, const area<f, s>& );
};


typedef enum { metropolis, swendsenwang, wolff, twolff} iterationtype;

template <class floatt, class sitetype> 
class isingmodel : public lattice<sitetype>
{
 

  
  floatt           m_K, m_one_minus_exp;
  floatt  m_B, m_Jx, m_exp_kx;
  //  m_Kx, m_Ky, m_B, m_one_minus_exp, m_beta;
  unsigned         m_nit; // number of iterations
  rnd::exponential * m_draw; 

  floatt m_t2p; // as in ani.f ! find a better name....

  static const unsigned m_ncalc = 9;
  static const unsigned m_npmax = 2000; // this was copied from ani.f, no idea why it should be this big.

  
  blitz::Array<floatt, 1> p;  // size are defined in the constructor.
  blitz::Array<floatt, 2> q;
  blitz::Array<floatt, 2> pp; 
  
  
  
  void (isingmodel<floatt, sitetype>::*m_iteration)(void); // function used for iterations.  , neat idea, but it's going to give to much trouble.
  
 public:

  void Iteration() { (this->*m_iteration)(); }

  void MetropolisIteration();         // how to do one simple Metropolis iteration
  void SwendsenWangIteration();       // Swendsen-Wang (not implemented)
  void WolffIteration();              // Simple Wolff.
 
  void TransverseWolffIteration();    // Wolff for the Transverse Model.
  inline area<floatt, sitetype> JumpsTransverseDirection(sites::c_location<floatt>&, typename sitetype::spintype::attype&);

  // constructor:
  isingmodel(iterationtype it, rnd::exponential * d, floatt K = 1, float t2p = 1) : lattice<sitetype>(),  
            p(m_ncalc), q(m_ncalc,m_ncalc), pp(m_ncalc,m_npmax) 
    {
      m_draw = d;
      m_K = K;
      m_nit = 0;
      m_B = 0;

      m_t2p = t2p;
      m_one_minus_exp = 1 - exp(- 2*m_K);
      
      switch(it)
	{	  
	case metropolis:   m_iteration = &MetropolisIteration; break;
	case swendsenwang: m_iteration = &SwendsenWangIteration; break;
	case wolff:        m_iteration = &WolffIteration; break;
	case twolff:       m_iteration = &TransverseWolffIteration; break;
	default:           m_iteration = 0;
	}
    }
  
  
  //  void Iterate() { m_iteration(); } // for the time being I resign about this.

  void     summarize(ostream& os = cout);
  unsigned getnit() { return m_nit;}
  
};

#endif // #define _ISING_H__
