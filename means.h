
#include "ising.h"
#include "blitz/array.h"

namespace means
{
  template <class floatt, class sitetype> 
    void init(isingmodel<floatt, sitetype>& model, int n)
    {
      for(int i = 0; i < n; i++)  
	{
	  //*debug*/ cout << "iteration : " << i << endl << *model;
	  model.Iteration();	
	}
    }


  template <class floatt, class sitetype>
    void calc(isingmodel<floatt, sitetype>& model, unsigned averageQ, unsigned averagei, ostream& os = cout)
    {
      
      blitz::Array<floatt, 1> sn(model.m_ncalc);
      
      floatt sumM2 = 0;  
      floatt sumM4 = 0;     
      floatt tsumM2 = 0; // totalsum
      floatt tsumM4 = 0;
      floatt M, M2, Q, sumQ = 0 , sumQ2 =0;
      unsigned nup;
      
      
      for(int k = 0; k < averageQ; k++)
	{      
	  sumM2 = 0; 
	  sumM4 = 0;
	  for(int j = 0; j < averagei; j++)
	    {
	      // loop to reach new equilibrium		 
	      //		 for(int i = 0; i < cmdl.iterations; i++) model.WolffIteration();	  
	      model.Iteration();	
	      /*
		// function 'core' ( or 'addpqr') of ani
		bsn = nit * wn;
		sps = sps * wn;
		st  = amt * wn;  //! magnetisation per `site': m
		 s2  = st  * st;//  ! m^2                                                  

		 sn(1) = s2;//            ! m^2
		 sn(2) = s2 * s2;//       ! m^4
		 sn(3) = bsn;//           ! ?
		 sn(4) = bsn * bsn;//     ! ?^2
		 sn(5) = sn(3) * sn(4);// ! ?^3
		 sn(6) = sn(4) * sn(4);// ! ?^4
		 sn(7) = sn(1) * sn(3);
		 sn(8) = sn(2) * sn(3);
		 sn(9) = sps;

		 */
		 // Magnetisation:
		 M = (floatt)2 * model.nup() - model.size();
		 M2 = M * M;
		 sumM2 += M2;
		 sumM4 += M2*M2;          
	       }
	     
	     tsumM2 += sumM2;
	     tsumM4 += sumM4;
	     Q = (floatt)sumM2*sumM2 / (averagei * sumM4); 
	     sumQ  += Q;
	     sumQ2 += Q*Q;
	   }
	 
	 //cout << "   <M^2> = " << (floatt)sumM2/cmdl.averagei 
	 //     << "\n <M^4> = " << (floatt)sumM4/cmdl.averagei
	 Q  = sumQ / averageQ; // help to calculte dQ
	 floatt dQ = sqrt(sumQ2 / (averageQ) - Q*Q); // SDn-1 
	 floatt realQ = (floatt)tsumM2*tsumM2 / (averagei * averageQ *tsumM4); 
	 os  << "    Q  = " << realQ << " +- " << dQ << endl;
	
    }
};
