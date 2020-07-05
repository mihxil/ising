/*
   february 1999 Michiel Meeuwissen


 */

#include "ising.h"
#include "cmdline.h"
#include "means.h"

#define COMPILE_TEMPLATES
#include "ising.cpp" // strange thing to do, but otherwise the template instantiations don't work.
#include "spins.cpp"

/* -------------------------------------------------------------------	\
   MAIN									\

*/

int main(int argc, char *argv[])
{
  
  c_cmdln cmdl(argc, argv, // reads the commandline parameters in 'cmdl'.
	       "Ising programma -- Michiel Meeuwissen, march 1999"); 
  
  // We only use one random-generator.
  // The several classes all have a pointer to it. 
  // further, we make sure that the seed is always postive...
  // Take care about that, because it depends from the generator (should make something better for that)
  rnd::exponential draw(abs(cmdl.seed) + 1 );
  
    // copy a pointer to this random number generator to the classes
  spins::spin::draw = &draw;
  sites::sitet::draw = &draw;
  
  // How to write floating point numbers:
  cout.setf(ios::fixed, ios::floatfield); // default precision 6 digits
  // Stroustrup sais ios -> ios_base, but this isn't recognized by the current state of this compiler.
  
   
  try
    {
      typedef double Float;
      typedef sites::transverse2D site;
      
      site::setmultfac(cmdl.multfac); // only has meaning in transversemodel.
      
      // set the size before we define our model (otherwise runtime error, because the model declares memory)
      site::setsizei(0, cmdl.sizex);
      if(site::dimension() > 1) site::setsizei(1, cmdl.sizey);	
      site::setsize(cmdl.sizea);
      
      // define model               algoritm random  K       t2p
      isingmodel<Float, site>  model(twolff, &draw, cmdl.K, cmdl.K / 2);
      
      // init model
      model.init(1);
      
      // show state of model before 
      model.summarize(cout);
      cout << "size: " << model.size().getvalue() << endl; 
      
      // init calculations:
      means::init(model, cmdl.ntoss);
      
      // show it after
      model.summarize(cout);
      // cout << model;
      
      /*
      // start the calculation:      
      for(int i = 0; i<10; i++)
	{	  
	  means::calc(model,cmdl.averageQ, cmdl.averagei, cout);
	  
	  model.summarize(cout);
	}
      */
      cout << "number of iterations done: " << model.getnit() << endl;


     
    }
  catch (bad_alloc e)
    {
      cerr << "ising: Not enough memory available\n" << e.what() << endl;
    }
  catch (ising_error e)
    {
      cerr << "ising error: " << e.what() << endl;
    }
  catch (exception e)
    {
      cerr << "ising: unanticipated exception:\n" << e.what() << endl;
    }
     
  return 0;
}
