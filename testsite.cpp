#include "sites.h"
#include <iostream>


main()
{

  sites::simple3D::setsize(10,20,30);
  
  sites::simple3D test (1899);

  cout << test() << "  " << test.x[0] << "   " << test.x[1] << "   " << test.x[2] << endl;

  /*
  for(int i = 0; i< test.maxnb(); i++)
    {
      cout << (test.nb(i))() << "  " << test.nb(i).x[0] << "   " << test.nb(i).x[1] << "   " << test.nb(i).x[2] << endl;
    }
  */
}
