#ifndef __ISING_ERROR_H__
#define __ISING_ERROR_H__


#include <g++-2/stdexcept>
#include <string>

  class ising_error : public exception
    { 
      string _what;
    public:
      ising_error(const string& what_arg) :_what(what_arg) {}  
      virtual const char* what() const { return _what.c_str();}
    };

#endif //  __ISING_ERROR_H__
