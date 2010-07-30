// 
// File    : Error.hh
// ------------------
// Created : Tue May 27 10:52:24 2008
// Authors : Rollin C. Thomas (RCT) - rcthomas@lbl.gov
// Purpose : Tired of not having errors show up right...
//  
// $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/Grid/Transfer/Error.hh,v 1.1 2010-07-30 16:45:52 xavier Exp $ 
// 

#ifndef __ERROR__
#define __ERROR__

#include <stdexcept>

class Error : public std::runtime_error
{

   public :

      Error( const std::string& message ) : std::runtime_error( message ) {}

};

#endif
