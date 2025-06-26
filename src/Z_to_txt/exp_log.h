// COPYRIGHT NOTICE
//  This file is part of isingZ, a program to compute partition
//      functions and free energies for domain walls in 2D Ising models.
//  Copyright 2012 by Creighton K. Thomas (creightonthomas@gmail.com),
//                    A. Alan Middleton (aam@syr.edu)
//  The development of this software was supported in part by the National Science
//  Foundation under grant DMR-1006731.
// 
//  The program isingZ is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//  
//  See http://www.gnu.org/licenses/gpl.txt for more details.
//  

// exp_log.h
//
// Hand-crafted exp and log intended for use with gmpxx.
//
// This code does not work with very very large numbers or for more than 4096 bits - take care to verify
// these routines for your usage.
//
// In the context of the isingZ code.
// These functions are not needed for the core calculations, but are needed to handle
// setting up and calculations to summarize the output.
// The exp is needed to create weights, given couplings: wt[bond] = exp(-2/T * J[bond])
// The log is needed for finding free energies via F = -kT log(Z), as the core computation
//  finds the partition function Z.
//
// This code has limitations, which can be fixed if desired, at a slight loss of efficiency
// (which will not be noticeable compared with the core Pfaffian computations). But they 
// have minimal impact also, unless one runs at very high precision, especially as these are
// used only in set up and result reporting.
//
// One of these limitations is defined by MAXPOW, which sets the precision of the exponential.
// There is also a hardcoded precision check in log to
// catch a demand for too high of a precision calculation in the find_log() method, given
// the hardcoded constants at the start of the source file exp_log.cc.
// 

#ifndef EXP_LOG
#define EXP_LOG

#include "dataType.h"

#define MAXPOW 300 // sets precision of exponential
// precision of calculations set in global

class exp_log {
    static void halve(dataType x, dataType &xh, int &n);
    static dataType exponentiate(dataType x);
    static void resquare(dataType &xhalved, int n);
  public:
    static dataType pi;
    static dataType icutoff;
    static dataType exp(const dataType &x) {
      dataType xhalved;
      int n;
      halve(x, xhalved, n);
      xhalved = exponentiate(xhalved);
      resquare(xhalved, n);
      return xhalved;
    }
    static dataType agm(const dataType &a, const dataType &b);
    static dataType find_log(const dataType &x);
    static bool close(const dataType &a, const dataType &b, int thresh);
};

#endif // EXP_LOG
