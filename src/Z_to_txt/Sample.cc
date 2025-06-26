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

// Sample.cc
//

#include "Sample.h"
#include <string_view>
#include <string>
#include <iostream>
#include "exp_log.h"

// Constructor takes a filename and temperature: reads in J_{ij}
// and computes bond weights. Input format for file is:
// First line: Lx, Ly
// Then Lx*Ly lines of x, y coords, direction of bond, J_{ij}
//Diagram of spin layout (x,y) values for spins =
//
// (   0,   0) (   1,   0) (   2,   0) ... (Lx-1,   0)
// (   0,   1) (   1,   1) (   2,   1) ... (Lx-1,   1)
// (   0,   2) (   1,   2) (   2,   2) ... (Lx-1,   2)
//      .           .           .               .
//      .           .           .               .
//      .           .           .               .
// (   0,Ly-1) (   1,Ly-1) (   2,Ly-1) ... (Lx-1,Ly-1)
//
//
// For this spin layout, directions =
//
//      N(0)
//
//  W(3)    E(1)
//
//      S(2)

Sample::Sample(std::string_view filename, dataType T)
{
  Z_prefactor = 1;
  std::ifstream infile(filename.data(), std::ifstream::in);
  infile >> Lx >> Ly;
  xbonds = new dataType*[Lx];
  ybonds = new dataType*[Lx];
  for (int i=0; i<Lx; i++)
  {
    xbonds[i] = new dataType[Ly];
    ybonds[i] = new dataType[Ly];
    for (int j=0; j<Ly; j++)
    {
      xbonds[i][j] = 0;
      ybonds[i][j] = 0;
    }
  }
  int nextx;
  int nexty;
  std::string direction;
  std::string Jchars;
  dataType J;
  while (infile >> nextx)
  {
    infile >> nexty >> direction >> Jchars;
    dataType J = ((dataType)Jchars.c_str());

    exp_log EL;
    Z_prefactor *= EL.exp(J/T);
    switch(direction[0])
    {
      case 'N':
      case '0':
	ybonds[ nextx         ][(nexty+Ly-1)%Ly] = EL.exp(-2*J/T);
	break;
      case 'E':
      case '1':
	xbonds[ nextx         ][ nexty         ] = EL.exp(-2*J/T);
	break;
      case 'S':
      case '2':
	ybonds[ nextx         ][ nexty         ] = EL.exp(-2*J/T);
	break;
      case 'W':
      case '3':
	xbonds[(nextx+Lx-1)%Lx][ nexty         ] = EL.exp(-2*J/T);
    }
  }
}

Sample::~Sample()
{
  for (int i=0; i<Lx; i++)
  {
    delete[] xbonds[i];
    delete[] ybonds[i];
  }
  delete[] xbonds;
  delete[] ybonds;
}

// spin (& plaquette) numbering runs left->right (W->E), and up->down (N->S)
// px,py are plaquette coords; dir: 0->N, 1->E, 2->S, 3->W
dataType Sample::get_p_bond(int px, int py, Dir dir)
{
  switch(dir)
  {
    case N:
      return -xbonds[px][py];
    case E:
      return ybonds[px+1][py];
    case S:
      return xbonds[px][py+1];
    case W:
      return -ybonds[px][py];
    default:
      return 0;
  }
}

int Sample::get_Lx()
{
  return Lx;
}

int Sample::get_Ly()
{
  return Ly;
}

dataType Sample::get_Z_prefactor()
{
  return Z_prefactor;
}

void Sample::printMe(dataType T) {
  std::cout << "#Sample of size " << Lx << " x " << Ly << " prefactor " << Z_prefactor << "\n";
  std::cout << Lx << " " << Ly << "\n";

  exp_log EL;
  for (int j = 0; j < Ly; ++j)
    for (int i = 0; i <Lx ; ++i) {
      dataType Jx = -EL.find_log(xbonds[i][j])*T/2;
      dataType Jy = -EL.find_log(ybonds[i][j])*T/2;
      std::cout << i << " " << j << " 1 " << Jx << "\n";
      std::cout << i << " " << j << " 2 " << Jy << "\n";
    }
}
