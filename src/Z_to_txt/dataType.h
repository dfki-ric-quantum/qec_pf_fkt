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
//

// dataType.h

#ifndef DATATYPE_H
#define DATATYPE_H

#include <gmpxx.h>
typedef mpf_class dataType;

enum Dir {N, E, S, W};

#endif // DATATYPE_H
