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

// FINDmatrix.cc

#include "FINDmatrix.h"
#include <iostream>
#include <cstdlib> // for exit


FINDmatrix::FINDmatrix(Sample* _S)
: offx(0), offy(0), S(_S)
{
  Lx = S->get_Lx();
  Ly = S->get_Ly();
  initialize();
}

/*
 * Copy constructor:
 * Only the matrix is copied, not the descendents (A and B are NULL here)
 */
FINDmatrix::FINDmatrix(FINDmatrix& other)
: Lx(other.Lx), Ly(other.Ly), offx(other.offx), offy(other.offy),
  mtx_L(other.mtx_L), S(other.S), A(NULL), B(NULL), prefactor(other.prefactor)
{
  copy_matrix(other.mat,&mat,mtx_L);
}

/*
 * FINDmatrix constructor with fixed boundaries
 */
FINDmatrix::FINDmatrix(int _Lx, int _Ly, int _offx, int _offy, Sample* _S)
: Lx(_Lx), Ly(_Ly), offx(_offx), offy(_offy), S(_S)
{
  initialize();
}

void FINDmatrix::initialize()
{
  if (Lx == 1 && Ly == 1)              // Base case: build a Kasteleyn city,
  {                                    // ..  K matrix    ->  Pfaffian storage
    A = NULL;                          // ..  0  1  1  1      1  1  1
    B = NULL;                          // .. -1  0  1  1  ->  1  1
    mtx_L = 4;                         // .. -1 -1  0  1      1
				       // .. -1 -1 -1  0
    allocate_matrix(&mat,mtx_L);

    mat[0][0]=1;                       // 0->1, N to E
    mat[0][1]=1;                       // 0->2, N to S
    mat[0][2]=1;                       // 0->3, N to W
    mat[1][0]=1;                       // 1->2, E to S
    mat[1][1]=1;                       // 1->3, E to W
    mat[2][0]=1;                       // 2->3, S to W
    prefactor = 1;
  }
  else if (Lx > Ly)                    // Recursion with a vertical separator
  {                                    // A=left sublattice, B=right sublattice
    A = new FINDmatrix(Lx/2,Ly,offx,offy,S);
    B = new FINDmatrix(Lx-Lx/2,Ly,offx+Lx/2,offy,S);
    prefactor = combine_vertical();
    delete A; A = NULL;
    delete B; B = NULL;
  }
  else                                 // Recursion with a horizontal separator
  {                                    // A=top sublattice, B=bottom sublattice
    A = new FINDmatrix(Lx,Ly/2,offx,offy,S);
    B = new FINDmatrix(Lx,Ly-Ly/2,offx,offy+Ly/2,S);
    prefactor = combine_horizontal();
    delete A; A = NULL;
    delete B; B = NULL;
  }
}

FINDmatrix::FINDmatrix(int _mtx_L, dataType** input_matrix)
{
  A = NULL;
  B = NULL;

  mtx_L = _mtx_L;
  mat = new dataType*[mtx_L-1];
  for (int i=0; i<mtx_L; i++)
    mat[i] = new dataType[mtx_L-1-i];
  for (int i=0; i<mtx_L; i++)
  {
    if (input_matrix[i][i] != 0)
    {
      std::cerr << "non-skew-symmetric matrix provided, aborting\n";
      exit(1);
    }
    for (int j=i+1; j<mtx_L; j++)
    {
      if (input_matrix[i][j] != -input_matrix[j][i])
      {
	std::cerr << "non-skew-symmetric matrix provided, aborting\n";
	exit(1);
      }
      mat[i][j-i-1] = input_matrix[i][j];
    }
  }

  prefactor = 1;
}

FINDmatrix::~FINDmatrix()
{
  if (A != NULL)
    delete A;
  A = NULL;
  if (B != NULL)
    delete B;
  B = NULL;
  if (mat != NULL)
    delete_matrix(&mat,mtx_L);
  mat = NULL;
}

void FINDmatrix::allocate_matrix(dataType*** mtx, int L)
{
  (*mtx) = new dataType*[L-1];
  for (int i=0; i<L-1; i++)
  {
    (*mtx)[i] = new dataType[L-1-i];
    for (int j=0; j<L-1-i; j++)
      (*mtx)[i][j] = 0;
  }
}

void FINDmatrix::copy_matrix(dataType** mtx1, dataType*** mtx2, int L)
{
  (*mtx2) = new dataType*[L-1];
  for (int i=0; i<L-1; i++)
  {
    (*mtx2)[i] = new dataType[L-1-i];
    for (int j=0; j<L-1-i; j++)
      (*mtx2)[i][j] = mtx1[i][j];
  }
}

void FINDmatrix::delete_matrix(dataType*** mtx, int L)
{
  for (int i=0; i<L-1; i++)
    delete[] (*mtx)[i];
  delete[] (*mtx);
}

void FINDmatrix::allocate_transpose_matrix(dataType*** mtx, int L)
{
  (*mtx) = new dataType*[L];
  for (int i=0; i<L; i++)
  {
    (*mtx)[i] = new dataType[1+i];
    for (int j=0; j<1+i; j++)
      (*mtx)[i][j] = (i==j);
  }
}

dataType FINDmatrix::Z()
{
  return prefactor * Pf_eliminate(mtx_L/2);
}

dataType FINDmatrix::Z(int vsep, int hsep)
{
  for (int i=0; i<Lx; i++)
    mat[i][2*Lx+Ly-2*i-2] += hsep*S->get_p_bond(offx+i,offy,N);
  for (int i=0; i<Ly; i++)
    mat[Lx+i][Lx+2*Ly-2*i-2] -= vsep*S->get_p_bond(offx,offy+i,W);
  //output();
  return prefactor * Pf_eliminate(mtx_L/2);
}

// presume wrapHorz already done
dataType FINDmatrix::Zvert(int vsep)
{
  for (int i=0; i<Ly; i++)
    mat[i][2*Ly-2*i-2] -= vsep*S->get_p_bond(offx,offy+i,W);
  return prefactor * Pf_eliminate(Ly);
}

dataType FINDmatrix::wrapHorz(int hsep)
{
  // add weights that wrap around the row at the bottom
  for (int i=0; i<Lx; i++)
    mat[i][2*Lx+Ly-2*i-2] += hsep*S->get_p_bond(offx+i,offy,N);
  // reorder self so that the horizontal bonds, those that
  // cross the vertical axis, are up front, to be eliminated.
  // (4 groups of bonds: bottom, right, top, left, traversed ccw.)
  // (reverse right group, top group, then right through top)
  int xchgfactor = 1;
  for (int i = 0; i < Ly/2; ++i) {
    swaprows(Lx+i, Lx+Ly-1-i);
    //std::cout << Lx+i << " <-> " << Lx+Ly-1-i-(Lx+i+1) << "\n";
    xchgfactor = -xchgfactor;
  }
  for (int i = 0; i < Lx/2; ++i) {
    swaprows(Lx+Ly+i, Lx+Ly+Lx-1-i);
    //std::cout << Lx+Ly+i << " <-> " << Lx+Ly+Lx-1-i-(Lx+Ly+i+1) << "\n";
    xchgfactor = -xchgfactor;
  }
  for (int i = 0; i < (Lx+Ly)/2; ++i) {
    swaprows(Lx+i, Lx+Ly+Lx-1-i);
    //std::cout << Lx+i << " <-> " << Lx+Ly+Lx-1-i-(Lx+i+1) << "\n";
    xchgfactor = -xchgfactor;
  }
  prefactor *= Pf_eliminate(Lx) * xchgfactor;
  return prefactor;
}

void FINDmatrix::output()
{
  std::cout << Lx << " " << Ly << "\n";
//  std::cout << Lx << " " << Lx << " " << offx << " " << offy << "\n";
  for (int i=0; i<mtx_L; i++)
  {
    for (int j=0; j<i; j++)
//      std::cout << std::setw(4) << -mat[j][i-j-1] << " ";
//setw() requires iomanip to be included
      std::cout << -mat[j][i-j-1] << " ";
    std::cout << "0";
    for (int j=0; j<mtx_L-1-i; j++)
//      std::cout << " " << std::setw(4) << mat[i][j];
      std::cout << " " << mat[i][j];
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout.flush();
}


dataType FINDmatrix::combine_vertical()
{
  mtx_L = A->mtx_L + B->mtx_L;
  allocate_matrix(&mat,mtx_L);

  int* Aordering = new int[A->mtx_L];
  int* Bordering = new int[B->mtx_L];
  int counter = 0;
  for (int i=0; i<Ly; i++)             // interleaving part
  {
    Bordering[2*B->Lx+2*Ly-1-i] = counter++;
    Aordering[A->Lx+i] = counter++;
    mat[counter-2][0] = -S->get_p_bond(B->offx,offy+i,W);
  }
  for (int i=0; i<A->Lx; i++)
    Aordering[i] = counter++;
  for (int i=0; i<2*B->Lx+Ly; i++)
    Bordering[i] = counter++;
  for (int i=0; i<A->Lx+Ly; i++)
    Aordering[A->Lx+Ly+i] = counter++;
  
  fill_mat(A,Aordering);
  fill_mat(B,Bordering);

  delete[] Aordering;
  delete[] Bordering;

  return A->prefactor * B->prefactor * Pf_eliminate(Ly);
}


dataType FINDmatrix::combine_horizontal()
{
  mtx_L = A->mtx_L + B->mtx_L;
  allocate_matrix(&mat,mtx_L);

  int* Aordering = new int[A->mtx_L];
  int* Bordering = new int[B->mtx_L];
  int counter = 0;
  for (int i=0; i<Lx; i++)             // interleaving part
  {
    Aordering[Lx+A->Ly+i] = counter++;
    Bordering[Lx-1-i] = counter++;
    mat[counter-2][0] = S->get_p_bond(offx+Lx-1-i,B->offy,N);
  }
  for (int i=0; i<Lx+A->Ly; i++)
    Aordering[i] = counter++;
  for (int i=0; i<Lx+2*B->Ly; i++)
    Bordering[Lx+i] = counter++;
  for (int i=0; i<A->Ly; i++)
    Aordering[2*Lx+A->Ly+i] = counter++;

  fill_mat(A,Aordering);
  fill_mat(B,Bordering);

  delete[] Aordering;
  delete[] Bordering;

  return A->prefactor * B->prefactor * Pf_eliminate(Lx);
}

void FINDmatrix::fill_mat(FINDmatrix* from, int* ordering)
{
  for (int i=0; i<from->mtx_L; i++)
  {
    int newi = ordering[i];
    for (int j=0; j<from->mtx_L-1-i; j++)
    {
      int newj = ordering[j+1+i];
      if (newi > newj)
        mat[newj][newi-newj-1] = -from->mat[i][j];
      else
        mat[newi][newj-newi-1] = from->mat[i][j];
    }
  }
}


// use a semi-pivot: allow pivoting only within the first 2*numEvenRows so the
// rest of the matrix is unchanged by this
dataType FINDmatrix::Pf_eliminate(int numEvenRows)
{
  int pivotfactor = 1;
  for (int i = 0; i < numEvenRows*2; i += 2)
  {
    dataType maxMag = 0;
    int pivotrow = 0;
//  for (int j = 0; j < numEvenRows*2-i; j += 2)
    for (int j = 0; j < numEvenRows*2-i-1; j++)
    {
      if (mat[i][j] > maxMag)
      {
	pivotrow = j;
	maxMag = mat[i][j];
      }
      if (mat[i][j] < -maxMag)
      {
	pivotrow = j;
	maxMag = -mat[i][j];
      }
    }
    if (pivotrow != 0)
    {
      pivotfactor = -pivotfactor;
      pivotrows(i, pivotrow);
    }
    if (mat[i][0] == 0)
    {
      std::cerr << "zero superdiag error\n";
      exit(1);
    }
    for (int j = 1; j < mtx_L - i - 1; j++)
    {				       // do the cross operation
      if (mat[i][j] != 0)
	crossOp(i, j);
    }
  }

  dataType superDiagProd = 1.;
  for (int i = 0; i < numEvenRows * 2; i += 2)
    superDiagProd *= mat[i][0];

  if (2*numEvenRows < mtx_L)
  {
    for (int i=0; i<2*numEvenRows; i++)
      delete[] mat[i];
    mtx_L -= 2*numEvenRows;
    dataType** newmat = new dataType*[mtx_L-1];
    for (int i=0; i<mtx_L-1; i++)
      newmat[i] = mat[i+2*numEvenRows];
    delete[] mat;
    mat = newmat;
  }
  return pivotfactor * superDiagProd;
}

// i and j are both row indices, not offsets
// does not assume zeros in previous rows, see pivotrows() for pivoting op
void FINDmatrix::swaprows(int i, int j) { // swap rows i and j; true row indices, not offset for j
        if (j < i) {int tmpr = i; i = j; j = tmpr;}  // order the two arguments so that i is smaller
        dataType tmp;
        int rowA, offA, rowB, offB, flag;
        mat[i][j-i-1] = -mat[i][j-i-1];
        for (int k = 0; k < mtx_L; ++k) {  // swap rows k
                if (k == i || k == j) continue; // skips two elements: zero (diag) and negated
                if (k < i) {
                        rowA = k;
                        offA = i - k - 1;
                        flag = 1;
                } else {
                        rowA = i;
                        offA = k - i - 1;
                        flag = -1;
                }
                if (k < j) {
                        rowB = k;
                        offB = j - k - 1;
                        //flag = flag;
                } else {
                        rowB = j;
                        offB = k - j - 1;
                        flag = -flag;
                }
                tmp = flag * mat[rowA][offA];
                mat[rowA][offA] = flag * mat[rowB][offB];
                mat[rowB][offB] = tmp;
        }
}

// i is native row - swap with row i+j (j is the offset)
// assumes zeros in previous rows, see swaprows() for full pre-elimination swap
void FINDmatrix::pivotrows(int i, int j)
{
  dataType tmp;
  // see diagram for cross op: i) swap x and y, ii) -1 with c, iii) : with -: iv) 2 with d
  // i) - swap x,y 
  tmp = mat[i][0];
  mat[i][0] = mat[i][j];
  mat[i][j] = tmp;

  // ii) - swap 1,c, with negation
  for (int k = 0; k < j-1; ++k)
  { // swap start of row i with column from [i][i+j] down
    tmp = -mat[i+1][k];
    mat[i+1][k] = -mat[i+2+k][j-k-2];
    mat[i+2+k][j-k-2] = tmp;
  }

  // iii) swap intersecting element
  mat[i+1][j-1] = -mat[i+1][j-1];  // negate intersection

  // iv) swap end of row i with end of row i+j
  for (int k = 0; j + k < mtx_L-i-2; k++)
  { 
    tmp = mat[i+1][j+k];
    mat[i+1][j+k] = mat[i+j+1][k];
    mat[i+j+1][k] = tmp;
  }
}

void FINDmatrix::crossOp(int i, int j)
{ // already tested that [i][0] != 0 and [i][j] != 0
  dataType scaleFactor = -mat[i][j]/mat[i][0]; // j >= 1
  mat[i][j] = 0;                      // zap [i][j] exactly
  for (int k = 0; k < j - 1; ++k)     // add column i+1 to col j - from -transpose:
    if (mat[i+1][k] != 0)           // items 1 to c in diagram above
      mat[i+2+k][j-k-2] -= scaleFactor * mat[i+1][k];
  for (int k = 0; j + k < mtx_L-i-2; k++) // add rows: adding 2 entries to d entries
    if (mat[i+1][j+k] != 0)
      mat[i+j+1][k] += scaleFactor * mat[i+1][j+k];
}
