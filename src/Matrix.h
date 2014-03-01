// -*- mode: c++; fill-column: 70; -*-

//////////////////////////////////////////////////////////////////////

// Copyright 2012 Jesse Windle - jwindle@ices.utexas.edu

// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////////////////////

// This class "implements" the MatrixFrame class.  You can think the
// Matrix class as a data container and the MatrixFrame class as a
// view of that container.  See MatrixFrame.h for documentation of
// MatrixFrame class.

// Everything here is COLUMN MAJOR for compatibility with Fortran.

// When compiling include -lblas -llapack.

//////////////////////////////////////////////////////////////////////

#ifndef __MATRIX__
#define __MATRIX__

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm> // For min,max.
#include <stdio.h>

#ifndef DISABLE_FIO
#include <fstream>
using std::ofstream;
using std::ifstream;
#endif

#include "MatrixFrame.h"

using std::vector;
using std::string;
using std::ostream;
using std::istream;
using std::stringstream;
using std::min;
using std::max;

//////////////////////////////////////////////////////////////////////
			     // BLOCK //
//////////////////////////////////////////////////////////////////////

// Inherentence with templated classes is slightly more tricky.
// http://stackoverflow.com/questions/1239908/why-doesnt-a-derived-template-class-have-access-to-a-base-template-class-ident

template<typename SCLR>
class Block : public Frame<SCLR>
{
 protected:
  std::vector<SCLR> v;   // vector to store data.

  using Frame<SCLR>::p;
  using Frame<SCLR>::nr;
  using Frame<SCLR>::nc;
  using Frame<SCLR>::nm;

 public:

  // Constructors.
  Block() : Frame<SCLR>(), v(1)
    { p = &v[0]; nr = 1; nc = 1; nm = 1; }
  Block(uint r, uint c=1, uint n=1) : Frame<SCLR>(), v(r*c*n)
    { p = &v[0]; nr = r; nc = c; nm = n; }
  Block( int r,  int c=1,  int n=1) : Frame<SCLR>(), v(r*c*n)
    { p = &v[0]; nr = (uint)r; nc = (uint)c; nm = (int)n; }
  Block(uint r, uint c, uint n, double f): Frame<SCLR>(), v(r*c*n, f)
    { p = &v[0]; nr = r; nc = c; nm = n; }
  Block(SCLR d, const Block<SCLR>& M) : Frame<SCLR>(), v(M.nr*M.nc*M.nm, d)
    { p = &v[0]; nr = M.nr; nc = M.nc; nm = M.nm; }
  Block(const Block<SCLR>& M) : Frame<SCLR>(), v(M.v)
    { p = &v[0]; nr = M.nr; nc = M.nc; nm = M.nm; }
  Block(const Frame<SCLR>& M, uint n=1) : Frame<SCLR>(), v(M.vol()*n)
    { p = &v[0]; nr = M.rows(); nc = M.cols(); nm = M.mats()*n;
      for(uint i = 0; i < nr*nc*nm; i++) v[i] = M(i % (nr*nc*M.mats()) ); }
  Block(SCLR d) : Frame<SCLR>(), v(1, d)
  { p = &v[0]; nr = 1; nc = 1; nm = 1; }

  Block(const double *ptr, uint r, uint c, uint n=1) : Frame<SCLR>(), v(r*c*n)
  { p = &v[0]; nr = r; nc = c; nm = n;
    for(uint i = 0; i < nr*nc*nm; i++) v[i] = ptr[i]; }
  Block(const double *ptr,  int r,  int c,  int n=1) : Frame<SCLR>(), v(r*c*n)
  { p = &v[0]; nr = (uint)r; nc = (uint)c; nm = (uint)n;
    for(uint i = 0; i < nr*nc*nm; i++) v[i] = ptr[i]; }

  // For predefined types of matrices.
  Block(const string& s, uint r, uint n=1);

  // ~Block() {};

  // Test equality and assign equality.
  Block<SCLR>& operator= (const Block<SCLR> &M);      // Makes copy.
  Block<SCLR>& operator= (const Frame<SCLR> &M); // Makes copy.
  // bool    operator==(const Block &M) const;

  // Iterators...
  typename vector<SCLR>::iterator begin()
  { return v.begin(); }
  typename vector<SCLR>::iterator end()   
  { return v.end(); }
  typename vector<SCLR>::const_iterator begin() const 
  { return v.begin(); }
  typename vector<SCLR>::const_iterator end() const 
  { return v.end(); }

  // Utility functions.
  void resize(uint r, uint c=1, uint n=1)
  { v.resize(r*c*n); p = &v[0]; nr = r; nc = c; nm = n; }
  void clone(const Frame<SCLR>& M);
  template<typename IDX> void clone(const Frame<SCLR>& M, const Frame<IDX>& rs, const Frame<IDX>& cs);
  template<typename IDX> void clone(const Frame<SCLR>& M, const Frame<IDX>& rs, uint c);
  template<typename IDX> void clone(const Frame<SCLR>& M, uint r, const Frame<IDX>& cs);
  // void copy(const Block& M);
  //void cbind(const Frame<SCLR>& M);
  //void rbind(const Frame<SCLR>& M);

  // Read //
  uint read(      istream&  is, bool header=0, bool binary=0);
  uint readstring(const string& s, bool header=0);
  #ifndef DISABLE_FIO
  uint read(const string& file, bool header=0, bool binary=0);
  #endif

  // Writing is taken care of in Frame<SCLR>.h.
  // Block operations are taken care of in Frame<SCLR>.h

  
}; // Block

//////////////////////////////////////////////////////////////////////
			    // TYPEDEF //
//////////////////////////////////////////////////////////////////////

#ifndef Matrix
typedef Block<double> Matrix;
#endif

#ifndef Mat
typedef Block<double> Mat;
#endif

#ifndef Matrices
typedef Block<double> Matrices;
#endif

//////////////////////////////////////////////////////////////////////
			  // DECLARATIONS //
//////////////////////////////////////////////////////////////////////

template<typename SCLR> 
void mult(Block<SCLR>& c, Frame<SCLR>& a, Frame<SCLR>& b, char ta='N', char tb='N', SCLR alpha=1.0, SCLR beta=0.0);

template<typename SCLR>
int cg(Frame<SCLR>& x, Frame<SCLR>& A, Frame<SCLR>& b, SCLR tol, int max_iter);

// ax = b -> b := x;
template<typename SCLR>
int symsolve(Frame<SCLR> a, Block<SCLR>& b, char uplo='L');
template<typename SCLR>
int syminv(Frame<SCLR> a, Block<SCLR>& ainv, char uplo='L');

// Get the Cholesky decomposition of a matrix a.
template<typename SCLR>
int chol(Block<SCLR>& c, Frame<SCLR> a, char uplo='L');

template<typename SCLR>
int symeigen(Block<SCLR>& evec, Block<SCLR>& eval, Frame<SCLR> symmat);
template<typename SCLR>
int symsqrt(Block<SCLR>& rt, Frame<SCLR> symmat);
template<typename SCLR>
int syminvsqrt(Block<SCLR>& rt, Frame<SCLR> symmat);
template<typename SCLR>
int svd2(Block<SCLR>& U, Block<SCLR>& S, Block<SCLR>& tV, Block<SCLR>& X);
template<typename SCLR>
int svd(Block<SCLR>& U, Block<SCLR>& S, Block<SCLR>& tV, Block<SCLR>& X, char jobz='A', bool pad_S=false);

//--------------------------------------------------------------------

#define BLASDEC(SCLR)							\
									\
  void rsyevd(char jobz, char uplo, int n, SCLR* a, int lda, SCLR* w, SCLR* work, int lwork, int* iwork, int liwork, int* info); \
  void rgesvd(char jobu, char jobvt, int m, int n, SCLR* a, int lda, SCLR* s, SCLR* u, int ldu, SCLR* vt, int ldvt, SCLR*work, int lwork, int* info); \
  void rgesdd(char jobz, int m, int n, SCLR* a, int lda, SCLR* s, SCLR* u, int ldu, SCLR* vt, int ldvt, SCLR* work, int lwork, int* iwork, int* info); 

BLASDEC(double)
BLASDEC(float)

#undef BLASDEC

//--------------------------------------------------------------------

extern "C" {

  // DOUBLE

  void dsyevd_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
  void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
  void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO); 

  // FLOAT

  void ssyevd_(char* JOBZ, char* UPLO, int* N, float* A, int* LDA, float* W, float* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
  void sgesvd_(char* JOBU, char* JOBVT, int* M, int* N, float* A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* INFO);
  void sgesdd_(char* JOBZ, int* M, int* N, float* A, int* LDA, float* S, float* U, int* LDU, float* VT, int* LDVT, float* WORK, int* LWORK, int* IWORK, int* INFO);

}

//////////////////////////////////////////////////////////////////////
			  // Constructors //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
Block<SCLR>::Block(const string& s, uint r, uint n) : Frame<SCLR>(), v(1)
{
  switch (s[0])
    {
    case 'I': // The identity matrix.
      resize(r, r, n);
      for(uint k = 0; k < nm; k++)
	for(uint i = 0; i < nr; i++)
	  Frame<SCLR>::operator()(i,i,k) = 1;
      break;
    case '1': // The "unity" column vectors.
      resize(r, 1, n);
      for(uint k = 0; k < nm; k++)
	for(uint i = 0; i < nr; i++)
	  Frame<SCLR>::operator()(i,0,k) = 1;
      break;
    case 'N': // The Natural numbers.
      resize(r, 1, 1);
      for(uint i = 0; i < nr; i++)
	Frame<SCLR>::operator()(i,0,0) = (SCLR)(i+1);
      break;
    case 'W': // The Whole numbers.
      resize(r, 1, 1);
      for(uint i = 0; i < nr; i++)
	Frame<SCLR>::operator()(i,0,0) = (SCLR)i;
      break;
    // case 'Z': // A sequence of integers.
    //   int diff = (int)n - (int)r;
    //   int dist = abs(diff);
    //   resize(dist+1, 1, 1);
    //   int sgn = diff / dist;
    //   int idx = 0;
    //   while(idx <= dist){
    // 	operator()(idx++) = r;
    // 	r = r + sgn;
    //   }
    //   break;
    default: // Set to scalar zero.
      resize(1,1,1);
    }
}

// template<typename SCLR>
// Block<SCLR>::Block(const Frame<SCLR>& a, const Frame<SCLR>& b, char ta, char tb, SCLR alpha) : Frame<SCLR>(), v(1)
// {
//   fprintf(stderr, "You should not be using this constructor unless type is double.\n");
//   resize(1,1,1);
// }

//////////////////////////////////////////////////////////////////////
		       // Utility Functions //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
void Block<SCLR>::clone(const Frame<SCLR>& M)
{
  resize(M.rows(), M.cols(), M.mats());
  Frame<SCLR>::copy(M);
} // copy

template<typename SCLR> template<typename IDX> 
void Block<SCLR>::clone(const Frame<SCLR>& M, const Frame<IDX>& rs, const Frame<IDX>& cs)
{
  resize(rs.area(), cs.area(), 1);
  Frame<SCLR>::copy(M, rs, cs);
}

template<typename SCLR> template<typename IDX> 
void Block<SCLR>::clone(const Frame<SCLR>& M, const Frame<IDX>& rs, uint c)
{
  resize(rs.area(), 1);
  Frame<SCLR>::copy(M, rs, c);
}

template<typename SCLR> template<typename IDX> 
void Block<SCLR>::clone(const Frame<SCLR>& M, uint r, const Frame<IDX>& cs)
{
  resize(1, cs.area());
  Frame<SCLR>::copy(M, r, cs);
}

// void Block::cbind(const Frame<SCLR>& M)
// {
//   sizecheck(mats()==M.mats() && rows()==M.rows());
//   Block temp(*this);
//   resize(rows(), cols() + M.cols(), mats());
//   for(uint m = 0; m < mats(); ++m){
//     copy(temp[m], 0, 0);
//     col(nc, M.cols()).copy(M[m], 0, 0);
//   }
// }

// Is it a bad idea to overload a function found in Frame<SCLR>?
// According to Effective C++ it is, but this makes things mroe
// intuitive.

// void Block::copy(const Block& M)
// {
//   resize(M.rows(), M.cols(), M.mats());
//   for(uint i = 0; i < vol(); i++) v[i] = M.vec(i);
// } // copy

//////////////////////////////////////////////////////////////////////
		  // Assgiment and Test Equality //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
Block<SCLR>& Block<SCLR>::operator= (const Block<SCLR> &M)
{
  clone(M);
  return *this;
} // operator=

// May not need both of these operators.

template<typename SCLR>
Block<SCLR>& Block<SCLR>::operator= (const Frame<SCLR> &M)
{
  clone(M);
  return *this;
} // operator=

// bool Block::operator==(const Block &M) const
// {
//   if(p==&M(0) && nr==M.rows() && nc==M.cols() && nm==M.mats()) return true;
//   if(vol() != M.vol()) return false;
//   for(uint i = 0; i < vol(); i++) if(vec(i)!=M.vec(i)) return false;
//   return true;
// }

//////////////////////////////////////////////////////////////////////
			      // READ //
//////////////////////////////////////////////////////////////////////

// Read in a matrix from a stream.  The header contains the
// dimensionality information, i.e. rows, cols, mats.  If header is
// set to true these values are read from the stream and used to set
// the dimensions of this matrix.

template<typename SCLR> 
uint Block<SCLR>::read( std::istream& is, bool header, bool binary)
{
  // Tell us if something is wrong.
  if (!is || is.eof())  return 0;

  uint n = 0; // The number of items read.

  // Read binary.
  if(binary){
    if(header){
      uint r,c,m;
      is.read((char*) &r, sizeof(nr));
      is.read((char*) &c, sizeof(nc));
      is.read((char*) &m, sizeof(nm));
      resize(r, c, m);
    }
    n = Frame<SCLR>::scan(is, false, true);
  }
  // Write human.
  if(!binary){
    if(header){
      uint r,c,m;
      is >> r;
      is >> c;
      is >> m;
      resize(r, c, m);
    }
    n = Frame<SCLR>::scan(is, false, false);
  }
  return n;
} // read

#ifndef DISABLE_FIO
template<typename SCLR>
uint Block<SCLR>::read(const string& file, bool header, bool binary)
{
  std::ifstream ifs(file.c_str());
  if(!ifs){
    fprintf(stderr, "Cannot read file %s.\n", file.c_str());
    return 0;
  }
  return read(ifs, header, binary);
} // read
#endif

template<typename SCLR>
uint Block<SCLR>::readstring(const string& s, bool header)
{
  stringstream ss(s);
  return read(ss, header, false);
} // readstring

//////////////////////////////////////////////////////////////////////
		      // END OF CLASS METHODS //
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
			    // CASTING //
//////////////////////////////////////////////////////////////////////

// Cast scalar to Block.
template<typename SCLR>
Block<SCLR> cast(SCLR d)
{
  return Block<SCLR>(1, 1, 1, d);
}

//////////////////////////////////////////////////////////////////////
			       // R //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
Block<SCLR> seq(SCLR start, SCLR end, SCLR delta=1)
{
  double sign = end - start < 0 ? -1.0 : 1.0;
  double new_delta = sign * fabs((double)delta);
  double N_double = (end-start) / new_delta;
  if (N_double < 0) fprintf(stderr, "Problem in seq: N_double < 0");
  int N = floor(N_double) + 1;

  Block<SCLR> a(N);
  a(0) = start;
  for (int i=1; i<N; i++)
    a(i) = a(i-1) + (SCLR)delta;

  return a;
}

template<typename SCLR>
Block<SCLR> rowSums(Frame<SCLR> M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Block<SCLR> a(nr);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(i) += M(i,j);

  return a;
}

template<typename SCLR>
Block<SCLR> rowMeans(Frame<SCLR> M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Block<SCLR> a(nr);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(i) += M(i,j);

  for (uint i=0; i<nr; i++)
    a(i) = a(i) / nc;

  return a;
}

template<typename SCLR>
Block<SCLR> rowSums(Block<SCLR>& a, Frame<SCLR> M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  a = M.col(0);

  for (uint j=1; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(i) += M(i,j);

  return a;
}

template<typename SCLR>
Block<SCLR> colSums(Frame<SCLR> M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Block<SCLR> a(nc);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(j) += M(i,j);

  for (uint j=0; j<nc; j++)
    a(j) = a(j) / nr;
  
  return a;
}

template<typename SCLR>
Block<SCLR> colMeans(Frame<SCLR> M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Block<SCLR> a(nc);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(j) += M(i,j) / nr;

  return a;
}

template<typename SCLR>
SCLR maxAll(Frame<SCLR> M)
{
  uint N = M.size();
  SCLR mx = M(0);
  for(uint i=1; i<N; i++) {
    mx = M(i) > mx ? M(i) : mx;
  }
  return mx;
}

template<typename SCLR>
SCLR minAll(Frame<SCLR> M)
{
  uint N = M.size();
  SCLR mx = M(0);
  for(uint i=1; i<N; i++) {
    mx = M(i) < mx ? M(i) : mx;
  }
  return mx;
}

//////////////////////////////////////////////////////////////////////
			   // transpose //
//////////////////////////////////////////////////////////////////////

// A = t(B)
template<typename SCLR>
void trans(Block<SCLR>& A, Block<SCLR>& B)
{
  uint T = B.mats();
  uint R = B.rows();
  uint C = B.cols();
  if(A.rows() != C || A.cols() != R || A.mats() != T) A.resize(C, R, T);
  for(uint t=0; t < T; ++t)
    for(uint i=0; i < R; ++i)
      for(uint j=0; j < C; ++j)
	A(j,i,t) = B(i,j,t);
}

//////////////////////////////////////////////////////////////////////
		   // SIMPLE CONJUGATE GRADIENTS //
//////////////////////////////////////////////////////////////////////

// http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// by Jonathan Richard Shewchuk

// Assumes A is positive definite.

template<typename SCLR>
int cg(Frame<SCLR>& x, Frame<SCLR>& A, Frame<SCLR>& b, SCLR tol, int max_iter)
{
  uint P = b.rows();

  int iter = 0;

  Block<SCLR> r(b);
  gemm(r, A, x, 'N', 'N', (SCLR)-1.0, (SCLR)1.0);

  Block<SCLR> d(r);

  SCLR delta_new = dot(r, r);
  // See note below...
  // SCLR delta_0   = delta_new;
  SCLR delta_old = 0;

  Block<SCLR> q(P);
  SCLR alpha;
  SCLR beta;

  //while(iter < max_iter && delta_new > tol*delta_0){
  //The above condition is suggested by Shewchuk, but we want an absolute tolerance.
  while(iter < max_iter && delta_new > tol){
    gemm(q, A, d);
    alpha = delta_new / dot(d, q);
    axpy(alpha, d, x);

    // You might want to discard this step.
    if (iter % 50 == 0) {
      r.clone(b);
      gemm(r, A, x, 'N', 'N', (SCLR)-1.0, (SCLR)1.0);
    }
    else{
      axpy((SCLR)-1.0 * alpha, q, r);
    }

    delta_old = delta_new;
    delta_new = dot(r, r);

    beta = delta_new / delta_old;
    hsum(d, d, r, beta, 0.0);
    //for(uint j = 0; j < P; j++)
    //  d(j) = beta * d(j) + r(j);

    iter++;
  }

  return iter;
}

//////////////////////////////////////////////////////////////////////
		     // BLAS / LAPACK WRAPPERS //
//////////////////////////////////////////////////////////////////////

// Solve a symmetric, positive definite system of equations.
// ax = b -> b := x;
// int symsolve(MF a, Block<SCLR>& b, char uplo='L')
template<typename SCLR>
int symsolve(Frame<SCLR> a, Block<SCLR>& b, char uplo)
{
  // b.fill(0.0);
  // for(uint i = 0; i < b.cols(); ++i)
  //   b(i,i) = 1.0;

  Block<SCLR> temp(a);
  int info = posv(temp, b, uplo);

  if (info) {
    fprintf(stderr, "Problem with symsolve: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);

    #ifndef NTHROW
    throw std::runtime_error("potrf failed\n");
    #endif
  }

  return info;
}

//------------------------------------------------------------------------------
// int syminv(MF a, Matrix& ainv, char uplo='L')
template<typename SCLR>
int syminv(Frame<SCLR> a, Block<SCLR>& ainv, char uplo)
{
  ainv.resize(a.rows(), a.cols());
  ainv.fill(0.0);
  for(uint i = 0; i < ainv.cols(); ++i)
     ainv(i,i) = 1.0;

  return symsolve(a, ainv, uplo);
}

//------------------------------------------------------------------------------
// Get the Cholesky decomposition of a matrix a.
// int chol(Matrix& c, MF a, char uplo='L')
template<typename SCLR>
int chol(Block<SCLR>& c, Frame<SCLR> a, char uplo)
{
  c.clone(a);
  int info = chol(c, uplo);

  // FIX FIX Set other entries to zero.
  // Do I want to have an option for letting a person not do this?

  if (uplo=='L') {
    for (uint j = 1; j < c.cols(); ++j)
      for (uint i = 0; i < j; ++i)
	c(i,j) = 0.0;
  }
  else {
    for (uint j = 0; j < c.cols(); ++j)
      for (uint i = j+1; i < c.cols(); ++i)
	c(i,j) = 0.0;
  }

  if (info) {
    fprintf(stderr, "Problem with chol: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);
    
    #ifndef NTHROW
    throw std::runtime_error("potrf failed\n");
    #endif
  }

  return info;
}

//------------------------------------------------------------------------------
template<typename SCLR>
int symeigen(Block<SCLR>& evec, Block<SCLR>& eval, Frame<SCLR> symmat)
{
  sizecheck(symmat.rows()==symmat.cols());
  int info = 0;
  int N = symmat.rows();
  
  evec.clone(symmat);
  eval.resize(N);

  // int lwork  = 1 + 6 * N + 2*N*N;
  // int liwork = 3 + 5*N;

  int lwork = -1;
  int liwork = -1;

  std::vector<SCLR> work(1);
  std::vector<int>    iwork(1);

  rsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);

  lwork = (int) work[0];
  liwork = (int) iwork[0];

  work.resize(lwork);
  iwork.resize(liwork);

  rsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);
  
  if (info != 0) {
    fprintf(stderr, "problem in symeigen; info=%i.\n", info);
    #ifndef NTHROW
    throw std::runtime_error("symeigen failed\n");
    #endif
  }

  return info;
}

//------------------------------------------------------------------------------
template<typename SCLR>
int symsqrt(Block<SCLR>& rt, Frame<SCLR> symmat)
{
  int N = symmat.rows();

  Block<SCLR> evec;
  Block<SCLR> eval;
  int info = symeigen(evec, eval, symmat);

  Block<SCLR> d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

//------------------------------------------------------------------------------
template<typename SCLR>
int syminvsqrt(Block<SCLR>& rt, Frame<SCLR> symmat)
{
  int N = symmat.rows();

  Block<SCLR> evec;
  Block<SCLR> eval;
  int info = symeigen(evec, eval, symmat);

  Block<SCLR> d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(1.0 / eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

//--------------------------------------------------------------------
template<typename SCLR>
int svd2(Block<SCLR>& U, Block<SCLR>& S, Block<SCLR>& tV, Block<SCLR>& X)
{
  // Always calculate a thinned U and a non-thinned V.  Pad out the
  // singular values by 0.0 when needed.
  Block<SCLR> A(X);

  int m = A.rows();
  int n = A.cols();

  tV.resize(n,n);
  S.resize(n); S.fill(0.0); // The singular values.
  // Should be of size min(m,n), but I want to pad things out.
  
  if (m >= n) U.resize(m,n);
  else U.resize(m,m);
  int ldu = U.rows();

  char jobu  = 'S';
  char jobvt = 'A';

  int maxmn = m > n ? m : n;
  int minmn = m < n ? m : n;

  vector<SCLR> work(1);
  int lwork = -1;

  int info;

  // Workspace query.
  rgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);
  // printf("lwork: %g\n", work[0]);
  
  // SVD.
  lwork = (int)work[0];
  lwork = lwork > 1 ? lwork : 5 * minmn + maxmn;

  work.resize(lwork);
  rgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);

  if (info != 0) {
    fprintf(stderr, "problem in svd; info=%i.\n", info);
    #ifndef NTHROW
    throw std::runtime_error("svd failed\n");
    #endif
  }

  return info;
}

//--------------------------------------------------------------------
template<typename SCLR>
int svd(Block<SCLR>& U, Block<SCLR>& S, Block<SCLR>& tV, Block<SCLR>& X, char jobz, bool pad_S)
{
  
  // Sigma is M x N.
  // U is M x M.
  // V is N x N.
  // Routine returns V^T.

  Block<SCLR> A(X);

  int m = A.rows();
  int n = A.cols();

  int maxmn = m > n ? m : n;
  int minmn = m < n ? m : n;

  // Dim
  if (jobz=='A') {
    U.resize(m, m);
    tV.resize(n, n);
  }
  else {
    // if (jobz=='S')
    U.resize(m, minmn);
    tV.resize(minmn, n);
  }

  int lenS = pad_S ? maxmn : minmn;
  S.resize(lenS); S.fill(0.0); // The singular values.

  int lda  = A.rows();
  int ldu  = U.rows();
  int ldvt = tV.rows();

  vector<int> iwork(8 * minmn);
  vector<SCLR> work(1);
  int lwork = -1;

  int info;

  // Workspace query.
  rgesdd(jobz, m, n, &A(0), lda, &S(0), &U(0), ldu, &tV(0), ldvt, &work[0], lwork, &iwork[0], &info); 
  // printf("lwork: %g\n", work[0]);
  
  // SVD.
  lwork = (int)work[0];
  lwork = lwork > 1 ? lwork : 3 * minmn + max(maxmn, 4 * minmn * minmn + 3*minmn + maxmn) ;

  work.resize(lwork);
  rgesdd(jobz, m, n, &A(0), lda, &S(0), &U(0), ldu, &tV(0), ldvt, &work[0], lwork, &iwork[0], &info); 

  if (info != 0) {
    fprintf(stderr, "problem in svd; info=%i.\n", info);
    #ifndef NTHROW
    throw std::runtime_error("svd failed\n");
    #endif
  }

  return info;
}

//------------------------------------------------------------------------------
// this = op(a) * op(b) + alpha * this.
template<typename SCLR>
void mult(Block<SCLR>& c, Frame<SCLR>& a, Frame<SCLR>& b, char ta, char tb, SCLR alpha, SCLR beta)
{
  uint opa_rows = ta=='T' ? a.cols() : a.rows();
  uint opb_cols = tb=='T' ? b.rows() : b.cols();
  c.resize(opa_rows, opb_cols, 1);
  
  gemm(c, a, b, ta, tb, alpha, beta);
}

//////////////////////////////////////////////////////////////////////
	  // Hadamard operations by OVERLOADED OPERATORS //
//////////////////////////////////////////////////////////////////////

// These form of multiplication, addition, division, and subtraction
// will not be as fast as hprodeq, hsumeq, etc since it involves
// returning a copy of a Block.  However, it will make it easier to
// read code.  There are matrix packages out there, e.g. eigen, that
// will cleverly reduce the amount of overhead for concatenated
// operations.  But I found it to be lacking in documentation given
// its complexity.

#define MHOP(NAME, OP)				\
  template<typename SCLR>					\
  Block<SCLR> operator OP(const Frame<SCLR>& a, const Frame<SCLR>& b)	\
  {								\
    Frame<SCLR> small = a.area() < b.area() ? a : b;		\
    Frame<SCLR> big   = a.area() < b.area() ? b : a;		\
    sizecheck(hconform(big,small));		\
    Block<SCLR> c(big);						\
    NAME(c, small, 1.0);					\
    return c;							\
  }								\
  template<typename SCLR>					\
  Block<SCLR> operator OP(const Frame<SCLR>& a, SCLR b)	\
  {								\
    Block<SCLR> c(a);						\
    NAME(c, b);							\
    return c;							\
  }								\
  template<typename SCLR>					\
  Block<SCLR> operator OP(SCLR b, const Frame<SCLR>& a)	\
  {								\
    Block<SCLR> c(a);						\
    NAME(c, b);							\
    return c;							\
  }								\

//MHOP(hprodeq, *, *=, Block) MHOP(hsumeq, +, +=, Block)
//MHOP(hdiveq,  /, /=, Block) MHOP(hsubeq, -, -=, Block)

// MHOP(hprodeq, *, *=, Frame<SCLR>) MHOP(hsumeq, +, +=, Frame<SCLR>)
// MHOP(hdiveq,  /, /=, Frame<SCLR>) MHOP(hsubeq, -, -=, Frame<SCLR>)

MHOP(hprodeq, *) MHOP(hsumeq, +)
MHOP(hdiveq,  /) MHOP(hsubeq, -)

#undef MHOP

//////////////////////////////////////////////////////////////////////
			// Basic Functions //
//////////////////////////////////////////////////////////////////////

#ifndef sq
template<typename RealType>
RealType sq(RealType x){return x * x;}
#endif

#define UNARY(FUNC, FUNCEQ)					\
  template<typename RealType>					\
  Frame<RealType> FUNC(Frame<RealType> a, Frame<RealType> b)	\
  {								\
    sizecheck(a.vol()==b.vol());				\
    for(uint l = 0; l < a.vol(); l++) a(l) = FUNC(b(l));	\
    return a;							\
  }								\
  template<typename RealType>					\
  Block<RealType> FUNC(Frame<RealType> a)			\
  {								\
    Block<RealType> c(a);					\
    FUNC(c, a);							\
    return c;							\
  }								\
  template<typename RealType>					\
  Frame<RealType> FUNCEQ(Frame<RealType> a)			\
  {								\
    for(uint l = 0; l < a.vol(); l++) a(l) = FUNC(a(l));	\
    return a;							\
  }								\

UNARY(log, logq)   UNARY(exp, expq)   UNARY(sqrt, sqrtq)
UNARY(sin, sinq)   UNARY(cos, cosq)   UNARY(tan , tanq)
UNARY(asin, asinq) UNARY(acos, acosq) UNARY(atan, atanq)
UNARY(sinh, sinhq) UNARY(cosh, coshq) UNARY(tanh, tanhq)
UNARY(fabs, fabsq) UNARY(ceil, ceilq) UNARY(floor, floorq)
UNARY(sq, sqq)     UNARY(log10, log10q)

#undef UNARY

// NEED TO FIX THIS.

#define BINARY(FUNC, NAME)						\
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> c, Frame<SCLR> a, Frame<SCLR> b)		\
  {									\
    sizecheck(c.vol() == a.vol() && a.vol()==b.vol());			\
    for(uint l = 0; l < a.vol(); l++) c(l) = FUNC(a(l), b(l));		\
    return c;								\
  }									\
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> c, double a, Frame<SCLR> b)		\
  {									\
    sizecheck(c.vol()==b.vol());					\
    for(uint l = 0; l < b.vol(); l++) c(l) = FUNC(a, b(l));		\
    return c;								\
  }									\
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> c, Frame<SCLR> a, double b)		\
  {									\
    sizecheck(c.vol()==a.vol());					\
    for(uint l = 0; l < a.vol(); l++) c(l) = FUNC(a(l), b);		\
    return c;								\
  }									\
  template<typename SCLR>						\
  Block<SCLR> NAME(Frame<SCLR> a, Frame<SCLR> b)			\
  {									\
    Block<SCLR> c(a);							\
    NAME(c, a, b);							\
    return c;								\
  }									\
  template<typename SCLR>						\
  Block<SCLR> NAME(double a, Frame<SCLR> b)				\
  {									\
    Block<SCLR> c(b);							\
    NAME(c, a, b);							\
    return c;								\
  }									\
  template<typename SCLR>						\
  Block<SCLR> NAME(Frame<SCLR> a, double b)				\
  {									\
    Block<SCLR> c(a);							\
    NAME(c, a, b);							\
    return c;								\
  }									\

BINARY(min, hmin) BINARY(max, hmax) BINARY(pow, pow)

#undef BINARY


//////////////////////////////////////////////////////////////////////
			  // END OF CLASS //
//////////////////////////////////////////////////////////////////////

#endif // MATRIX

/*
  class block_iterator
  {

  private:

    int row;
    int col;
    int stride;

    SCLR *begin;
    double *end;
    int     offset;

  public:

    block_iterator(double *p, int r, int c)
      : row(r)
      , col(c)
      , stride(0)
      , begin(p)
      , end(p + row*col)
      , offset(0) {};

    virtual inline double& operator()(int i, int j){
      offset = (row + stride) * j + i;
      return *(begin + offset);
      // ptr = begin + (row + stride) * j + i;
      // return *ptr;
    }

    virtual inline double& operator[](int i){
      offset = (row + stride) * (i / col) + i % col;
      return *(begin + offset);
      // ptr = begin + (row + stride) * (i / col) + i % col;
      // return *ptr;
    }

    virtual inline double& operator++(){
      offset += offset % stride != 0 ? 1 : stride;
      return *(begin + offset);
      // ptr += (ptr - begin) % stride != 0 ? 1 : stride;
      // return *ptr;
    }

  };
 */
