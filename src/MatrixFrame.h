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
// See appendix at end of document for additional information.
//////////////////////////////////////////////////////////////////////

#ifndef __MATRIX_FRAME__
#define __MATRIX_FRAME__

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdio.h>

#ifndef DISABLE_FIO
#include <fstream>
using std::ofstream;
using std::ifstream;
#endif

#ifdef USE_R
#include <R_ext/Utils.h>
#include <R_ext/PrtUtil.h>
#endif

using std::vector;
using std::string;
using std::ostream;
using std::istream;
using std::stringstream;
using std::ostream;

#ifndef uint
typedef unsigned int uint;
#endif

//////////////////////////////////////////////////////////////////////
                          // MatrixFrame //
//////////////////////////////////////////////////////////////////////

inline bool idxcheck(bool b){
#ifndef NDEBUG
  if (!b) throw std::runtime_error("Index out of bounds.\n");
  return b;
#else
  return true;
#endif
}

inline bool sizecheck(bool b){
#ifndef NDEBUG
  if (!b) throw std::runtime_error("Incompatible dimension.\n");
  return b;
#else
  return true;
#endif
}

inline bool memcheck(bool b){
#ifndef NDEBUG
  if (!b) throw std::runtime_error("Memory overlap.\n");
  return b;
#else
  return true;
#endif
}

template<typename SCLR>
class Frame
{
 protected:
  // I do not want to use *const p since I want to be able to resize my matrices.
  SCLR *p;  // Pointer to the array of matrices.
  uint nr;  // Number of rows.
  uint nc;  // Number of columns.
  uint nm;  // Number of matrices.

  // Matrix Check
  bool allow(uint l) const
    { return (p!=NULL && l < (nr*nc)); }
  bool allow(int  l) const
    { return (p!=NULL && l>=0 && l < (int)(nr*nc)); }
  bool allow(uint r, uint c) const
    { return (p!=NULL && r < nr && c < nc); }
  bool allow(int  r,  int c) const
    { return (p!=NULL && r>=0 && c>=0 && r < (int)nr && c < (int)nc); }

  // Array of Matrix Check
  bool indexok(uint i) const
  { return i < nm; }
  bool indexok( int i) const
  { return (0 <= i && i < (int)nm); }

 public:
  // Constructors, etc.  Do not make a constructor which copies.
  Frame()
    { p = NULL; nr = 0; nc = 0; nm = 0; }
  Frame(SCLR *ptr, uint r=1, uint c=1, uint m=1)
    { p = ptr; nr = r; nc = c; nm = m; }
  Frame(SCLR *ptr,  int r=1,  int c=1,  int m=1)
    { p = ptr; nr = (uint)r; nc = (uint)c; nm = (uint)m; }
  Frame(const Frame<SCLR> & M)
    { p = M.p; nr = M.nr; nc = M.nc; nm = M.nm; }

  ~Frame()
    { p = NULL; }

  // Test equality and assign equality.
  Frame<SCLR>& operator= (const Frame<SCLR> &M); // Makes copy.
  bool         operator==(const Frame<SCLR> &M) const;
  bool         sameframe (const Frame<SCLR> &M) const;

  // rows(), cols(), mats, area(), vol();
  inline uint rows() const { return nr; }
  inline uint cols() const { return nc; }
  inline uint mats() const { return nm; }
  inline uint area() const { return nr * nc; }
  inline uint vol()  const { return nr * nc * nm; }
  inline uint size() const { return nr * nc * nm; }

  // Matrix Access

  // Returns the (r,c) element of matrix.
  SCLR& operator()(uint r, uint c)
    { 
    #ifndef NDEBUG
        idxcheck(allow(r, c)); 
    #endif      
        return p[c * nr + r]; 
    }
  const SCLR& operator()(uint r, uint c) const
    { 
    #ifndef NDEBUG
        idxcheck(allow(r, c)); 
    #endif
        return p[c * nr + r]; 
    }

  // Array of Matrix Access

  // Returns the lth element of array of matrices.
  // I debate whether this is confusing notation.
  const SCLR& operator()(uint l) const
    { 
    #ifndef NDEBUG
        idxcheck(l < nr*nc*nm);    
    #endif    
        return p[l]; 
    }
    
  SCLR& operator()(uint l)
    { 
    #ifndef NDEBUG
        idxcheck(l < nr*nc*nm);   
    #endif
        return p[l]; 
    }

  // Returns the (r,c) element of matrix[t].
  SCLR& operator()(uint r, uint c, uint t)
    {   
    #ifndef NDEBUG  
        idxcheck(indexok(t) && allow(r,c)); 
    #endif
        return p[t * nr*nc + c * nr + r]; 
    }
    
  const SCLR& operator()(uint r, uint c, uint t) const
    { 
    #ifndef NDEBUG
        idxcheck(indexok(t) && allow(r,c)); 
    #endif
        return p[t * nr*nc + c * nr + r]; 
    }    
    
  SCLR& get(uint r, uint c=0, uint t=0)
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nr*nc + c * nr + r]; }
  const SCLR& get(uint r, uint c=0, uint t=0) const
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nr*nc + c * nr + r]; }

  // Returns the ith element of the array p.
  SCLR& vec(uint i)
  { idxcheck(i < nr*nc*nm); return p[i]; }
  const SCLR& vec(uint i) const
  { idxcheck(i < nr*nc*nm); return p[i]; }

  // Returns a Frame pointing to the ith matrix or,
  // if there is one matrix, to the ith column.
  Frame<SCLR> operator[](uint i)
  { idxcheck(indexok(i)); return Frame<SCLR>(&p[0+i*area()], nr, nc); }

  // Get the pointer.  Be wary.
  // const double* const getp()
  // { return p; }
  inline SCLR* getp()
  { return p; }
  inline void setp(SCLR *p_)
  { p = p_; }

  // Array of Matrix Functions.

  void copy(const Frame<SCLR>& M);     // Copy values.
  // void thincopy(Frame& M);          // Copy pointer and dimensions.
  Frame<SCLR> fill(const SCLR& x);            // Fill with value.
  Frame<SCLR> col(uint c, uint num=1); // The c-th to c+num-1th col.
  Frame<SCLR> dim(uint r, uint c, uint m=1); // Return a MF with different, compatible dim.

  // Read / Write.

  bool write(      ostream&  os, bool header=0, bool binary=0);
  uint  scan(      istream&  is, bool header=0, bool binary=0);

  #ifndef DISABLE_FIO
  bool write(const string& file, bool header=0, bool binary=0);
  uint  scan(const string& file, bool header=0, bool binary=0);
  #endif
  // bool  readstring(const string& s, bool header=0);

  // Matrix Functions.

  // Fill this matrix from M starting at (r,c).
  void copy(const Frame<SCLR>& M, uint r, uint c);
  // Copy the matrix M along rows rs and columns cs.
  template<typename IDX> void copy(const Frame<SCLR>& M, const Frame<IDX>& rs, const Frame<IDX>& cs);
  template<typename IDX> void copy(const Frame<SCLR>& M, const Frame<IDX>& rs, uint c);
  template<typename IDX> void copy(const Frame<SCLR>& M, uint r, const Frame<IDX>& cs);
  void copy_transpose(const Frame<SCLR>& M);

  // Set the elements in rows rs and columns cs using the matrix M.
  template<typename IDX> void set(const Frame<IDX>& rs, const Frame<IDX>& cs, const Frame<SCLR>& M);
  template<typename IDX> void set(const Frame<IDX>& rs, uint c, const Frame<SCLR>& M);
  template<typename IDX> void set(uint r, const Frame<IDX>& cs, const Frame<SCLR>& M);

}; // MatrixFrame

//////////////////////////////////////////////////////////////////////
			    // TYPEDEF //
//////////////////////////////////////////////////////////////////////

#ifndef MF
typedef Frame<double> MF;
#endif

#ifndef MatrixFrame
typedef Frame<double> MatrixFrame;
#endif

//////////////////////////////////////////////////////////////////////
		       // WRAPPER TO FORTRAN //
//////////////////////////////////////////////////////////////////////

// y = alpha x + y.
template<typename SCLR>
void axpy(SCLR alpha, Frame<SCLR> x, Frame<SCLR> y); 

 // x'y
template<typename SCLR>
SCLR dot(Frame<SCLR> x, Frame<SCLR> y);

// c = alpha op(a) * op(b) + beta c.
template<typename SCLR>
void gemm(Frame<SCLR> c, Frame<SCLR> a, Frame<SCLR> b, char ta='N', char tb='N', SCLR alpha=1.0, SCLR beta=0.0); 

// b = alpha op(a) * b  OR  b = alpha b * op(a) where a is triangular.
template<typename SCLR>
void trmm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side='L', char ta='N', char diag='N', SCLR alpha=1.0);

// Solve x:  op(a) x = alpha b  OR  x op(a) = alpha b, a triangular.
// i.e: x = alpha inv(op(a)) b  OR  x = alpha b inv(op(a)).
// The solution is overwriten into b.
template<typename SCLR>
void trsm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side='L', char ta='N', char diag='N', SCLR alpha=1.0);

// c = alpha a' a + beta c, ta='T'
// c = alpha a a' + beta c, ta='N'
template<typename SCLR>
void syrk(Frame<SCLR> c, Frame<SCLR> a, char ta='N', SCLR alpha=1.0, SCLR beta=0.0);

// Solve a general linear system, ax = b for x.
template<typename SCLR>
int gesv(Frame<SCLR> a, Frame<SCLR> b);

// Solves ax = b for x where a is sym. pos. def.  Note: the lower (or
// upper) portion of A is overwritten with the Cholesky decomposition.
template<typename SCLR>
int posv(Frame<SCLR> a, Frame<SCLR> b, char uplo);

// Cholesky Decomposition (in place)
template<typename SCLR>
int potrf(Frame<SCLR> a, char uplo);

template<typename SCLR>
int chol(Frame<SCLR> a, char uplo='L');

//--------------------------------------------------------------------
  
// BLAS Level 1
// BLAS Level 3
// LAPACK

#define BLASDEC(TYPE) \
  void raxpy(int n, TYPE da, TYPE* dx, int incx, TYPE* dy, int incy);	\
  TYPE rdot(int n, TYPE* dx, int incx, TYPE* dy, int incy);		\
									\
  void rgemm(char transa, char transb, int m, int n, int k, TYPE alpha, TYPE* a, int lda, TYPE* b, int ldb, TYPE beta, TYPE* c, int ldc); \
  void rsyr2k(char uplo, char trans, int n, int k, TYPE alpha, TYPE* a, int lda, TYPE* b, int ldb, TYPE beta, TYPE* c, int ldc); \
  void rtrmm(char side, char uplo, char transa, char diag, int m, int n, TYPE alpha, TYPE* a, int lda, TYPE* b, int ldb); \
  void rtrsm(char side, char uplo, char transa, char diag, int m, int n, TYPE alpha, TYPE* a, int lda, TYPE* b, int ldb); \
  void rsyrk(char uplo, char trans, int n, int k, TYPE alpha, TYPE* a, int lda, TYPE beta, TYPE* c, int ldc); \
									\
  void rgesv(int n, int nrhs, TYPE* a, int lda, int* ipiv, TYPE* b, int ldb, int& info); \
  void rposv(char uplo, int n, int nrhs, TYPE* a, int lda, TYPE* b, int ldb, int& info); \
  void rpotrf(char uplo, int n, TYPE* a, int lda, int& info);

BLASDEC(double)
BLASDEC(float)

#undef BLASDEC

//--------------------------------------------------------------------
extern "C" {

  // BLAS Level 1
  // BLAS Level 3
  // LAPACK

  // double

  void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);
  double ddot_(int* N, double* DX, int* INCX, double* DY, int* INCY);

  void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
  void dsyr2k_(char* UPLO, char* TRANS, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

  void dtrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);
  void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);
  void dsyrk_(char* UPLO, char* TRANS, int* N, int* K, double* ALPHA, double* A, int* LDA, double* BETA, double* C, int* LDC);

  void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
  void dposv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, int* INFO);
  void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);

  // float
 
  void saxpy_(int* N, float* DA, float* DX, int* INCX, float* DY, int* INCY);
  float sdot_(int* N, float* DX, int* INCX, float* DY, int* INCY);

  void sgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, float* ALPHA, float* A, int* LDA, float* B, int* LDB, float* BETA, float* C, int* LDC);
  void ssyr2k_(char* UPLO, char* TRANS, int* N, int* K, float* ALPHA, float* A, int* LDA, float* B, int* LDB, float* BETA, float* C, int* LDC);

  void strmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* ALPHA, float* A, int* LDA, float* B, int* LDB);
  void strsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* ALPHA, float* A, int* LDA, float* B, int* LDB);
  void ssyrk_(char* UPLO, char* TRANS, int* N, int* K, float* ALPHA, float* A, int* LDA, float* BETA, float* C, int* LDC);

  void sgesv_(int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO);
  void sposv_(char* UPLO, int* N, int* NRHS, float* A, int* LDA, float* B, int* LDB, int* INFO);
  void spotrf_(char* UPLO, int* N, float* A, int* LDA, int* INFO);

}

//////////////////////////////////////////////////////////////////////
			 // IMPLEMENTATION //
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
                       // Utility Functions //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
void Frame<SCLR>::copy(const Frame<SCLR>& M)
{
  if (this==&M) return;
  idxcheck(nr==M.rows() && nc==M.cols() && nm==M.mats());
  for(uint i = 0; i < vol(); i++) p[i] = M.vec(i);
} // copy

// I am now forcing a person to do this only in a constructor.
// void Frame::thincopy(Frame& M)
// {
//   if (this==&M) return;
//   p  = &M(0);
//   nr = M.rows();
//   nc = M.cols();
//   nm = M.mats();
// } // thincopy

template<typename SCLR>
Frame<SCLR> Frame<SCLR>::col(uint c, uint num)
{
  // Check that these are valid parameters.
  idxcheck(allow((uint)0,c));
  idxcheck(allow((uint)0,c+num-1));
  double *ptr = &operator()(0,c);
  return Frame<SCLR>(ptr, nr, num);
}

template<typename SCLR>
Frame<SCLR> Frame<SCLR>::fill(const SCLR& d)
{
  for(uint i = 0; i < vol(); i++) p[i] = d;
  return *this;
}

template<typename SCLR>
Frame<SCLR> Frame<SCLR>::dim(uint r, uint c, uint m)
{
  sizecheck (r*c*m==nr*nc*nm);
  return Frame<SCLR>(p, r, c, m);
}

//////////////////////////////////////////////////////////////////////
                  // Assigment and Test equality //
//////////////////////////////////////////////////////////////////////

/*
  The copy constructor and assignment operator work in two different
  ways.  The copy constructor simply copies the pointer p, and
  dimension info nr and nc into a new object.  These quantities are
  all small and hence there isn't much overhead to this operation.

  The assignment operator works differently.  It copies the contents
  of one MatrixFrame into the contents of another MatrixFrame.  I
  don't want to do this accidentally so I've included a warning.  In
  general, you should not use assignement and instead you should use
  the copy _function_.

  My concern is that one could be using the derived class Matrix and
  write something like c[i] = a[i], which seems to have the more
  intuitive notion of copying the conents of matrix a[i] into matrix
  c[i], but which is still ambiguous since one could just want to copy
  an address.

  To recap, the behavior you should expect:
  shallow copy:
    MatrixFrame mf1(mf2);
    MatrixFrame mf1 = mf2;
  hard copy:
    MatrixFrame mf1, mf2;
    mf1 = mf2;

  Here is something I need to think about.
  Matrix a, b;
  a.col(i) = b.col(j);
  or
  Matrix a = b.col(j);

 */

// Assignment.  See discussion above.
template<typename SCLR>
Frame<SCLR>& Frame<SCLR>::operator= (const Frame<SCLR>& M)
{
  //cerr << "Warning: be careful with the assignment operator.\n"
  //     << "MatrixFrame::operator= makes a deep copy.\n";
  if (this==&M) return *this;
  copy(M);
  return *this;
} // operator=

// Test equality.
template<typename SCLR>
bool Frame<SCLR>::operator==(const Frame<SCLR>& M) const
{
  if(sameframe(M)) return true;
  if(vol() != M.vol()) return false;
  for(uint i = 0; i < vol(); i++) if(M(i)!=operator()(i)) return false;
  return true;
} // operator==

template<typename SCLR>
bool Frame<SCLR>::sameframe(const Frame<SCLR>& M) const
{
  return (p==&M(0) && nr==M.rows() && nc==M.cols() && nm==M.mats());
}

//////////////////////////////////////////////////////////////////////
			   // Comparison //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
Frame<SCLR> lt(Frame<SCLR> c, const Frame<SCLR>& a, const Frame<SCLR>& b)
{
  sizecheck(c.vol()==a.vol() && a.vol()==b.vol());
  if(a.sameframe(b)) c.fill(1.0);
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) (a(i) < b(i));
  return c;
}

template<typename SCLR>
Frame<SCLR> lteq(Frame<SCLR> c, const Frame<SCLR>& a, const Frame<SCLR>& b)
{
  sizecheck(c.vol()==a.vol() && ( a.vol()%b.vol() )==0 );
  if(a.sameframe(b)) c.fill(0.0);
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) (a(i) <= b(i));
  return c;
}

template<typename SCLR>
Frame<SCLR> between(Frame<SCLR> c, const Frame<SCLR> a, const Frame<SCLR> lower, const Frame<SCLR> upper)
{
  sizecheck(c.vol()==a.vol() && a.vol()==lower.vol() && lower.vol()==upper.vol());
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) ( lower(i) <= a(i) && a(i) <= upper(i) );
  return c;
}

// Frame<SCLR> within(Frame<SCLR> c, const Frame<SCLR> a, double lower, double upper)
// {
//   sizecheck(c.vol()==a.vol() && a.vol()==lower.vol() && lower.vol()==upper.vol());
//   for(uint i = 0; i < a.vol(); i++)
//     c(i) = (double) ( lower <= a(i) && a(i) <= upper );
//   return c;
// }

//////////////////////////////////////////////////////////////////////
		    // Input / Output Operators //
//////////////////////////////////////////////////////////////////////

// Output matrix to stream.  This is done in a column major way so
// that a human reading the output will see an array of transposed
// matrices.

template<typename SCLR>
ostream& operator<<(ostream& os, Frame<SCLR> M)
{
  M.write(os, false, false);
  return os;
}

// Read in data from a string using scan.

template<typename SCLR>
Frame<SCLR>& operator<<(Frame<SCLR>& M, const string& s)
{
  stringstream ss(s);
  M.scan(ss, false, false);
  return M;
}

//////////////////////////////////////////////////////////////////////
			  // Read / Write //
//////////////////////////////////////////////////////////////////////

// Writes a matrix.  You may chose to include a header, which are the
// dimensions of the array of matrices.

template<typename SCLR>
bool Frame<SCLR>::write(std::ostream& os, bool header, bool binary)
{
  if (!os) return false;
  // Write binary.
  if(binary){
    if(header){
      os.write((char*) &nr, sizeof(nr));
      os.write((char*) &nc, sizeof(nc));
      os.write((char*) &nm, sizeof(nm));
    }
    for (uint i = 0; i < vol(); i++)
      os.write((char*) &p[i], sizeof(SCLR));
  }
  // Write human.
  if(!binary){
    if(header)
      os << nr << " " << nc << " " << nm << "\n";
    for(uint k = 0; k < nm; k++){
      for(uint j = 0; j < nc; j++){
	for(uint i = 0; i < nr; i++){
	  os << operator()(i,j,k) << " ";
	}
	os << "\n";
      }
      if ((k+1) != nm) os << "\n";
    }
  }
  os.flush();
  return true;
} // write

#ifndef DISABLE_FIO

template<typename SCLR>
bool Frame<SCLR>::write(const string& file, bool header, bool binary)
{
  std::ofstream ofs(file.c_str());
  if (!ofs) return false;
  return write(ofs, header, binary);
} // write

#endif

// Reads in data from a string of values until the end of the stream
// or the end of the array of matrices is reached.  You are alerted if
// you do not read in enough data to fill the array.

template<typename SCLR>
uint Frame<SCLR>::scan( std::istream& is, bool header, bool binary)
{
  // Tell us if something is wrong.
  if (!is || is.eof())  return 0;

  uint i = 0; // The nubmer of items read.

  // Read binary.
  if(binary){
    if(header){
      uint r,c,m;
      is.read((char*) &r, sizeof(nr));
      is.read((char*) &c, sizeof(nc));
      is.read((char*) &m, sizeof(nm));
      // sizecheck(vol() == r*c*m); // A somewhat strict condition.
    }
    while(!is.eof() && i < vol())
      is.read((char*) &p[i++], sizeof(SCLR));
  }
  // Write human.
  if(!binary){
    if(header){
      uint r,c,m;
      is >> r >> c >> m;
      // sizecheck(vol() == r*c*m); // A somewhat strict condition.
    }
    while(!is.eof() && i < vol())
      { is >> p[i++]; ws(is); } // ws extracts intermediate white space.
                                // Needed in case the stream is padded by white space.
  }
  // Warnings:
  if (i != vol())
    fprintf(stderr, "In scan: Number of items read (%i) different \
                     than number of elements in Matrix.\n", i);
  if (!is.eof())
    fprintf(stderr, "In scan: Did not reach end of file.\n");

  return i;
} // scan

#ifndef DISABLE_FIO

template<typename SCLR>
uint Frame<SCLR>::scan(const string& file, bool header, bool binary)
{
  std::ifstream ifs(file.c_str());
  if(!ifs){
    fprintf(stderr, "Cannot read file %s.\n", file.c_str());
    return 0;
  }
  return scan(ifs, header, binary);
} // read

#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                                                                  //
//     MATRIX OPERATIONS - IE THINGS THAT ONLY USE FIRST MATRIX     //
//                                                                  //
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
                  // Conformity and Overlap Check //
//////////////////////////////////////////////////////////////////////

// Check if MF a and MF b overlap in memory.

template<typename SCLR>
bool overlap(const Frame<SCLR>& a, const Frame<SCLR>& b)
{
  // Get the low and the high pointer.
  const Frame<SCLR> low  = &a(0) < &b(0) ? a : b;
  const Frame<SCLR> high = &a(0) < &b(0) ? b : a;
  // Check if there is overlap.  Do we need to subtract 1?  Yes.
  return (&low(0) + low.area() - 1) < &high(0) ? false : true;
} // overlap


// Hadamard conform.  This is not commutative.  Checks if the size of
// a is equally divided by the size of b.

template<typename SCLR>
bool hconform(const Frame<SCLR>& a, const Frame<SCLR>& b)
{
  return ( a.area() % b.area() )==0;
} // hconform

// Matrix product conform.  Checks if c = op(a) * op(b) is valid.
// Returns 0 if invalid and the matching dimension if valid.
template<typename SCLR>
uint pconform(const Frame<SCLR>& c, const Frame<SCLR>& a, const Frame<SCLR>& b, char transa='N', char transb='N')
{
  uint opa_rows = transa=='T' ? a.cols() : a.rows();
  uint opa_cols = transa=='T' ? a.rows() : a.cols();
  uint opb_rows = transb=='T' ? b.cols() : b.rows();
  uint opb_cols = transb=='T' ? b.rows() : b.cols();
  bool conform  = (opa_cols==opb_rows) && (c.rows()==opa_rows) && (c.cols()==opb_cols);
  if (!conform) {
    printf("c_rows: %u\n", c.rows());
    printf("c_cols: %u\n", c.cols());
    printf("opa_rows: %u\n", opa_rows);
    printf("opa_cols: %u\n", opa_cols);
    printf("opb_rows: %u\n", opb_rows);
    printf("opb_cols: %u\n", opb_cols);
  }
  return conform ? opa_cols : 0;
} // pconform

// Matrices are the same dimension.
template<typename SCLR>
bool dconform(const Frame<SCLR>& a, const Frame<SCLR>& b)
{
  return ( a.rows()==b.rows() && a.cols()==b.cols() );
}

//////////////////////////////////////////////////////////////////////
			    // COPYING //
//////////////////////////////////////////////////////////////////////

// Need to check for alisasing in memory!!!

// This matrix is smaller or the same area as M.
// Fill this matrix from M starting at (r,c).
template<typename SCLR>
void Frame<SCLR>::copy(const Frame<SCLR>& M, uint r, uint c)
{
  sizecheck( (r + nr) <= M.rows() && (c + nc) <= M.cols() );
  for(uint j = 0; j < nc; j++)
    for(uint i = 0; i < nr; i++)
      operator()(i,j) = M(r+i, c+j);
} // copy

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::copy(const Frame<SCLR>& M, const Frame<IDX>& rs, const Frame<IDX>& cs)
{
  sizecheck(rs.area()==rows() && cs.area()==cols());
  // Should check min and max of rs and cs to make sure you are in bounds.
  // I suppose this is checked by indexing.
  for(uint j = 0; j < cs.area(); j++){
    for(uint i = 0; i < rs.area(); i++){
      operator()(i,j) = M(rs(i), cs(j));
    }
  }
} // copy

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::copy(const Frame<SCLR>& M, const Frame<IDX>& rs, uint c)
{
  sizecheck(rs.area()==rows() && 1==cols());
  for(uint i = 0; i < rs.area(); i++){
    operator()(i,0) = M(rs(i), c);
  }
} // copy

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::copy(const Frame<SCLR>& M, uint r, const Frame<IDX>& cs)
{
  sizecheck(1==rows() && cs.area()==cols());
  for(uint j = 0; j < cs.area(); j++){
    operator()(0,j) = M(r, cs(j));
  }
} // copy

template<typename SCLR>
void Frame<SCLR>::copy_transpose(const Frame<SCLR>& M)
{
  sizecheck(nr==M.cols() && nc==M.rows() && nm==M.mats());
  for(uint k = 0; k < nm; ++k){
    for(uint j = 0; j < nc; ++j){
      for(uint i = 0; i < nr; ++i){
	get(i,j,k) = M.get(j,i,k);
      }
    }
  }
}

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::set(const Frame<IDX>& rs, const Frame<IDX>& cs, const Frame<SCLR>& M)
{
  sizecheck(rs.area()==M.rows() && cs.area()==M.cols());
  for(uint j = 0; j < cs.area(); j++){
    for(uint i = 0; i < rs.area(); i++){
      operator()(rs(i),cs(j)) = M(i, j);
    }
  }
}

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::set(const Frame<IDX>& rs, uint c, const Frame<SCLR>& M)
{
  sizecheck(rs.area()==M.area());
  for(uint i = 0; i < rs.area(); i++){
    operator()(rs(i), c) = M(i);
  }
}

template<typename SCLR> template<typename IDX>
void Frame<SCLR>::set(uint r, const Frame<IDX>& cs, const Frame<SCLR>& M)
{
  sizecheck(cs.area()==M.area());
  for(uint j = 0; j < cs.area(); j++){
    operator()(r, cs(j)) = M(j);
  }
}

//////////////////////////////////////////////////////////////////////
                      // Operation along rows //
//////////////////////////////////////////////////////////////////////

// Sometime we want to take an operation "along rows" i.e. the
// operation a(i,j) op= b(j % b.area()) where op may be *,+,/,-.  The
// functions to do this are <op>onrow where <op> is prod, sum, div, or
// sub.  We require that a.cols() to be a multiple of area(b).

#define ROWOP(NAME, OPEQ)			 \
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> a, Frame<SCLR>& b)	\
  {						 \
    memcheck(!overlap(a, b));			 \
    sizecheck(a.cols()%b.area()==0);		 \
    uint barea = b.area();			 \
    for(uint j = 0; j < a.cols(); j++)		 \
      for(uint i = 0; i < a.rows(); i++)	 \
        a(i,j) OPEQ b(j % barea);		 \
    return a;					 \
  }						 \

ROWOP(prodonrow, *=) ROWOP(sumonrow, +=)
ROWOP(divonrow,  /=) ROWOP(subonrow, -=)

#undef ROWOP

//////////////////////////////////////////////////////////////////////
			// Statistics //
//////////////////////////////////////////////////////////////////////

template<typename SCLR>
SCLR sum(const Frame<SCLR>& a)
{
  double total = 0.0;
  for(uint i = 0; i < a.vol(); i++)
    total += a(i);
  return total;
}

template<typename SCLR>
SCLR mean(const Frame<SCLR>& a)
{
  return sum(a) / a.vol();
}

//////////////////////////////////////////////////////////////////////
			 // BLAS / LAPACK //
//////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------//
// y = alpha x + y.
template<typename SCLR>
void axpy(SCLR alpha, Frame<SCLR> x, Frame<SCLR> y)
{
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  raxpy((int)x.rows(), alpha, &x(0), 1, &y(0), 1);
}

//------------------------------------------------------------------//
// x'y
template<typename SCLR>
SCLR dot(Frame<SCLR> x, Frame<SCLR> y)
{
  #ifndef NDEBUG
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  #endif
  return rdot(x.rows(), x.getp(), 1, y.getp(), 1);
}

//------------------------------------------------------------------//
// c = alpha op(a) * op(b) + beta c.
// void gemm(Frame<SCLR> c, Frame<SCLR> a, Frame<SCLR> b, char ta='N', char tb='N', double alpha=1.0, double beta=0.0)
template<typename SCLR>
void gemm(Frame<SCLR> c, Frame<SCLR> a, Frame<SCLR> b, char ta, char tb, SCLR alpha, SCLR beta)
{
  #ifndef NDEBUG
  memcheck(!overlap(c,a) && !overlap(c,b));
  #endif
  // Get the dimensionality information we need.
  int cnr = (int)c.rows(); int cnc = (int)c.cols();
  int anr = (int)a.rows(); int bnr = (int)b.rows();
  int k   = (int)pconform(c, a, b, ta, tb);
  // Make sure things conform.
  #ifndef NDEBUG
  sizecheck(k!=0);
  #endif
  rgemm(ta, tb, cnr, cnc, k, alpha, &a(0), anr, &b(0), bnr, beta, &c(0), cnr);
} // gemm

//------------------------------------------------------------------//
// b = alpha op(a) * b  OR  b = alpha b * op(a) where a is triangular.

// void trmm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
template<typename SCLR>
void trmm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side, char ta, char diag, SCLR alpha)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='L' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  rtrmm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trmm

//------------------------------------------------------------------//
// Solve x:  op(a) x = alpha b  OR  x op(a) = alpha b, a triangular.
// i.e: x = alpha inv(op(a)) b  OR  x = alpha b inv(op(a)).
// The solution is overwriten into b.

// void trsm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
template<typename SCLR>
void trsm(Frame<SCLR> a, Frame<SCLR> b, char uplo, char side, char ta, char diag, SCLR alpha)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='L' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  rtrsm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trsm

//------------------------------------------------------------------------------
// C := alpha*A*A**T + beta*C, ta='N'
// C := alpha*A**T*A + beta*C. ta='T'

template<typename SCLR>
void syrk(Frame<SCLR> c, Frame<SCLR> a, char ta, SCLR alpha, SCLR beta)
{
  memcheck(!overlap(c,a));
  char tb = ta=='N' ? 'T' : 'N';
  pconform(c, a, a, ta, tb);
  int k = ta=='N' ? a.cols() : a.rows();
  
  rsyrk('U', ta, c.rows(), k, alpha, &a(0), a.rows(), beta, &c(0), c.rows());

  // Better way?
  for(uint j=0; j<c.cols()-1; j++)
    for(uint i=j+1; i<c.rows(); i++)
      c(i,j) = c(j,i);
}

//////////////////////////////////////////////////////////////////////
		  // MATRIX FRAME LAPACK WRAPPER //
//////////////////////////////////////////////////////////////////////

// Solve a general linear system, ax = b for x.

template<typename SCLR>
int gesv(Frame<SCLR> a, Frame<SCLR> b)
{
  memcheck(!overlap(a, b));       // No overlap in memory.
  sizecheck(pconform(b, a, b)!=0); // a is square and b conforms.
  int info;
  std::vector<int> ipiv(a.rows());
  rgesv(a.rows(), b.cols(), &a(0), a.rows(), &ipiv[0], &b(0), b.rows(), info);
  return info;
}

// Shorthand.
template<typename SCLR>
int solve(Frame<SCLR> a, Frame<SCLR> b)
{
  return gesv(a, b);
}

//------------------------------------------------------------------//
// Solves ax = b for x where a is sym. pos. def.  Note: the lower (or
// upper) portion of A is overwritten with the Cholesky decomposition.

template<typename SCLR>
int posv(Frame<SCLR> a, Frame<SCLR> b, char uplo)
{
  memcheck(!overlap(a,b));
  sizecheck(pconform(b, a, b)!=0);
  int info;
  rposv(uplo, a.rows(), b.cols(), &a(0), a.rows(), &b(0), b.rows(), info);

  if (info != 0) {
    printf("Error in posv: info = %i\n", info);
    #ifndef NTHROW
    throw std::runtime_error("aborted in posv\n");
    #endif
  }

  return info;
}

//------------------------------------------------------------------//
// Cholesky Decomposition

template<typename SCLR>
int potrf(Frame<SCLR> a, char uplo)
{
  sizecheck(a.rows()==a.cols());
  int info = 0;
  rpotrf(uplo, a.rows(), &a(0), a.rows(), info);
  return info;
}

// int chol(Frame<SCLR> a, char uplo='L')
template<typename SCLR>
int chol(Frame<SCLR> a, char uplo)
{
  return potrf(a, uplo);
}

//------------------------------------------------------------------//

//////////////////////////////////////////////////////////////////////
                      // Hadamard Operations //
//////////////////////////////////////////////////////////////////////

// A Hadamard operation is an operation done element wise.
// Notationally, given op in {*,+,/,-} these functions perform
// something like

//    a[i] op= b[i % b.area()]
//    c[i] = alpha * a[i] op b[i % b.area()] + beta c[i].

// Thus the elements of b are recycled when area(b) < area(a).  We
// require that area(b) be a multiple of area(a).

// The functions available are h<op>eq(a, b, sc) and h<op>(c, a, b,
// alpha=0.0, beta=1.0) where <op> is prod, sum, div, or sub.

// Hadamard Operation Equals (HOPEQ) a op= sc * b
#define HOPEQ(NAME, OPEQ)						\
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> a, const Frame<SCLR>& b, SCLR sc=1.0) \
  {                                                                     \
    memcheck(!overlap(a,b));                                              \
    sizecheck(hconform(a,b));                                              \
    uint barea = b.area();						\
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ sc * b(i % barea);					\
    return a;								\
  }                                                                     \
  template<typename SCLR>						\
  Frame<SCLR> NAME(Frame<SCLR> a, SCLR b)			\
  {                                                                     \
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ b;							\
    return a;								\
  }                                                                     \
  template<typename SCLR>						\
  Frame<SCLR>& operator OPEQ(Frame<SCLR>& a, const Frame<SCLR>& b) \
  {                                                                     \
    memcheck(!overlap(a,b));                                              \
    sizecheck(hconform(a,b));                                              \
    uint barea = b.area();						\
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ b(i % barea);						\
    return a;								\
  }									\

HOPEQ(hprodeq, *=) HOPEQ(hsumeq, +=)
HOPEQ(hdiveq,  /=) HOPEQ(hsubeq, -=)

#undef HOPEQ

// Hadamard Operation (HOP) c = a op sc * b
#define HOP(NAME, OP)                                                   \
  template<typename SCLR>						\
  void NAME(Frame<SCLR> c, const Frame<SCLR>& a, const Frame<SCLR>& b, double alpha=1.0, double beta=0.0)	\
  {                                                                     \
    bool okay = (!overlap(c,b) && !overlap(a,b)) ||			\
      (!overlap(c,a) && !overlap(c,b));					\
    memcheck(okay);							\
    sizecheck(hconform(a,b));                                              \
    sizecheck(c.area()==a.area());						\
    uint barea = b.area();						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = alpha * a(i) OP b(i % barea) + beta * c(i);		\
  }                                                                     \
  template<typename SCLR>						\
  void NAME(Frame<SCLR> c, const Frame<SCLR>& a, double b)                                      \
  {                                                                     \
    sizecheck(c.area()==a.area());						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = a(i) OP b;                                                 \
  }                                                                     \
  template<typename SCLR>						\
  void NAME(Frame<SCLR> c, double b, const Frame<SCLR>& a)					\
  {                                                                     \
    sizecheck(c.area()==a.area());						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = b OP a(i);                                                 \
  }                                                                     \

HOP(hprod, *) HOP(hsum, +)
HOP(hdiv,  /) HOP(hsub, -)

#undef HOP

//////////////////////////////////////////////////////////////////////
			  // END OF CODE //
//////////////////////////////////////////////////////////////////////

#endif // MATRIX_FRAME

//////////////////////////////////////////////////////////////////////
			    // APPENDIX //
//////////////////////////////////////////////////////////////////////

/*********************************************************************

  The goal of the matrixMatrixFrame class and the Matrix class is to
  represent arrays of matrices.  An array with only one matrix can be
  thought of simply as a matrix.  We had three goals in mind when
  creating this clss: 1) keep it simple/transparent, 2) make it easy
  to use in MCMC, and 3) make it easy to talk to R.

  Regarding (1), there is a tradeoff between ease of understanding the
  code and ease of calculations.  Eigen is a great Matrix package, but
  it is hard to understand what exactly is going on.  I don't
  understand expression templates, but they make calculations nice and
  easy.  Since this will be used in MCMC simulations we wanted to know
  exactly what is going on under the hood, which means that this code
  is more understandable but that you will pay a cost when expressing
  your computations.

  Regarding (2), we hope that an array of Matrices will be sufficient
  for most MCMC algorithms.  It is possible that one may need an array
  of an array of matrices, such as an Matrix Normal DLM within a Gibbs
  sampler, but we haven't ecountered that too often in our MCMC
  experience.

  Regarding (3), you should be able to pass the pointer to some data
  from R directly to the MatrixFrame class to set up your data
  structure.

  You can think of the MatrixFrame class as simply a wrapper to
  BLAS/LAPACK, and you can think of the Matrix class as a simple data
  container.  That's all there is too it.  In that sense, MatrixFrame
  is an "interface".  The matrix operations are ALWAYS performed on
  the first matrix in the array.  If you want to perform an operation
  on different matrix in the array use the [] operator.

  Below you will find that the class is split into two categories:
  those parts of the class that have to do with the array of matrices
  and those parts of the class that have to do with the first matrix
  in the array, i.e. the matrix on which we are operating.

  In general, we need to be careful with this class.  There are lots
  of opportunities to break things since the class possess a pointer
  to some data, which the class itself did not create.

  Everything is COLUMN MAJOR for compatibility with Fortran.

  Idea: You could template this to the dimension of the array.

*********************************************************************/
