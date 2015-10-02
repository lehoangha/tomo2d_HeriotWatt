/*
 * sparse_rect.h - sparse rectangular matrix interface
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_SPARSE_RECT_H_
#define _TOMO_SPARSE_RECT_H_

#include <map>
#include <array.h> // from mconv

struct SparseMatAux {
    SparseMatAux(double d=1.0, const Array1d<int> *r=0,
		 const Array1d<int> *c=0)
	: scale_factor(d), p_row(r), p_col(c) {}
    
    double scale_factor;
    const Array1d<int> *p_row;
    const Array1d<int> *p_col;
};
    
class SparseRectangular {
public:
    typedef Array1d< map<int,double> > sparseMap;
    
    SparseRectangular(){};
    SparseRectangular(const sparseMap& A, int ncol=0);
    SparseRectangular(const sparseMap& A, const Array1d<int>& r,
		      const Array1d<int>& c, int ncol=0);
    SparseRectangular(const Array1d<const sparseMap*>& As,
		      const Array1d<SparseMatAux>& spec, int ncol=0);

    int nRow() const { return nrow; }
    int nCol() const { return ncol; }
    void Ax(const Array1d<double>&, Array1d<double>&, bool init=true) const;
    void Atx(const Array1d<double>&, Array1d<double>&, bool init=true) const;
    void dump(const char*) const;
    
private:
    int nrow, ncol;
    Array1d<double> val;
    Array1d<int> index_j;
    Array1d<int> first_k;
};

#endif /* _TOMO_SPARSE_RECT_H_ */
