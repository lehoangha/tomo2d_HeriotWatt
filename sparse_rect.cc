/*
 * sparse_rect.cc - sparse rectangular matrix implementation
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "sparse_rect.h"
#include "error.h"
#include <iostream>
#include <fstream>

SparseRectangular::SparseRectangular(const sparseMap& A, int n)
    : ncol(n)
{
    typedef map<int,double>::const_iterator mapBrowser;

    nrow = A.size();
    int ntotal=0;
    for (int i=1; i<=nrow; i++){
	ntotal += A(i).size();
    }

    val.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);

    int k=1;
    for (int i=1; i<=nrow; i++){
	first_k(i) = k;
	for (mapBrowser p=A(i).begin(); p!=A(i).end(); p++){
	    index_j(k) = p->first;
	    val(k) = p->second;
	    k++;
	}
    }
    first_k(nrow+1) = k;

    if (ncol==0){
	// scan index_j to get the largest colmun index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k(i)-1;
	    int jmax = index_j(last_k);
	    if (jmax>ncol) ncol = jmax;
	}
    }
}

SparseRectangular::SparseRectangular(const sparseMap& A,
				     const Array1d<int>& row_index,
				     const Array1d<int>& col_index, int n)
    : ncol(n)
{
    typedef map<int,double>::const_iterator mapBrowser;

    if (A.size() != row_index.size())
	error("SparseRectangular::size mismatch detected.");

    nrow = 0;
    int ntotal=0;
    if (A.size() != row_index.size())
	error("SparseRectangular::size mismatch for row_index.");
    for (int i=1; i<=A.size(); i++){
	if (row_index(i)>0){
	    nrow++;
	    for (mapBrowser p=A(i).begin(); p!=A(i).end(); p++){
		int j = p->first;
		if (col_index.size() < j)
		    error("SparseRectangular::size mismatch for col_index.");
		if (col_index(j)>0) ntotal++;
	    }
	}
    }
    
    val.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);

    int k=1, ii=1;
    for (int i=1; i<=A.size(); i++){
	if (row_index(i)>0){
	    first_k(ii++) = k;
	    for (mapBrowser p=A(i).begin(); p!=A(i).end(); p++){
		int j = p->first;
		if (col_index(j)>0){
		    index_j(k) = col_index(j);
		    val(k) = p->second;
		    k++;
		}
	    }
	}
    }
    first_k(nrow+1) = k;

    if (ncol==0){
	// scan index_j to get the largest colmun index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k(i)-1;
	    int jmax = index_j(last_k);
	    if (jmax>ncol) ncol = jmax;
	}
    }
}

SparseRectangular::SparseRectangular(const Array1d<const sparseMap*>& As,
				     const Array1d<SparseMatAux>& spec, int n)
    : ncol(n)
{
    typedef map<int,double>::const_iterator mapBrowser;

    if (As.size() != spec.size())
	error("SparseRectangular::size mismatch detected.");

    nrow = 0;
    int ntotal=0;
    for (int iA=1; iA<=As.size(); iA++){
	if (spec(iA).p_row->size() != As(iA)->size())
	    error("SparseRectangular::size mismatch for row_index.");
	
	for (int i=1; i<=As(iA)->size(); i++){
	    if ((*spec(iA).p_row)(i)>0){
		nrow++;
		for (mapBrowser p=(*As(iA))(i).begin(); p!=(*As(iA))(i).end(); p++){
		    int j = p->first;
		    if (spec(iA).p_col->size() < j)
			error("SparseRectangular::size mismatch for col_index.");
		    if ((*spec(iA).p_col)(j)>0) ntotal++;
		}
	    }
	}
    }

    val.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);

    int k=1, ii=1;
    for (int iA=1; iA<=As.size(); iA++){
	double factor=spec(iA).scale_factor;
//	cerr << "sparse_rect " << iA << " " << factor << '\n';
	for (int i=1; i<=As(iA)->size(); i++){
	    if ((*spec(iA).p_row)(i)>0){
		first_k(ii++) = k;
		for (mapBrowser p=(*As(iA))(i).begin(); p!=(*As(iA))(i).end(); p++){
		    int j = p->first;
		    if ((*spec(iA).p_col)(j)>0){
			index_j(k) = (*spec(iA).p_col)(j);
			val(k) = factor*p->second;
//			cerr << ii << " " << index_j(k) << " " << val(k)
//			     << " " << factor << '\n';
			k++;
		    }
		}
	    }
	}
    }
    first_k(nrow+1) = k;

    if (ncol==0){
	// scan index_j to get the largest colmun index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k(i)-1;
	    int jmax = index_j(last_k);
	    if (jmax>ncol) ncol = jmax;
	}
    }
}

// calculate Ax = b
void SparseRectangular::Ax(const Array1d<double>& x, Array1d<double>& b,
			   bool init) const
{
    if (x.size() != ncol ||
	b.size() != nrow) error("SparseRectangular::Ax - size mismatch");

    if (init==true){
	for (int i=1; i<=nrow; i++) b(i) = 0.0;
    }
    for (int i=1; i<=nrow; i++){
	for (int k=first_k(i); k<first_k(i+1); k++){
	    b(i) += val(k)*x(index_j(k));
	}
    }
}

// calculate At*x = b
void SparseRectangular::Atx(const Array1d<double>& x, Array1d<double>& b,
			    bool init) const
{
    if (x.size() != nrow ||
	b.size() != ncol) error("SparseRectangular::Atx - size mismatch");

    if (init==true){
	for (int i=1; i<=ncol; i++) b(i) = 0.0;
    }
    for (int i=1; i<=nrow; i++){
	for (int k=first_k(i); k<first_k(i+1); k++){
	    int j=index_j(k);
	    b(j) += val(k)*x(i);
	}
    }
}

// dump out the contents (for debugging)
void SparseRectangular::dump(const char* fn) const
{
    ofstream os(fn);

    os << "SparseRectangular dump\n";

    os << "index_j=\n";
    for (int i=1; i<=index_j.size(); i++){
	os << index_j(i) << " ";
	if (i%10==0) os << '\n';
    }
    os << '\n';

    os << "first_k=\n";
    for (int i=1; i<=first_k.size(); i++){
	os << first_k(i) << " ";
	if (i%10==0) os << '\n';
    }
    os << '\n';

    os << "val=\n";
    for (int i=1; i<=val.size(); i++){
	os << val(i) << " ";
	if (i%10==0) os << '\n';
    }
    os << '\n';
}


