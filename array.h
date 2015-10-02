/*
 * array.h
 */

#ifndef _JK_NEW_ARRAY_H_
#define _JK_NEW_ARRAY_H_

#include <vector>
#include <cstddef>
#include <iostream>
#include "error.h"

//using namespace std;

template<class T>
class Array1d : public vector<T> {
public:
    Array1d(size_t n=0) : vector<T>(n) {}
    Array1d(const T*, size_t);	// array initialization
    
    T& operator()(size_t i) {
//	if (i>size()){cerr << "Array1d OoR " << i << " " << size() << '\n';}
	return operator[](i); }
    const T& operator()(size_t i) const {
//	if (i>size()){cerr << "Array1d OoR " << i << " " << size() << '\n';}
	return operator[](i); }

    size_t size() const { return vector<T>::size(); }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back()))+1; }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back()))+1; }


    void resize(size_t n, const T& val = T()) { vector<T>::resize(n, val); }
    void push_back (const T& x) { vector<T>::push_back(x); }

    // for Numerical recipes' vector
    T*       toRecipe() { return begin()-1; }
    
    // unary operators
    Array1d<T> operator-();
    
    // binary operators
    Array1d<T>& operator=(const T& val);
    Array1d<T>& operator+=(const Array1d<T>&);
    Array1d<T>& operator+=(const T& val);
    Array1d<T>& operator-=(const Array1d<T>&);
    Array1d<T>& operator-=(const T& val);
    Array1d<T>& operator*=(const T&);
    Array1d<T>& operator/=(const T&);
    
protected:
          T& operator[](size_t i) { return vector<T>::operator[](i-1); }
    const T& operator[](size_t i) const { return vector<T>::operator[](i-1); }
};

template<class T>
inline
Array1d<T>::Array1d(const T* init, size_t n)
    : vector<T>(n)
{
    T* dest = begin()+n;
    const T* src = init+n;
    while (dest > begin()) *--dest = *--src;
}
    
template<class T>
inline
Array1d<T> Array1d<T>::operator-()
{
    Array1d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array1d<T>& Array1d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}


// nonmember functions
template<class T>
inline
Array1d<T> operator+(const Array1d<T>& a, const Array1d<T>& b)
{
    Array1d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array1d<T> operator-(const Array1d<T>& a, const Array1d<T>& b)
{
    Array1d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array1d<T> operator*(const Array1d<T>& a, const T& val)
{
    Array1d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array1d<T> operator/(const Array1d<T>& a, const T& val)
{
    Array1d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array1d<T>& a)
{
    for (int i=1; i<a.size(); i++){
	s << a(i) << ", ";
    }
    s << a(a.size()) << '\n';
    return s;
}

// Array2d
template<class T>
class Array2d : public vector<T> {
public:
    Array2d() : vector<T>(0) {}
    Array2d(size_t n1, size_t n2)
	: vector<T>(n1*n2), nrow(n1), ncol(n2) {}
    ~Array2d();
    
    T&       operator()(size_t i, size_t j)
	{ return vector<T>::operator[](offset(i,j)); }
    const T& operator()(size_t i, size_t j) const
	{ return vector<T>::operator[](offset(i,j)); }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back())+1); }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back())+1); }
    
    void resize(size_t n1, size_t n2, const T& val = T())
	{ vector<T>::resize(n1*n2, val); nrow=n1, ncol=n2; }

    size_t nCol() const { return ncol; }
    size_t nRow() const { return nrow; }

    T** toRecipe();

    // unary operators
    Array2d<T> operator-();
    
    // binary operators
    Array2d<T>& operator=(const T& val);
    Array2d<T>& operator+=(const Array2d<T>&);
    Array2d<T>& operator+=(const T& val);
    Array2d<T>& operator-=(const Array2d<T>&);
    Array2d<T>& operator-=(const T& val);
    Array2d<T>& operator*=(const T&);
    Array2d<T>& operator/=(const T&);
    
private:
    size_t size() const { return vector<T>::size(); }

    size_t offset(size_t i, size_t j) const;

    size_t ncol;
    size_t nrow;

    vector<T*> datapp;			// for toRecipe();
};

template<class T> 
inline size_t Array2d<T>::offset(size_t i, size_t j) const
{
//    if ((i-1)*ncol+j-1 >=size()){
//	cerr << "Array2d OoR " << i << " " << j << " " << nrow << " " << ncol << '\n';
//    }
    return (i-1)*ncol+j-1;
}

template<class T>
Array2d<T>::~Array2d()
{
}

template<class T>
T** Array2d<T>::toRecipe()
{
    datapp.resize(nrow+1);
    datapp[1] = begin()-1;
    for (int i=2; i<=nrow; i++) datapp[i] = datapp[i-1]+ncol;
    return &(datapp.front());
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array2d<T> Array2d<T>::operator-()
{
    Array2d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array2d<T>& Array2d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}

// nonmember functions
template<class T>
inline
Array2d<T> operator+(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array2d<T> operator-(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array2d<T> operator*(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array2d<T> operator*(const T& val, const Array2d<T>& a)
{
    return operator*(a,val);
}

template<class T>
inline
Array2d<T> operator/(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array2d<T>& a)
{
    for (int i=1; i<=a.nRow(); i++){
	for (int j=1; j<a.nCol(); j++){
	    s << a(i,j) << ", ";
	}
	s << a(i,a.nCol()) << '\n';
    }
    return s;
}

// multiplication
// note: this can be more optimized by directly accessing member data
//       if these are friend of Array class.
//       How to deal with offset() is the key?
//
// matrix * vector
template<class T>
Array1d<T> operator*(const Array2d<T>& A, const Array1d<T>& x)
{
    int m = A.nRow();
    int n = x.size();
    if (A.nCol() != n) error("matrix*vector::size mismatch");

    Array1d<T> b(m);
    for (int i=1; i<=m; i++){
	double val=0;
	for (int j=1; j<=n; j++){
	    val += A(i,j)*x(j);
	}
	b(i) = val;
    }

    return b;
}

// matrix * matrix
template<class T>
Array2d<T> operator*(const Array2d<T>& A, const Array2d<T>& B)
{
    int l = A.nRow();
    int m = A.nCol();
    int n = B.nCol();
    if (B.nRow() != m) error("matrix*matrix::size mismatch");

    Array2d<T> C(l,n);
    for (int i=1; i<=l; i++){
	for (int j=1; j<=n; j++){
	    double val=0;
	    for (int k=1; k<=m; k++){
		val += A(i,k)*B(k,j);
	    }
	    C(i,j) = val;
	}
    }

    return C;
}

// Array3d
template<class T>
class Array3d : public vector<T> {
public:
    Array3d() : vector<T>(0) {}
    Array3d(size_t n1, size_t n2, size_t n3)
	: vector<T>(n1*n2*n3), nn1(n1), nn2(n2), nn3(n3) {}
    ~Array3d();
    
    T&       operator()(size_t i, size_t j, size_t k)
	{ return vector<T>::operator[](offset(i,j,k)); }
    const T& operator()(size_t i, size_t j, size_t k) const
	{ return vector<T>::operator[](offset(i,j,k)); }
    T* begin() { return vector<T>::begin(); }
    T* end() { return vector<T>::end(); }
    const T* begin() const { return vector<T>::begin(); }
    const T* end() const { return vector<T>::end(); }
    
    void resize(size_t n1, size_t n2, size_t n3, const T& val = T())
	{
	    vector<T>::resize(n1*n2*n3, val);
	    nn1=n1, nn2=n2; nn3=n3;
	}

    size_t nSize1() const { return nn1; }
    size_t nSize2() const { return nn2; }
    size_t nSize3() const { return nn3; }

    // unary operators
    Array3d<T> operator-();
    
    // binary operators
    Array3d<T>& operator=(const T& val);
    Array3d<T>& operator+=(const Array3d<T>&);
    Array3d<T>& operator+=(const T& val);
    Array3d<T>& operator-=(const Array3d<T>&);
    Array3d<T>& operator-=(const T& val);
    Array3d<T>& operator*=(const T&);
    Array3d<T>& operator/=(const T&);
    
private:
    size_t size() const { return vector<T>::size(); }
    size_t offset(size_t i, size_t j, size_t k) const;

    size_t nn1, nn2, nn3;
};

template<class T> 
inline size_t Array3d<T>::offset(size_t i, size_t j, size_t k) const
{
//    if (((i-1)*nn2+(j-1))*nn3+(k-1) >=size()){
//	cerr << "Array3d OoR " << i << " " << j << " " << k << " "
//	     << nn1 << " " << nn2 << " " << nn3 << '\n';
//    }
    return ((i-1)*nn2+(j-1))*nn3+(k-1);
}

template<class T>
Array3d<T>::~Array3d()
{
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array3d<T> Array3d<T>::operator-()
{
    Array3d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator+=(const Array3d<T>& a)
{
    if (size() != a.size()) error("Array3d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator-=(const Array3d<T>& a)
{
    if (size() != a.size()) error("Array3d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array3d<T>& Array3d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}

// nonmember functions
template<class T>
inline
Array3d<T> operator+(const Array3d<T>& a, const Array3d<T>& b)
{
    Array3d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array3d<T> operator-(const Array3d<T>& a, const Array3d<T>& b)
{
    Array3d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array3d<T> operator*(const Array3d<T>& a, const T& val)
{
    Array3d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array3d<T> operator*(const T& val, const Array3d<T>& a)
{
    return operator*(a,val);
}

template<class T>
inline
Array3d<T> operator/(const Array3d<T>& a, const T& val)
{
    Array3d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array3d<T>& a)
{
    for (int i=1; i<=a.nSize1(); i++){
	for (int j=1; j<=a.nSize2(); j++){
	    for (int k=1; k<=a.nSize3(); k++){
		s << a(i,j,k) << ", ";
	    }
	    s << '\n';
	}
	s << "\n\n";
    }
    return s;
}

#endif /* _JK_NEW_ARRAY_H_ */
