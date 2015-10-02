/*
 * heap_deque.h - heap-structured deque class
 *                (allow iteration over elements unlike priority_queue)
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_HEAP_DEQUE_H_
#define _TOMO_HEAP_DEQUE_H_

#include <deque>
#include <algorithm>
#include <functional>
#include <iostream>

template<class T, class Cmp>
class heap_deque : public deque<T> {
public:
    typedef typename heap_deque<T,Cmp>::const_iterator const_iterator;

    heap_deque(const Cmp& cmp) : _cmp(cmp) { make_heap(deque<T>::begin(),deque<T>::end(),_cmp); }

	 void push(const T& a){ push_back(a); push_heap(deque<T>::begin(),deque<T>::end(),_cmp); }
	 void pop() { pop_heap(deque<T>::begin(),deque<T>::end(),_cmp); deque<T>::pop_back(); }
	 const T& top() const { return deque<T>::front(); }

    void promote(const_iterator);
    void demote(const_iterator){ cerr << "heap_deque.demote(): not implemented yet\n"; };
    void push_only(const T& a){ push_back(a); }
    void reheap(){ make_heap(deque<T>::begin(),deque<T>::end(),_cmp); }
    
private:
    const Cmp& _cmp;
};

template<class T, class Cmp>
inline
void heap_deque<T, Cmp>::promote(const_iterator p)
{
     int orig_index = p-deque<T>::begin()+1;
    if (orig_index==1) return; // no need to be promoted
    int sup_index  = orig_index/2;
    while (_cmp(*(deque<T>::begin()+sup_index-1),*(deque<T>::begin()+orig_index-1))){
	T tmp = *(deque<T>::begin()+orig_index-1);
	*(deque<T>::begin()+orig_index-1) = *(deque<T>::begin()+sup_index-1);
	*(deque<T>::begin()+sup_index-1) = tmp;

	if (sup_index==1) break;
	orig_index = sup_index;
	sup_index = orig_index/2;
    }
}

#endif /* _TOMO_HEAP_DEQUE_H_ */
