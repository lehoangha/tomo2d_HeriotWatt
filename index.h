/*
 * index.h
 *
 * Jun Korenaga, MIT/WHOI
 * January 19999
 */

#ifndef _TOMO_INDEX_H_
#define _TOMO_INDEX_H_

class Index2d {
public:
    Index2d(){}
    Index2d(int i1, int i2){ i_[0]=i1; i_[1]=i2; }

    int i() const { return i_[0]; }
    int k() const { return i_[1]; }
    void set(int i1, int i2) { i_[0]=i1; i_[1]=i2; }
    
private:
    int i_[2];
};


#endif /* _TOMO_INDEX_H_ */
