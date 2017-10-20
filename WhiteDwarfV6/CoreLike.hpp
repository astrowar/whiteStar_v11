//
//  CoreLike.hpp
//  whiteStar
//
//  Created by Eraldo Rangel on 26/12/15.
//  Copyright © 2015 Eraldo Rangel. All rights reserved.
//

#ifndef CoreLike_hpp
#define CoreLike_hpp

#include <stdio.h>
#include <vector>
#include <cstdlib>

class argumento {
    
};

class vrange : argumento {
    
    
public:
     float xa,xb,dx;
    vrange(  float xinitial, float xfinal, float _dx)
    {
        xa= xinitial;
        xb = xfinal;
        dx =_dx;
    }
    
    int nn(){ return (xb- xa)/dx;  }
};

class nrange :public  vrange {
    
    int id;
public:
    nrange(int _id,float xinitial, float xfinal, float _dx  ): vrange( xinitial,xfinal,_dx)
    {
        id =_id;
    }
    
};

int totalSize ( std::vector< vrange > v ,int i = 0){
    if(i< v.size()-1  )
    {
        return  v[i].nn()*totalSize(v,i+1);
    }
    return v[i].nn();
}
int totalSize ( std::vector< nrange > v ,int i = 0){
    if(i< v.size()-1  )
    {
        return  v[i].nn()*totalSize(v,i+1);
    }
    return v[i].nn();
}

class xarg {
    int id;
    float x;
public:
    xarg(int _id, float _x)
    {
        id = _id;
        x = _x;
    }
};

class lmatrix {
    std::vector< nrange> ranges;
    std::vector< size_t> stride;
    float *data;
    lmatrix( std::vector< nrange >  _ranges ){
        ranges = _ranges;
         int  mm = sizeof( float) * totalSize(ranges);
         data = (float*) malloc(mm);
        size_t n = ranges.size();
        stride.resize(n,1);
        size_t st = 1;
        
        for(int i = n-1 ; i >= 0;--i)
        {
            stride[i - (n-1)] = st; // ultimo sempre eh 1
            st *= ranges[i].nn();
        }
   }
    
    
};
//template< unsigned int R,  typename... Ts >
//void setRange(tMatrix<  R  > &m , int n, trange range1 , Ts... range2 ) {

float getLM_j(lmatrix *m , nrange xa ){
    return 0.0;
}

float getLM_j(lmatrix *m , xarg xa ){
    return 0.0;
}

template< typename T, typename...  Ts  >
float getLM( lmatrix *m , nrange xa , T xb, Ts... args )
{
    return getLM(m, xb, args...);
}

template< typename T, typename...  Ts  >
float getLM(   lmatrix *m , xarg xa ,T xb, Ts... args )
{
    return getLM(m, xb, args...);
}

template<typename T >
float getLM( lmatrix *m , nrange xa , T xb )
{
    return getLM_j(m,  xb);
}

template<typename T >
float getLM( lmatrix *m , xarg xa , T xb )
{
    return getLM_j(m,  xb);
}



#endif /* CoreLike_hpp */
