
#include "tnt/tnt.h"

void getrowdiffs(TNT::Array2D<int> &m, int *res, int *len) {
  // assumes that res is already set to zero
  for (int i=1;i<m.dim1();i++) {
    for (int j=0;j<i;j++) {
      int count=0;
      for (int k=0;k<m.dim2();k++) {
	        count += (m[i][k]!=m[j][k]);
      }
      if (count>=*len) count=*len-1;
      res[count]+=1;
    }
  }
} 

void getcoldiffmat(TNT::Array2D<int> &m, int *res) {
  // assumes that res is already set to zero
  int count=0;
  for (int i=0;i<m.dim2();i++) {
    for (int j=i+1;j<m.dim2();j++) {
      int diff=0;
      for (int k=0;k<m.dim1();k++) {
	diff += (m[k][i]!=m[k][j]);
      }
      res[count++]=diff;
    }
  }
} 

extern "C" {
  void rowdiffs(int *data, int *rows, int *cols, int *res, int *len) {    
    TNT::Array2D<int> d(*rows,*cols,data);
    getrowdiffs(d,res,len);
  }
  void coldiffmat(int *data, int *rows, int *cols, int *res) {    
    TNT::Array2D<int> d(*rows,*cols,data);
    getcoldiffmat(d,res);
  }

}
