#include <Rcpp.h>
using namespace Rcpp;





// --------------------Function of summing all elements in a matrix
// [[Rcpp::export]]
double summ(NumericMatrix a){
  double col=a.ncol();
  NumericVector tt(col);
  
  for(int i=0; i<col; ++i){
    tt[i]=sum(a(_,i));
  }
  
  return sum(tt);
  
}
// 
// /*** R
// a<-matrix(1:9, 3,3)
// summ(a)
//   */
// 




// // [[Rcpp::export]]
// Function of sum of outer product on a vector and a matrix
// the elements of the vector is either 0 or 1
// double outerm(NumericMatrix a, NumericMatrix b){
//   double asize=a.size();
//   NumericVector out(asize);
//
//   for(int i=0; i<asize; ++i){
//    if(a[i]==0){out[i]=0;}else{out[i]=summ(b);}
//   }
//   return sum(out);
//
// }


//-------------Function of sum of outer product on a vector and a matrix
// [[Rcpp::export]]
double outerm(NumericMatrix a, NumericMatrix b){

  double aa=sum(a);

  return summ(b*aa);

}
//
// /*** R
// aa<-matrix(c(1,0,1), 3,1)
// bb<-matrix(1:9, 3,3)
// outerm(aa,bb)
//   */



// ------------------Submatrix by performs similar to matrix a[-c(low:up), -b] in R 
// [[Rcpp::export]]
NumericMatrix nsubmat(NumericMatrix a, int low, int up, int column){
  double row=a.nrow();
  double col=a.ncol();
  NumericMatrix out(row-up+low-1, col);
  NumericMatrix outf(row-up+low-1, col-1);
  int counteri=0;
  int counterj=0;

  for(int i=0; i<row; ++i){
    if(( (i<low) & (i<up))|( (i>low) & (i>up))){out(counteri,_)=a(i,_);++counteri;}
    else{continue;}
  }

  for(int j=0; j<col; ++j){
    if(j==column){continue;}
    else{outf(_,counterj)=out(_,j);++counterj;}
  }

  return outf;
}
//
// /*** R
// bb<-matrix(1:9, 3,3)
// nsubmat(bb, 1,2,0)
// */
//



// //----------------- Re-generate the R function: xi12GTy12
// // [[Rcpp::export]]
// double xi12GTy12_Rcpp(NumericMatrix win, int nf, int nv, NumericVector njf, NumericVector njv){
// 
//   NumericVector ef=cumsum(njf);
//   NumericVector sf=ef-njf;
//   NumericVector ev=cumsum(njv);
//   NumericVector sv=ev-njv;
//   NumericVector s3(nf);
// 
//   for(int i=0; i<nf; ++i){
//     NumericMatrix w = win(_, Range(sf[i], ef[i]-1) );
//     int ncol=w.ncol();
//     NumericVector s2(ncol);
// 
//     for(int j=0; j<ncol; ++j){
//       NumericVector s1(nv);
//       for(int k=0; k<nv; ++k){
//         NumericMatrix w1=w(Range(sv[k], ev[k]-1), Range(j,j));
//         NumericMatrix w2=nsubmat(w, sv[k], ev[k]-1, j);
//         s1[k]=outerm(w1, w2);
//       }
//       s2[j]=sum(s1);
//     }
//     s3[i]=sum(s2);
//   }
// 
//   return sum(s3);
// 
// }
// 
// //
// // /**** R
// // xi12GTy12_Rcpp(t(win_t),m, n, Nj_t, Nj_c)
// // */
// //





//---------------------- Loop of finding sum per column by row clusters from win matrix.
// [[Rcpp::export]]
NumericMatrix loopsum(NumericMatrix a, int nv, NumericVector njv){
  int col=a.ncol();
  NumericMatrix b(nv, col);
  NumericVector ev=cumsum(njv);
  NumericVector sv=ev-njv;

  for(int i=0; i<col; ++i){
    for(int j=0; j<nv; ++j){
      NumericMatrix w1=a(Range(sv[j], ev[j]-1), Range(i,i));
      b(j,i)=sum(w1);
    }
  }

  return b;
}
// 
// /****R
// a<-matrix(c(1:(20000*20000)), 20000,20000)
// nv<-100
// njv<-rep(200,100)
// out<-loopsum(a, nv, njv)
// */
// 








//---------------- Submatrix by performs similar to matrix a[-b, -d]in R
// [[Rcpp::export]]
NumericMatrix offsubmat(NumericMatrix a, int b, int d){
  double row=a.nrow();
  double col=a.ncol();
  NumericMatrix out(row-1, col);
  NumericMatrix outf(row-1, col-1);
  int counteri=0;
  int counterj=0;

  for(int i=0; i<row; ++i){
    if(i==b){continue;}
    else{out(counteri,_)=a(i,_);++counteri;}
  }

  for(int j=0; j<col; ++j){
    if(j==d){continue;}
    else{outf(_,counterj)=out(_,j);++counterj;}
  }

  return outf;
}
// 
// /****R
// a<-matrix(c(1:9), 3,3)
// b<-2
// d<-1
// offsubmat(a, b, d)
// */
// 




//---------------- Submatrix by performs similar to matrix a[-b, fixed column]in R
// [[Rcpp::export]]
NumericMatrix offsubmat_col(NumericMatrix a, int b){
  double row=a.nrow();
  double col=a.ncol();
  NumericMatrix out(row-1, col);
  int counteri=0;
  
  for(int i=0; i<row; ++i){
    if(i==b){continue;}
    else{out(counteri,_)=a(i,_);++counteri;}
  }
  
  return out;
}
// 
// /****R
// a<-matrix(c(1:9), 3,3)
// b<-2
// offsubmat_col(a, b)
// */
// 







//------------------ Submatrix by performs similar to matrix a[-r, -c] in R for each element and sum the submatrix.
// [[Rcpp::export]]
NumericMatrix sumsubmat(NumericMatrix a){
  double row=a.nrow();
  double col=a.ncol();
  NumericMatrix outi(row-1, col-1);
  NumericMatrix out(row,col);

  for(int i=0; i<row; ++i){
    for(int j=0; j<col; ++j){
      outi=offsubmat(a, i, j);
      out(i,j)=summ(outi);
    }
  }
  return out;
}
// 
// 
// /******R
// a<-matrix(c(1:(100*20000)), 100,20000)
// test<-sumsubmat(a)
// dim(test)
// */






//---------------------------------- Element-wise matrix multiplication
// [[Rcpp::export]]
NumericVector multMat(NumericMatrix m1, NumericMatrix m2) {
  NumericVector multMatrix = m1 * m2;
  multMatrix.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  return multMatrix;
}
// 
// /********R
// multMat( matrix(1:9, nrow=3), matrix(1:9, nrow=3) )
// */
// 



//------------------ Submatrix by performs similar to matrix a[-r, fixed column] in R for each element and sum the submatrix.
// [[Rcpp::export]]
NumericMatrix sumsubmat_col(NumericMatrix a){
  double row=a.nrow();
  double col=a.ncol();
  NumericMatrix outi(row-1, col);
  NumericMatrix out(row, col);
  
  for(int j=0; j<col; ++j){
    NumericMatrix temp=a(_,Range(j,j));
    for(int i=0; i<row; ++i){
      outi=offsubmat_col(temp, i);
      out(i,j)=summ(outi);
    }
  }
  return out;
}

// 
// /******R
// a<-matrix(c(1:(3*3)), 3,3)
// test<-sumsubmat_col(a)
// dim(test)
// */





//---------------------MAIN FUNCTION---------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------

//---------------- Within pair comparisons---------------------------

//---------Treatment individual wins:
// condition (1)
// [[Rcpp::export]]
NumericMatrix cond1(NumericVector cat_c,
                    NumericVector y2_c,
                    NumericVector y2_t
                    ){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<col; i++){
    for(int j=0; j<row; j++){
      if((cat_c[i]==1) | (cat_c[i]==2)){
        cond(j,i)=y2_t[j]>=y2_c[i];
      }else{cond(j,i)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond1(cat_c=data_c$cat, y2_c=data_c$y2, y2_t=data_t$y2)
// sum(test)
// */
// 


// condition (3)
// (3)_1
// [[Rcpp::export]]
NumericMatrix cond3_1(NumericVector cat_c,
                    NumericVector y2_c,
                    NumericVector y2_t,
                    NumericVector y1_c,
                    NumericVector y1_t
){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<col; i++){
    for(int j=0; j<row; j++){
      if(cat_c[i]==4){
        cond(j,i)=(y2_t[j]>=y2_c[i])&(y1_t[j]>=y1_c[i]);
      }else{cond(j,i)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond3_1(cat_c=data_c$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
// sum(test)
// */
// 

// (3)_2
// [[Rcpp::export]]
NumericMatrix cond3_2(NumericVector cat_c,
                      NumericVector cat_t,
                      NumericVector y2_c,
                      NumericVector y2_t,
                      NumericVector y1_c,
                      NumericVector y1_t
){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<col; i++){
    for(int j=0; j<row; j++){
      if(((cat_c[i]==1)|(cat_c[i]==4))&((cat_t[j]==3)|(cat_t[j]==4))){
        cond(j,i)=(y2_t[j]<y2_c[i])&(y1_t[j]>=y1_c[i]);
      }else{cond(j,i)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond3_2(cat_c=data_c$cat, cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
// sum(test)
// */



//---------Control individual wins:
// condition (2)
// [[Rcpp::export]]
NumericMatrix cond2(NumericVector cat_t,
                    NumericVector y2_c,
                    NumericVector y2_t
){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      if((cat_t[i]==1) | (cat_t[i]==2)){
        cond(i,j)=y2_t[i]<y2_c[j];
      }else{cond(i,j)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond2(cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2)
// sum(test)
// */



// condition (4)
// (4)_1
// [[Rcpp::export]]
NumericMatrix cond4_1(NumericVector cat_t,
                      NumericVector y2_c,
                      NumericVector y2_t,
                      NumericVector y1_c,
                      NumericVector y1_t
){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      if(cat_t[i]==4){
        cond(i,j)=(y2_t[i]<y2_c[j])&(y1_t[i]<y1_c[j]);
      }else{cond(i,j)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond4_1(cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
// sum(test)
// */
// 


// (4)_2
// [[Rcpp::export]]
NumericMatrix cond4_2(NumericVector cat_c,
                      NumericVector cat_t,
                      NumericVector y2_c,
                      NumericVector y2_t,
                      NumericVector y1_c,
                      NumericVector y1_t
){
  int row=y2_t.size();
  int col=y2_c.size();
  NumericMatrix cond(row,col);
  
  for(int i=0; i<row; i++){
    for(int j=0; j<col; j++){
      if(((cat_t[i]==1)|(cat_t[i]==4))&((cat_c[j]==3)|(cat_c[j]==4))){
        cond(i,j)=(y2_t[i]>=y2_c[j])&(y1_t[i]<y1_c[j]);
      }else{cond(i,j)=0;}
    }
  }
  return cond;
}
// 
// /*******R
// test<-cond4_2(cat_c=data_c$cat, cat_t=data_t$cat, y2_c=data_c$y2, y2_t=data_t$y2, y1_c=data_c$y1, y1_t=data_t$y1)
// sum(test)
// */
// 
// 





//------------------ Re-generate R function: com_less
// [[Rcpp::export]]
NumericVector com_less_Rcpp(NumericVector a,
                            NumericVector b, 
                            NumericVector c, 
                            NumericVector d){
  int len0=a.size();
  int len1=c.size();
  NumericVector temp(3);
  NumericMatrix com(len1,len0);
  
  for(int i=0; i<len0; ++i){
    for(int j=0; j<len1; ++j){
      temp[0]=b[i];
      temp[1]=c[j];
      temp[2]=d[j];
      if(a[i]<min(temp)){com(j,i)=1;}
      else{com(j,i)=0;}
    }
  }
  return com;
}
// 
// /****R
// test<-com_less_Rcpp(data_c$time_Fatal, data_c$time_censor, data_t$time_censor, data_t$time_Fatal)
// sum(test)
// */
// 


//------------------ Re-generate R function: com_more
// [[Rcpp::export]]
NumericMatrix com_more_Rcpp(NumericVector a, NumericVector b, NumericVector c, NumericVector d){
  int len1=a.size();
  int len0=c.size();
  NumericVector temp(3);
  NumericMatrix com(len1,len0);
  
  for(int j=0; j<len1; ++j){
    for(int i=0; i<len0; ++i){
      temp[0]=b[j];
      temp[1]=c[i];
      temp[2]=d[i];
      if(a[j]>=min(temp)){com(j,i)=1;}
      else{com(j,i)=0;}
    }
  }
  return com;
}
// 
// /****R
// test2<-com_more_Rcpp(data_t$time_Fatal, data_t$time_censor, data_c$time_censor, data_c$time_Fatal)
// sum(test2)
// */
// 






//------------- Re-generate R function: x12GTy12
// [[Rcpp::export]]
double x12GTy12_Rcpp_test(NumericMatrix win, 
                           int nf, 
                           int nv, 
                           NumericVector njf, 
                           NumericVector njv){
    
  NumericVector ef=cumsum(njf);
  NumericVector sf=ef-njf;
  NumericVector s3(nf);
  int col=win.ncol();
  NumericMatrix out1(nv,col);
  NumericMatrix out2(nv,col);
  
  for(int i=0; i<nf; ++i){
    NumericMatrix w = win(_, Range(sf[i], ef[i]-1) );
    
    out1=loopsum(w, nv, njv);
    out2=sumsubmat_col(out1);
    
    s3[i]=sum(multMat(out1,out2));
  }
  
  return sum(s3);
  
}
// 
// /**** R
// x12GTy12_Rcpp_test(t(win_t),m, n, Nj_t, Nj_c)
// */
// 




//----------------- Re-generate the R function: xi12GTy12
// [[Rcpp::export]]
double xi12GTy12_Rcpp_test(NumericMatrix win, 
                           int nf, 
                           int nv, 
                           NumericVector njf, 
                           NumericVector njv){

  NumericVector ef=cumsum(njf);
  NumericVector sf=ef-njf;
  NumericVector s3(nf);
  int col=win.ncol();
  NumericMatrix out1(nv,col);
  NumericMatrix out2(nv,col);

  for(int i=0; i<nf; ++i){
    NumericMatrix w = win(_, Range(sf[i], ef[i]-1) );

    out1=loopsum(w, nv, njv);
    out2=sumsubmat(out1);

    s3[i]=sum(multMat(out1,out2));
  }

  return sum(s3);

}
// 
// /**** R
// xi12GTy12_Rcpp_test(t(win_t),m, n, Nj_t, Nj_c)
// */
// 




//----------------- Regenerate R function: covx12GTy12
// [[Rcpp::export]]
double covx12GTy12_Rcpp_test(NumericMatrix win1,
                             NumericMatrix win2,
                             int nf, 
                             int nv, 
                             NumericVector njf, 
                             NumericVector njv){
  
  NumericVector ef=cumsum(njf);
  NumericVector sf=ef-njf;
  NumericVector s3(nf);
  int col=win1.ncol();
  NumericMatrix out1(nv,col);
  NumericMatrix out2(nv,col);
  NumericMatrix out3(nv,col);
  
  for(int i=0; i<nf; ++i){
    NumericMatrix w1 = win1(_, Range(sf[i], ef[i]-1) );
    NumericMatrix w2 = win2(_, Range(sf[i], ef[i]-1) );
    
    out1=loopsum(w1, nv, njv);
    out2=loopsum(w2, nv, njv);
    out3=sumsubmat_col(out2);
    
    s3[i]=sum(multMat(out1,out3));
  }
  
  return sum(s3);
  
}
// 
// /**** R
// covx12GTy12_Rcpp_test(t(win_t), t(win_c),m, n, Nj_t, Nj_c)
// */
// 




//------------------ Re-generate R function: covxi12GTy12
// [[Rcpp::export]]
double covxi12GTy12_Rcpp_test(NumericMatrix win1,
                              NumericMatrix win2,
                              int nf, 
                              int nv, 
                              NumericVector njf, 
                              NumericVector njv){
  
  NumericVector ef=cumsum(njf);
  NumericVector sf=ef-njf;
  NumericVector s3(nf);
  int col=win1.ncol();
  NumericMatrix out1(nv,col);
  NumericMatrix out2(nv,col);
  NumericMatrix out3(nv,col);
  
  for(int i=0; i<nf; ++i){
    NumericMatrix w1 = win1(_, Range(sf[i], ef[i]-1) );
    NumericMatrix w2 = win2(_, Range(sf[i], ef[i]-1) );
    
    out1=loopsum(w1, nv, njv);
    out2=loopsum(w2, nv, njv);
    out3=sumsubmat(out2);
    
    s3[i]=sum(multMat(out1,out3));
  }
  
  return sum(s3);
  
}
// 
// /**** R
// covxi12GTy12_Rcpp_test(t(win_t), t(win_c),m, n, Nj_t, Nj_c)
// */













