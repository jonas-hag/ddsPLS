#include <Rcpp.h>
using namespace Rcpp;

void prodMat(const NumericMatrix a,const  NumericMatrix b,NumericMatrix out,
             const  int n,const  int p,const  int q) {
  for (int i = 0; i < p; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < q; j++)
        out(i,j) += a(i,k) * b(k,j);
}

void prodMatVector(const NumericMatrix a,const  NumericVector b,NumericVector out,
                   const int q,const  int p,const  bool normalize=true) {
  for (int i = 0; i < q; i++)
    for (int k = 0; k < p; k++)
      out(i) += a(i,k) * b(k);
  if (normalize==true){
    out = out/sqrt(sum(out*out));
  }
}

void crossprodMatVector(const NumericMatrix a,const  NumericVector b, NumericVector out,
                        const int p,const  int q,
                        const  bool normalize=true,const  bool normalize2=false,
                        const double norm2t=0.0) {
  for (int i = 0; i < p; i++)
    for (int k = 0; k < q; k++){
      if (normalize2==true){
        out(i) += a(k,i) * b(k)/norm2t;
      } else {
        out(i) += a(k,i) * b(k);
      }
    }
  if (normalize==true){
    if (normalize2==false){
      out = out/sqrt(sum(out*out));
    }
  }
}

void prodVectVec(const NumericVector t, const  NumericVector u, NumericMatrix out,
                 const int n,const  int p) {
  for (int i = 0; i < n; i++)
    for (int k = 0; k < p; k++)
      out(i,k) += t(i) * u(k);
}

void crossprodMat(const NumericMatrix a,const NumericMatrix b,NumericMatrix out,
                  int n, int p, int q) {
  for (int i = 0; i < p; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < q; j++)
        out(i,j) += a(k,i) * b(k,j);
}

double covMat_uncenter(NumericMatrix out,
                       const NumericMatrix a,const NumericMatrix b,
                       const int n,const int q,const int p) {
  double coeff = 1.0/double(n-1);
  double outij = 0.0;
  double maxCOV = 0.0;
  for (int i = 0; i < q; i++){
    for (int k = 0; k < n; k++){
      for (int j = 0; j < p; j++){
        outij = a(k,i) * b(k,j);
        outij *= coeff;
        out(i,j) += outij;
      }
    }
  }
  for (int i = 0; i < q; i++){
    for (int j = 0; j < p; j++){
      outij = abs(out(i,j));
      if(outij > maxCOV){
        maxCOV = outij;
      }
    }
  }
  return maxCOV;
}

void tcrossprodMat(NumericMatrix a, NumericMatrix b,NumericMatrix out, int n, int p, int q) {
  for (int i = 0; i < p; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < q; j++)
        out(i,j) += a(i,k) * b(j,k);
}

double detCpp(const NumericMatrix m, const int n) {
  NumericMatrix subM(n-1,n-1);
  double det = 0.0;
  double deti = 0.0;
  if (n == 2)
    det = (m(0,0) * m(1,1)) - (m(1,0) * m(0,1));
  else {
    for (int k = 0; k < n; k++) {
      int subi = 0;
      for (int i = 1; i < n; i++) {
        int subj = 0;
        for (int j = 0; j < n; j++) {
          if (j == k)
            continue;
          subM(subi,subj) = m(i,j);
          subj++;
        }
        subi++;
      }
      deti = detCpp(subM,n-1);
      deti *= pow(-1, k) * m(0,k);
      det += deti;
    }
  }
  return det;
}

void inverseM(const NumericMatrix M, NumericMatrix invM , const int n){
  if ( n>2 ){
    int count1=0,count2=0;
    double signCo=0.0;
    NumericMatrix subM(n-1,n-1);
    double detCurr=0.0;
    double detM = detCpp(M,n);
    double coeff = 1.0/detM;
    for (int I = 0; I < n; I++) {
      for (int J = 0; J < n; J++) {
        count1 = 0;
        count2 = 0;
        for (int ii = 0; ii < n; ii++) {
          if(ii != I){
            for (int jj = 0; jj < n; jj++) {
              if(jj != J){
                subM(count1,count2) = M(ii,jj);
                count2 += 1;
              }
            }
            count1 += 1;
            count2 = 0;
          }
        }
        detCurr = detCpp(subM,n-1);
        signCo = pow(-1.0,I+J);
        invM(J,I) = detCurr*coeff*signCo;
      }
    }
  }
  else if (n == 1) {
    invM(0,0) = 1/M(0,0);
  }
  else {
    double detM = detCpp(M,n);
    double coeff = 1.0/detM;
    invM(0,0) = coeff*M(1,1);
    invM(0,1) = (-1.0)*coeff*M(0,1);
    invM(1,1) = coeff*M(0,0);
    invM(1,0) = (-1.0)*coeff*M(1,0);
  }
}


// [[Rcpp::export]]
List do_one_componentCpp(const NumericMatrix x0,const NumericMatrix y0,
                         const NumericMatrix COV,
                         double lam,double errorMin=1e-9){
  int n = x0.nrow();
  int p = x0.ncol();
  int q = y0.ncol();
  NumericVector max_cov_y(q);
  NumericVector max_cov_x(p);
  LogicalVector id_y_high(q);
  int countNoNullY=0;
  LogicalVector id_x_high(p);
  int countNoNullX=0;
  // NumericVector roro(p);
  // NumericVector coco(q);
  // NumericVector::iterator it;

  // Get max values per column and per row
  for (int i = 0u; i < q; ++i) {
    // roro = abs(COV.row(i));
    // it = std::max_element(roro.begin(), roro.end());
    max_cov_y(i) = max(abs(COV.row(i)));//*it;
    if (max_cov_y(i)>lam){
      id_y_high(i) = true;
      countNoNullY += 1;
    }else{
      id_y_high(i) = false;
    }
  }
  for (int j = 0u; j < p; ++j) {
    // coco = abs(COV.column(j));
    // it = std::max_element(coco.begin(), coco.end());
    max_cov_x(j) = max(abs(COV.column(j)));
    if (max_cov_x(j)>lam){
      id_x_high(j) = true;
      countNoNullX += 1;
    }else{
      id_x_high(j) = false;
    }
  }
  // Get reduced covariance matrix
  NumericMatrix COV_high(countNoNullY,countNoNullX);
  NumericVector U0(p);
  NumericVector V0(q);
  NumericVector V_svd(q);
  NumericVector t(n);
  int countY = 0u;
  int countX = 0u;
  if (countNoNullY>0){
    for (int i = 0u; i < q; ++i) {
      countX = 0u;
      if (id_y_high(i)==true){
        countX = 0u;
        for (int j = 0u; j < p; ++j) {
          if (id_x_high(j)==true){
            double coefIJ = COV(i,j);
            double value = abs(coefIJ)-lam;
            if (value>0) {
              if (coefIJ>0){
                COV_high(countY,countX) = value;
              } else {
                COV_high(countY,countX) = value*(-1.0);
              }
            }
            else {
              COV_high(countY,countX) = 0;
            }
            countX += 1;
          }
        }
        countY += 1;
      }
    }
    double error = 2;
    NumericVector u0In = rnorm(countNoNullX);
    NumericVector v0In(countNoNullY);
    u0In = u0In/sqrt(sum(u0In*u0In));
    NumericVector u0In2(countNoNullX);
    while(error>errorMin){
      // for (int i = 0u; i < countNoNullY; ++i) {
      //   v0In(i) = sum(COV_high.row(i)*u0In);
      // }
      prodMatVector(COV_high,u0In,v0In,countNoNullY,countNoNullX,false);
      u0In2 = NumericVector(countNoNullX);
      crossprodMatVector(COV_high,v0In,u0In2,countNoNullX,countNoNullY,true);
      // for (int j = 0u; j < countNoNullX; ++j) {
      //   u0In2(j) = sum(COV_high.column(j)*v0In);
      // }
      // u0In2 = u0In2/sqrt(sum(u0In2*u0In2));
      error = sum((u0In2-u0In)*(u0In2-u0In));
      u0In = Rcpp::clone(u0In2);
    }
    v0In = v0In/sqrt(sum(v0In*v0In));
    countY = 0u;
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==true){
        V0(i) = v0In(countY);
        countY += 1;
      }
    }
    countX = 0u;
    for (int j = 0u; j < p; ++j) {
      if (id_x_high(j)==true){
        U0(j) = u0In(countX);
        countX += 1;
      }
    }
    // Build score
    // for (int k = 0u; k < n; ++k) {
    //   t(k) = sum(x0.row(k)*U0);
    // }
    t = NumericVector(n);
    prodMatVector(x0,U0,t,n,p,false);
    // Build y0 masked
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==true){
        V_svd(i) = sum(y0.column(i)*t);
      }
    }
    V_svd = V_svd/sqrt(sum(V_svd*V_svd));
  }
  // Generate output
  List out;
  out["t"] = t;
  out["U0"] = U0;
  out["V_svd"] = V_svd;
  out["V0"] = V0;
  return(out);
}

NumericMatrix modelddsPLSCpp(NumericMatrix U_out,const NumericMatrix x,const NumericMatrix y,
                             const int n,const int p,const int q,
                             double lam=0,int R=1,
                             double errorMin=1e-9){
  // NumericMatrix y_init(n,q);
  NumericMatrix muX(n,p);
  NumericMatrix sdXInvMat(p,q);
  NumericMatrix muY(n,q);
  NumericMatrix sdYMat(p,q);
  NumericMatrix sdYXInvMat(p,q);
  // double muCurrent;
  double RSS0=0;
  NumericVector vectHere(n);
  NumericMatrix y_plus_un(n,q);
  NumericMatrix x_plus_un(n,p);
  // NumericMatrix U_out(p,R);
  NumericMatrix U_star(p,R);
  NumericMatrix V_out(q,R);
  NumericMatrix bXr(R,p);
  NumericMatrix bYr(R,q);
  // NumericMatrix B(p,q);
  List B_r;
  NumericMatrix y_est(n,q);
  NumericVector var_expl(R);
  NumericVector t_r(n);
  NumericVector U0(p);
  NumericVector V_svd(q);
  NumericMatrix V0(q,R);
  NumericVector V0_r(q);
  NumericVector bt(p);
  NumericMatrix B_r_0(p,q);
  NumericVector deltaU(p);
  // NumericMatrix B_r_r(p,q);
  double normU02;
  double RSSr;
  NumericMatrix COV(q,p);
  // NumericVector coli(n);
  // NumericVector colj(n);
  double maxCOV=0.0;
  // double myCOV=0.0;
  // Initialize matrices
  NumericMatrix y0(Rcpp::clone(y));
  NumericMatrix x0(Rcpp::clone(x));
  // for (int k = 0u; k < n; ++k) {
  //   for (int i = 0u; i < q; ++i) {
  //     muCurrent = mean(y.column(i));
  //     muY(k,i) = muCurrent;
  //     y_init(k,i) = y(k,i)-muCurrent;
  //     y0(k,i) = y_init(k,i);
  //   }
  //   for (int j = 0u; j < p; ++j) {
  //     muCurrent = mean(x.column(j));
  //     muX(k,j) = muCurrent;
  //     x0(k,j) = x(k,j)-muCurrent;
  //   }
  // }
  // Compute initial residual sum of squares
  for (int i = 0u; i < q; ++i) {
    vectHere = y0.column(i);
    RSS0 += sum(vectHere*vectHere);
  }

  // Begin to build subspaces
  for (int r = 0u; r < R; ++r) {
    // Build empirical covariance matrix
    // for (int i = 0u; i < q; ++i) {
    //   coli = y0.column(i);
    //   for (int j = 0u; j < p; ++j) {
    //     colj = x0.column(j);
    //     myCOV = 0.0;
    //     for (int k = 0u; k < n; ++k) {
    //       myCOV += (x0(k,j)*y0(k,i))/double(n-1);
    //     }
    //     COV(i,j) = myCOV;
    //     if(abs(myCOV)>maxCOV){
    //       maxCOV = abs(myCOV);
    //     }
    //   }
    // }
    if(r>0){
      COV = NumericMatrix(q,p);
    }
    maxCOV = covMat_uncenter(COV,y0,x0,n,q,p);
    // covMat(y0, x0,COV, n, q, p);
    if(maxCOV>lam){
      List c_h = do_one_componentCpp(x0,y0,COV,lam,errorMin);
      t_r = c_h["t"];
      U0 = c_h["U0"];
      U_out(_,r) = Rcpp::clone(U0);
      V_svd = c_h["V_svd"];
      V0_r = c_h["V0"];
      V0(_,r) = Rcpp::clone(V0_r);
      // Build regression matrices
      normU02 = sum(U0*U0);
      if(normU02>errorMin){
        // for (int j = 0u; j < p; ++j){
        //   bt(j) = sum(t_r*x0.column(j))/sum(t_r*t_r);
        //   x_plus_un(_,j) = t_r*double(bt(j));
        // }
        bt = NumericVector(p);
        crossprodMatVector(x0,t_r,bt,p,n,true,true,normU02);
        prodVectVec(t_r,bt,x_plus_un,n,p);
        U_star(_,r) = U0;
        V_out(_,r) = V_svd;
        prodVectVec(U0,V_svd,B_r_0,p,q);
        prodVectVec(t_r,V_svd,y_plus_un,n,q);
        y_est += y_plus_un;
        // for (int i = 0u; i < q; ++i){
        //   // for (int j = 0u; j < p; ++j){
        //   //   B_r_0(j,i) = U0(j)*V_svd(i);
        //   // }
        //   // y_plus_un(_,i) = t_r*double(V_svd(i));
        //   // y_est(_,i) = y_est(_,i) + y_plus_un.column(i);
        // }
        bXr(r,_) = bt;

        if(r>0){
          for (int s_r = r-1; s_r > 0; --s_r){
            // for(s_r in (r-1):1){
            deltaU = U_star(_,s_r)*sum(bXr(s_r,_)*U_star(_,s_r));
            U_star(_,r) = U_star(_,r) - deltaU;
          }
        }
        // for (int i = 0u; i < q; ++i){
        //   for (int j = 0u; j < p; ++j){
        //     B_r_r(j,i) = U_star(j,r)*V_svd(i);
        //     B(j,i) += B_r_r(j,i);
        //   }
        // }
        // Computation of explained variance
        RSSr = 0;
        for (int i = 0u; i < q; ++i){
          vectHere = y.column(i)-y_plus_un.column(i);
          RSSr += sum(vectHere*vectHere);
        }
        var_expl(r) = 1-RSSr/RSS0;
        for (int k = 0u; k < n; ++k){
          y0(k,_) = y0.row(k) - y_plus_un.row(k);
          x0(k,_) = x0.row(k) - x_plus_un.row(k);
        }
      }
    }
  }
  //
  //   List out;
  //   out["U_out"] = U_out;
  //   out["U_star"] = U_star;
  //   out["V_out"] = V_out;
  //   out["V_optim"] = V0;
  //   out["B"] = B;
  //   out["B_r"] = B_r;
  //   out["var_expl"] = var_expl;
  //   out["y_est"] = y_est;
  //   out["bXr"] = bXr;
  //   out["bYr"] = bYr;
  //   out["e_x"] = x0;
  //   out["e_y"] = y0;
  //   return(out);
  return V0;
}

// [[Rcpp::export]]
NumericMatrix bootstrap_pls_CT_Cpp(const NumericMatrix X_init,const NumericMatrix Y_init,
                                   const NumericVector lambdas,const NumericVector lambda_prev,
                                   NumericMatrix uIN, NumericMatrix vIN,
                                   const int h=1, const double errorMin=1.0e-9){
  int N_lambdas = lambdas.length();
  int n = X_init.nrow();
  int p = X_init.ncol();
  int q = Y_init.ncol();
  NumericVector ids(n);
  NumericVector idIB(n);
  NumericVector idOOB;
  NumericMatrix X_train(n,p);
  NumericMatrix Y_train(n,q);
  NumericMatrix y_train_pred(n,q);
  NumericMatrix y_train_pred_next(n,q);
  int countOOB = 0;
  bool test = true;
  NumericVector colD(n);
  // NumericVector mu_k(p);
  // NumericVector sd_k(p);
  // NumericVector mu_y(q);
  // NumericVector sd_y(q);
  double mui=0.0,sdi=0.0;
  NumericMatrix u(p,h);
  NumericMatrix V_reconstruct(q,h);
  NumericMatrix t_all(n,h);
  NumericMatrix X_r(n,p);
  NumericMatrix Y_r(n,q);
  NumericMatrix Y_r_mask(n,q);
  NumericMatrix x_r(n,p);
  NumericMatrix y_r(n,q);
  NumericMatrix B_youyou(p,q);
  int r=0;
  // Build bootstrap indexes, IB and OOB
  for (int k = 0u; k < n; ++k) {
    ids(k) = k;
  }
  while (test){
    idIB = sample(ids,n,true);
    for (int k = 0u; k < n; ++k) {
      if (std::find(idIB.begin(), idIB.end(), k)==idIB.end()){
        idOOB.push_back( k );
        countOOB += 1;
      }
    }
    if(countOOB>0) {
      test = false;
    }
  }
  // Build train matrices
  for (int k = 0u; k < n; ++k) {
    X_train(k,_) = X_init.row(idIB(k));
    Y_train(k,_) = Y_init.row(idIB(k));
  }
  // Build test matrices
  NumericMatrix X_test_normalize(countOOB,p);
  NumericMatrix Y_test_normalize(countOOB,q);
  NumericMatrix y_test_pred(countOOB,q);
  NumericMatrix y_test_pred_RSS(countOOB,q);
  for (int k = 0u; k < countOOB; ++k) {
    X_test_normalize(k,_) = X_init.row(idOOB(k));
    Y_test_normalize(k,_) = Y_init.row(idOOB(k));
  }
  // Get mean and sd for Y and X on train dataset
  for (int i = 0u; i < q; ++i) {
    colD = Y_train.column(i);
    // mu_y(i) = mean(colD);
    // sd_y(i) = sd(colD);
    mui = mean(colD);
    sdi = sd(colD);
    for(int k = 0u; k < n; ++k) {
      Y_train(k,i) -= mui;
      if(sdi>0){
        Y_train(k,i) /= sdi;
      }
      Y_r(k,i) = Y_train(k,i);
    }
    for(int k = 0u; k < countOOB; ++k) {
      Y_test_normalize(k,i) -= mui;
      if(sdi>0){
        Y_test_normalize(k,i) /= sdi;
      }
    }
  }
  for (int j = 0u; j < p; ++j) {
    colD = X_train.column(j);
    mui = mean(colD);
    sdi = sd(colD);
    for(int k = 0u; k < n; ++k) {
      X_train(k,j) -= mui;
      if(sdi>0){
        X_train(k,j) /= sdi;
      }
      X_r(k,j) = X_train(k,j);
    }
    for(int k = 0u; k < countOOB; ++k) {
      X_test_normalize(k,j) -= mui;
      if(sdi>0){
        X_test_normalize(k,j) /= sdi;
      }
    }
  }
  // Build model on past components if needed
  NumericMatrix u_r(p,1);
  NumericMatrix v_r(q,1);
  NumericVector t_r(n);
  NumericVector bt(p);
  NumericMatrix C(p,h);
  NumericMatrix D(q,h);
  double norm2t=0.0;
  bool test_previous_ok = true;
  bool test_no_null=true;
  // Need to build the h-1 components for the current bootstrap dataset
  if (h>1){
    for (int s = 0u; s < h-1; ++s) {
      u(_,s) = uIN(_,s);
    }
    r=0;
    test=true;
    while (test){
      NumericMatrix onlyToInv(r+1,r+1);
      NumericMatrix InversedMat(r+1,r+1);
      NumericMatrix U_star_cl(p,r);
      NumericMatrix v_r = modelddsPLSCpp(u_r,X_r,Y_r,n,p,q,lambda_prev(r));
      // u_r = m_gogo["U_out"];
      // v_r = m_gogo["V_optim"];
      norm2t = 0.0;
      for (int j = 0u; j < p; ++j){
        u(j,r) = u_r(j,0);
      }
      // for (int k = 0u; k < n; ++k) {
      //   t_r_k = sum(X_r.row(k)*u_r(r));
      //   t_r(k) = t_r_k;
      //   norm2t += t_r_k*t_r_k;
      // }
      t_r = NumericVector(n);
      prodMatVector(X_r,u_r.column(0),t_r,n,p,false);
      norm2t = sum(t_r*t_r);
      if (norm2t>1e-9) {
        // Build X regressors and estimated matrix
        // for (int j = 0u; j < p; ++j){
        //   bt(j) = sum(t_r*X_r.column(j))/norm2t;
        //   C(j,r) = bt(j);
        //   x_r(_,j) = t_r*double(bt(j));
        // }
        bt = NumericVector(p);
        crossprodMatVector(X_r,t_r,bt,p,n,true,true,norm2t);
        C(_,r) = Rcpp::clone(bt);
        prodVectVec(t_r,bt,x_r,n,p);
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          if (abs(v_r(i,0))<errorMin){
            D(i,r) = 0.0;
          } else {
            D(i,r) = sum(Y_r.column(i)*t_r)/norm2t;
          }
          y_r(_,i) = t_r*double(D(i,r));
        }

        // Only matrix inversion of the problem
        // for (int s = 0u; s < r; ++s) {
        //   for (int t = 0u; t < r; ++t) {
        //     onlyToInv(s,t) = sum(C.column(s)*u.column(t));
        //   }
        // }
        onlyToInv = NumericMatrix(r+1,r+1);
        crossprodMat(C,u,onlyToInv,p,r+1,r+1);
        inverseM(onlyToInv,InversedMat,r+1);
        // Matrix regression building : B_youyou
        // for (int s = 0u; s < p; ++s) {
        //   for (int t = 0u; t < r; ++t) {
        //     U_star_cl(s,t) = 0.0;
        //     for (int g = 0u; g < r; g++){
        //       U_star_cl(s,t) += u(s,g)*InversedMat(g,t);
        //     }
        //   }
        // }
        U_star_cl = NumericMatrix(p,r+1);
        prodMat(u,InversedMat,U_star_cl,r,p,r+1);
        // for (int s = 0u; s < p; ++s) {
        //   for (int t = 0u; t < q; ++t) {
        //     B_youyou(s,t) = 0.0;
        //     for (int g = 0u; g < r; g++){
        //       B_youyou(s,t) += U_star_cl(s,g)*D(t,g);
        //     }
        //   }
        // }
        B_youyou = NumericMatrix(p,q);
        tcrossprodMat(U_star_cl,D,B_youyou,r+1,p,q);
        // Do deflations
        for (int k = 0u; k < n; ++k) {
          for (int s = 0u; s < p; ++s) {
            X_r(k,s) -= x_r(k,s);
          }
          for (int i = 0u; i < q; ++i) {
            Y_r(k,i) -= y_r(k,i);
          }
        }
      }
      else {
        test_no_null = false;
        test_previous_ok = false;
      }
      if (r==h-2) {
        test = false;
      }
      // Update r
      r += 1;
    }
  } else {
    norm2t = 1.0;
    test_previous_ok = true;
  }
  // Begin to test each lambda
  r = h-1;
  NumericMatrix onlyToInv(h,h);
  NumericMatrix InversedMat(h,h);
  NumericMatrix B_all(p,q);
  NumericMatrix U_star_cl(p,h);
  NumericVector diff_B(p);
  NumericMatrix u_il(p,1);//,V_il(q);
  double coeffTest = 0.0;
  NumericVector vars_expl(N_lambdas), vars_expl_h(N_lambdas), Q2(N_lambdas), Q2_all(N_lambdas);
  for (int iLam = 0u; iLam < N_lambdas; ++iLam){
    if (test_previous_ok==true) {
      NumericMatrix V_il = modelddsPLSCpp(u_il,X_r,Y_r,n,p,q,lambdas(iLam));
      // u_il = m_1["U_out"];
      u(_,r) = u_il.column(0);
      // V_il = m_1["V_optim"];
      norm2t = 0;
      // for (int k = 0u; k < n; ++k) {
      //   t_r_k = sum(X_r.row(k)*u_il);
      //   t_r(k,1) = t_r_k;
      //   norm2t += t_r_k*t_r_k;
      // }
      t_r = NumericVector(n);
      prodMatVector(X_r,u_il.column(0),t_r,n,p,false);
      norm2t = sum(t_r*t_r);
      if (norm2t>1e-9) {
        // Build X regressors and estimated matrix
        // for (int j = 0u; j < p; ++j){
        //   bt(j,1) = sum(t_r(0)*X_r.column(j))/norm2t;
        //   C(j,r) = bt(j,1);
        //   // x_r(_,j) = t_r*double(bt(j));
        // }
        bt = NumericVector(p);
        crossprodMatVector(X_r,t_r,bt,p,n,true,true,norm2t);
        C(_,r) = bt;
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          coeffTest = abs(V_il(i,0));
          if (coeffTest<1.0e-9){
            D(i,r) = 0.0;
          } else {
            D(i,r) = sum(Y_r.column(i)*t_r)/norm2t;
          }
          // y_r(_,i) = t_r*double(D(i,r));
        }
        // Only matrix inversion of the problem
        // for (int s = 0u; s < h; ++s) {
        //   for (int t = 0u; t < h; ++t) {
        //     onlyToInv(s,t) = sum(C.column(s)*u.column(t));
        //   }
        // }
        // inverseM(onlyToInv,InversedMat,h);
        // // Matrix regression building : B_youyou
        // for (int s = 0u; s < p; ++s) {
        //   for (int t = 0u; t < h; ++t) {
        //     U_star_cl(s,t) = 0;
        //     for (int g = 0u; g < h; g++){
        //       U_star_cl(s,t) += u(s,g)*InversedMat(g,t);
        //     }
        //   }
        // }
        // for (int s = 0u; s < p; ++s) {
        //   for (int t = 0u; t < q; ++t) {
        //     B_all(s,t) = 0;
        //     for (int g = 0u; g < h; g++){
        //       B_all(s,t) += U_star_cl(s,g)*D(t,g);
        //     }
        //   }
        // }
        onlyToInv = NumericMatrix(h,h);
        crossprodMat(C,u,onlyToInv,p,h,h);
        inverseM(onlyToInv,InversedMat,h);
        U_star_cl = NumericMatrix(p,h);
        prodMat(u,InversedMat,U_star_cl,h,p,h);
        B_all = NumericMatrix(p,q);
        tcrossprodMat(U_star_cl,D,B_all,h,p,q);
      }
    }
    // Create the prediction matrices and metrics
    double err=0.0;
    double n_t_p_i=0.0,d_t_i=0.0;
    double n_t_p_n_i=0.0;
    double n_Q2_i=0.0,d_Q2_i=0.0;
    double d_Q2_a_i=0.0;
    for (int i = 0u; i < q; ++i) {
      // n_t_p_i=0.0, d_t_i=0.0, n_t_p_n_i=0.0;
      for (int k = 0u; k < n; ++k) {
        y_train_pred(k,i) = sum(X_train.row(k)*B_all.column(i));
        err = y_train_pred(k,i)-Y_train(k,i);
        // numer_train_pred += err*err;
        n_t_p_i += err*err;
        diff_B = B_all.column(i)-B_youyou.column(i);
        y_train_pred_next(k,i) = sum(X_train.row(k)*diff_B);
        err = y_train_pred_next(k,i)-Y_train(k,i);
        // numer_train_pred_next += err*err;
        n_t_p_n_i += err*err;
        err = Y_train(k,i);
        // denom_train += err*err;
        d_t_i += err*err;
      }
      // vars_expl(iLam) += (1-n_t_p_i/d_t_i)/double(q);
      // vars_expl_h(iLam) += (1-n_t_p_n_i/d_t_i)/double(q);
      // n_Q2_i=0.0, d_Q2_i=0.0; d_Q2_a_i=0.0;
      for (int k = 0u; k < countOOB; ++k) {
        y_test_pred(k,i) = sum(X_test_normalize.row(k)*B_all.column(i));
        err = y_test_pred(k,i)-Y_test_normalize(k,i);
        // numer_Q2 += err*err;
        n_Q2_i += err*err;
        err = y_test_pred_RSS(k,i)-Y_test_normalize(k,i);
        // denom_Q2 += err*err;
        d_Q2_i += err*err;
        y_test_pred_RSS(k,i) = sum(X_test_normalize.row(k)*B_youyou.column(i));
        err = Y_test_normalize(k,i);
        // denom_Q2_all += err*err;
        d_Q2_a_i += err*err;
      }
      // Q2(iLam) += (1-n_Q2_i/d_Q2_i)/double(q);
      // Q2_all(iLam) += (1-n_Q2_i/d_Q2_a_i)/double(q);
    }
    vars_expl(iLam) = 1.0-n_t_p_i/d_t_i;
    vars_expl_h(iLam) = 1.0-n_t_p_n_i/d_t_i;
    Q2(iLam) = 1.0-n_Q2_i/d_Q2_i;
    Q2_all(iLam) = 1.0-n_Q2_i/d_Q2_a_i;
  }
  // Prepare output
  // List out;
  // out["vars_expl"] = vars_expl;
  // out["vars_expl_h"] = vars_expl_h;
  // out["Q2"] = Q2;
  // out["Q2_all"] = Q2_all;
  // out["idOOB"] = idOOB;
  NumericMatrix out(4,N_lambdas);
  out(0,_) = vars_expl; // R^2
  out(1,_) = vars_expl_h; // R^2_h
  out(2,_) = Q2_all; // Q^2
  out(3,_) = Q2; // Q^2_h
  return out;
}
