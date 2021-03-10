#include "iostream"
#include "Eigen/Dense"
#include "Eigen/SVD"
#include <random>
 
using namespace Eigen;
using namespace std;

class oneComponent {
  public:
    VectorXd t;
    VectorXd U0;
    VectorXd V0;
    VectorXd V_svd;
};

class multiComponent {
public:
    MatrixXd U_out;
    MatrixXd V0;
};


class ddsPLSCpp {
  public:
    VectorXd R2;
    VectorXd R2h;
    VectorXd Q2;
    VectorXd Q2h;
    VectorXd idIB;
    VectorXd idOOB;
};

oneComponent do_one_componentCpp(const MatrixXd x0,const MatrixXd y0,
                         const MatrixXd COV,
                         const int n,const int p,const int q,
                         double lam,double errorMin=1e-9){
  VectorXd max_cov_y = VectorXd::Zero(q);
  VectorXd max_cov_x = VectorXd::Zero(p);
  VectorXi id_y_high = VectorXi::Zero(q);
  int countNoNullY=0;
  VectorXi id_x_high = VectorXi::Zero(p);
  int countNoNullX=0;
  // Get max values per column and per row
  for (int i = 0u; i < q; ++i) {
    max_cov_y(i) = COV.block(i,0,1,p).array().abs().maxCoeff();
    if (max_cov_y(i)>lam){
      id_y_high(i) = 1;
      countNoNullY += 1;
    }
  }
  for (int j = 0u; j < p; ++j) {
    max_cov_x(j) = COV.block(0,j,q,1).array().abs().maxCoeff();
    if (max_cov_x(j)>lam){
      id_x_high(j) = 1;
      countNoNullX += 1;
    }
  }
  // Get reduced covariance matrix
  MatrixXd COV_high = MatrixXd::Zero(countNoNullY,countNoNullX);
  VectorXd U0 = VectorXd::Zero(p);
  VectorXd V0 = VectorXd::Zero(q);
  VectorXd V_svd = VectorXd::Zero(q);
  VectorXd t = VectorXd::Zero(n);
  int countY = 0u;
  int countX = 0u;
  if (countNoNullY>0){
    for (int i = 0u; i < q; ++i) {
      countX = 0u;
      if (id_y_high(i)==1){
        countX = 0u;
        for (int j = 0u; j < p; ++j) {
          if (id_x_high(j)==1){
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
    double error = 2.0;
    VectorXd u0In;
    u0In.setRandom(countNoNullX);
    VectorXd v0In(countNoNullY);
    u0In = u0In/sqrt(u0In.squaredNorm());
    VectorXd u0In2(countNoNullX);
    while(error>errorMin){
      v0In = COV_high*u0In;
      u0In2 = COV_high.transpose()*v0In;
      u0In2 /= sqrt(u0In2.squaredNorm());
      error = (u0In2-u0In).squaredNorm();
      u0In = u0In2;
    }
    v0In  /= sqrt(v0In.squaredNorm());
    countY = 0u;
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==1){
        V0(i) = v0In(countY);
        countY += 1;
      }
    }
    countX = 0u;
    for (int j = 0u; j < p; ++j) {
      if (id_x_high(j)==1){
        U0(j) = u0In(countX);
        countX += 1;
      }
    }
    // Build score
    t = x0*U0;
    // Build y0 masked
    for (int i = 0u; i < q; ++i) {
      if (id_y_high(i)==1){
        V_svd(i) = (y0.block(0,i,n,1).transpose()*t).sum();
      }
    }
    V_svd /= sqrt(V_svd.squaredNorm());
  }
  // Generate output
  oneComponent out;
  out.t = t;
  out.U0 = U0;
  out.V_svd = V_svd;
  out.V0 = V0;
  return(out);
}


multiComponent modelddsPLSCpp(MatrixXd U_out,MatrixXd V0,const MatrixXd x,
  const MatrixXd y,const int n,const int p,const int q,double lam=0.0,int R=1,
  double errorMin=1e-9){
  double RSS0=0.0;
  double normU02;
  double RSSr;
  double maxCOV=0.0;
  MatrixXd muX = MatrixXd::Zero(n,p);
  MatrixXd sdXInvMat = MatrixXd::Zero(p,q);
  MatrixXd muY = MatrixXd::Zero(n,q);
  MatrixXd sdYMat = MatrixXd::Zero(p,q);
  MatrixXd sdYXInvMat = MatrixXd::Zero(p,q);
  MatrixXd y_plus_un = MatrixXd::Zero(n,q);
  MatrixXd x_plus_un = MatrixXd::Zero(n,p);
  MatrixXd U_star = MatrixXd::Zero(p,R);
  MatrixXd V_out = MatrixXd::Zero(q,R);
  MatrixXd bXr = MatrixXd::Zero(R,p);
  MatrixXd bYr = MatrixXd::Zero(R,q);
  MatrixXd y_est = MatrixXd::Zero(n,q);
  MatrixXd B_r_0 = MatrixXd::Zero(p,q);
  MatrixXd y0 = y;
  MatrixXd x0 = x;
  MatrixXd COV = MatrixXd::Zero(q,p);
  VectorXd vectHere = VectorXd::Zero(n);
  VectorXd var_expl = VectorXd::Zero(R);
  VectorXd t_r = VectorXd::Zero(n);
  VectorXd U0 = VectorXd::Zero(p);
  VectorXd V_svd = VectorXd::Zero(q);
  VectorXd V0_r = VectorXd::Zero(q);
  VectorXd bt = VectorXd::Zero(p);
  VectorXd deltaU = VectorXd::Zero(p);
  // Compute initial residual sum of squares
  RSS0 = y0.squaredNorm();
  // Begin to build subspaces
  for (int r = 0u; r < R; ++r) {
    // Build empirical covariance matrix
    COV = y0.transpose()*x0/double(n-1.0);
    maxCOV = COV.array().abs().maxCoeff();
    if(maxCOV>lam){
      oneComponent c_h = do_one_componentCpp(x0,y0,COV,n,p,q,lam,errorMin);
      t_r = c_h.t;
      U0 = c_h.U0;
      U_out.block(0,r,p,1) = U0;
      V_svd = c_h.V_svd;
      V0_r = c_h.V0;
      V0.block(0,r,q,1) = V0_r;
      // Build regression matrices
      normU02 = U0.squaredNorm();
      if(normU02>errorMin){
        bt = x0.transpose()*t_r/normU02;
        x_plus_un = t_r*bt.transpose();
        U_star.block(0,r,p,1) = U0;
        V_out.block(0,r,q,1) = V_svd;
        B_r_0 = U0*V_svd.transpose();
        y_plus_un = t_r*V_svd.transpose();
        y_est += y_plus_un;
        bXr.block(r,0,1,p) = bt.transpose();
        if(r>0){
          for (int s_r = r-1; s_r > 0; --s_r){
            deltaU = U_star.block(0,s_r,p,1)*(bXr.block(s_r,0,1,p)*U_star.block(0,s_r,p,1)).sum();
            U_star.block(0,r,p,1) = U_star.block(0,r,p,1) - deltaU;
          }
        }
        // Computation of explained variance
        RSSr = (y-y_plus_un).squaredNorm();
        var_expl(r) = 1.0-RSSr/RSS0;
        y0 -= y_plus_un;
        x0 -= x_plus_un;
      }
    }
  }
  multiComponent out;
  out.U_out = U_out;
  out.V0 = V0;
  return out;
}

ddsPLSCpp bootstrap_pls_CT_Cpp(const MatrixXd X_init,const MatrixXd Y_init,
                                   const VectorXd lambdas,const VectorXd lambda_prev,
                                   MatrixXd uIN, MatrixXd vIN,
                                   const int h=1, const double errorMin=1.0e-9){
  int N_lambdas = lambdas.size();
  int n = X_init.rows();
  int p = X_init.cols();
  int q = Y_init.cols();
  VectorXd idIB = VectorXd::Zero(n);
  VectorXd idIBChosen = VectorXd::Zero(n);
  VectorXd idOOBChosen = VectorXd::Zero(n);
  MatrixXd X_train = MatrixXd::Zero(n,p);
  MatrixXd Y_train = MatrixXd::Zero(n,q);
  MatrixXd y_train_pred = MatrixXd::Zero(n,q);
  MatrixXd y_train_pred_next = MatrixXd::Zero(n,q);
  int countOOB = 0;
  int nbOOB=0;
  bool test = true;
  VectorXd colD = VectorXd::Zero(n);
  VectorXd MU = VectorXd::Zero(n);
  double mui=0.0,sdi=0.0;
  MatrixXd u = MatrixXd::Zero(p,h);
  MatrixXd V_reconstruct = MatrixXd::Zero(q,h);
  MatrixXd t_all = MatrixXd::Zero(n,h);
  MatrixXd X_r = MatrixXd::Zero(n,p);
  MatrixXd Y_r = MatrixXd::Zero(n,q);
  MatrixXd Y_r_mask = MatrixXd::Zero(n,q);
  MatrixXd x_r = MatrixXd::Zero(n,p);
  MatrixXd y_r = MatrixXd::Zero(n,q);
  MatrixXd B_youyou = MatrixXd::Zero(p,q);
  int r=0;
  // Build bootstrap indexes, IB and OOB
  while (test){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,n-1);
    for (int i=0; i<n; ++i) {
      int number = distribution(generator);
      ++idIB[number];
      idIBChosen[number] = 1;
      idOOBChosen[number] = 0;
    }
    nbOOB = n - idIBChosen.sum();
    if(nbOOB>0){
      test = false;
    }
  }
  VectorXd idOOB(nbOOB);
  countOOB = 0;
  for (int k = 0u; k < n; ++k) {
    if (idIB(k)==1){
      idOOB(countOOB) = k;
      countOOB += 1;
    }
  }
  // Build train matrices
  for (int k = 0u; k < n; ++k) {
    X_train.block(k,0,1,p) = X_init.block(idIB(k),0,1,p);
    Y_train.block(k,0,1,q) = Y_init.block(idIB(k),0,1,q);
  }
  // Build test matrices
  MatrixXd X_test_normalize(countOOB,p);
  MatrixXd Y_test_normalize(countOOB,q);
  MatrixXd y_test_pred(countOOB,q);
  MatrixXd y_test_pred_RSS(countOOB,q);
  for (int k = 0u; k < countOOB; ++k) {
    X_test_normalize.block(k,0,1,p) = X_init.block(idOOB(k),0,1,p);
    Y_test_normalize.block(k,0,1,q) = Y_init.block(idOOB(k),0,1,q);
  }
  // Get mean and sd for Y and X on train dataset
  for (int i = 0u; i < q; ++i) {
    colD = Y_train.block(0,i,n,1);
    // mu_y(i) = mean(colD);
    // sd_y(i) = sd(colD);
    mui = colD.mean();
    for(int k = 0; k < n; k++){
      MU(k)     = mui;
    }
    sdi = sqrt((colD - MU).squaredNorm()/(n-1));
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
    colD = X_train.block(0,j,n,1);
    mui = colD.mean();
    for(int k = 0; k < n; k++){
      MU(k)     = mui;
    }
    sdi = sqrt((colD - MU).squaredNorm()/(n-1));
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
  MatrixXd u_r = MatrixXd::Zero(p,1);
  MatrixXd v_r = MatrixXd::Zero(q,1);
  VectorXd t_r = VectorXd::Zero(n);
  VectorXd bt = VectorXd::Zero(p);
  MatrixXd C = MatrixXd::Zero(p,h);
  MatrixXd D = MatrixXd::Zero(q,h);
  double norm2t=0.0;
  bool test_previous_ok = true;
  bool test_no_null=true;
  // Need to build the h-1 components for the current bootstrap dataset
  if (h>1){
    for (int s = 0u; s < h-1; ++s) {
      u.block(0,s,p,1) = uIN.block(0,s,p,1);
    }
    r=0;
    test=true;
    while (test){
      MatrixXd onlyToInv(r+1,r+1);
      MatrixXd InversedMat(r+1,r+1);
      MatrixXd U_star_cl(p,r);
      multiComponent res1 = modelddsPLSCpp(u_r,v_r,X_r,Y_r,n,p,q,lambda_prev(r),r);
      v_r = res1.V0;
      u_r = res1.U_out;
      for (int j = 0u; j < p; ++j){
        u(j,r) = u_r(j,0);
      }
      t_r = X_r*u_r;
      norm2t = (t_r.transpose()*t_r).sum();
      if (norm2t>errorMin) {
        // Build X regressors and estimated matrix
        bt = VectorXd(p);
        bt = X_r.transpose()*t_r/norm2t;
        C.block(0,r,p,1) = bt;
        x_r = t_r*bt.transpose();
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          if (abs(v_r(i,0))<errorMin){
            D(i,r) = 0.0;
          } else {
            D(i,r) = (Y_r.block(0,i,n,1).transpose()*t_r).sum()/norm2t;
          }
          y_r.block(0,i,n,1) = t_r*D(i,r);
        }
        // Only matrix inversion of the problem
        onlyToInv = C.transpose()*u;
        InversedMat = onlyToInv.inverse();
        U_star_cl = u*InversedMat;
        B_youyou = U_star_cl*D.transpose();
        X_r -= x_r;
        Y_r -= y_r;
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
  MatrixXd onlyToInv(h,h);
  MatrixXd InversedMat(h,h);
  MatrixXd B_all(p,q);
  MatrixXd U_star_cl(p,h);
  VectorXd diff_B(p);
  MatrixXd u_il(p,1);
  MatrixXd V_il(q,1);
  double coeffTest = 0.0;
  VectorXd vars_expl(N_lambdas), vars_expl_h(N_lambdas), Q2(N_lambdas), Q2_all(N_lambdas);
  for (int iLam = 0u; iLam < N_lambdas; ++iLam){
    if (test_previous_ok==true) {
      multiComponent res1 = modelddsPLSCpp(u_il,V_il,X_r,Y_r,n,p,q,lambdas(iLam),r);
      V_il = res1.V0;
      u_il = res1.U_out;
      u.block(0,r,p,1) = u_il.block(0,r,p,1);
      t_r = X_r*u_il.block(0,r,p,1);
      norm2t = (t_r.transpose()*t_r).sum();
      if (norm2t>errorMin) {
        // Build X regressors and estimated matrix
        bt = X_r.transpose()*t_r/norm2t;
        C.block(0,r,p,1) = bt;
        x_r = t_r*bt.transpose();
        // Build Y regressors and estimated matrix
        for (int i = 0u; i < q; ++i) {
          coeffTest = abs(V_il(i,0));
          if (coeffTest<errorMin){
            D(i,r) = 0.0;
          } else {
            D(i,r) = (Y_r.block(0,i,n,1).transpose()*t_r).sum()/norm2t;
          }
        }
        onlyToInv = C.transpose()*u;
        InversedMat = onlyToInv.inverse();
        U_star_cl = u*InversedMat;
        B_all = U_star_cl*D.transpose();
      }
    }
    // Create the prediction matrices and metrics
    MatrixXd diff_B = B_all - B_youyou;
    MatrixXd y_train_pred_next = X_train*diff_B;
    MatrixXd y_test_pred = X_test_normalize*B_all;
    MatrixXd y_test_pred_RSS = X_test_normalize*B_youyou;
    double n_t_p_i = (y_train_pred-Y_train).squaredNorm();
    double n_t_p_n_i = (y_train_pred_next-Y_train).squaredNorm();
    double d_t_i = (Y_train).squaredNorm();
    vars_expl(iLam) = 1.0-n_t_p_i/d_t_i;
    vars_expl_h(iLam) = 1.0-n_t_p_n_i/d_t_i;
    double n_Q2_i = (y_test_pred-Y_test_normalize).squaredNorm();
    double d_Q2_i = (y_test_pred_RSS-Y_test_normalize).squaredNorm();
    double d_Q2_a_i = (Y_test_normalize).squaredNorm();
    Q2(iLam) = 1.0-n_Q2_i/d_Q2_i;
    Q2_all(iLam) = 1.0-n_Q2_i/d_Q2_a_i;
  }
  ddsPLSCpp out;
  out.R2 = vars_expl; // R^2
  out.R2h = vars_expl_h; // R^2_h
  out.Q2 = Q2_all; // Q^2
  out.Q2h = Q2; // Q^2_h
  return out;
}
 
int main()
{
int n = 50;
    int p = 5;
    int q = 4;
    MatrixXd x;
    x.setRandom(n,p);
    MatrixXd y;
    y.setRandom(n,q);
    MatrixXd COV(q,p);
    COV = y.transpose()*x/double(n-1.0);
    double lam=0.01;
    int R=1;
    MatrixXd U_out;
    U_out.setRandom(p,1);
    U_out /= sqrt(U_out.squaredNorm());
    MatrixXd V0;
    V0.setRandom(q,1);
    V0 /= sqrt(V0.squaredNorm());
    // ddsPLSCpp res = modelddsPLSCpp( U_out, V0, x,y, n, p, q,lam,R);
    Vector3d lambdas(0.0,0.1,0.2);
    Vector2d lambda_prev(0.0,0.0);
    ddsPLSCpp res = bootstrap_pls_CT_Cpp(x,y,lambdas,lambda_prev,U_out,V0,2);

    cout << "------------------------" << endl;
    cout << res.R2 << endl;
    cout << "------------------------" << endl;
    cout << res.R2h << endl;
    cout << "------------------------" << endl;
    cout << res.Q2 << endl;
    cout << "------------------------" << endl;
    cout << res.Q2h << endl;


// cout << res.t << endl;
// cout << "------------------------" << endl;
// cout << res.U_out << endl;
// cout << "------------------------" << endl;
// cout << res.V_svd << endl;
// cout << "------------------------" << endl;
// cout << res.V0 << endl;

//     MatrixXd bt;
//     bt = a*b.transpose();
//     cout << bt << endl;
//     MatrixXd out2 = MatrixXd::Zero(p,3);
// out2.block(0,0,p,1) = a;
// cout << a << endl;
// cout << out2 << endl;
// cout << "____________________________________" << endl;
//     MatrixXd xxx;
//     xxx.setRandom(n,p);
// cout << xxx.block(2,0,1,p).array() << endl;
// cout << xxx.block(2,0,1,p).array().abs().maxCoeff() << endl;
// cout << "____________________________________" << endl;
// double lam = 2.0;
// cout << abs(xxx(0,0))-lam << endl;
// cout << "____________________________________" << endl;

//     VectorXd u0In;
//     u0In.setRandom(p);
//     cout << u0In << endl;
//     u0In /= sqrt(u0In.squaredNorm());
// cout << "____________________________________" << endl;
//     cout << u0In << endl;


   // VectorXd DV = VectorXd::Zero(30);
 //cout << DV << " \n";   

    return 0;
}