// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

// ===========================
// nllfun_cpp
// ===========================

// [[Rcpp::export]]
double nllfun_cpp(
    const Eigen::VectorXd& obs_E,
    const Eigen::MatrixXd& obs_V,
    const Eigen::VectorXd& pred_E,
    const Eigen::MatrixXd& pred_V,
    int n = 1
) {
  Eigen::MatrixXd inv_pred_V = pred_V.inverse();
  Eigen::VectorXd resids = obs_E - pred_E;
  double log_det = std::log(pred_V.determinant());
  double trace_term = (obs_V * inv_pred_V).trace();
  double quad_term = resids.transpose() * inv_pred_V * resids;
  return 0.5 * n * (log_det + trace_term + quad_term);
}

// ===========================
// nllfun_var_cpp
// ===========================

// [[Rcpp::export]]
double nllfun_var_cpp(
    const Eigen::VectorXd& obs_mean,
    const Eigen::VectorXd& obs_var,
    const Eigen::VectorXd& pred_mean,
    const Eigen::VectorXd& pred_var,
    int n = 1
) {
  Eigen::VectorXd ll_components =
    (obs_var.array() / pred_var.array()) +
    ((obs_mean - pred_mean).array().square() / pred_var.array()) +
    pred_var.array().log();
  double total_ll = -n * ll_components.sum();
  return -total_ll;
}

// ===========================
// logdens2wt_cpp
// ===========================

// [[Rcpp::export]]
Eigen::VectorXd logdens2wt_cpp(const Eigen::VectorXd& x) {
  double max_x = x.maxCoeff();
  Eigen::VectorXd x2 = (x.array() - max_x).exp();

  double max_finite = 0.0;
  for (int i = 0; i < x2.size(); i++)
    if (std::isfinite(x2[i]) && x2[i] > max_finite)
      max_finite = x2[i];

    for (int i = 0; i < x2.size(); i++)
      if (!std::isfinite(x2[i]))
        x2[i] = max_finite;

      double sum_x2 = x2.sum();
      if (sum_x2 > 0) x2 /= sum_x2;

      return x2;
}

// ===========================
// compute_variance_cpp
// ===========================

// [[Rcpp::export]]
Eigen::MatrixXd compute_variance_cpp(
    const Eigen::MatrixXd& d_f_d_bi,
    const Eigen::MatrixXd& Omega
) {
  return d_f_d_bi * Omega * d_f_d_bi.transpose();
}

// ===========================
// gen_bi2_cpp
// ===========================

// [[Rcpp::export]]
Eigen::MatrixXd gen_bi2_cpp(
    const Eigen::MatrixXd& Omega,
    const Eigen::MatrixXd& biseq
) {
  Eigen::MatrixXd Omega_adj = Omega;
  for (int i = 0; i < Omega_adj.rows(); i++)
    if (std::abs(Omega_adj(i, i)) < 1e-10)
      Omega_adj(i, i) = 1e-10;

    Eigen::LLT<Eigen::MatrixXd> llt(Omega_adj);
    if (llt.info() != Eigen::Success)
      stop("Cholesky decomposition failed");

    Eigen::MatrixXd L = llt.matrixL();
    Eigen::MatrixXd U = L.transpose(); // Upper like R
    return biseq * U;
}

// ===========================
// meancov_cpp
// ===========================

// [[Rcpp::export]]
List meancov_cpp(const Eigen::MatrixXd& M, const Eigen::VectorXd& w) {
  double ws = w.sum();
  Eigen::VectorXd mean = (M.array().colwise() * w.array()).colwise().sum() / ws;
  Eigen::MatrixXd centered = M.rowwise() - mean.transpose();
  Eigen::MatrixXd cov =
    (centered.array().colwise() * w.array()).matrix().transpose() * centered / ws;
  return List::create(Named("E") = mean, Named("V") = cov);
}

// ===========================
// gen_bi_batch_cpp
// ===========================

// [[Rcpp::export]]
Eigen::MatrixXd gen_bi_batch_cpp(const Eigen::MatrixXd& Omega, int n_samples) {
  Eigen::LLT<Eigen::MatrixXd> llt(Omega);
  Eigen::MatrixXd L = llt.matrixL();
  int d = Omega.rows();
  Eigen::MatrixXd out(n_samples, d);

  for (int i = 0; i < n_samples; ++i) {
    Eigen::VectorXd z = Eigen::VectorXd::NullaryExpr(d, []() { return R::rnorm(0, 1); });
    out.row(i) = (L * z).transpose();
  }
  return out;
}

// ===========================
// jacobiann_vec_fast_cpp
// ===========================

// [[Rcpp::export]]
NumericMatrix jacobiann_vec_fast_cpp(Function f, NumericVector bi, double eps = 1e-8) {
  int n = bi.size();
  NumericVector f0 = f(bi);
  int m = f0.size();
  NumericMatrix J(m, n);
  NumericVector bi_pert = clone(bi);

  for (int j = 0; j < n; ++j) {
    bi_pert = clone(bi);
    bi_pert[j] += eps;
    NumericVector f1 = f(bi_pert);
    for (int i = 0; i < m; ++i)
      J(i, j) = (f1[i] - f0[i]) / eps;
  }
  return J;
}

// ===========================
// MCapprEV_cpp
// ===========================

// [[Rcpp::export]]
List MCapprEV_cpp(NumericMatrix m, NumericVector logwt) {
  int N = m.nrow();
  int p = m.ncol();
  NumericVector wt = exp(logwt - max(logwt)); // stabilize
  double sumw = sum(wt);

  NumericVector mean(p);
  for (int j = 0; j < p; ++j)
    mean[j] = sum(m(_, j) * wt) / sumw;

  NumericMatrix cov(p, p);
  for (int i = 0; i < N; ++i) {
    NumericVector diff = m(i, _) - mean;
    for (int a = 0; a < p; ++a)
      for (int b = 0; b < p; ++b)
        cov(a, b) += wt[i] * diff[a] * diff[b];
  }
  cov = cov / sumw;

  return List::create(_["E"] = mean, _["V"] = cov);
}

// ===========================
// samplogdensfun_cpp
// ===========================

// [[Rcpp::export]]
Eigen::VectorXd samplogdensfun_cpp(
    const Eigen::MatrixXd& bi,
    const Eigen::MatrixXd& Omega,
    double omega_expansion
) {
  int N = bi.rows();
  int d = bi.cols();
  Eigen::VectorXd out(N);

  Eigen::MatrixXd OmExp = Omega * omega_expansion;

  Eigen::LLT<Eigen::MatrixXd> lltOm(Omega);
  Eigen::LLT<Eigen::MatrixXd> lltOmExp(OmExp);
  if (lltOm.info() != Eigen::Success || lltOmExp.info() != Eigen::Success)
    stop("Cholesky decomposition failed");

  Eigen::MatrixXd L = lltOm.matrixL();
  Eigen::MatrixXd Lexp = lltOmExp.matrixL();

  double logdetOm = 2.0 * L.diagonal().array().log().sum();
  double logdetOmExp = 2.0 * Lexp.diagonal().array().log().sum();
  double normConst = -0.5 * d * std::log(2 * M_PI);

  for (int i = 0; i < N; ++i) {
    Eigen::VectorXd x = bi.row(i);
    double qOm = (L.triangularView<Lower>().solve(x)).squaredNorm();
    double qOmExp = (Lexp.triangularView<Lower>().solve(x)).squaredNorm();
    out[i] = -0.5 * (logdetOm + qOm) + 0.5 * (logdetOmExp + qOmExp);
  }
  return out;
}

// ===========================
// g_iter_generic_cpp
// ===========================

// [[Rcpp::export]]
NumericMatrix g_iter_generic_cpp(Function gfun,
                                 NumericVector beta,
                                 NumericMatrix bi,
                                 Nullable<NumericVector> ai = R_NilValue) {
  int n = bi.nrow();
  int p = bi.ncol();
  NumericMatrix theta_i(n, p);

  NumericVector ai_vec;
  if (ai.isNull()) {
    ai_vec = NumericVector(n, 1.0);  // default ai = 1 for all
  } else {
    ai_vec = as<NumericVector>(ai);
    if (ai_vec.size() == 0)
      ai_vec = NumericVector(n, 1.0);
    else if (ai_vec.size() == 1)
      ai_vec = NumericVector(n, ai_vec[0]);
  }

  for (int i = 0; i < n; ++i) {
    NumericVector bi_row = bi(i, _);
    NumericVector res = gfun(beta, bi_row, ai_vec[i]);
    for (int j = 0; j < p; ++j)
      theta_i(i, j) = res[j];
  }

  return theta_i;
}

