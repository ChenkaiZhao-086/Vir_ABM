#include <Rcpp.h>
#include <cmath>
#include <string>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp17)]]

// Linear interpolation for daily series:
// x[0] corresponds to day = start_day, x[1] to start_day+1, ...
static inline double interp_daily(const NumericVector &x, double time,
                                  double start_day) {
  const R_xlen_t L = x.size();
  if (L <= 0)
    return 0.0;
  if (L == 1)
    return x[0];

  double u = time - start_day; // "day index" from 0
  if (u <= 0.0)
    return x[0];

  double max_u = static_cast<double>(L - 1);
  if (u >= max_u)
    return x[L - 1];

  int k = static_cast<int>(std::floor(u));
  double frac = u - static_cast<double>(k);

  return (1.0 - frac) * x[k] + frac * x[k + 1];
}

// [[Rcpp::export]]
List ParmInferenceCpp(double time, const NumericVector &state,
                      const List &parms) {
  const double N = as<double>(parms["num_of_agent"]);
  const double beta_seasonal = as<double>(parms["beta_seasonal"]);
  const double phi = as<double>(parms["phi"]);
  const double beta_amplify = as<double>(parms["beta_amplify"]);
  const double added_cases = as<double>(parms["added_cases"]);
  const double Penalty = as<double>(parms["Penalty"]);
  const double Restrain = as<double>(parms["Restrain"]);
  const NumericVector PartialProtection = parms["PartialProtection"];

  const NumericVector beta0 = parms["beta0"];
  const NumericVector gamma = parms["gamma"];
  const NumericVector omega = parms["omega"];
  const NumericVector comp = parms["comp"];
  const NumericVector I_ref = parms["I_ref"];

  const int n = beta0.size(); // number of viruses

  // ---- NPI
  const bool use_npi =
      parms.containsElementNamed("NPI") ? as<bool>(parms["NPI"]) : false;

  double npi_intensity = 0.0;
  NumericVector NPIses(n);

  if (use_npi) {
    NPIses = parms["NPISes"];
    const NumericVector NPI_value = parms["NPI_value"];
    const double NPI_start = as<double>(parms["NPI_start"]);
    const double NPI_end = as<double>(parms["NPI_end"]);

    if (time >= NPI_start && time <= NPI_end) {
      npi_intensity = interp_daily(NPI_value, time, NPI_start);
    }
  }

  // ---- Seasonal force
  const double seasonal_force =
      1.0 + beta_seasonal *
                cos((2.0 * M_PI * time / 365.0) - (2.0 * M_PI / 365.0) * phi);

  // ---- Precompute sums for competition logic:
  // stress_i = sum_{j != i} ((I_j + X_j) / I_ref_j) * (comp_j - comp_i)
  //          = A - comp_i * B
  // where
  //   A = sum_j ((I_j + X_j) / I_ref_j) * comp_j
  //   B = sum_j ((I_j + X_j) / I_ref_j)
  double A = 0.0;
  double B = 0.0;

  // For Penalty: only infectious I from other viruses
  double sumI_j = 0.0;
  for (int j = 0; j < n; ++j) {
    const int idxS = 4 * j;
    const int idxI = idxS + 1;
    const int idxX = idxS + 2;

    const double Ij = state[idxI];
    const double Xj = state[idxX];

    const double wj = (Ij + Xj) / I_ref[j];
    A += wj * comp[j];
    B += wj;

    sumI_j += Ij + Xj;
  }

  // ---- Derivatives + incidence output
  NumericVector dydt(4 * n);
  NumericVector out(n);
  CharacterVector out_names(n);

  for (int i = 0; i < n; ++i) {
    const int idxS = 4 * i;
    const int idxI = idxS + 1;
    const int idxX = idxS + 2;
    const int idxR = idxS + 3;

    const double Si = state[idxS];
    const double Ii = state[idxI];
    const double Xi = state[idxX];
    const double Ri = state[idxR];
    const double pRi = PartialProtection[i];

    // Mechanism 1: virus competition
    const double stress_i = A - comp[i] * B;
    const double S_eff = Si * exp(-stress_i);

    // beta with seasonality and NPI
    double beta_i = beta0[i] * seasonal_force;
    if (use_npi) {
      beta_i *= (1.0 - NPIses[i] * npi_intensity);
      // NPI: beta * (1 - sens * intensity)
    }

    // lambda_i(t) = beta_i * I_i / N
    const double lambda_i = beta_i * Ii * beta_amplify / N;

    // Mechanism 2: coinfection suppression (uniform)
    double I_other_prev = (sumI_j - Ii - Xi) / N;
    if (I_other_prev < 0.0)
      I_other_prev = 0.0;
    const double coinf_factor = std::exp(-Penalty * I_other_prev);

    // Incidence (cases per day)
    const double incS_i = lambda_i * S_eff * coinf_factor;
    const double incR_i = lambda_i * (pRi * Ri) * coinf_factor;
    const double inc_i = incS_i + incR_i;

    // SIR with waning
    dydt[idxS] = -incS_i - added_cases + omega[i] * Ri;
    dydt[idxI] = inc_i + added_cases - gamma[i] * Ii;
    dydt[idxX] = gamma[i] * Ii - Restrain * Xi;
    dydt[idxR] = gamma[i] * Ii - omega[i] * Ri - incR_i;

    out[i] = inc_i + added_cases;
    out_names[i] = "Inc_" + to_string(i + 1);
  }

  out.attr("names") = out_names;

  return List::create(dydt, out);
}
