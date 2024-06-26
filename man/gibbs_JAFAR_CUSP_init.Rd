% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAFAR_D-CUSP_initialization.R
\name{gibbs_JAFAR_CUSP_init}
\alias{gibbs_JAFAR_CUSP_init}
\title{Initialize all unknown/latent random variables in the Gibbs sampler for JAFAR under the D-CUSP prior}
\usage{
gibbs_JAFAR_CUSP_init(
  n,
  M,
  p_m,
  K0,
  K0_m,
  a_sig,
  b_sig,
  a_theta,
  b_theta,
  var_spike_theta,
  a_xi,
  b_xi,
  a_m,
  b_m,
  prec0,
  prec0m,
  var_spike,
  a_chi,
  b_chi,
  alpha,
  alpha_loc,
  seed
)
}
\arguments{
\item{n}{number of observations}

\item{M}{number of views}

\item{p_m}{views dimensions (length M)}

\item{K0}{Initial number of latent factors in shared components}

\item{K0_m}{Initial number of latent factors in view-specific components (length M)}

\item{a_sig}{Shape parameter of inverse-gamma prior on response noise}

\item{b_sig}{Scale parameter of inverse-gamma prior on response noise}

\item{a_theta}{Shape parameter in slab element of prior for response loadings}

\item{b_theta}{Scale parameter in slab element of prior for response loadings}

\item{var_spike_theta}{Variance parameter in spike element of prior for response-loadings}

\item{a_xi}{Shape1 parameters in beta of prior on response-loadings spike and slab weights}

\item{b_xi}{Shape2 parameters in beta of prior on response-loadings spike and slab weights}

\item{a_m}{Shape parameters of inverse-gamma prior on predictors idiosyncratic components}

\item{b_m}{Scale parameters of inverse-gamma prior on predictors idiosyncratic components}

\item{prec0}{Prior precision for response intercept}

\item{prec0m}{Prior precisions for predictors intercepts}

\item{var_spike}{Variance parameter in spike element of prior for predictors-loadings}

\item{a_chi}{Shape parameters in slab element of prior fo shared-component loadings}

\item{b_chi}{Scale parameters in slab element of prior for shared-component loadings}

\item{alpha}{Stick-breaking parameter in shared-component loadings}

\item{alpha_loc}{Stick-breaking parameter in view-specific loadings}

\item{seed}{See for random number generation}
}
\value{
A list of initial values for all unknown/latent random variables in the Gibbs sampler for JAFAR under the D-CUSP prior
}
\description{
Initialize all unknown/latent random variables in the Gibbs sampler for JAFAR under the D-CUSP prior
}
