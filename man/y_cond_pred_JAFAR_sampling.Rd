% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JAFAR_D-CUSP_predictions.R
\name{y_cond_pred_JAFAR_sampling}
\alias{y_cond_pred_JAFAR_sampling}
\title{Compute the linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration via sampling of the latent factors}
\usage{
y_cond_pred_JAFAR_sampling(
  Xpred,
  nPred,
  M,
  p_m,
  K,
  K_m,
  Theta,
  Lambda_m,
  Gamma_m,
  mu_y,
  s2_inv_y,
  mu_m,
  s2_inv_m,
  rescale_pred = FALSE
)
}
\arguments{
\item{Xpred}{Multi-view features in out-of-sample observations}

\item{nPred}{Number of out-of-sample observations}

\item{M}{number of views}

\item{p_m}{views dimensions (length M)}

\item{K}{Number of latent factors in shared components}

\item{K_m}{Number of latent factors in view-specific components (length M)}

\item{Theta}{Response loadings}

\item{Lambda_m}{Shared-component loadings}

\item{Gamma_m}{View-specific loadings}

\item{mu_y}{Response Intercepts}

\item{s2_inv_y}{Response noise precision}

\item{mu_m}{Intercepts of idiosyncratic components in multi-view features}

\item{s2_inv_m}{Precisions of idiosyncratic components in multi-view features}

\item{rescale_pred}{Rescale loadings to induce correlation matrix}
}
\value{
The values of linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration
}
\description{
Compute the linear predictor of y|X=x in out-of-sample observations for a single MCMC iteration via sampling of the latent factors
}
