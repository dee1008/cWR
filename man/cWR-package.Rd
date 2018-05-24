\name{cWR-package}
\alias{cWR-package}
\alias{cWR}
\title{
  Composite endpoint analysis for time-to-event data using Win Ratio
}
\description{
  This package uses Win Ratio (WR) as a summary statistic to compare the composite     endpoints of time-to-event data between two groups. Options of clustered or          independent time-to-event data can be specified.
}
\usage{
cWR(treatment, cluster, 
    y1, y2, 
    delta1, delta2, 
    null.WR=1, alpha.sig=0.05)
}
\arguments{
	\item{treatment}{An integer vector with code 0 as control group and 1 as treatment group for each subject}
	\item{cluster}{An integer vector with unique cluster ID for each cluster. When subjects are independent, the cluster ID is unique for each subject.}
	\item{y1}{Let T_H, T_D and T_C be time to non-fatal event, time to fatal event and censoring time, respectively. 
		y1 is a numeric vector with min(T_H, T_D, T_C) for each subject.}
	\item{y2}{A numeric vector with min(T_D, T_C) for each subject}
	\item{delta1}{An integer vector with code 1 indicating that T_H is observed, 0 otherwise}
	\item{delta2}{An integer vector with code 1 indicating that T_D is observed, 0 otherwise}
	\item{null.WR}{The null hypothesis of the WR statistic. The default is H0: WR=1 or log(WR)=0.}
	\item{alpha.sig}{The significance level, with default value 0.05}
}
\value{
	\item{name}{The test name}
	\item{U1}{First estimated clustered U-statistic}
	\item{U2}{Second estimated clustered U-statistic}
	\item{logWR}{Estimated WR on log scale}
	\item{se}{Estimated standard error of the WR on log scale}
	\item{z}{Test statistic}
	\item{ci}{100(1-alpha.sig)\% confidence interval}
	\item{p}{P-value of the significance testing}
	\item{var_cov}{Variance and covariance matrix of the first and second clustered U-statistics}
}

\details{
The function "cWR" performs significance testing of comparing two composite time-to-event outcomes between groups. The Win Ratio summary statistic is built on the "unmatched" approach described by Pocock et al. (2011). We assume that the composite endpoints can be formulated as semi-competing risk data. Each individual in the study is measured on time to non-fatal (non-terminal) event (e.g. hospitalization) and time to fatal (terminal) event (e.g. death). Specifically, the fatal event is considered clinically more important compared to the non-fatal event. Censoring is allowed, but time to censor needs to be observed. 

This function can handle independent data, as well as clustered data. The inference of clustered data is based on the generalized bivariate clustered U-statistics proposed by Zhang and Jeong (2018). This clustered U-statistic accounts for the potential correlations among subjects within a cluster. When the cluster size is 1, it's the independent setting and the inference is the same as the method proposed by Bebu and Lachin (2015). 
   
Note: The option "treatment", "cluster", "y1", "y2", "delta1", "delta2" are required and no defaults are provided. These options have to be vectors with the same length. No missing values are allowed. 
}
\author{
Di Zhang <diz11@pitt.edu>;
Jong H. Jeong <jjeong@pitt.edu>

Maintainer: Di Zhang <diz11@pitt.edu>
}
\references{
Pocock, S. J., Ariti, C. A., Collier, T. J., andWang, D. (2011). The win ratio: a new approach
to the analysis of composite endpoints in clinical trials based on clinical priorities.
European heart journal 33, 176-182.

Bebu, I. and Lachin, J. M. (2015). Large sample inference for a win ratio analysis of a
composite outcome based on prioritized components. Biostatistics 17, 178-187.

Zhang, D. and Jeong, H. J. Inference on the Win Ratio for Clustered Semi-Competing Risk Data. (Drafted, 2018)
}
\keyword{ Clinical Trials; Cluster Randomization; Cramer-Wold Theorem; Pragmatic Trials; Semi-Competing
Risks; Suvival; Time-to-event. }
\seealso{

}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
