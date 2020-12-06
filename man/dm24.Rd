\name{dm24}
\alias{dm24}
\docType{data}
\title{
Diabetes Mellitus and COVID-19 severity data
%%   ~~ data name/kind ... ~~
}
\description{
Meta-analysis data set including 24 observational studies for exploring the relationship between underlying diabetes mellitus
and the severe course of COVID-19 disease. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(dm24)}
\format{
  A data frame with 24 observations on the following 11 variables:
  \describe{
    \item{\code{studyName}}{character string, label of the study}
    \item{\code{y}}{numeric, log odds ratio estimate within a study }
    \item{\code{sigma}}{numeric, standard error of the log odds ratio estimate within a study}
    \item{\code{labels}}{integer, the index of the study}
    \item{\code{n_total}}{integer, total number of patients within a study}
    \item{\code{age_mean}}{numeric, mean age within a study}
    \item{\code{age_sd}}{numeric, standard deviation of the age of patients}
    \item{\code{n_males}}{integer, number of male patients of patients within a study}
    \item{\code{n_severe}}{integer, number of patients with a severe course of COVID19 disease.}
    \item{\code{reason}}{character string, the reason of categorisation}
    \item{\code{n_diabetes}}{integer, number of patients with diabetes mellitus within a study}
  }
}
\details{
The data set was collected and analyzed by Kumar et al. (2020). Overall, 13954 patients were included in 24 studies. Most studies were from mainland China. Diabetes was present in 1405 patients. The largest study (number 20) based on registry data from 6637 patients was only classified as a fair quality study. Kumar et al. (2020) omitted this study and conducted a second meta-analysis based on 23 studies. 


%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Kumar, A., Arora, A., Sharma, P., Anikhindi, S. A., Bansal, N., Singla, V., Khare, S., Srivastava, A. (2020). Is diabetes mellitus associated with mortality and severity of COVID-19? A meta-analysis. 
\emph{Diabetes & Metabolic Syndrome: Clinical Research & Reviews} \bold{14}, 535-545.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(dm24)
str(dm24)

}
\keyword{datasets}
