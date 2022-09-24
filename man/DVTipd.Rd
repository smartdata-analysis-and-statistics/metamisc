\name{DVTipd}
\alias{DVTipd}
\docType{data}
\title{Hypothetical dataset for diagnosis of Deep Vein Thrombosis (DVT)
}
\description{
A hypothetical dataset with 500 subjects suspected of having deep vein thrombosis (DVT).
}
\usage{data(DVTipd)}
\format{
  A data frame with 500 observations of 16 variables.
  \describe{
    \item{\code{sex}}{gender (0=female, 1=male)}
    \item{\code{malign}}{active malignancy (0=no active malignancy, 1=active malignancy)}
    \item{\code{par}}{paresis (0=no paresis, 1=paresis)}
    \item{\code{surg}}{recent surgery or bedridden}
    \item{\code{tend}}{tenderness venous system}
    \item{\code{oachst}}{oral contraceptives or hst}
    \item{\code{leg}}{entire leg swollen}
    \item{\code{notraum}}{absence of leg trauma}
    \item{\code{calfdif3}}{calf difference >= 3 cm}
    \item{\code{pit}}{pitting edema}
    \item{\code{vein}}{vein distension}
    \item{\code{altdiagn}}{alternative diagnosis present}
    \item{\code{histdvt}}{history of previous DVT}
    \item{\code{ddimdich}}{dichotimized D-dimer value}
    \item{\code{dvt}}{final diagnosis of DVT}
    \item{\code{study}}{study indicator}
  }
}
\details{
Hypothetical dataset derived from the Individual Participant Data Meta-Analysis from Geersing \emph{et al} (2014). The dataset consists of consecutive outpatients with suspected deep vein thrombosis, with documented information on the presence or absence of proximal deep vein thrombosis (\code{dvt}) by an acceptable reference test. Acceptable such tests were either compression ultrasonography or venography at initial presentation, or, if venous imaging was not performed, an uneventful follow-up for at least three months.
}
\source{
Geersing GJ, Zuithoff NPA, Kearon C, Anderson DR, Ten Cate-Hoek AJ, Elf JL, et al. Exclusion of deep vein thrombosis using the Wells rule in clinically important subgroups: individual patient data meta-analysis. \emph{BMJ}. 2014;348:g1340. 
}
\examples{
data(DVTipd)
str(DVTipd) 
summary(apply(DVTipd,2,as.factor))

## Develop a prediction model to predict presence of DVT
model.dvt <- glm("dvt~sex+oachst+malign+surg+notraum+vein+calfdif3+ddimdich", 
                  family=binomial, data=DVTipd)
summary(model.dvt)
}
\keyword{datasets}
