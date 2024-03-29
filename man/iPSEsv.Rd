% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sPSE.R
\name{iPSEsv}
\alias{iPSEsv}
\title{General approach of causal mediation analysis for survival outcome under sequential mediators.}
\usage{
iPSEsv(data, exposure, base.conf, time.conf, mediators, DAGs,
  suv.model = "Aalen", med.model = "null", time.conf.model = "null",
  a1 = 1, a0 = 0)
}
\arguments{
\item{data}{A data frame, where the column is the variable and row is the sample.}

\item{exposure}{The name of exposure.}

\item{base.conf}{The names of baseline confounders.}

\item{time.conf}{The names of time-varying confounders.}

\item{mediators}{The name of mediators.}

\item{DAGs}{The ordered variable names in the DAG.}

\item{suv.model}{The survival model "Aalen" and "Cox". Default is "Aalen".}
}
\value{
A list of iPSEsv and partPSEsv
}
\description{
The main function of iPSEsv.
}
\examples{
data <- data.frame(time=c(4,3,1,1,2,2,3,3,5,10,2,5,1,7),
status=c(1,1,1,0,1,1,0,1,1,1,0,1,1,0),
x1=rnorm(14,0,1),
x2=rnorm(14,0,1),
x3=rnorm(14,0,1),
c1=rnorm(14,0,1),
c2=rnorm(14,0,1),
c3=rnorm(14,0,1),
sex=c(0,0,1,0,1,1,1,1,1,0,0,0,0,1),
age=c(10,4,2,19,22,31,18,21,41,22,31,29,11,32),
E=sample(c(0,1),size = 14,replace = T))

base.conf <- c("sex","age")
time.conf <- c("c1","c2","c3")
mediators <- c("x1","x2","x3")
DAGs <- c("c1","x1","c2","x2","c3","x3")
exposure <- "E"

output <- iPSEsv(data=data,exposure=exposure,base.conf=base.conf,time.conf=time.conf,
                mediators=mediators,DAGs=DAGs,suv.model="Aalen")
}
\author{
An-Shun Tai \email{daansh13@gmail.com}, Pei-Hsuan Lin \email{a52012232@gmail.com}, and Sheng-Hsuan Lin \email{shenglin@nctu.edu.tw}
}
