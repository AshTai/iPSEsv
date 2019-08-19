# iPSEsv
General approach of causal mediation analysis for survival outcome under sequential mediators

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/iPSEsv")


- Required packages:

1. require(timereg)
2. require(survival)


## Example
```R
library(BayICE)
data <- data.frame(time=c(4,3,1,1,2,2,3,3,5,10,2,5,1,7)
status=c(1,1,1,0,1,1,0,1,1,1,0,1,1,0)
x1=rnorm(14,0,1)
x2=rnorm(14,0,1)
x3=rnorm(14,0,1)
c1=rnorm(14,0,1)
c2=rnorm(14,0,1)
c3=rnorm(14,0,1)
sex=c(0,0,1,0,1,1,1,1,1,0,0,0,0,1)
age=c(10,4,2,19,22,31,18,21,41,22,31,29,11,32)
E=sample(c(0,1),size = 14,replace = T))

base.conf <- c("sex","age")
time.conf <- c("c1","c2","c3")
mediators <- c("x1","x2","x3")
DAGs <- c("c1","x1","c2","x2","c3","x3")
exposure <- "E"

output <- iPSEsv(data=data,exposure=exposure,base.conf=base.conf,time.conf=time.conf,
                 mediators=mediators,DAGs=DAGs,suv.model="Aalen")
```

## Contact information
An-Shun Tai ([daansh13@gmail.com](mailto:daansh13@gmail.com)) https://anshuntai.weebly.com
Pei-Hsuan Lin ([a52012232@gmail.com](mailto:a52012232@gmail.com))
Sheng-Hsuan Lin ([shenglin@nctu.edu.tw](mailto:shenglin@nctu.edu.tw))
