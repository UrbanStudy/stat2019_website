data {
  int n;
  int N;
  real aT;
  real bT;
  real aS;
  real aC;
  real bS;
  real bC;
}

parameters {
  real<lower=0,upper=1> pT; //true prevalence
  real<lower=0,upper=1> S; //sensitivity P(test=1 | T=1)
  real<lower=0,upper=1> C; //specificity P(test=0 | T=0)
}

transformed parameters {
  real<lower=0,upper=1> Aprev;
  
  Aprev = pT*S+(1-pT)*(1-C);
}

model {
  pT ~ beta(aT,bT);
  
  S  ~ beta(aS,bS);
  C  ~ beta(aC,bC);
  
  n ~ binomial(N,Aprev);
}

