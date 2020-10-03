data {
  int Nvec[4];
  int N;
  real aT;
  real bT;
  real aS1;
  real bS1;
  real aC1;
  real bC1;
  real aS2p;
  real bS2p;
  real aC2p;
  real bC2p;
  real aS2n;
  real bS2n;
  real aC2n;
  real bC2n;
}

parameters {
  real<lower=0,upper=1> pT; //true prevalence
  real<lower=0,upper=1> S1; //sensitivity P(test1=1 | T=1)
  real<lower=0,upper=1> C1; //specificity P(test1=0 | T=0)
  real<lower=0,upper=1> S2_t1p; //cond sensitivity P(test2=1 | T=1, test1=1)
  real<lower=0,upper=1> S2_t1n; //cond sensitivity P(test2=1 | T=1, test1=0)
  real<lower=0,upper=1> C2_t1p; //cond specificity P(test2=0 | T=0, test1=1)
  real<lower=0,upper=1> C2_t1n; //cond specificity P(test2=0 | T=0, test1=0)
}

transformed parameters {
  real<lower=0,upper=1> S2;
  real<lower=0,upper=1> C2;
  real<lower=0,upper=1> rho12s;
  real<lower=0,upper=1> rho12c;
  vector[4] multp;
  
  S2 = S1*S2_t1p + (1-S1)*S2_t1n;
  C2 = (1-C1)*C2_t1p + C1*C2_t1n;
  
  //multinomial probablities
  multp[1] = pT*S1*S2_t1p+(1-pT)*(1-C1)*(1-C2_t1p); //P(test1=1,test2=1)
  multp[2] = pT*S1*(1-S2_t1p)+(1-pT)*(1-C1)*C2_t1p; //P(test1=1,test2=0)
  multp[3] = pT*(1-S1)*S2_t1n + (1-pT)*C1*(1-C2_t1n); //P(test1=0,test2=1)
  multp[4] = pT*(1-S1)*(1-S2_t1n) + (1-pT)*C1*C2_t1n; //P(test1=0,test2=0)
  
  //correlations
  rho12s = S1*(S2_t1p - S2)/sqrt(S1*(1-S1)*S2*(1-S2)); 
  rho12c = C1*(C2_t1n - C2)/sqrt(C1*(1-C1)*C2*(1-C2));
}

model {
  pT ~ beta(aT,bT);
  
  S1  ~ beta(aS1,bS1);
  C1  ~ beta(aC1,bC1);
  
  S2_t1p ~ beta(aS2p,bS2p);
  C2_t1p ~ beta(aC2p,bC2p);
  S2_t1n ~ beta(aS2n,bS2n);
  C2_t1n ~ beta(aC2n,bC2n);
  
  Nvec ~ multinomial(multp);
}

