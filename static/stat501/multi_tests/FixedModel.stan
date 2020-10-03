data {
  int Nvec[4];
  real aT;
  real bT;
  real aS1;
  real aC1;
  real bS1;
  real bC1;
  real aS2;
  real aC2;
  real bS2;
  real bC2;
  real acovs;
  real bcovs;
  real acovc;
  real bcovc;
}

parameters {
  real<lower=0,upper=1> pT; //true prevalence
  vector<lower=0,upper=1>[2] S; //sensitivity P(test1=1 | T=1), P(test2=1 | T=1)
  vector<lower=0,upper=1>[2] C; //specificity P(test1=0 | T=0), P(test2=0 | T=0)
  vector<lower=0,upper=1>[2] cov12base; //covs=cov(T1,T2 | T=1), covc=cov(T1,T2 | T=0)
}

transformed parameters {
  real ulcovs = min(S)-S[1]*S[2]; 
  real ulcovc = min(C)-C[1]*C[2]; 
  vector<lower=0,upper=1>[4] multp;
  vector<lower=0,upper=1>[2] cov12;
  
  cov12[1] = cov12base[1]*ulcovs; //cov12s
  cov12[2] = cov12base[2]*ulcovc; //cov12c
  {
    multp[1] = pT*(prod(S) + cov12[1]) + (1-pT)*(prod(1-C) + cov12[2]); // P(T1=1, T2=1)
    multp[2] = pT*(S[1]*(1-S[2]) - cov12[1]) + (1-pT)*((1-C[1])*C[2] - cov12[2]); // P(T1=1, T2=0)
    multp[3] = pT*((1-S[1])*S[2] - cov12[1]) + (1-pT)*(C[1]*(1-C[2]) - cov12[2]); // P(T1=0, T2=1)
    multp[4] = pT*(prod(1-S) + cov12[1]) + (1-pT)*(prod(C) + cov12[2]); // P(T1=0, T2=0)
  }
}

model {
  pT ~ beta(aT,bT);
  
  cov12base[1] ~ beta(acovs,bcovs);
  cov12base[2] ~ beta(acovc,bcovc);
  
  S[1]  ~ beta(aS1,bS1);
  C[1]  ~ beta(aC1,bC1);
  S[2]  ~ beta(aS2,bS2);
  C[2]  ~ beta(aC2,bC2);
  
  Nvec ~ multinomial(multp);
}

