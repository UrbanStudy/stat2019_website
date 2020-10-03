data {
  int N;            //number of subjects
  int T1vec[N];     //response vector test1 
  int T2vec[N];     //response vector test2 
  real aT;          //prior parameter for true prevalence pT
  real bT;          //prior parameter for true prevalence pT
  real aS1;         //prior mean for a11 (sensitivity par test1)
  real bS1;         //prior sd for a11 (sensitivity par test1)
  real aC1;         //prior mean for a10 (specificity par test1)
  real bC1;         //prior sd for a10 (specificity par test1)
  real aS2;         //prior mean for a21 (sensitivity par test2)
  real bS2;         //prior sd for a21 (sensitivity par test2)
  real aC2;         //prior mean for a20 (specificity par test2)
  real bC2;         //prior sd for a20 (specificity par test2)
  real b1mu;        //prior mean coeff random effect b1 (shared for both tests)
  real b1sd;        //prior sd coeff random effect b1 (shared for both tests)
  real b0mu;        //prior mean coeff random effect b0 (shared for both tests)
  real b0sd;        //prior sd coeff random effect b0 (shared for both tests)
}
parameters {
  real<lower=0,upper=1> pT; //true prevalence
  vector[N] Dlatent; // latent (continuous) true response
  real a11; //
  real a21; //
  real a10; //
  real a20; //
  real b11; //
  real b10; //
  // real b21; //
  // real b20; //
  vector[N] ivec; //
}

transformed parameters {
  matrix<lower=0,upper=1>[N,4] multp; //col1: TP1, col2: FP1, col3: TP2, col4 FP2
  real zpT; //quantile such that P(Z<zpT) = pT
  real<lower=0,upper=1> S1; //global sensitivity test 1
  real<lower=0,upper=1> C1; //global specificity test 1
  real<lower=0,upper=1> S2; //global sensitivity test 2
  real<lower=0,upper=1> C2; //global specificity test 2
  
  
  for(i in 1:N){
    multp[i,1] = normal_cdf(a11+b11*ivec[i],0,1); //S1: sensitivity test 1
    multp[i,2] = normal_cdf(-a10-b10*ivec[i],0,1); //1-C1: 1-specificity test 1
    multp[i,3] = normal_cdf(a21+b11*ivec[i],0,1); //S2: sensitivity test 2
    multp[i,4] = normal_cdf(-a20-b10*ivec[i],0,1); //1-C2: 1-specificity test 2
    // multp[i,3] = normal_cdf(a21+b21*ivec[i],0,1); //S2: sensitivity test 2
    // multp[i,4] = normal_cdf(-a20-b20*ivec[i],0,1); //1-C2: 1-specificity test 2
  }
  
  zpT = inv_Phi(pT);
  
  S1 = normal_cdf((a11)/sqrt(1+b11*b11),0,1);
  C1 = normal_cdf((a10)/sqrt(1+b10*b10),0,1);
  S2 = normal_cdf((a21)/sqrt(1+b11*b11),0,1);
  C2 = normal_cdf((a20)/sqrt(1+b10*b10),0,1);
}

model {
  
  pT ~ beta(aT,bT); //tue prevalence
  
  a11~normal(aS1,bS1); //sensitivity par test1
  a21~normal(aS2,bS2); //sensitivity par test2
  a10~normal(aC1,bC1); //specificity par test1
  a20~normal(aC2,bC2); //specificity par test2
  b11~normal(b1mu,b1sd); //sensitivity slope rnd effect both tests
  b10~normal(b0mu,b0sd); //specificity slope rnd effect both tests
  
  // specifying the likelihood using a latent random effect (Dlatent)
  // if Dlatent > 0 --> D = 1, if Dlatent <= 0 --> D = 0 (D is true but unobserved status)
  for(i in 1:N){
    Dlatent[i] ~ normal(zpT,1);  
    ivec[i] ~ normal(0,1);
    
    if(Dlatent[i]>0){
      T1vec[i] ~ bernoulli(multp[i,1]); //likelihood for subject i with T1 when D=1
      T2vec[i] ~ bernoulli(multp[i,3]); //likelihood for subject i with  T2 when D=1
    }else{
      T1vec[i] ~ bernoulli(multp[i,2]); //likelihood for subject i with  T1 when D=0
      T2vec[i] ~ bernoulli(multp[i,4]); //likelihood for subject i with  T2 when D=0  
    }
  }
}


