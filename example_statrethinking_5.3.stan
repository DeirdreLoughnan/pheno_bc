//This is code adapted from the statistical rethinking example in section 5.3

data {
  int<lower = 0 > N;
  int<lower = 0 > n_sp;
  int<lower = 1, upper= n_sp> sp[N]; 
  vector[N] chill; 
  vector[N] bb;
  
  int site[N];
}

parameters {
  real mu_a;
  real mu_chill;
  
  vector[n_sp] a_sp;
  vector[n_sp] b_chill;
    
  real<lower=0> sigma_chill;
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  vector[4] a_site;
}

transformed parameters {
  vector[N] mu;
  vector[N] y_hat;

  mu = a_site[site];
  
  for(i in 1:N){
	y_hat[i] = a_sp[sp[i]]  +
		b_chill[sp[i]] * chill[i] 
		;  }
}
model {
  mu_chill ~ normal(0, 35);
  sigma_chill ~ normal(0, 30);
  b_chill ~ normal(mu_chill, sigma_chill);
  
  bb ~ normal(y_hat, sigma_y);
  sigma_y ~ uniform(0, 50);
  a_site ~ normal(0, .5);
  a_sp ~ normal(mu_a,sigma_a);
}
