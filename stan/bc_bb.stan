//November 2020 started by D. Loughnan 

// This Stan program aims to test for differences in budburst across woody species that dominate BC's forests 

// simple linar model do.bb = a + chill + forcing + photo
data {
  int<lower=0> N;
  int<lower=0> n_sp;
  int<lower=0> n_site;
  int<lower=1, upper=n_sp> sp[N]; // not sure what this is doing
  vector[N] chill; 
  vector[N] site;
  vector[N] force;
  vector[N] photo;
  vector[N] bb; //response var
}

// The parameters accepted by the model. 

parameters {
  real mu_a;
  real mu_force; 
  real mu_chill;
  real mu_photo;
  real mu_site;
  
  vector[n_sp] a_sp;
  vector[n_sp] b_force;
  vector[n_sp] b_photo;
  vector[n_sp] b_chill;
  vector[n_sp] b_site;
  
  real<lower=0> sigma_a;
  real<lower=0> sigma_force;
  real<lower=0> sigma_photo;
  real<lower=0> sigma_chill;
  real<lower=0> sigma_site;
 
  real<lower=0> sigma_y; 
}

transformed parameters{
  vector[N] y_hat;

  for(i in 1:N){
		y_hat[i] = a_sp[sp[i]] + 
		b_site[sp[i]] * site[i] + 
		b_force[sp[i]] * force[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill[sp[i]] * chill[i];
  }
}

model {
  // Priors. Make them flat
	mu_force ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
	mu_photo ~ normal(0, 40);
	mu_chill ~ normal(0, 35);
	mu_site ~ normal(0, 35);
	
	sigma_force ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	sigma_photo ~ normal(0, 10); 
	sigma_chill ~ normal(0, 10);
	sigma_site ~ normal(0, 10);
	
	b_force ~ normal(mu_force, sigma_force);
	b_photo ~ normal(mu_photo, sigma_photo);
	b_chill ~ normal(mu_chill, sigma_chill);
	b_site ~ normal(mu_site, sigma_site);
	
	a_sp ~ normal(mu_a,sigma_a);
  bb ~ normal(y_hat, sigma_y);
}

