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

transformed data { 
// 9 interaction terms
  vector[N] inter_fp;       // vector[N] inter_fp  = force .* photo;
  vector[N] inter_fs;       // force*site         
  vector[N] inter_ps;       // photoperiod by site   
  vector[N] inter_fc;       // force*chill            
  vector[N] inter_pc;       // photo*chill            
  vector[N] inter_sc;       // chill * site           
          
//what is the dot doing? bc vectors of length n?
  inter_fp    = force .* photo; 
  inter_fs    = force .* site;  
  inter_ps    = photo .* site;  
  inter_fc    = force .* chill;   
  inter_pc    = photo .* chill;  
  inter_sc    = site .* chill;   
}

// The parameters accepted by the model. 

parameters {
  real mu_a;
  real mu_force; 
  real mu_chill;
  real mu_photo;
  real mu_site;
  real mu_inter_fp;
  real mu_inter_fs;
  real mu_inter_ps;
  real mu_inter_fc;
  real mu_inter_pc;
  real mu_inter_sc;
  
  vector[n_sp] a_sp;
  vector[n_sp] b_force;
  vector[n_sp] b_chill;
  vector[n_sp] b_site;
  vector[n_sp] b_photo_ncp; //
  vector[n_sp] b_inter_fp_ncp;
  vector[n_sp] b_inter_fs_ncp;
  vector[n_sp] b_inter_ps_ncp;
  vector[n_sp] b_inter_fc_ncp;
  vector[n_sp] b_inter_pc_ncp;
  vector[n_sp] b_inter_sc_ncp;
  
  real<lower=0> sigma_a;
  real<lower=0> sigma_force;
  real<lower=0> sigma_chill;
  real<lower=0> sigma_site;
  real<lower=0> sigma_photo;
  real<lower=0> sigma_b_inter_fp;
  real<lower=0> sigma_b_inter_fs;
  real<lower=0> sigma_b_inter_ps;
  real<lower=0> sigma_b_inter_fc;
  real<lower=0> sigma_b_inter_pc;
  real<lower=0> sigma_b_inter_sc;
 
  real<lower=0> sigma_y; 
}

transformed parameters{
  vector[n_sp] b_photo;
  
  vector[n_sp] b_inter_fp;
  vector[n_sp] b_inter_fs;
  vector[n_sp] b_inter_ps;
  vector[n_sp] b_inter_fc;
  vector[n_sp] b_inter_pc;
  vector[n_sp] b_inter_sc;
  
  vector[N] y_hat;
  
  b_photo = mu_photo + sigma_photo * b_photo_ncp;
  
  b_inter_fp = mu_inter_fp + sigma_b_inter_fp*b_inter_fp_ncp;
  b_inter_fs = mu_inter_fs + sigma_b_inter_fs*b_inter_fs_ncp;
  b_inter_ps = mu_inter_ps + sigma_b_inter_ps*b_inter_ps_ncp;
  b_inter_fc = mu_inter_fc + sigma_b_inter_fc*b_inter_fc_ncp;
  b_inter_pc = mu_inter_pc + sigma_b_inter_pc*b_inter_pc_ncp;
  b_inter_sc = mu_inter_sc + sigma_b_inter_sc*b_inter_sc_ncp;

  for(i in 1:N){
		y_hat[i] = a_sp[sp[i]] + 
		b_site[sp[i]] * site[i] + 
		b_force[sp[i]] * force[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill[sp[i]] * chill[i] +
		b_inter_fp[sp[i]] * inter_fp[i] +
		b_inter_fs[sp[i]] * inter_fs[i] +
		b_inter_ps[sp[i]] * inter_ps[i] +
		b_inter_fc[sp[i]] * inter_fc[i] +
		b_inter_pc[sp[i]] * inter_pc[i] +
		b_inter_sc[sp[i]] * inter_sc[i] 
		;
  }
}

model {
  // Priors. Make them flat
	mu_force ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
	mu_photo ~ normal(0, 35);
	mu_chill ~ normal(0, 35);
	mu_site ~ normal(0, 35);
	
	mu_inter_fp ~ normal(0,35);
	mu_inter_fc ~ normal(0,35);
	mu_inter_pc ~ normal(0,35);
	mu_inter_fs ~ normal(0,35);
	mu_inter_ps ~ normal(0,35);
	mu_inter_sc ~ normal(0,35);
	
	sigma_force ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	sigma_photo ~ normal(0, 10); 
	sigma_chill ~ normal(0, 10);
	sigma_site ~ normal(0, 10);
	sigma_b_inter_fp ~ normal(0, 10);
	sigma_b_inter_fs ~ normal(0, 10);
	sigma_b_inter_ps ~ normal(0, 10);
	sigma_b_inter_fc ~ normal(0, 10);	
	sigma_b_inter_pc ~ normal(0, 10);	
	sigma_b_inter_sc ~ normal(0, 10);
	
	b_photo_ncp ~normal(0,1);
	
	b_force ~ normal(mu_force, sigma_force);
//  b_photo ~ normal(mu_photo, sigma_photo); // bc still need this info
	b_chill ~ normal(mu_chill, sigma_chill);
	b_site ~ normal(mu_site, sigma_site);
	b_inter_fp_ncp ~ normal(0, 1); 
  b_inter_fs_ncp ~ normal(0, 1);
	b_inter_ps_ncp ~ normal(0, 1);		
	b_inter_fc_ncp ~ normal(0, 1);
	b_inter_pc_ncp ~ normal(0, 1);		
	b_inter_sc_ncp ~ normal(0, 1);	
	
// 	b_inter_fp ~ normal(0, 35); 
//   b_inter_fs ~ normal(0, 35);
// 	b_inter_ps ~ normal(0, 35);		
// 	b_inter_fc ~ normal(0, 35);
// 	b_inter_pc ~ normal(0, 35);		
// 	b_inter_sc ~ normal(0, 35);
	
	a_sp ~ normal(mu_a,sigma_a);
  bb ~ normal(y_hat, sigma_y);
}

generated quantities{
   real ypred_new[N];
   
   for (i in 1:N)
   ypred_new[i] = normal_rng(y_hat[i], sigma_y);
}

