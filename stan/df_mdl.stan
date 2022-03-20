// Stan model following Meetup 2016-03-28
// Now with chilling levels as dummy variables for each level, two levels
// Including species as intercept
// Leafout day as a function of species and site of origin as modeled group level factors, and site, temperature, photoperiod, and chilling as unmodeled factors (experimental manipulation)

// Note by Lizzie: This model includes pooled intercepts
// Adding non-centered parameterization (NCP or ncp) on b_inter_sp to start 




data {
  int<lower=0> N;
  int<lower=0> n_sp;
  int<lower=0> n_site;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] lday;
  vector[N] warm;
  vector[N] photo;
  vector[N] chill1; 
  // vector[N] chill2;
  vector[N] site;
  // vector <lower=0, upper = 1>[N] site2;
  // vector <lower=0, upper = 1>[N] site3;
}

transformed data { 	      // 9 interaction terms
  vector[N] inter_wp;         // vector[N] inter_wp  = warm .* photo;
  vector[N] inter_wc1;           
  // vector[N] inter_wc2;           
  vector[N] inter_pc1;           
  // vector[N] inter_pc2;           
  vector[N] inter_sc1;
  vector[N] inter_ws;           
  vector[N] inter_ps; 
  // vector[N] inter_sc2;       
  // vector[N] inter_s2c1;
  // vector[N] inter_ws2;           
  // vector[N] inter_ps2; 

  inter_wp    = warm .* photo; 
  inter_ws    = warm .* site;
  inter_ps    = photo .* site;
  // inter_ws    = warm .* site2;  
  // inter_ws    = warm .* site3; 
  // inter_ps    = photo .* site2;  
  // inter_ps    = photo .* site3;  
  inter_wc1    = warm .* chill1;  
  // inter_wc2    = warm .* chill2;  
  inter_pc1    = photo .* chill1;  
  // inter_pc2    = photo .* chill2;  
  // inter_sc1    = site2 .* chill1;  
  // inter_sc1    = site3 .* chill1; 
  inter_sc1    = site .* chill1;  

}

parameters {
  vector[n_sp] a_sp;
  vector[n_sp] b_warm;
  vector[n_sp] b_photo;
  vector[n_sp] b_chill1;
  // vector[n_sp] b_chill2;
  vector[n_sp] b_site;
  // real b_site2;
  // real b_site3;

  vector[n_sp] b_inter_wp_ncp;
  vector[n_sp] b_inter_wc1_ncp;
  // vector[n_sp] b_inter_wc2_ncp;
  vector[n_sp] b_inter_pc1_ncp;
  // vector[n_sp] b_inter_pc2_ncp;
  vector[n_sp] b_inter_sc1_ncp;
  vector[n_sp] b_inter_ws_ncp;
  vector[n_sp] b_inter_ps_ncp;
  // vector[n_sp] b_inter_sc2_ncp;
  // vector[n_sp] b_inter_s2c1_ncp;
  // vector[n_sp] b_inter_ws2_ncp;
  // vector[n_sp] b_inter_ps2_ncp;

  real mu_a; 
  real mu_b_warm; 
  real mu_b_chill1;
  // real mu_b_chill2;
  real mu_b_photo;
  real mu_b_site;

  real mu_b_inter_wp;
  real mu_b_inter_wc1;
  // real mu_b_inter_wc2;
  real mu_b_inter_pc1;
  // real mu_b_inter_pc2;
  real mu_b_inter_sc1;
  real mu_b_inter_ws;
  real mu_b_inter_ps;
  // real mu_b_inter_sc2;
  // real mu_b_inter_s2c1;
  // real mu_b_inter_ws2;
  // real mu_b_inter_ps2;

  real<lower=0> sigma_b_warm;
  real<lower=0> sigma_b_photo;
  real<lower=0> sigma_b_chill1;
  // real<lower=0> sigma_b_chill2;
  real<lower=0> sigma_b_site;

  real<lower=0> sigma_a;

  real<lower=0> sigma_b_inter_wp;
  real<lower=0> sigma_b_inter_wc1;
  // real<lower=0> sigma_b_inter_wc2;
  real<lower=0> sigma_b_inter_pc1;
  // real<lower=0> sigma_b_inter_pc2;
  real<lower=0> sigma_b_inter_sc1;
  real<lower=0> sigma_b_inter_ws;
  real<lower=0> sigma_b_inter_ps;
  // real<lower=0> sigma_b_inter_sc2;
  // real<lower=0> sigma_b_inter_s2c1;
  // real<lower=0> sigma_b_inter_ws2;
  // real<lower=0> sigma_b_inter_ps2;
  
  real<lower=0> sigma_y; 
  }


transformed parameters { // Vectorize: Won't save time probably here (no scalar x vector)
  vector[n_sp] b_inter_wp;
  vector[n_sp] b_inter_wc1;
  // vector[n_sp] b_inter_wc2;
  vector[n_sp] b_inter_pc1;
  // vector[n_sp] b_inter_pc2;
  vector[n_sp] b_inter_sc1;
   vector[n_sp] b_inter_ws;
  vector[n_sp] b_inter_ps;
  // vector[n_sp] b_inter_sc2;
  //   vector[n_sp] b_inter_s2c1;
  //  vector[n_sp] b_inter_ws2;
  // vector[n_sp] b_inter_ps2;
  
  vector[N] y_hat; // Note to self: all these declarations must happen together!

  b_inter_wp = mu_b_inter_wp + sigma_b_inter_wp*b_inter_wp_ncp;
  b_inter_wc1 = mu_b_inter_wc1 + sigma_b_inter_wc1*b_inter_wc1_ncp;
  // b_inter_wc2 = mu_b_inter_wc2 + sigma_b_inter_wc2*b_inter_wc2_ncp;
  b_inter_pc1 = mu_b_inter_pc1 + sigma_b_inter_pc1*b_inter_pc1_ncp;
  // b_inter_pc2 = mu_b_inter_pc2 + sigma_b_inter_pc2*b_inter_pc2_ncp;
  b_inter_sc1 = mu_b_inter_sc1 + sigma_b_inter_sc1*b_inter_sc1_ncp;
  b_inter_ws = mu_b_inter_ws + sigma_b_inter_ws*b_inter_ws_ncp;
  b_inter_ps = mu_b_inter_ps + sigma_b_inter_ps*b_inter_ps_ncp;

  // b_inter_sc2 = mu_b_inter_sc2 + sigma_b_inter_sc2*b_inter_sc2_ncp;
// 	b_inter_s2c1 = mu_b_inter_s2c1 + sigma_b_inter_s2c1*b_inter_s2c1_ncp;
//   b_inter_ws2 = mu_b_inter_ws + sigma_b_inter_ws2*b_inter_ws2_ncp;
//   b_inter_ps2 = mu_b_inter_ps + sigma_b_inter_ps2*b_inter_ps2_ncp;
  
	for(i in 1:N){
		y_hat[i] = a_sp[sp[i]] + 
		b_site[sp[i]] * site[i] + 
//    b_site2 * site2[i] +
// 	  b_site3 * site3[i] +
		b_warm[sp[i]] * warm[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill1[sp[i]] * chill1[i] + 
		// b_chill2[sp[i]] * chill2[i] +
		b_inter_wp[sp[i]] * inter_wp[i] +
		b_inter_wc1[sp[i]] * inter_wc1[i] +
		// b_inter_wc2[sp[i]] * inter_wc2[i] +
		b_inter_pc1[sp[i]] * inter_pc1[i] +
		// b_inter_pc2[sp[i]] * inter_pc2[i] +
		b_inter_sc1[sp[i]] * inter_sc1[i] +
		b_inter_ws[sp[i]] * inter_ws[i] +
		b_inter_ps[sp[i]] * inter_ps[i] 
		// b_inter_sc2[sp[i]] * inter_sc2[i] 
		// b_inter_s2c1[sp[i]] * inter_s2c1[i] +
		// b_inter_ws2[sp[i]] * inter_ws2[i] +
		// b_inter_ps2[sp[i]] * inter_ps2[i] 
		;
				
		}
	
}

model {
	// Priors //
	mu_b_warm ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
	mu_b_photo ~ normal(0, 35);
	mu_b_chill1 ~ normal(0, 35);
	// mu_b_chill2 ~ normal(0, 35);
	mu_b_site ~ normal(0, 35);
  // b_site2 ~ normal(0,5);
  // b_site3 ~ normal(0,5);

/*	mu_b_inter_wp ~ normal(0, 35); // Delete because all in NCP now
	mu_b_inter_ws ~ normal(0, 35);
	mu_b_inter_ps ~ normal(0, 35);
	mu_b_inter_wc1 ~ normal(0, 35);	
	mu_b_inter_wc2 ~ normal(0, 35);	
	mu_b_inter_pc1 ~ normal(0, 35);	
	mu_b_inter_pc2 ~ normal(0, 35);	
	mu_b_inter_sc1 ~ normal(0, 35);	
	mu_b_inter_sc2 ~ normal(0, 35);	*/
	
	sigma_b_warm ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	sigma_b_photo ~ normal(0, 10); 
	sigma_b_chill1 ~ normal(0, 10);
	// sigma_b_chill2 ~ normal(0, 10);
	sigma_b_site ~ normal(0, 10);

	sigma_b_inter_wp ~ normal(0, 10);
	sigma_b_inter_wc1 ~ normal(0, 10);	
	// sigma_b_inter_wc2 ~ normal(0, 10);	
	sigma_b_inter_pc1 ~ normal(0, 10);	
	// sigma_b_inter_pc2 ~ normal(0, 10);	
	sigma_b_inter_sc1 ~ normal(0, 10);	
	sigma_b_inter_ws ~ normal(0, 10);
	sigma_b_inter_ps ~ normal(0, 10);
	// sigma_b_inter_s2c1 ~ normal(0, 10);	
	// sigma_b_inter_ws2 ~ normal(0, 10);
	// sigma_b_inter_ps2 ~ normal(0, 10);
	// sigma_b_inter_sc2 ~ normal(0, 10);

	a_sp ~ normal(mu_a, sigma_a);  
	
	b_warm ~ normal(mu_b_warm, sigma_b_warm);
	b_photo ~ normal(mu_b_photo, sigma_b_photo);
	b_chill1 ~ normal(mu_b_chill1, sigma_b_chill1);
	// b_chill2 ~ normal(mu_b_chill2, sigma_b_chill2);
	b_site ~ normal(mu_b_site, sigma_b_site);

/*	b_inter_wp ~ normal(mu_b_inter_wp, sigma_b_inter_wp); // Delete because all in NCP now
	b_inter_ws ~ normal(mu_b_inter_ws, sigma_b_inter_ws);
	b_inter_ps ~ normal(mu_b_inter_ps, sigma_b_inter_ps);
	b_inter_wc1 ~ normal(mu_b_inter_wc1, sigma_b_inter_wc1);		
	b_inter_wc2 ~ normal(mu_b_inter_wc2, sigma_b_inter_wc2);
	b_inter_pc1 ~ normal(mu_b_inter_pc1, sigma_b_inter_pc1);		
	b_inter_pc2 ~ normal(mu_b_inter_pc2, sigma_b_inter_pc2);
	b_inter_sc1 ~ normal(mu_b_inter_sc1, sigma_b_inter_sc1);		
	b_inter_sc2 ~ normal(mu_b_inter_sc2, sigma_b_inter_sc2); */
  b_inter_wp_ncp ~ normal(0, 35); 
	b_inter_wc1_ncp ~ normal(0, 35);		
	// b_inter_wc2_ncp ~ normal(0, 35);
	b_inter_pc1_ncp ~ normal(0, 35);		
	// b_inter_pc2_ncp ~ normal(0, 35);
	b_inter_sc1_ncp ~ normal(0, 35);	
	b_inter_ws_ncp ~ normal(0, 35);
	b_inter_ps_ncp ~ normal(0, 35);
	// b_inter_sc2_ncp ~ normal(0, 35);
	// b_inter_s2c1_ncp ~ normal(0, 35);	
	// b_inter_ws2_ncp ~ normal(0, 35);
	// b_inter_ps2_ncp ~ normal(0, 35);

	
	lday ~ normal(y_hat, sigma_y);

}

