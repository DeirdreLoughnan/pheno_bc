// Stan model following Meetup 2016-03-28
// Now with chilling levels as dummy variables for each level, two levels
// Including species as intercept
// Leafout day as a function of species and site of origin as modeled group level factors, and site, temperature, photoperiod, and chilling as unmodeled factors (experimental manipulation)

functions {
  matrix lambda_vcv(matrix vcv, real lambda, real sigma){
    matrix[rows(vcv),cols(vcv)] local_vcv;
    matrix[rows(vcv),cols(vcv)] sigma_mat;  
    local_vcv = vcv * lambda;
    for(i in 1:rows(local_vcv))
      local_vcv[i,i] = vcv[i,i];
    sigma_mat = diag_matrix(rep_vector(sigma, rows(vcv)));
    return(sigma_mat * local_vcv * sigma_mat);
  }
}

data {
  int<lower=0> N;
  int<lower=0> n_sp;
  //int<lower=0> n_site;
  array[N] int<lower=1, upper=n_sp> sp;
  vector[N] lday;
  vector[N] warm;
  vector[N] photo;
  vector[N] chill1; 
  vector [N] site2;
  vector [N] site3;
  vector [N] site4;
  matrix[n_sp,n_sp]Vphy;  
}

 transformed data { 	      // 9 interaction terms
   vector[N] inter_wp;         // vector[N] inter_wp  = warm .* photo;
   vector[N] inter_wc1;           
   vector[N] inter_pc1;
  vector[N] inter_s2c1;
  vector[N] inter_ws2;
  vector[N] inter_ps2;
  
   vector[N] inter_s3c1;
  vector[N] inter_ws3;
  vector[N] inter_ps3;
  // 
   vector[N] inter_s4c1;
  vector[N] inter_ws4;
  vector[N] inter_ps4;

 
 
  inter_wp    = warm .* photo;
  inter_wc1    = warm .* chill1;
  inter_pc1    = photo .* chill1;

  inter_ws2    = warm .* site2;
  inter_ws3    = warm .* site3;
  inter_ws4    = warm .* site4;
  inter_ps2    = photo .* site2;
  inter_ps3    = photo .* site3;
  inter_ps4    = photo .* site4;
  inter_s2c1    = chill1 .* site2;
  inter_s3c1    = chill1 .* site3;
  inter_s4c1    = chill1 .* site4;

 }

parameters {
  
  //real mu_a; 
  real mu_b_warm; 
  real mu_b_chill1;
  real mu_b_photo;

  real mu_b_inter_wp;
  real mu_b_inter_wc1;
  real mu_b_inter_pc1;
  real mu_b_inter_s2c1;
  real mu_b_inter_ws2;
  real mu_b_inter_ps2;
  real mu_b_inter_s3c1;
  real mu_b_inter_ws3;
  real mu_b_inter_ps3;
  real mu_b_inter_s4c1;
  real mu_b_inter_ws4;
  real mu_b_inter_ps4;

  real a_z;
  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;

  vector[n_sp] a_sp;
  //vector[n_sp] b_warm;
  //vector[n_sp] b_photo;
  vector[n_sp] b_chill1;
  
  // vector[n_sp] b_inter_wp;
  // vector[n_sp] b_inter_wc1;
  // vector[n_sp] b_inter_pc1;
  // 
  // vector[n_sp] b_inter_s2c1;
  // vector[n_sp] b_inter_ws2;
  // vector[n_sp] b_inter_ps2;
  // //
  // vector[n_sp] b_inter_s3c1;
  // vector[n_sp] b_inter_ws3;
  // vector[n_sp] b_inter_ps3;
  // //
  // vector[n_sp] b_inter_s4c1;
  // vector[n_sp] b_inter_ws4;
  // vector[n_sp] b_inter_ps4;
  vector[n_sp] b_force_ncp;
  vector[n_sp] b_photo_ncp;
  vector[n_sp] b_inter_wp_ncp;
  vector[n_sp] b_inter_wc1_ncp;
  vector[n_sp] b_inter_pc1_ncp;

  vector[n_sp] b_inter_s2c1_ncp;
  vector[n_sp] b_inter_ws2_ncp;
  vector[n_sp] b_inter_ps2_ncp;

  vector[n_sp] b_inter_s3c1_ncp;
  vector[n_sp] b_inter_ws3_ncp;
  vector[n_sp] b_inter_ps3_ncp;

  vector[n_sp] b_inter_s4c1_ncp;
  vector[n_sp] b_inter_ws4_ncp;
  vector[n_sp] b_inter_ps4_ncp;
  
  real b_site2;
  real b_site3;
  real b_site4;

  //real<lower=0> sigma_a;
  real<lower=0> sigma_b_warm;
  real<lower=0> sigma_b_photo;
  real<lower=0> sigma_b_chill1;

  real<lower=0> sigma_b_inter_wp;
  real<lower=0> sigma_b_inter_wc1;
  real<lower=0> sigma_b_inter_pc1;
 
  real<lower=0> sigma_b_inter_s2c1;
  real<lower=0> sigma_b_inter_ws2;
  real<lower=0> sigma_b_inter_ps2;
  // 
  real<lower=0> sigma_b_inter_s3c1;
  real<lower=0> sigma_b_inter_ws3;
  real<lower=0> sigma_b_inter_ps3;
  // 
  real<lower=0> sigma_b_inter_s4c1;
  real<lower=0> sigma_b_inter_ws4;
  real<lower=0> sigma_b_inter_ps4;
  
  real<lower=0> sigma_y; 
  }

transformed parameters { // Vectorize: Won't save time probably here (no scalar x vector)
  vector[n_sp] b_warm;
  vector[n_sp] b_photo;
  vector[n_sp] b_inter_wp;
  vector[n_sp] b_inter_wc1;
  vector[n_sp] b_inter_pc1;
  vector[n_sp] b_inter_s2c1;
  vector[n_sp] b_inter_ws2;
  vector[n_sp] b_inter_ps2;

  vector[n_sp] b_inter_s3c1;
  vector[n_sp] b_inter_ws3;
  vector[n_sp] b_inter_ps3;

  vector[n_sp] b_inter_s4c1;
  vector[n_sp] b_inter_ws4;
  vector[n_sp] b_inter_ps4;
  
  vector[N] y_hat; // Note to self: all these declarations must happen together!
  b_warm = mu_b_warm + sigma_b_warm*b_force_ncp;
  b_photo = mu_b_photo + sigma_b_photo*b_photo_ncp;

  b_inter_wp = mu_b_inter_wp + sigma_b_inter_wp*b_inter_wp_ncp;
  b_inter_wc1 = mu_b_inter_wc1 + sigma_b_inter_wc1*b_inter_wc1_ncp;
  b_inter_pc1 = mu_b_inter_pc1 + sigma_b_inter_pc1*b_inter_pc1_ncp;
  b_inter_s2c1 = mu_b_inter_s2c1 + sigma_b_inter_s2c1*b_inter_s2c1_ncp;
  b_inter_ws2 = mu_b_inter_ws2 + sigma_b_inter_ws2*b_inter_ws2_ncp;
  b_inter_ps2 = mu_b_inter_ps2 + sigma_b_inter_ps2*b_inter_ps2_ncp;

  b_inter_s3c1 = mu_b_inter_s3c1 + sigma_b_inter_s3c1*b_inter_s3c1_ncp;
  b_inter_ws3 = mu_b_inter_ws3 + sigma_b_inter_ws3*b_inter_ws3_ncp;
  b_inter_ps3 = mu_b_inter_ps3 + sigma_b_inter_ps3*b_inter_ps3_ncp;

  b_inter_s4c1 = mu_b_inter_s4c1 + sigma_b_inter_s4c1*b_inter_s4c1_ncp;
  b_inter_ws4 = mu_b_inter_ws4 + sigma_b_inter_ws4*b_inter_ws4_ncp;
  b_inter_ps4 = mu_b_inter_ps4 + sigma_b_inter_ps4*b_inter_ps4_ncp;

	for(i in 1:N){
		y_hat[i] = a_sp[sp[i]] + 
    b_site2 * site2[i] +
	  b_site3 * site3[i] +
	  b_site4 * site4[i] +
		b_warm[sp[i]] * warm[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill1[sp[i]] * chill1[i]  
		+ b_inter_wp[sp[i]] * inter_wp[i] +
		b_inter_wc1[sp[i]] * inter_wc1[i] +
		b_inter_pc1[sp[i]] * inter_pc1[i]  
		+ b_inter_s2c1[sp[i]] * inter_s2c1[i] +
		b_inter_ws2[sp[i]] * inter_ws2[i] +
		b_inter_ps2[sp[i]] * inter_ps2[i] +
		b_inter_s3c1[sp[i]] * inter_s3c1[i] +
		b_inter_ws3[sp[i]] * inter_ws3[i] +
		b_inter_ps3[sp[i]] * inter_ps3[i]+
		b_inter_s4c1[sp[i]] * inter_s4c1[i] 
		+ b_inter_ws4[sp[i]] * inter_ws4[i] 
		+ b_inter_ps4[sp[i]] * inter_ps4[i]
		;
		}
	
}

model {
	// Priors //
	a_sp ~ multi_normal(rep_vector(a_z,n_sp), lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa)); 
	a_z ~ normal(5,3); //(5,10)
	lam_interceptsa ~ beta(1,1);
  sigma_interceptsa ~ normal(1,15); //(1,10)

//	mu_b_warm ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
	// mu_b_photo ~ normal(0, 35);
	mu_b_chill1 ~ normal(0, 35);
	
	// mu_b_inter_wp ~ normal(0, 35);
	// mu_b_inter_wc1 ~ normal(0, 35);	
	// mu_b_inter_pc1 ~ normal(0, 35);	
	
  b_site2 ~ normal(0,5);
  b_site3 ~ normal(0,5);
  b_site4 ~ normal(0,5);
	
	sigma_b_warm ~ normal(0, 20); // Start big at 10, go smaller if introduces problems
	sigma_b_photo ~ normal(0, 10); 
	sigma_b_chill1 ~ normal(0, 10);

	sigma_b_inter_wp ~ normal(0, 10);
	sigma_b_inter_wc1 ~ normal(0, 10);
	sigma_b_inter_pc1 ~ normal(0, 10);
	sigma_b_inter_s2c1 ~ normal(0, 10);
	sigma_b_inter_ws2 ~ normal(0, 10);
	sigma_b_inter_ps2 ~ normal(0, 10);

	sigma_b_inter_s3c1 ~ normal(0, 10);
	sigma_b_inter_ws3 ~ normal(0, 10);
	sigma_b_inter_ps3 ~ normal(0, 10);

	sigma_b_inter_s4c1 ~ normal(0, 10);
	sigma_b_inter_ws4 ~ normal(0, 10);
	sigma_b_inter_ps4 ~ normal(0, 10);

	//a_sp ~ normal(mu_a, sigma_a);  
	
	//b_warm ~ normal(mu_b_warm, sigma_b_warm);
	// b_photo ~ normal(mu_b_photo, sigma_b_photo);
	b_chill1 ~ normal(mu_b_chill1, sigma_b_chill1);

//   b_inter_wp ~ normal(mu_b_inter_wp, sigma_b_inter_wp); // Delete because all in NCP now
//   b_inter_wc1 ~ normal(mu_b_inter_wc1, sigma_b_inter_wc1);		
// 	b_inter_pc1 ~ normal(mu_b_inter_pc1, sigma_b_inter_pc1);
// 	
// 	b_inter_ws2 ~ normal(mu_b_inter_ws2, sigma_b_inter_ws2);
// 	b_inter_ps2 ~ normal(mu_b_inter_ps2, sigma_b_inter_ps2);
// 	b_inter_s2c1 ~ normal(mu_b_inter_s2c1, sigma_b_inter_s2c1);
// 	// 
// 	b_inter_ws3 ~ normal(mu_b_inter_ws3, sigma_b_inter_ws3);
// 	b_inter_ps3 ~ normal(mu_b_inter_ps3, sigma_b_inter_ps3);
// 	b_inter_s3c1 ~ normal(mu_b_inter_s3c1, sigma_b_inter_s3c1);
// 	// 
// 	b_inter_ws4 ~ normal(mu_b_inter_ws4, sigma_b_inter_ws4);
// 	b_inter_ps4 ~ normal(mu_b_inter_ps4, sigma_b_inter_ps4);
// 	b_inter_s4c1 ~ normal(mu_b_inter_s4c1, sigma_b_inter_s4c1);
  b_force_ncp ~ normal(0, 35); 
	b_photo_ncp ~ normal(0, 35); 

  b_inter_wp_ncp ~ normal(0, 35);
	b_inter_wc1_ncp ~ normal(0, 35);
	b_inter_pc1_ncp ~ normal(0, 35);
	b_inter_s2c1_ncp ~ normal(0, 35);
	b_inter_ws2_ncp ~ normal(0, 35);
	b_inter_ps2_ncp ~ normal(0, 35);

	b_inter_s3c1_ncp ~ normal(0, 35);
	b_inter_ws3_ncp ~ normal(0, 35);
	b_inter_ps3_ncp ~ normal(0, 35);

	b_inter_s4c1_ncp ~ normal(0, 35);
	b_inter_ws4_ncp ~ normal(0, 35);
	b_inter_ps4_ncp ~ normal(0, 35);
	
	lday ~ normal(y_hat, sigma_y);

}

