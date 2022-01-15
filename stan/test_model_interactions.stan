//November 2020 started by D. Loughnan 

// This Stan program aims to test for differences in budburst across woody species that dominate BC's forests 

// simple linar model do.bb = a + chill + forcing + photo

// Oct 12, 2021: I was not properly incorporating site as a dummy variable. Statistical Rethinking has categorical variables included as both dummy variables and indexing. Here I try the indexing approach

// Jan 7: trying to add interactions for both cues and for site

data {
  int<lower = 0 > N;
  int<lower = 0 > n_sp;
  //int<lower = 0 > n_site;
  int<lower = 1, upper=n_sp> sp[N];
  //int<lower = 1, upper = n_site>site[N];
  vector[N] chill; 
  vector <lower=0, upper = 1>[N] force;
  vector <lower=0, upper = 1>[N] photo;
  vector[N] bb; //response var
  // added upper and lower bounds to site dummy variable
  // vector <lower=0, upper = 1>[N] site2;
  // vector <lower=0, upper = 1>[N] site3;
  // vector <lower=0, upper = 1>[N] site4;
}

parameters {
  real mu_grand;
  
  real mu_force; 
  real mu_chill;
  real mu_photo;
  
  real mu_fc;
  real mu_cp;
  real mu_fp;

  vector[n_sp] a_sp;
  vector[n_sp] b_force;
  vector[n_sp] b_photo;
  vector[n_sp] b_chill;
  
  vector[n_sp] b_fc;
  vector[n_sp] b_cp;
  vector[n_sp] b_fp;
  
  real <lower=0> sigma_a;
  real <lower=0> sigma_force;
  real <lower=0> sigma_chill;
  real <lower=0> sigma_photo;
  
  real <lower=0> sigma_fc;
  real <lower=0> sigma_cp;
  real <lower=0> sigma_fp;
  real<lower=0> sigma_y; 
  
}

transformed parameters{

  real y_hat[N];
  for(i in 1:N){
		y_hat[i] = mu_grand + a_sp[sp[i]] + 
b_force[sp[i]] * force[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill[sp[i]] * chill[i] + 
		b_fc[sp[i]] * (force[i] * chill[i]) +
		b_cp[sp[i]] * (photo[i] * chill[i]) +
		b_fp[sp[i]] * (force[i] * photo[i]) ;
  }
}

model {
  // Priors. Make them flat
  mu_grand ~ normal(50, 5);
	mu_force ~ normal(0, 35); 
	mu_photo ~ normal(0, 35);
	mu_chill ~ normal(0, 35);
	// b_site2 ~ normal(0,5);
	// b_site3 ~ normal(0,5);
	// b_site4 ~ normal(0,5);

	mu_fp ~ normal(0,35);
	mu_fc ~ normal(0,35);
	mu_cp ~ normal(0,35);

	sigma_force ~ normal(1, 5); 
	sigma_photo ~ normal(1, 5); 
	sigma_chill ~ normal(1, 5);
	sigma_a ~ normal(0,5);
	sigma_y ~ normal(0,5);
	sigma_fp ~ normal(1, 10);
	sigma_fc ~ normal(1, 10);
	sigma_cp ~ normal(1, 10);
	
	//b_photo_ncp ~normal(0,1);
	b_force ~ normal(mu_force, sigma_force);
  b_photo ~ normal(mu_photo, sigma_photo); // bc still need this info
	b_chill ~ normal(mu_chill, sigma_chill);
	
	b_fc ~ normal(mu_fc, sigma_fc);
  b_fp ~ normal(mu_fp, sigma_fp); // bc still need this info
	b_cp ~ normal(mu_cp, sigma_cp);

	
	a_sp ~ normal(0,sigma_a);
  bb ~ normal(y_hat, sigma_y);
}

generated quantities{
   real ypred_new[N];
   
   for (i in 1:N)
   ypred_new[i] = normal_rng(y_hat[i], sigma_y);
}


	// sigma_b_inter_cs2 ~ normal(0, 10);
	// sigma_b_inter_fs2 ~ normal(0, 10);
	// sigma_b_inter_ps2 ~ normal(0, 10);
	// 
	// sigma_b_inter_cs3 ~ normal(0, 10);
	// sigma_b_inter_fs3~ normal(0, 10);
	// sigma_b_inter_ps3 ~ normal(0, 10);
	// 
	// sigma_b_inter_cs4 ~ normal(0, 10);
	// sigma_b_inter_fs4 ~ normal(0, 10);
	// sigma_b_inter_ps4 ~ normal(0, 10);

//   b_inter_fp ~ normal(mu_inter_fp, sigma_b_inter_fp);
//    b_inter_fc ~ normal(mu_inter_fc, sigma_b_inter_fc);
// 	 b_inter_pc ~ normal(mu_inter_pc, sigma_b_inter_pc);
// 
// 	 b_inter_cs2 ~ normal(mu_inter_cs2, sigma_b_inter_cs2);
//   b_inter_fs2 ~ normal(mu_inter_fs2, sigma_b_inter_fs2);
// 	 b_inter_ps2 ~ normal(mu_inter_ps2, sigma_b_inter_ps2);
// 	 
// 	 	 b_inter_cs3 ~ normal(mu_inter_cs3, sigma_b_inter_cs3);
//   b_inter_fs3 ~ normal(mu_inter_fs3, sigma_b_inter_fs3);
// 	 b_inter_ps3 ~ normal(mu_inter_ps3, sigma_b_inter_ps3);
// 	 
// 	 	 b_inter_cs4 ~ normal(mu_inter_cs4, sigma_b_inter_cs4);
//   b_inter_fs4 ~ normal(mu_inter_fs4, sigma_b_inter_fs4);
// 	 b_inter_ps4 ~ normal(mu_inter_ps4, sigma_b_inter_ps4);
	 
//    b_inter_fp_ncp ~ normal(0, 1);
//    b_inter_fc_ncp ~ normal(0, 1);
// 	 b_inter_pc_ncp ~ normal(0, 1);
// 
// 	 b_inter_cs2_ncp ~ normal(0, 1);
//   b_inter_fs2_ncp ~ normal(0, 1);
// 	 b_inter_ps2_ncp ~ normal(0, 1);
// 	 
// 	 	 b_inter_cs3_ncp ~ normal(0, 1);
//   b_inter_fs3_ncp ~ normal(0, 1);
// 	 b_inter_ps3_ncp ~ normal(0, 1);
// 	 
// 	 	 b_inter_cs4_ncp ~ normal(0, 1);
//   b_inter_fs4_ncp ~ normal(0, 1);
// 	 b_inter_ps4_ncp ~ normal(0, 1);
	


