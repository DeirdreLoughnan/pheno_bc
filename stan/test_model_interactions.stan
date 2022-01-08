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
  vector <lower=0, upper = 1>[N] site2;
  vector <lower=0, upper = 1>[N] site3;
  vector <lower=0, upper = 1>[N] site4;
}

transformed data { 
  // standardizing the forcing, photo, and chill portions
  // vector[N] force_std;
  // vector[N] photo_std;
  // vector[N] chill_std;

// 9 interaction terms
  vector[N] inter_fp;       // vector[N] inter_fp  = force .* photo;
  vector[N] inter_fc;       // force*chill
  vector[N] inter_pc;       // photo*chill
 
  vector[N] inter_cs2;       // chill * site
  vector[N] inter_fs2;       // force*site
  vector[N] inter_ps2;       // photoperiod by site

  vector[N] inter_cs3;       // chill * site
  vector[N] inter_fs3;       // force*site
  vector[N] inter_ps3;       // photoperiod by site

  vector[N] inter_cs4;       // chill * site
  vector[N] inter_fs4;       // force*site
  vector[N] inter_ps4;       // photoperiod by site

  // force_std = (force-mean(force))/(2*sd(force));
  // photo_std = (photo-mean(photo))/(2*sd(photo));
  // chill_std = (chill-mean(chill))/(2*sd(chill));

//what is the dot doing? bc vectors of length n?
  inter_fp    = force .* photo;
  inter_fc    = force .* chill;
  inter_pc    = photo .* chill;
  
  inter_fs2    = force .* site2;
  inter_ps2   = photo .* site2;
  inter_cs2    = chill .* site2;
  
  inter_fs3    = force .* site3;
  inter_ps3    = photo .* site3;
  inter_cs3    = chill .* site3;
  
  inter_fs4    = force .* site4;
  inter_ps4    = photo .* site4;
  inter_cs4    = chill .* site4;
}

// The parameters accepted by the model. 

parameters {
  real mu_grand;
  //real mu_a;
  real mu_force; 
  real mu_chill;
  real mu_photo;
  real mu_inter_fp;
  real mu_inter_fc;
  real mu_inter_pc;
  
  real mu_inter_fs2;
  real mu_inter_ps2;
  real mu_inter_cs2;
  
  real mu_inter_fs3;
  real mu_inter_ps3;
  real mu_inter_cs3;
  
  real mu_inter_fs4;
  real mu_inter_ps4;
  real mu_inter_cs4;
  
  vector[n_sp] a_sp;
  vector[n_sp] b_force;
  //vector[n_sp] b_photo;
  vector[n_sp] b_chill;
  
  
  vector[n_sp] b_inter_fp;
  vector[n_sp] b_inter_fc;
  vector[n_sp] b_inter_pc;

  vector[n_sp] b_inter_cs2;
  vector[n_sp] b_inter_fs2;
  vector[n_sp] b_inter_ps2;

  vector[n_sp] b_inter_cs3;
  vector[n_sp] b_inter_fs3;
  vector[n_sp] b_inter_ps3;

  vector[n_sp] b_inter_cs4;
  vector[n_sp] b_inter_fs4;
  vector[n_sp] b_inter_ps4;

  
  real b_site2;
  real b_site3;
  real b_site4;
  
  vector[n_sp] b_photo_ncp; //
  
  //start without the interactions being ncp
  // vector[n_sp] b_inter_fp_ncp;
  // vector[n_sp] b_inter_fc_ncp;
  // vector[n_sp] b_inter_pc_ncp;
  // 
  // vector[n_sp] b_inter_cs2_ncp;
  // vector[n_sp] b_inter_fs2_ncp;
  // vector[n_sp] b_inter_ps2_ncp;
  // 
  // vector[n_sp] b_inter_cs3_ncp;
  // vector[n_sp] b_inter_fs3_ncp;
  // vector[n_sp] b_inter_ps3_ncp;
  // 
  // vector[n_sp] b_inter_cs4_ncp;
  // vector[n_sp] b_inter_fs4_ncp;
  // vector[n_sp] b_inter_ps4_ncp;
  
  real <lower=0> sigma_a;
  real <lower=0> sigma_force;
  real <lower=0> sigma_chill;
  real <lower=0> sigma_photo;
  real<lower=0> sigma_b_inter_fp;
  real<lower=0> sigma_b_inter_fc;
  real<lower=0> sigma_b_inter_pc;

  real<lower=0> sigma_b_inter_cs2;
  real<lower=0> sigma_b_inter_fs2;
  real<lower=0> sigma_b_inter_ps2;

  real<lower=0> sigma_b_inter_cs3;
  real<lower=0> sigma_b_inter_fs3;
  real<lower=0> sigma_b_inter_ps3;

  real<lower=0> sigma_b_inter_cs4;
  real<lower=0> sigma_b_inter_fs4;
  real<lower=0> sigma_b_inter_ps4;

  real<lower=0> sigma_y; 
}

transformed parameters{
   vector[n_sp] b_photo;
// 
//   vector[n_sp] b_inter_fp;
//   vector[n_sp] b_inter_fc;
//   vector[n_sp] b_inter_pc;
// 
//   vector[n_sp] b_inter_cs2;
//   vector[n_sp] b_inter_fs2;
//   vector[n_sp] b_inter_ps2;
// 
//   vector[n_sp] b_inter_cs3;
//   vector[n_sp] b_inter_fs3;
//   vector[n_sp] b_inter_ps3;
// 
//   vector[n_sp] b_inter_cs4;
//   vector[n_sp] b_inter_fs4;
//   vector[n_sp] b_inter_ps4;

  vector[N] y_hat;
  
  b_photo = mu_photo + sigma_photo * b_photo_ncp;

  // b_inter_fp = mu_inter_fp + sigma_b_inter_fp*b_inter_fp_ncp;
  // b_inter_fc = mu_inter_fc + sigma_b_inter_fc*b_inter_fc_ncp;
  // b_inter_pc = mu_inter_pc + sigma_b_inter_pc*b_inter_pc_ncp;
  // 
  // b_inter_cs2 = mu_inter_cs2 + sigma_b_inter_cs2*b_inter_cs2_ncp;
  // b_inter_fs2 = mu_inter_fs2 + sigma_b_inter_fs2*b_inter_fs2_ncp;
  // b_inter_ps2 = mu_inter_ps2 + sigma_b_inter_ps2*b_inter_ps2_ncp;
  // 
  //   b_inter_cs3 = mu_inter_cs3 + sigma_b_inter_cs3*b_inter_cs3_ncp;
  // b_inter_fs3 = mu_inter_fs3 + sigma_b_inter_fs3*b_inter_fs3_ncp;
  // b_inter_ps3 = mu_inter_ps3 + sigma_b_inter_ps3*b_inter_ps3_ncp;
  // 
  //   b_inter_cs4 = mu_inter_cs4 + sigma_b_inter_cs4*b_inter_cs4_ncp;
  // b_inter_fs4 = mu_inter_fs4 + sigma_b_inter_fs4*b_inter_fs4_ncp;
  // b_inter_ps4 = mu_inter_ps4 + sigma_b_inter_ps4*b_inter_ps4_ncp;

  for(i in 1:N){
		y_hat[i] = mu_grand + a_sp[sp[i]] + 
	  b_site2 * site2[i] +
	  b_site3 * site3[i] +
	  b_site4 * site4[i] + 
		b_force[sp[i]] * force[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill[sp[i]] * chill[i] + 
		b_inter_fp[sp[i]] * inter_fp[i] +
		b_inter_fc[sp[i]] * inter_fc[i] +
		b_inter_pc[sp[i]] * inter_pc[i]
		+ b_inter_cs2[sp[i]] * inter_cs2[i] +
		b_inter_fs2[sp[i]] * inter_fs2[i] +
		b_inter_ps2[sp[i]] * inter_ps2[i] +
		b_inter_cs3[sp[i]] * inter_cs3[i] +
		b_inter_fs3[sp[i]] * inter_fs3[i] +
		b_inter_ps3[sp[i]] * inter_ps3[i] +
		b_inter_cs4[sp[i]] * inter_cs4[i] +
		b_inter_fs4[sp[i]] * inter_fs4[i] +
		b_inter_ps4[sp[i]] * inter_ps4[i]
		;
  }
}

model {
  // Priors. Make them flat
  mu_grand ~ normal(50, 5);
	mu_force ~ normal(0, 35); 
	mu_photo ~ normal(0, 35);
	mu_chill ~ normal(0, 35);
	b_site2 ~ normal(0,5);
	b_site3 ~ normal(0,5);
	b_site4 ~ normal(0,5);

	mu_inter_fp ~ normal(0,35);
	mu_inter_fc ~ normal(0,35);
	mu_inter_pc ~ normal(0,35);
	
	mu_inter_fs2 ~ normal(0,35);
	mu_inter_ps2 ~ normal(0,35);
	mu_inter_cs2 ~ normal(0,35);
	mu_inter_fs3 ~ normal(0,35);
	mu_inter_ps3 ~ normal(0,35);
	mu_inter_cs2 ~ normal(0,35);
	mu_inter_fs4 ~ normal(0,35);
	mu_inter_ps4 ~ normal(0,35);
	mu_inter_cs4 ~ normal(0,35);
	
	sigma_force ~ normal(1, 5); 
	sigma_photo ~ normal(1, 5); 
	sigma_chill ~ normal(1, 5);
	sigma_a ~ normal(0,5);
	sigma_y ~ normal(0,5);
	sigma_b_inter_fp ~ normal(1, 10);
	sigma_b_inter_fc ~ normal(1, 10);
	sigma_b_inter_pc ~ normal(1, 10);

	sigma_b_inter_cs2 ~ normal(0, 10);
	sigma_b_inter_fs2 ~ normal(0, 10);
	sigma_b_inter_ps2 ~ normal(0, 10);
	
	sigma_b_inter_cs3 ~ normal(0, 10);
	sigma_b_inter_fs3~ normal(0, 10);
	sigma_b_inter_ps3 ~ normal(0, 10);
	
	sigma_b_inter_cs4 ~ normal(0, 10);
	sigma_b_inter_fs4 ~ normal(0, 10);
	sigma_b_inter_ps4 ~ normal(0, 10);
	
	
	b_photo_ncp ~normal(0,1);
	b_force ~ normal(mu_force, sigma_force);
  //b_photo ~ normal(mu_photo, sigma_photo); // bc still need this info
	b_chill ~ normal(mu_chill, sigma_chill);

  b_inter_fp ~ normal(mu_inter_fp, sigma_b_inter_fp);
   b_inter_fc ~ normal(mu_inter_fc, sigma_b_inter_fc);
	 b_inter_pc ~ normal(mu_inter_pc, sigma_b_inter_pc);

	 b_inter_cs2 ~ normal(mu_inter_cs2, sigma_b_inter_cs2);
  b_inter_fs2 ~ normal(mu_inter_fs2, sigma_b_inter_fs2);
	 b_inter_ps2 ~ normal(mu_inter_ps2, sigma_b_inter_ps2);
	 
	 	 b_inter_cs3 ~ normal(mu_inter_cs3, sigma_b_inter_cs3);
  b_inter_fs3 ~ normal(mu_inter_fs3, sigma_b_inter_fs3);
	 b_inter_ps3 ~ normal(mu_inter_ps3, sigma_b_inter_ps3);
	 
	 	 b_inter_cs4 ~ normal(mu_inter_cs4, sigma_b_inter_cs4);
  b_inter_fs4 ~ normal(mu_inter_fs4, sigma_b_inter_fs4);
	 b_inter_ps4 ~ normal(mu_inter_ps4, sigma_b_inter_ps4);
	 
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
	
	
	a_sp ~ normal(0,sigma_a);
  bb ~ normal(y_hat, sigma_y);
}

generated quantities{
   real ypred_new[N];
   
   for (i in 1:N)
   ypred_new[i] = normal_rng(y_hat[i], sigma_y);
}



