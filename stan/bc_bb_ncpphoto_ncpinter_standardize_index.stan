//This is code adapted from the statistical rethinking example in section 5.3

data {
  int<lower = 0 > N;
  int<lower = 0 > n_sp;
  int<lower = 0 > n_site;
  int<lower = 1, upper= n_sp> sp[N]; 
  int<lower = 1, upper = n_site> site[N];
  vector[N] chill; 
  //vector[N] site;
  vector[N] force;
  vector[N] photo;
  vector[N] bb;
  
  //int site[N];
}

transformed data { 
  // standardizing the forcing, photo, and chill portions
  vector[N] force_std;
  vector[N] photo_std;
  vector[N] chill_std;
  //vector[N] site_std;
  
// 9 interaction terms
  vector[N] inter_fp;       // vector[N] inter_fp  = force .* photo;
  //vector[N] inter_fs;       // force*site         
  //vector[N] inter_ps;       // photoperiod by site   
  vector[N] inter_fc;       // force*chill            
  vector[N] inter_pc;       // photo*chill            
  //vector[N] inter_sc;       // chill * site    
  
  force_std = (force-mean(force))/(2*sd(force));
  photo_std = (photo-mean(photo))/(2*sd(photo));
  chill_std = (chill-mean(chill))/(2*sd(chill));   
  //site_std = (site-mean(site))/(2*sd(site));   

  inter_fp    = force_std .* photo_std;
  inter_fc    = force_std .* chill_std;
  inter_pc    = photo_std .* chill_std;
  //inter_fs    = force_std .* site_std;
  //inter_ps    = photo_std .* site_std;
  //inter_cs    = chill_std .* site_std;  
}

parameters {
  real mu_grand;
  real mu_a;
  real mu_force; 
  real mu_chill;
  real mu_photo;
  
  real mu_inter_fp;
  real mu_inter_fc;
  real mu_inter_pc;
  
  vector[n_sp] a_sp;
  
  vector[n_sp] b_force;
  vector[n_sp] b_chill;
  vector[n_sp] b_photo_ncp; //
  vector[n_sp] b_inter_fp_ncp;
  //vector[n_sp] b_inter_fs_ncp;
  //vector[n_sp] b_inter_ps_ncp;
  vector[n_sp] b_inter_fc_ncp;
  vector[n_sp] b_inter_pc_ncp;
  //vector[n_sp] b_inter_sc_ncp;
    
  real<lower=0> sigma_force;
  real<lower=0> sigma_chill;
  real<lower=0> sigma_photo;
  real<lower=0> sigma_b_inter_fp;
  //real<lower=0> sigma_b_inter_fs;
  //real<lower=0> sigma_b_inter_ps;
  real<lower=0> sigma_b_inter_fc;
  real<lower=0> sigma_b_inter_pc;
  
  real<lower=0> sigma_a;
  real<lower=0> sigma_y;
  vector[n_site] a_site;
}

transformed parameters {
  vector[n_sp] b_photo;
  
  vector[n_sp] b_inter_fp;
  //vector[n_sp] b_inter_fs;
  //vector[n_sp] b_inter_ps;
  vector[n_sp] b_inter_fc;
  vector[n_sp] b_inter_pc;
  //vector[n_sp] b_inter_sc;
  
  vector[N] mu;
  vector[N] y_hat;

  mu = a_site[site];
  
  b_photo = mu_photo + sigma_photo * b_photo_ncp;
  
  b_inter_fp = mu_inter_fp + sigma_b_inter_fp*b_inter_fp_ncp;
  //b_inter_fs = mu_inter_fs + sigma_b_inter_fs*b_inter_fs_ncp;
 // b_inter_ps = mu_inter_ps + sigma_b_inter_ps*b_inter_ps_ncp;
  b_inter_fc = mu_inter_fc + sigma_b_inter_fc*b_inter_fc_ncp;
  b_inter_pc = mu_inter_pc + sigma_b_inter_pc*b_inter_pc_ncp;
  //b_inter_sc = mu_inter_sc + sigma_b_inter_sc*b_inter_sc_ncp;
  
  for(i in 1:N){
	y_hat[i] = mu_grand + a_sp[sp[i]]  + a_site[site[i]]+
    b_force[sp[i]] * force_std[i] + 
		b_photo[sp[i]] * photo_std[i] + 
		b_chill[sp[i]] * chill_std[i] +
		b_inter_fp[sp[i]] * inter_fp[i] +
		//b_inter_fs[sp[i]] * inter_fs[i] + 
		//b_inter_ps[sp[i]] * inter_ps[i] +
		b_inter_fc[sp[i]] * inter_fc[i] +
		b_inter_pc[sp[i]] * inter_pc[i] 
		;  }
}
model {
  mu_grand ~ normal(50,5);
  mu_force ~ normal(0, 50); // 100 = 3 months on either side. Narrow down to 35
	mu_photo ~ normal(0, 35);
	mu_chill ~ normal(0, 35);
	
	mu_inter_fp ~ normal(0,35);
	mu_inter_fc ~ normal(0,35);
	mu_inter_pc ~ normal(0,35);
	
  sigma_force ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	sigma_photo ~ normal(0, 10); 
	sigma_chill ~ normal(0, 30);
  
  sigma_b_inter_fp ~ normal(0, 10);
	//sigma_b_inter_fs ~ normal(0, 10);
	//sigma_b_inter_ps ~ normal(0, 10);
	sigma_b_inter_fc ~ normal(0, 10);	
	sigma_b_inter_pc ~ normal(0, 10);	
	
  b_photo_ncp ~normal(0,1);
	b_force ~ normal(mu_force, sigma_force);
  //b_photo ~ normal(mu_photo, sigma_photo); // bc still need this info
	b_chill ~ normal(mu_chill, sigma_chill);
  
  b_inter_fp_ncp ~ normal(0, 1); 
  //b_inter_fs_ncp ~ normal(0, 1);
	//b_inter_ps_ncp ~ normal(0, 1);		
	b_inter_fc_ncp ~ normal(0, 1);
	b_inter_pc_ncp ~ normal(0, 1);	
	
  bb ~ normal(y_hat, sigma_y);
  sigma_y ~ uniform(0, 50);
  a_site ~ normal(0, 0.5);
  a_sp ~ normal(mu_a,sigma_a);
}
