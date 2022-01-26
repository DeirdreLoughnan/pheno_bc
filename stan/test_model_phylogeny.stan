
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
  int<lower=1> N;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  vector[N] bb; 		// response
  vector[N] force; 	// predictor
  vector[N] chill; 	// predictor
  vector[N] photo; 	// predictor
  matrix[n_sp,n_sp]Vphy;     // phylogeny
  
  // Priors
  real a_z_prior_mu;
  real a_z_prior_sigma;
  real lam_interceptsa_prior_alpha;
  real lam_interceptsa_prior_beta;
  real sigma_interceptsa_prior_mu;
  real sigma_interceptsa_prior_sigma;
  // real b_zf_prior_mu;
  // real b_zf_prior_sigma;
  // real lam_interceptsbf_prior_alpha;
  // real lam_interceptsbf_prior_beta;
  // real sigma_interceptsbf_prior_mu;
  // real sigma_interceptsbf_prior_sigma;
  real mu_grand_prior_mu;
  real mu_grand_prior_sigma;
  real mu_force_prior_mu;
  real mu_force_prior_sigma;
  real mu_chill_prior_mu;
  real mu_chill_prior_sigma;
  real mu_photo_prior_mu;
  real mu_photo_prior_sigma;
  real sigma_force_prior_mu;
  real sigma_force_prior_sigma;
  real sigma_chill_prior_mu;
  real sigma_chill_prior_sigma;
  real sigma_photo_prior_mu;
  real sigma_photo_prior_sigma;
  real sigma_a_prior_mu;
  real sigma_a_prior_sigma;
  real sigma_y_mu_prior;
  real sigma_y_mu_sigma;  
}

parameters {
  real mu_grand;
  real mu_force;
  real mu_chill;
  real mu_photo;
  
  vector[n_sp] a;
  vector[n_sp] b_force;
  vector[n_sp] b_chill;
  vector[n_sp] b_photo;
  
  real a_z;
  
  real<lower=0> sigma_y;    
  real <lower=0> sigma_force;
  real <lower=0> sigma_chill;
  real <lower=0> sigma_photo;

  real<lower=0, upper=1> lam_interceptsa;       
  real<lower=0> sigma_interceptsa;
  // real<lower=0, upper=1> lam_interceptsbf;       
  // real<lower=0> sigma_interceptsbf;    
	}

model {
       real yhat[N];
       	for(i in 1:N){
		yhat[i] = mu_grand + a[sp[i]] + 
    b_force[sp[i]] * force[i] + 
		b_photo[sp[i]] * photo[i] + 
		b_chill[sp[i]] * chill[i];
			     	}
  a ~ multi_normal(rep_vector(a_z,n_sp), lambda_vcv(Vphy, lam_interceptsa, sigma_interceptsa)); 
  
  mu_grand ~ normal(mu_grand_prior_mu, mu_grand_prior_sigma);
  b_force ~ normal(mu_force_prior_mu, mu_force_prior_sigma); 
  b_chill ~ normal(mu_chill_prior_mu, mu_chill_prior_sigma); 
  b_photo ~ normal(mu_photo_prior_mu, mu_photo_prior_sigma); 
  bb ~ normal(yhat, sigma_y);
  
  sigma_force ~ normal(sigma_force_prior_mu, sigma_force_prior_sigma);
  sigma_chill ~ normal(sigma_chill_prior_mu, sigma_chill_prior_sigma);
  sigma_photo ~ normal(sigma_photo_prior_mu, sigma_photo_prior_sigma);
  // Priors
  a_z ~ normal(a_z_prior_mu, a_z_prior_sigma);
  lam_interceptsa ~ beta(lam_interceptsa_prior_alpha, lam_interceptsa_prior_beta);
  sigma_interceptsa ~ normal(sigma_interceptsa_prior_mu, sigma_interceptsa_prior_sigma);
  //b_zf ~ normal(b_zf_prior_mu, b_zf_prior_sigma);
  //lam_interceptsbf ~ beta(lam_interceptsbf_prior_alpha, lam_interceptsbf_prior_beta);
  //sigma_interceptsbf ~ normal(sigma_interceptsbf_prior_mu, sigma_interceptsbf_prior_sigma);
  sigma_y ~ normal(sigma_y_mu_prior, sigma_y_mu_sigma);

}

