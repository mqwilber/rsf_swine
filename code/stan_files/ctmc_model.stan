data {

	int<lower=0> N; // Number of data points
	int<lower=0> P; // Number of unique pigs
	int<lower=0> C; // Number of covariates (not including intercept)

	int z[N]; // Response variable for Poisson analysis
	int pigID[N]; // Pig IDs
	matrix[N, C] X; // Design matrix
	matrix[N, P] Z; // Random effect matrix
	vector[N] offset;

} parameters{
	
	real beta0;
	vector[C] betas;
	vector[P] alphas; // Random effects of pig on intercept

	real<lower=0> sigma_alpha;
	
} model {

	sigma_alpha ~ cauchy(0, 1);

	// beta0 ~ normal(0, 5); 
	// betas ~ normal(0, 5);
	alphas ~ normal(0, sigma_alpha);

	z ~ poisson(exp(beta0 + X*betas + Z*alphas + offset));
	
} 

// generated quantities {
// 	vector[C] betas;
// 	betas = R_ast_inverse * thetas;
// }