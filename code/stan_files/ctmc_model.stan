data {

	int<lower=0> N; // Number of data points
	int<lower=0> P; // Number of unique pigs
	int<lower=0> C; // Number of covariates (not including intercept)

	real z[N]; // Response variable for Poisson analysis
	int pigID[N]; // Pig IDs
	real X[N, C]; // Design matrix
	real Z[N, P]; // Random effect matrix
	real offset[N];

} parameters{
	
	real beta0;
	real betas[C];
	real alphas[P]; // Random effects of pig

	real<lower=0> sigma_alpha;
	
} model {

	sigma_alpha ~ cauchy(0, 1);

	beta0 ~ normal(0, 5); 
	betas ~ normal(0, 5);
	alphas ~ normal(0, sigma_alpha);

	z ~ poisson(exp(beta0 + X*betas + Z*alphas));
	
}