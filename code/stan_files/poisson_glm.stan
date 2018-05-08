data {
	int<lower=0> N;
	int<lower=0> p;
	matrix[N, p] X; // Design matrix
	int z[N]; // Response variable
	vector[N] tau; // Weighting time

} parameters {

	vector[p] beta; // Design matrix without an intercept
	real b0; // Intercept

} model{

	vector[N] lambda;
	beta ~ normal(0, 5);

	lambda = exp(log(tau) + X*beta);
	z ~ poisson(lambda);

}