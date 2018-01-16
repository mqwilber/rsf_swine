data{
    int<lower=0> N; // Number of data points ordered
    matrix[N, 4] obs_vals; // Hold the observation vector with two locations
    real deltas[N - 1]; // Hold the time steps between points.  N - 1 one of these
    
} parameters{
    
    // Model parameters
    real<lower=0> beta; // Autocorrelation parameter in the random walk
    real<lower=0> sigma; // Process error

} model{
    
    row_vector[4] mean_next; // Hold the mean 
    matrix[4, 4] Q;// Holds the covariance matrix
    matrix[4, 4] L; // Cholesky decomposition
    real zeta_var;
    real xi_var;
    real vl_covar;

    // Priors on the parameters
    beta ~ cauchy(0, 1); // Autocorrelation parameter
    sigma ~ cauchy(0, 1);


    // Loop through the data and get likelihood for each point 
    for(i in 2:N){

        row_vector[4] prev_vect;
        real delta;
        delta = deltas[i - 1];

        prev_vect = obs_vals[i - 1];

        // Calculate the mean vector
        mean_next[1] = prev_vect[1] + prev_vect[2]*(1 - exp(-beta*delta)) / beta;
        mean_next[2] = prev_vect[2]*exp(-beta*delta);
        mean_next[3] = prev_vect[3] + prev_vect[4]*(1 - exp(-beta*delta)) / beta;
        mean_next[4] = prev_vect[4]*exp(-beta*delta);

        // Calculate the covariance matrix Q: 4 X 4
        zeta_var = sigma^2*(1 - exp(-2*beta*delta)) / 2*beta;

        xi_var = (sigma^2 / beta^2)*(delta - (2/beta)*(1 - exp(-beta*delta)) + (1 / (2*beta))*(1 - exp(-2*beta*delta)));

        vl_covar = (sigma^2 / (2*beta^2))*(1 - 2*exp(-beta*delta) + exp(-2*beta*delta));

        
        Q[1, 1] = xi_var; Q[1, 2] = vl_covar; Q[1, 3] = 0; Q[1, 4] = 0;
        Q[2, 1] = vl_covar; Q[2, 2] = zeta_var; Q[2, 3] = 0; Q[2, 4] = 0;
        Q[3, 1] = 0; Q[3, 2] = 0; Q[3, 3] = xi_var; Q[3, 4] = vl_covar;
        Q[4, 1] = 0; Q[4, 2] = 0; Q[4, 3] = vl_covar; Q[4, 4] = zeta_var;

        L = cholesky_decompose(Q); // This makes the sampling more efficient...thought not much more lol.

        // Get the likelihood
        obs_vals[i] ~ multi_normal_cholesky(mean_next, L);

    } // End for

}
