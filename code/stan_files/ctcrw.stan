data{
    int<lower=0> N; // Number of data points ordered
    matrix[N, 2] obs_vals; // Hold the observation vector with two locations
    real deltas[N - 1]; // Hold the time steps between points.  N - 1 one of these

} parameters{
    
    // Model parameters
    real<lower=0> beta; // Autocorrelation parameter in the random walk
    real<lower=0> sigma; // Process error
    real<lower=0> sigma_obs; // Observation error
    real true_velx[N]; // The unobserved, true velocities in the x direction
    real true_vely[N]; // The unobserved, true velocities in the y direction

    real true_locx[N]; // The unobserved, true locations in the x direction
    real true_locy[N]; // The unobserved, true locations in the y direction


} model{
    
    // Define helper variables
    row_vector[2] mean_nextx; // Hold the mean 
    row_vector[2] mean_nexty; // Hold the mean 
    row_vector[2] full_trajx; // Holds the sampled trajectory
    row_vector[2] full_trajy; // Holds the sampled trajectory
    row_vector[2] prev_vectx; // Holds previous x loc, vel
    row_vector[2] prev_vecty; // Holds previous y loc, vel
    matrix[2, 2] Q;// Holds the covariance matrix
    matrix[2, 2] L; // Cholesky decomposition

    real zeta_var;
    real xi_var;
    real vl_covar;
    real delta;

    // Priors on the parameters
    beta ~ cauchy(0, 1); // Autocorrelation parameter...a finicky parameter
    sigma ~ cauchy(0, 1); // Process error
    sigma_obs ~ cauchy(0, 1); // Process error

    // Priors on initial velocities and locations
    true_velx[1] ~ normal(0, 1); // Tight priors to start with
    true_vely[1] ~ normal(0, 1);
    true_locx[1] ~ normal(obs_vals[1, 1], sigma_obs);
    true_locy[1] ~ normal(obs_vals[1, 2], sigma_obs);

    // Loop through the data and get likelihood for each point 
    for(i in 2:N){

        delta = deltas[i - 1];

        // Set observed locations
        prev_vectx[1] = true_locx[i - 1]; // x location
        prev_vecty[1] = true_locy[i - 1]; // y location

        // Set velocities
        prev_vectx[2] = true_velx[i - 1]; // x velocity
        prev_vecty[2] = true_vely[i - 1]; // y velocity

        // Calculate the mean vector
        mean_nextx[1] = prev_vectx[1] + prev_vectx[2]*(1 - exp(-beta*delta)) / beta;
        mean_nextx[2] = prev_vectx[2]*exp(-beta*delta);
        mean_nexty[1] = prev_vecty[1] + prev_vecty[2]*(1 - exp(-beta*delta)) / beta;
        mean_nexty[2] = prev_vecty[2]*exp(-beta*delta);

        // Calculate the covariance matrix Q: 4 X 4
        zeta_var = sigma^2*(1 - exp(-2*beta*delta)) / 2*beta;

        xi_var = (sigma^2 / beta^2)*(delta - (2/beta)*(1 - exp(-beta*delta)) + (1 / (2*beta))*(1 - exp(-2*beta*delta)));

        vl_covar = (sigma^2 / (2*beta^2))*(1 - 2*exp(-beta*delta) + exp(-2*beta*delta));

        
        Q[1, 1] = xi_var; Q[1, 2] = vl_covar;
        Q[2, 1] = vl_covar; Q[2, 2] = zeta_var;

        L = cholesky_decompose(Q); // This makes the sampling more efficient.
        // Pretty cool, to sample from a multivariable normal with MVM(mu, SIGMA)
        // You can write mu + LZ where Sigma = LL^T and Z are independent standard
        // normals

        full_trajx[1] = true_locx[i]; // x location
        full_trajx[2] = true_velx[i]; // Unobserved x vel

        full_trajy[1] = true_locy[i]; // observed y location
        full_trajy[2] = true_vely[i]; // unobserved y vel

        // Set "priors" on the true locations and velocities
        full_trajx ~ multi_normal_cholesky(mean_nextx, L);
        full_trajy ~ multi_normal_cholesky(mean_nexty, L);

        // Evaluate the likelihood for the observed locations (they are independent)
        obs_vals[i, 1] ~ normal(full_trajx[1], sigma_obs);
        obs_vals[i, 2] ~ normal(full_trajy[1], sigma_obs);


    } // End for

}
