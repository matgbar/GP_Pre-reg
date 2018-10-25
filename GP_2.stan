functions{
	//covariance function for main portion of the model
	matrix main_GP(
		int Nx,
		vector x,
		int Ny,
		vector y, 
		real alpha1,
		real alpha2,
		real alpha3,
		real rho1,
		real rho2,
		real rho3,
		real rho4,
		real rho5,
		real HR_f,
		real R_f){
					matrix[Nx, Ny] K1;
					matrix[Nx, Ny] K2;
					matrix[Nx, Ny] K3;
					matrix[Nx, Ny] Sigma;
	
					//periodic covariance that does not decay
					for(i in 1:Nx){
						for (j in 1:Ny){
							K1[i, j] = alpha1*exp(-square(x[i]-y[j])/2/square(rho1));
						}
					}
					
					//specifying first quasi-periodic process that incorporates heart rate
					for(i in 1:Nx){
						for(j in 1:Ny){
							K2[i, j] = alpha2*exp(-2*square(sin(pi()*fabs(x[i]-y[j])*HR_f))/square(rho2))*
							exp(-square(x[i]-y[j])/2/square(rho3));
						}
					}
					
					//specifying second quasi-periodic process that incorporates heart rate
					for(i in 1:Nx){
						for(j in 1:Ny){
							K3[i, j] = alpha3*exp(-2*square(sin(pi()*fabs(x[i]-y[j])*HR_f))/square(rho4))*
							exp(-2*square(sin(pi()*fabs(x[i]-y[j])*R_f))/square(rho5));
						}
					}
					Sigma = K1+K2+K3;
					return Sigma;
				}
	//function for posterior calculations
	vector post_pred_rng(
		real a1,
		real a2,
		real a3,
		real r1, 
		real r2,
		real r3,
		real r4,
		real r5,
		real HR,
		real R,
		real sn,
		int No,
		vector xo,
		int Np, 
		vector xp,
		vector yobs){
				matrix[No,No] Ko;
				matrix[Np,Np] Kp;
				matrix[No,Np] Kop;
				matrix[Np,No] Ko_inv_t;
				vector[Np] mu_p;
				matrix[Np,Np] Tau;
				matrix[Np,Np] L2;
				vector[Np] yp;
	
	//--------------------------------------------------------------------
	//Kernel Multiple GPs for observed data
	Ko = main_GP(No, xo, No, xo, a1, a2, a3, r1, r2, r3, r4, r5, HR, R);
	for(n in 1:No) Ko[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for predicted data
	Kp = main_GP(Np, xp, Np, xp, a1, a2, a3, r1, r2, r3, r4, r5, HR, R);
	for(n in 1:Np) Kp[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for observed and predicted cross 
	Kop = main_GP(No, xo, Np, xp, a1, a2, a3, r1, r2, r3, r4, r5, HR, R);
	
	//--------------------------------------------------------------------
	//Algorithm 2.1 of Rassmussen and Williams... 
	Ko_inv_t = Kop'/Ko;
	mu_p = Ko_inv_t*yobs;
	Tau=Kp-Ko_inv_t*Kop;
	L2 = cholesky_decompose(Tau);
	yp = mu_p + L2*rep_vector(normal_rng(0,1), Np);
	return yp;
	}
}

data { 
	int<lower=1> N1;
	int<lower=1> N2;
	vector[N1] X; 
	vector[N1] Y;
	vector[N2] Xp;
	real<lower=0> mu_HR;
	real<lower=0> mu_R;
	real<lower=0> sigma_HR;
	real<lower=0> sigma_R;
	real<lower=0> a1;
	real<lower=0> a2;
	real<lower=0> a3;
	real<lower=0> r1;
	real<lower=0> r2;
	real<lower=0> r3;
	real<lower=0> r4;
	real<lower=0> r5;
}

transformed data { 
	vector[N1] mu;
	for(n in 1:N1) mu[n] = 0;
}

parameters {
	real<lower = 0.8333, upper = 3.3333> HR;
	real<lower = 0.1667, upper = 0.5000> R;
	real<lower=1e-15> sigma_sq; //helps keep the matrices positive definite by ensuring sigma_sq never reaches 0
}

model{ 
	matrix[N1,N1] Sigma;
	matrix[N1,N1] L_S;
	
	//using GP function from above 
	Sigma = main_GP(N1, X, N1, X, a1, a2, a3, r1, r2, r3, r4, r5, HR, R);
	for(n in 1:N1) Sigma[n,n] += sigma_sq;
	
	L_S = cholesky_decompose(Sigma);
	Y ~ multi_normal_cholesky(mu, L_S);
	
	//priors for parameters
	sigma_sq ~ normal(0,2);
	HR ~ normal(mu_HR,sigma_HR);
	R ~ normal(mu_R, sigma_R);
}

generated quantities {
	vector[N2] Ypred = post_pred_rng(a1, a2, a3, r1, r2, r3, r4, r5, HR, R, sigma_sq, N1, X, N2, Xp, Y);
}

