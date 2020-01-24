// =================================================================================
//Add gaussian noise generation
// =================================================================================

/*
#include <cmath>
#include <limits>

//double generateGaussianNoise(double mu, double sigma, double seedVal)
double generateGaussianNoise(double mu, double sigma)


{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;	
	
	double u1, u2;
	do
	 {
	   //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	   //  std::srand(seedVal);
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

*/


// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
  // Variable 0
  set_variable_name				(0,"u");
  set_variable_type				(0,SCALAR);
  set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

	
  set_dependencies_value_term_RHS(0, "u,mu,phi,grad(phi)");	
  set_dependencies_gradient_term_RHS(0, "grad(u)");

	

  // Variable 1
  set_variable_name				(1,"phi");
  set_variable_type				(1,SCALAR);
  set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "phi,mu");
  set_dependencies_gradient_term_RHS(1, "");

	
  // Variable 2
  set_variable_name				(2,"mu");
  set_variable_type				(2,SCALAR);
  set_variable_equation_type		(2,AUXILIARY);

  set_dependencies_value_term_RHS(2, "phi,u,grad(phi)");
  //  set_dependencies_value_term_RHS(2, "phi,u");
  set_dependencies_gradient_term_RHS(2, "grad(phi)");
	
  
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
						dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

  // --- Getting the values and derivatives of the model variables ---

  // The temperature and its derivatives
  scalarvalueType u = variable_list.get_scalar_value(0);
  scalargradType ux = variable_list.get_scalar_gradient(0);

  // The order parameter and its derivatives
  scalarvalueType phi = variable_list.get_scalar_value(1);
  scalargradType phix = variable_list.get_scalar_gradient(1);

  // The order parameter chemical potential and its derivatives
  scalarvalueType mu = variable_list.get_scalar_value(2);

  // --- Setting the expressions for the terms in the governing equations ---

  //prec-conditioned PF
  scalarvalueType TANPF ; 

  for (unsigned i=0; i< phi.n_array_elements;i++) {
    TANPF[i]= std::tanh(phi[i]/std::sqrt(2.0)) ;
  }
 
  // The azimuthal angle
  scalarvalueType theta,gradMAG;
  for (unsigned i=0; i< phi.n_array_elements;i++) {
    theta[i] = std::atan2(phix[1][i],phix[0][i]);
    gradMAG[i]= std::sqrt(phix[0][i]*phix[0][i] + phix[1][i]*phix[1][i]);
  }
  
 // Anisotropic gradient energy coefficient, its derivative and square. W is non-dimenisonalized so W0 is 1 and not used here
  scalarvalueType W ;
  for (unsigned i=0; i< phi.n_array_elements;i++){
    if (q_point_loc[0][i]<PART) {
	W[i] = ((1.0)+(epsilonM)*std::cos((mult)*(theta[i]-(theta0))));
      }
    else
      {
	W[i] =	((1.0)+(epsilonM)*std::cos((mult)*(theta[i]-(theta01))));
      }
  }
  
  

  scalarvalueType noise1;
  for (unsigned i=0; i< phi.n_array_elements;i++) {
    //square root
    std::srand(std::sqrt(q_point_loc[0][i]*q_point_loc[0][i] + q_point_loc[1][i]*q_point_loc[1][i]));
    //square
    //  std::srand((q_point_loc[0][i]*q_point_loc[0][i] + q_point_loc[1][i]*q_point_loc[1][i]));    
    noise1[i]=(0.01)*(1-TANPF[i]*TANPF[i])*((double)(std::rand() % 100)/100.0 - 0.5 )*std::sqrt(userInputs.dtValue) ;
  }

   double TIME= this->currentTime;
  
  //quadrature point along directional solidfication
  scalarvalueType TAU,QQ,UTAU;

  for (unsigned i=0; i< phi.n_array_elements;i++){
    TAU[i]= (W[i]*W[i])*( 1.0 - ((1-ke)/LT)*(q_point_loc[1][i]-VV*TIME)) ;
     if (TAU[i]==0) { TAU[i]=ke*W[i]*W[i];}
    //UTAU[i] = 1.0+ke-(1.0-ke)*phi[i];
    UTAU[i] = 1.0+ke-(1.0-ke)*TANPF[i];
  }

 for (unsigned i=0; i< phi.n_array_elements;i++) {
   if (1/*1-TANPF[i]>=1.0e-06*/)
     QQ[i] = 1.0-TANPF[i]; 
   else {
     QQ[i] = 0; 
   }
     
 }

  
 
  scalarvalueType tau = TAU;
  //Define Antiflux current, gaussian current
  scalargradType JATF ;
  for (unsigned i=0; i< phi.n_array_elements;i++) {

    if (gradMAG[i]==0) {JATF[0][i]=0;JATF[1][i]=0;}
    else {
      if (1 /*std::abs(1-TANPF[i]*TANPF[i]) >=1.0e-06*/ ) {
	JATF[0][i] = (1+(1-ke)*u[i])*(0.5*(1-TANPF[i]*TANPF[i]))*(mu[i]/tau[i])*(phix[0][i]/gradMAG[i]);   
	JATF[1][i] = (1+(1-ke)*u[i])*(0.5*(1-TANPF[i]*TANPF[i]))*(mu[i]/tau[i])*(phix[1][i]/gradMAG[i]);

	//JATF[0][i] = (1+(1-ke)*u[i])*(0.5*(1-TANPF[i]*TANPF[i]))*(mu[i]/tau[i])*(1.0/gradMAG[i]);   
	//JATF[1][i] = (1+(1-ke)*u[i])*(0.5*(1-TANPF[i]*TANPF[i]))*(mu[i]/tau[i])*(1.0/gradMAG[i]);

      }
      else {
	JATF[0][i]=0;
	JATF[1][i]=0;
      }
    }
  }


    // Define required equations

  scalarvalueType eq_u = (u + constV(1.0/std::sqrt(2.0))*(constV(1.0)-TANPF*TANPF)*(constV(1.0)+constV(1.0-ke)*u)*(mu/tau)*(constV(userInputs.dtValue)/UTAU)+ noise1*constV(userInputs.dtValue)/UTAU )  ;
  //scalarvalueType eq_u = (u+ (constV(1.0)+constV(1.0-ke)*u)*(mu/tau/UTAU)*constV(userInputs.dtValue) )  ;
  
  scalargradType eqx_u = (constV(-D_tild*userInputs.dtValue)*QQ*ux/UTAU) + (constV(-1.0)*JATF)*(constV(userInputs.dtValue)/UTAU)  ;
  //scalargradType eqx_u = ((constV(-D_tild)*QQ*ux*(constV(userInputs.dtValue)/UTAU)) + (constV(-1.0)*JATF*phix)*(constV(userInputs.dtValue)/UTAU))  ;
			  
  //scalarvalueType eq_phi = (phi+constV(userInputs.dtValue)*mu/tau);
  scalarvalueType eq_phi = (phi+ constV(userInputs.dtValue)*mu/tau);
 
  // --- Submitting the terms for the governing equations ---

  // Terms for the equation to evolve the concentration
  variable_list.set_scalar_value_term_RHS(0,eq_u);
  variable_list.set_scalar_gradient_term_RHS(0,eqx_u);

  // Terms for the equation to evolve the order parameter
  variable_list.set_scalar_value_term_RHS(1,eq_phi);


}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
						   dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

  // --- Getting the values and derivatives of the model variables ---

  // The concentration  and its derivatives
  scalarvalueType u = variable_list.get_scalar_value(0);

  // The order parameter and its derivatives
  scalarvalueType phi = variable_list.get_scalar_value(1);
  scalargradType phix = variable_list.get_scalar_gradient(1);
  double TIME= this->currentTime;


  //prec-conditioned PF
  scalarvalueType TANPF ;
    for (unsigned i=0; i< phi.n_array_elements;i++) {
    TANPF[i]= std::tanh(phi[i]/std::sqrt(2.0) ) ;
  }
 
  scalarvalueType THERM;
   for (unsigned i=0; i< u.n_array_elements;i++) {
     THERM[i]= (1.0/LT)*(q_point_loc[1][i]-VV*TIME)  ;
  }

  // Derivative of the free energy density with respect to phi
   scalarvalueType f_phi= -(TANPF*constV(std::sqrt(2.0))-constV(std::sqrt(2.0))*constV(lam)*(constV(1.0)-TANPF*TANPF)*(u+THERM));

   
  // The azimuthal angle
   scalarvalueType theta,gradMAG;
  for (unsigned i=0; i< phi.n_array_elements;i++) {
    theta[i] = std::atan2(phix[1][i],phix[0][i]);
    gradMAG[i]= std::sqrt(phix[0][i]*phix[0][i] + phix[1][i]*phix[1][i]);
  }

  // Anisotropic gradient energy coefficient, its derivative and square
  // scalarvalueType W = constV(1.0)*(constV(1.0)+constV(epsilonM)*std::cos(constV(mult)*(theta-constV(theta0))));
  //scalarvalueType W_theta = constV(-1.0)*(constV(epsilonM)*constV(mult)*std::sin(constV(mult)*(theta-constV(theta0))));

 scalarvalueType W ;
 scalarvalueType W_theta;
  for (unsigned i=0; i< phi.n_array_elements;i++) {
    if (q_point_loc[0][i]<PART) {
      W[i] = ((1.0)+(epsilonM)*std::cos((mult)*(theta[i]-(theta0))));
      W_theta[i] = (-1.0)*((epsilonM)*(mult)*std::sin((mult)*(theta[i]-(theta0))));
    }
    else {
      W[i] = ((1.0)+(epsilonM)*std::cos((mult)*(theta[i]-(theta01))));
      W_theta[i] = (-1.0)*((epsilonM)*(mult)*std::sin((mult)*(theta[i]-(theta01))));
    }
  }

  //extra term to incorporate in pde for phi
  scalarvalueType extraTerm = W*W*constV(std::sqrt(2.0))*TANPF*gradMAG*gradMAG;
  
  // The anisotropy term that enters in to the  equation for mu
  scalargradType aniso;

  for (unsigned i=0; i< phi.n_array_elements;i++) {
    if (1 /*std::abs(1-TANPF[i]*TANPF[i]) >=1.0e-06*/) {
       aniso[0][i] = W[i]*W[i]*phix[0][i]-W[i]*W_theta[i]*phix[1][i];
       aniso[1][i] = W[i]*W[i]*phix[1][i]+W[i]*W_theta[i]*phix[0][i];       
     }
     else {
       aniso[0][i] = 0;
       aniso[1][i] = 0;
     }
  }

  // Define the terms in the equations
  scalarvalueType eq_mu = (-f_phi-extraTerm);
  scalargradType eqx_mu = (-aniso);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_RHS(2,eq_mu);
  variable_list.set_scalar_gradient_term_RHS(2,eqx_mu);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
