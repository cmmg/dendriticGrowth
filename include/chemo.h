//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include <math.h>

// for normal_distribution
#include <random>
#include <chrono>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <limits>

double generateGaussianNoise(double mu, double sigma, double seedVal)
//double generateGaussianNoise(double mu, double sigma)


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
	   std::srand(seedVal);
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

double generateGaussianNoise1(double mu, double sigma, double seedVal)

//double generateGaussianNoise1(double mu, double sigma)
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
	   // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	   std::srand(seedVal);
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}



//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  dealii::Table<1,double> c_conv(n_q_points);
  dealii::Table<1,double> phi_conv(n_q_points);
	
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), phi(n_q_points) , mu(n_q_points); 
  dealii::Table<2,Sacado::Fad::DFad<double> > c_j(n_q_points, dim), phi_j(n_q_points, dim) ;
  dealii::Table<2,double> c_conv_j(n_q_points, dim), mu_conv_j(n_q_points,dim);

  for (unsigned int q=0; q<n_q_points; ++q) {
    c[q]=0.0; c_conv[q]=0.0; phi[q]=0.0; phi_conv[q]=0.0 ;  mu[q]=0.0;  
    for (unsigned int j=0; j<dim; j++) {c_j[q][j]=0.0;  phi_j[q][j]=0.0; }
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ck==0) { c[q]+=fe_values.shape_value(i, q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i]; }
      else if (ck==1){ phi[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];  }

      for (unsigned int j=0; j<dim; j++) {
	if (ck==0) {
	  c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}
	else if (ck==1)
	  { phi_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];	   
	  }

	
      }
    }
        
  }

  
  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;

      
    for (unsigned int q=0; q<n_q_points; ++q) {

      //normal distribution randomness
           
      //   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      //std::default_random_engine generator (seed);

      double fact = pow(2.11,-7);
      //   Sacado::Fad::DFad<double> std_dev = std::sqrt(2*( 1-phi[q] + ke*(1+phi[q])*(DS/DL) )*(1+(1-ke)*c[q])*(fact)) ;
      double std_dev = std::sqrt(abs(2*( 1-phi_conv[q] + ke*(1+phi_conv[q])*(DS/DL) )*(1+(1-ke)*c_conv[q])*(fact))) ;

      //double std_dev = std::sqrt(2*(1-0.5 + ke*(1+0.5)*(DS/DL)  )*(1+(1-ke)*c[q].val())*(fact)) ;
      // std::cout<<"std_dev is  "<<std_dev <<std::endl;
      
     
      Point<dim> quadPoint=fe_values.quadrature_point(q);   
      
      //noise term
      std::srand(5);       
      double noise1 = OMEGA1*( 2.0*(double)(std::rand() % 100)/100.0 -1.0 ) ;
           
      Sacado::Fad::DFad<double>  coupling = lam*(1-phi[q]*phi[q])*(1-phi[q]*phi[q]);
      Sacado::Fad::DFad<double>  dfdc = phi[q]- (phi[q]*phi[q]*phi[q]);
      dfdc+= -coupling*(c[q] +(1.0/thermal)*(quadPoint[1]-VV*currentTime)) ;
      dfdc+= (1-phi[q]*phi[q])*noise1;

    
      double  theta = std::atan2(phi_j[q][1].val(), phi_j[q][0].val());
      //Sacado::Fad::DFad<double> angle = theta;

      double angle = theta;
      //  angle = 0.125;     
      
      Sacado::Fad::DFad<double> E = (1+ ep*std::cos(mm*(angle-theta0))) ;
      Sacado::Fad::DFad<double> Eprime = -ep*mm*std::sin(mm*(angle-theta0)) ;

      dealii::Table<1,Sacado::Fad::DFad<double> > BIGTERM(dim), DIFF(dim), JATF(dim), normal(dim) ;
      BIGTERM[0]=E*E*phi_j[q][0] - E*Eprime*phi_j[q][1] ;
      BIGTERM[1]=E*E*phi_j[q][1] + E*Eprime*phi_j[q][0] ;
      
      //std::cout<<"number1 is "<<generateGaussianNoise(0, 1) <<std::endl;
      //normal[0]=generateGaussianNoise(0, std_dev.val(), quadPoint.square()); normal[1]=generateGaussianNoise1(0, std_dev.val(), std::sqrt(quadPoint.square()));
      //normal[0]=generateGaussianNoise(0, std_dev, quadPoint.square()); normal[1]=generateGaussianNoise1(0, std_dev, std::sqrt(quadPoint.square()));
      normal[0]=generateGaussianNoise(0, std_dev, quadPoint[0]); normal[1]=generateGaussianNoise1(0, std_dev, quadPoint[1]);
      
      Sacado::Fad::DFad<double> QQ = 1-phi[q] + ke*(1+phi[q])*(DS/DL);
      double root2 = pow(2.0,-0.5);
      
      DIFF[0] = D*QQ*c_j[q][0];
      DIFF[1] = D*QQ*c_j[q][1];

      
      Sacado::Fad::DFad<double> gradMAG = phi_j[q][0]*phi_j[q][0] + phi_j[q][1]*phi_j[q][1] ;
      gradMAG = std::sqrt(gradMAG);

      if (gradMAG == 0) {JATF[0]=0;JATF[1]=0; }
      else {
	JATF[0]=(1-ke*(DS/DL))*root2*(1+(1-ke)*c[q])*(1.0/dt)*(phi[q]-phi_conv[q])*(phi_j[q][0]/gradMAG);
	JATF[1]=(1-ke*(DS/DL))*root2*(1+(1-ke)*c[q])*(1.0/dt)*(phi[q]-phi_conv[q])*(phi_j[q][1]/gradMAG);
      }

      Sacado::Fad::DFad<double> crossTerm = (1+(1-ke)*c[q])*(1.0/dt)*(phi[q]-phi_conv[q]);

      //std::srand(6);
      //double noise2 = OMEGA2*(2.0*(double)(std::rand() % 100)/100.0 - 1.0 ) ;
            
      if (ck==0) {
	
	Sacado::Fad::DFad<double> multiply0 = (1+ke-(1-ke)*phi[q]);
	R[i] += (multiply0)*(1.0/dt)*fe_values.shape_value(i, q)*(c[q]-c_conv[q])*fe_values.JxW(q);  	
	R[i]+=-fe_values.shape_value(i, q)*(crossTerm)*fe_values.JxW(q);  
	//R[i]+=-fe_values.shape_value(i, q)*((1-phi[q])*noise2)*fe_values.JxW(q); 
	
	
	for (unsigned int j = 0; j < dim; j++){	
	  R[i] += fe_values.shape_grad(i, q)[j]*DIFF[j]*fe_values.JxW(q);
	  R[i] += fe_values.shape_grad(i, q)[j]*JATF[j]*fe_values.JxW(q);
	  R[i] += -fe_values.shape_grad(i, q)[j]*normal[j]*fe_values.JxW(q);
	}
	
      }
      
      else if(ck==1) {
	
	Sacado::Fad::DFad<double> multiply1 = (1 + (1-ke)*(quadPoint[1]-VV*currentTime)/(-thermal))*E*E;
	//Sacado::Fad::DFad<double> multiply1 = E*E;
	R[i] += (multiply1)*(1.0/dt)*fe_values.shape_value(i, q)*(phi[q] - phi_conv[q])*fe_values.JxW(q);  //d is added
	R[i] += -fe_values.shape_value(i, q)*(dfdc)*fe_values.JxW(q);

	for (unsigned int j = 0; j < dim; j++) {
	  //R[i]+= (L)*fe_values.shape_grad(i, q)[j]*BIGTERM[j]*fe_values.JxW(q);
	  R[i]+= fe_values.shape_grad(i, q)[j]*BIGTERM[j]*fe_values.JxW(q);	  
	  
	}

	
      }

           
    }
  } 
  
}

#endif /* CHEMO_H_ */
