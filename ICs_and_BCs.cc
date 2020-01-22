// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

	  //double center[1][3] = {{1.0/2.0,1.0/2.0,1.0/2.0}};
  /*  double center1[1][3] = {{104,0,0}};
	  double center2[1][3] = {{204,0,0}};
	  double center3[1][3] = {{304,0,0}};
	  
	  double rad[1] = {6.0};
  */
  double centrex1 = 50;
  double centrex2 =170;
  double centrex3 = 290;
  double centrey = 0;
  
  double rad =1.0;
  scalar_IC = 0;

  // Initial condition for the concentration field
  if (index == 0){
    scalar_IC = -0.3 ; //-delta;
  }

  else if (index==1) {
	double dist1 = (p[0]-centrex1)*(p[0]-centrex1) + (p[1]-centrey)*(p[1]-centrey) ;
	dist1 = std::sqrt(dist1) ;

	double dist2 = (p[0]-centrex2)*(p[0]-centrex2) + (p[1]-centrey)*(p[1]-centrey) ;
	dist2 = std::sqrt(dist2) ;

	double dist3 = (p[0]-centrex3)*(p[0]-centrex3) + (p[1]-centrey)*(p[1]-centrey) ;
	dist3 = std::sqrt(dist3) ;


	if (dist1> rad && dist2 >rad && dist3 >rad) {
          if (dist1 > dist2 && dist1> dist3) scalar_IC = -std::tanh((dist1-rad)/sqrt(2.0));
          else if (dist2 > dist1 && dist2> dist3) scalar_IC = -std::tanh((dist2-rad)/sqrt(2.0));
          else if (dist3 > dist2 && dist3> dist1) scalar_IC = -std::tanh((dist3-rad)/sqrt(2.0));
          else if (dist1 == dist3) scalar_IC = -std::tanh((dist1-rad)/sqrt(2.0));

        }
        else if (dist1 <= rad && dist2 >rad && dist3 >rad) {
          scalar_IC = -std::tanh((dist1-rad)/sqrt(2.0));

        }
        else if (dist1 > rad && dist2 <=rad && dist3 >rad) {
	  scalar_IC = -std::tanh((dist2-rad)/sqrt(2.0));
        }

        else if (dist1 > rad && dist2 > rad && dist3 <=rad) {
	  scalar_IC = -std::tanh((dist3-rad)/sqrt(2.0));
        }
	
	 
  }

  
  
	  

    /*
  // Initial condition for the order parameter field
  else if (index == 1) {
    // Initial condition for the order parameter field 
    for (unsigned int i=0; i<1; i++) {
      dist1 = 0.0;dist2=0.0; dist3=0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
	dist1 += (p[dir]-center1[i][dir])*(p[dir]-center1[i][dir]);
	dist2 += (p[dir]-center2[i][dir])*(p[dir]-center2[i][dir]);
	dist3 += (p[dir]-center3[i][dir])*(p[dir]-center3[i][dir]);				  
      }
      dist1 = std::sqrt(dist1);
      dist2 = std::sqrt(dist2);
      dist3 = std::sqrt(dist3);
			  
      if (dist1<=rad[i]) { scalar_IC = 1.0; }
      else if (dist1>rad[i]) { scalar_IC = -1.0; }
      if (dist2<=rad[i]) { scalar_IC = 1.0; }
      else if (dist2>rad[i]) { scalar_IC = -1.0; }
      if (dist3<=rad[i]) { scalar_IC = 1.0; }
      else if (dist3>rad[i]) { scalar_IC = -1.0; }
			  
    }
		  
  }

  */

  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
