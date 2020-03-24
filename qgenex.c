/**********************************************************************
* File name: qgen.c
* Date start/finish: 16-05-1996/17-01-1997 tested
* By: Jason Thorne
* Name: Q element generator.
* Project: Matrix analytic methods in Quasi-stationary birth-death
*      processes Honours project 1996-1997
* Description: Generates the block elements for the infinitesimal 
  generator matrix.
* Updates:
* JX.X:DD-MM-YYYY|XXX|Description
* Version J1.0
***********************************************************************/
#include"matutil.h"
/* external variables whose scope is this source code and matmain.c */
    extern float lambda;
    extern float theta;
    extern float mu;
extern float delta;
/*************************************************************************/
/* J1.8 q_gen follows my model described in Example 2.4.1 */
matrix_ptr q_gen(int dim, int k)
{
    matrix_ptr Q = salloc(3);
    matrix_ptr A, B;
    float elem, temp_factor, muk, lak, mukj, lakj;
    int i, j, p;
/* begin */
    muk = (k+1)*mu;
    lak = (k+1)*lambda;
    /* Q0 */
    Q->row = dim;
    Q->col = dim;
    Q->element = (float *) calloc(dim*dim, sizeof(float));
   
    /* Q1 */
    (Q+1)->row = dim;
    (Q+1)->col = dim;
    (Q+1)->element = (float *) calloc(dim*dim, sizeof(float));
    if (k==0) 
      (Q+2)->element = (float *) NULL ;
    else {
          /* Q2 */
      (Q+2)->row = dim;
      (Q+2)->col = dim;
      (Q+2)->element = (float *) calloc(dim*dim, sizeof(float));
    }
      
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++) {
            p = i*dim + j;
            if (i==j) {
	     /* the temperature influence on the rates */
	      temp_factor = 1.0 / (float) (j+1);
	      lakj = temp_factor*lak;
	      mukj = temp_factor*muk;
              *(Q->element+p) = lakj;
	    if (i==0) {
		*((Q+1)->element+p) = -(mukj + lakj + delta);
		*((Q+1)->element+(p+1)) = delta;
	      } else if (i==dim-1) {
	      	*((Q+1)->element+p) = -(mukj + lakj + delta);
		*((Q+1)->element+(p-1)) = delta;
	      } else {
		*((Q+1)->element+(p-1)) = delta;
		*((Q+1)->element+p) = -(mukj + lakj + 2.0*delta);
		*((Q+1)->element+(p+1)) = delta;
	      } /* end case */
	      if ((Q+2)->element != (float *) NULL)
		*((Q+2)->element+p) = mukj;
            } else {
	      *(Q->element+p) = 0.0;
              if (j!=i+1) /* q1 is tri-diagonal */
		*((Q+1)->element+p) = 0.0;
	      if ((Q+2)->element != (float *) NULL)
		*((Q+2)->element+p) = 0.0;
           } /* end if */
	  }/* next j */
  
  
    return Q;
} /* end q_gen */
