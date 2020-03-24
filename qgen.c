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
/* J1.8 q_gen follows model on page 97 */
matrix_ptr q_gen(int dim, int k)
{
    matrix_ptr Q = salloc(3);
    matrix_ptr A, B;
    float elem;
    int i, j, p;
/* begin */
    /* Q0 */
    A = identity_gen(dim);
    B = matscalmult(lambda, *A);
    Q->row = dim;
    Q->col = dim;
    Q->element = (float *) calloc(dim*dim, sizeof(float));
    mat_copy(Q, B);
    free_matrix(A);
    free_matrix(B);
    /* Q1 */
    (Q+1)->row = dim;
    (Q+1)->col = dim;
    (Q+1)->element = (float *) calloc(dim*dim, sizeof(float));
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++) {
            p = i*dim + j;
            if (i==j)
	      if (i==dim-1)
                *((Q+1)->element+p) = - i*mu - lambda;
	      else
                *((Q+1)->element+p) = -k*theta - i*mu - lambda;		
            else if (j==(i-1))
                *((Q+1)->element+p) = i*mu;
            else
                *((Q+1)->element+p) = 0.0;
            /* end switch */
       } /* next j */
    /* Q2 */
    (Q+2)->row = dim;
    (Q+2)->col = dim;
    (Q+2)->element = (float *) calloc(dim*dim, sizeof(float));
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++) {
            p = i*dim + j;
            if (j==(i+1))
                *((Q+2)->element+p) = k*theta;
            else
                *((Q+2)->element+p) = 0.0;
        }   /* next j */
    return Q;
} /* end q_gen */
