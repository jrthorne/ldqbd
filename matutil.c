/**********************************************************************
* File name: matutil.c
* Date start/finish: 21-12-1996/dd-mm-yyyy
* By: Jason Thorne
* Name: Math Utilities
* Project: Matrix analytic methods in Quasi-stationary birth-death
*      processes Honours project 1996-1997
* Description: This is the implementation file for the mathematic utilities
* needed to implement the methods as described in Les Bright's thesis
*****************************************************************
* Updates: NOTE ALL VALUES SHOULD BE MATTYPE
* JX.X|DD-MM-YYYY|XXX|Description
* J1.1|12-01-1997|JRT|Made routines print a message to
*                stderr if there is an error
*                also, updated write matrix.
* J1.2|15-01-1997|JRT|Changes at uni. Free return memory if error. Also
*                     added the salloc routine fot allocating and initialising
* J1.3|16-01-1997|JRT|Changes made at uni.
* J1.4|29-01-1997|JRT|free_matrix_array added
* J1.5|04-02-1997|JRT|added the external variable mat_null;
* J1.6|28-03-1997|JRT|Made write_matrix function take file pionter from
*                     argument list.
*     |29-03-1997|JRT|Wrote myabs because the C fabs routine is a load of shit.
* J1.7|03-04-1997|JRT|Found bug in code that showed up when multiplying
*             matrix*vector.
* J1.8|15-05-1997|JRT|Changed matabs to return the maximum element of the
*            matrix, to comply with the rules of an infinity norm.
* J1.9|16-05-1997|JRT|Changed write matrix to output fixed dec places.
* J1.10|17-05-1997|JRT|Changed matrix definitions to be MATTYPE precision
*            and all (matrix_ptr) mallocs to salloc(1)
* J1.11|03-08-2001|JRT| Removed inverse matrix routine, as it uses 
*            libraries from adelaide uni. Now running at CSU. I guess
*            I should write my own inverse routine.
*        NOTE: I turned comments to /% %/ in commented code
* J1.12|15-08-2001|JRT| Put in code to calculate the inverse of a matrix
*            using gaussian elimenation. This involved writing two new 
*            routines and rewriting the inverse function.
*            Changed all MATTYPE values to float. Marked each change
*            reversed all float values back to MATTYPE.
* J1.13|16-05-2001|JRT| Development of mat_invert routine. got working
*            with no check for pivot row in 2001081701
* J1.14|23-08-2001|JRT| Swapping the rows was kind of a good idea,
*            but the rows have to be swapped back again for the
*            inverse to be strictly correct. What needs to be done
*            is ponters to the rows to be swapped. A much better idea.
* J2.0|24-05-2002|JRT| Cocoa on OSX
* J2.1|01-11-2002|JRT| On code warrior this time. Found that the swap I was
*            doing in 1.14 did not work. When I got it to work, I then found
*            that if I swap rows in the augmented matrix, I need to swap 
*            columns in the resultant inverse.
*            use defined MATTYPE instead of float or MATTYPE. defined in header.
* J2.2|02-11-2002|JRT| Found a bug. No break in salloc case statement
*****************************************************************
* Version:		J2.2
* SerialNo:		20021102
***********************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include "matutil.h"
/***************************************************************/
/* hidden typedefs */
/* J1.14 */
typedef struct swaptype *swaptypeptr;
struct swaptype {
	int from;
	int to;
	};
/***********************************************************************/
    const struct matrix_type null_mat = {0, 0, (MATTYPE *) NULL}; /* J1.5 */
/* J1.12 float then MATTYPE| J1.13 made a constant */
/* this dummy routine has to be included so that I can use the fourtran
   F77 library */
//void MAIN_() {}
/***********************************************************************/
/* Implementation */
/* This routine multiplies two matricies together and the resultant matrix
   is returned as a pointer of type structure matrix */
matrix_ptr matmult(struct matrix_type left, struct matrix_type right)
{
    int i, j, k, no_elem;
    MATTYPE sum; /* J1.12 */
    matrix_ptr result = salloc(1);
/* begin */
    if (left.col != right.row){
        fprintf(stderr, "Error in matmult\n");  /* J1.1 */
        free(result);            /* J1.2 */
        result = (matrix_ptr) NULL;    /* J1.3 */
    } else {
        result->row = left.row;
        result->col = right.col;
        no_elem = left.row*right.col;
        result->element = (MATTYPE *) calloc(no_elem, sizeof(MATTYPE));
	/* J1.12 */
        for (i=0; i<left.row; i++)
            for (j=0; j<right.col; j++) {
                sum = 0.0;
                for (k=0; k<left.col; k++) /* J1.7 */
                    sum += *(left.element + i*left.col + k) *
                           *(right.element + k*right.col + j);
                *(result->element + i*result->col + j) = sum;
            }
        return result;
    }
} /* end matmult */
/**********************************************************************/
/* add matricies together */
matrix_ptr matadd(struct matrix_type a, struct matrix_type b)
{
    matrix_ptr result = salloc(1);
    int i, no_elem;
    MATTYPE *elem, *aelem, *belem;
    /* J1.12 */
    /* a and b have to be of the same dimention */
    if (a.row != b.row || a.col != b.col){
        fprintf(stderr, "Error in matadd\n");  /* J1.1 */
        free(result);
        result = (matrix_ptr) NULL;
    } else {
        no_elem = a.row * a.col;
        aelem = a.element;
        belem = b.element;
        result->row = a.row;
        result->col = a.col;
        result->element = elem = (MATTYPE *) calloc(no_elem, sizeof(MATTYPE));
	/* J1.12 */
        for (i=0; i<no_elem; i++) {
            *elem = *(aelem + i) + *(belem + i);
            elem++;
        }
    }
    return result;
} /* end matadd */
/************************************************************************/
/* take matricies from each other */
matrix_ptr mattake(struct matrix_type a, struct matrix_type b)
{
    matrix_ptr result = salloc(1);
    int i, no_elem;
    MATTYPE *elem, *aelem, *belem;
    /* J1.12 */
    /* a and b have to be of the same dimention */
    if (a.row != b.row || a.col != b.col) {
        fprintf(stderr, "Error in mattake\n");  /* J1.1 */
        free(result);
        result = (matrix_ptr) NULL;
    } else {
        no_elem = a.row * a.col;
        aelem = a.element;
        belem = b.element;
        result->row = a.row;
        result->col = a.col;
        result->element = elem = (MATTYPE *) calloc(no_elem, sizeof(MATTYPE));
	/* J1.12 */
        for (i=0; i<no_elem; i++) {
            *elem = *(aelem + i) - *(belem + i);
            elem++;
        }
    }
    return result;
} /* end mattake */
/**********************************************************************/
/* This routine writes a matrix to a file. J1.6 the file pointer is not
  global */
void write_matrix(FILE *fp, matrix_ptr a)
{
    va_list ap;
    int i,j;
    int row = a->row;
    int col = a->col;
    MATTYPE *elem = a->element;
    /* J1.12 */
    for (i=0; i<row; i++) {
        for (j=0; j<col; j++)
           fprintf(fp, "%11.8lf ", *(elem + i*col + j)); /* J1.9 */
        fprintf(fp, "\n");
    }
    return;
} /* end write_matrix */
/**********************************************************************/
/* This routine is the assignment operator for the matrix structure
   the number of rows followed by the number of columns is sent followed
   by a variable length argument list of the elements. These elements are
   expected to be real */
matrix_ptr mat_assign(int row, int col, ...)
{
    va_list ap;
    int i, j;
    int no_elem = row*col;
    matrix_ptr result = salloc(1);
    MATTYPE *elem;
    /* J1.12 */
    va_start(ap, col);
    result->row = row;
    result->col = col;
    result->element = elem = (MATTYPE *) calloc(no_elem, sizeof(MATTYPE));
                           /* J1.12 allocate memory for the return value */
    for (i=0; i<row; i++)
        for (j=0; j<col; j++) {
            *elem = (MATTYPE) va_arg(ap, MATTYPE); /* J1.12 */
            elem++;
        }
    elem = (MATTYPE *) NULL;  /* J1.12 Null terminate */
    va_end(ap);
    return result;
} /* end mat_assign */
/*********************************************************************/
void free_matrix(matrix_ptr a)
{
    MATTYPE *p; /* J1.12 */
    if (a !=(matrix_ptr) NULL) {
        p = a->element;
        free(p);
        free(a);
        a = (matrix_ptr) NULL;
    }
} /* end free_matrix */
/***********************************************************************/
/* scalar multiply */
matrix_ptr matscalmult(MATTYPE x, struct matrix_type a) /* J1.12 */
{
    matrix_ptr result = salloc(1);
    int i;
    int no_elem = a.row * a.col;
    MATTYPE *elem, *aelem = a.element; /* J1.12 */
    result->row = a.row;
    result->col = a.col;
    result->element = elem = (MATTYPE *) calloc(no_elem, sizeof(MATTYPE));
              /* J1.12 */
    for (i=0; i<no_elem; i++)
        *(elem+i) = *(aelem + i) * x;
    return result;
} /* end matscalmult */
/*************************************************************************/
/* This routine reads a matrix from a file. The structure used in the file
   is first the number of rows. followed by the number of columns, then the
   matrix.
*/
matrix_ptr mat_read(FILE *fp)
{
   matrix_ptr result = salloc(1);
   int i;
   int no_row, no_col;
   MATTYPE a, *elem; /* J1.12 */
   char c; /* to get white spaces */
    /* read in the dimentions of the matrix first */
   fscanf(fp, "%i %i%c",  &no_row, &no_col, &c);
   result->row = no_row;
   result->col = no_col;
   result->element = elem = (MATTYPE *) calloc(no_row*no_col, sizeof(MATTYPE));
   /* J1.12 */
   for(i=0; i<no_row*no_col && (fscanf(fp, "%lf%c", &a, &c)) != EOF; i++)
        *(elem + i) = a;
   if (i<no_row*no_col) {
        fprintf(stderr, "Error in matread\n");  /* J1.1 */
        free(result->element);
        free(result);
        return (matrix_ptr) NULL;
   } else
       return result;
}
/********************************************************************/
/* This function generates the identity matrix (allways square) */
matrix_ptr identity_gen(int dim)
{
   matrix_ptr result = salloc(1);
   int i, size;
   int row, col;
   MATTYPE *elem; /* J1.12 */
/* begin */
   size = dim*dim;
   result->row = result->col = dim;
   result->element = elem = (MATTYPE *) calloc(size, sizeof(MATTYPE));
   /* J1.12 */
   for (i=0; i<size; i++) {
       row = i/dim;
       col = i - row*dim;
       if (row == col)
           *(elem+i) = 1.0;
       else
           *(elem+i) = 0.0;
   }
   return result;
} /* end identity_gen */
/****************************************************************
* Mat invert uses Gaussian elimination and choosing of a pivot
* row both forward to resolve to echelon form, then backward
* to diagonalise. This is rarely the best way to find an
* inverse.
* If the inverse can not be calculated, it returns null
****************************************************************/
matrix_ptr mat_invert(matrix_ptr A)
{
     int Arow 	= A->row;
     int Acol 	= A->col;
     int AugRow, AugCol; /* the row and columns of the augmented matrix */
     int lhs, rhs;
     MATTYPE coeff, denominator;
     matrix_ptr result, Aug; /* the result and Augmented matrix */
     int pivot; /* The pivot row */
     /* J1.14 to store the swaps. Swaptype has a from and to integer field,
        and there is one swaptype for each row of the matrix to be inverted
        twoint defined in matutil.h */
     swaptypeptr swapstore = calloc(Arow+1, sizeof(struct swaptype));
     swaptypeptr swapstoreIndex = swapstore;
     swaptypeptr i; /* used in loop to swap back */
/* begin */
     /* can only calculate the inverse of a square matrix */
     if (Arow == Acol && Arow > 0) 
     {
	 	/* Augment A with I */
     	Aug = augment_with_I(A);
	 	AugRow = Aug->row;
	 	AugCol = Aug->col;
        /* perform gaussian elimination to get A in echelon form */
    	for (rhs=0; rhs<(AugRow-1); rhs++)
		{
	    	/* Choose a pivot row from those below, if the current row */
	    	pivot = 0;
	    	do {
	        	denominator = *(Aug->element + (rhs+pivot)*AugCol + rhs);
	    	}while (++pivot < (Arow-rhs) && myabs(denominator) <= CLOSEZERO);
	    	/* if there is some appropriate row then swap it for rhs
	       	else we can not solve for the inverse */
        	if (myabs(denominator) > CLOSEZERO)
	    	{
	    		if (--pivot > 0) 
				{
		   			swap_elem_row_op(Aug, (pivot+rhs), rhs);
		   			swapstoreIndex++; /* J2.1 can't start at the root, must be one further for 
		   			                     reverse swap of columns at end of algorithm */
		   			swapstoreIndex->from = pivot+rhs;
		   			swapstoreIndex->to = rhs;
		   			
        		} /* end if */
        	} else /* can not solve for this inverse */
            	return (matrix_ptr) NULL;
        	/* end if */
	    	for (lhs=rhs+1; lhs<AugRow; lhs++)
	    	{   
				coeff = -( *(Aug->element + lhs*AugCol + rhs) /
			       denominator);
                if (myabs(coeff) > CLOSEZERO) /* if value effectively != 0 */
                    arith_elem_row_op(Aug, lhs, rhs, coeff);
                /* end if */
            } /* end for lhs */
         } /* end for rhs */
		/* now use Gaussian substitution backwards to diagonalise
	   	the matrix */
        for (rhs=(AugRow-1); rhs>=1; rhs--)
	    for (lhs=rhs-1; lhs>=0; lhs--)
	    {   /* if the pivot is zero, we can not choose this row
		   to orthogonalise the matrix. */
			if (myabs(*(Aug->element + rhs*AugCol + rhs)) > CLOSEZERO )
		    	coeff = -( *(Aug->element + lhs*AugCol + rhs) /
		               *(Aug->element + rhs*AugCol + rhs) );
        	else 
	        	return (matrix_ptr) NULL; 
        	/* end if */
        	if (myabs(coeff) > CLOSEZERO) /* if value effectively != 0 */
            	arith_elem_row_op(Aug, lhs, rhs, coeff);
         	/* end if */ /* If coeff is zero, this is not a problem */
         } /* end for lhs */
         /* end for rhs */
         /* now normalise the rows */
	 	 for (lhs=0; lhs<AugRow; lhs++)
	 	 {
	     	 if (myabs(*(Aug->element + lhs*AugCol + lhs)) > CLOSEZERO )
	         	coeff = 1.0 / *(Aug->element + lhs*AugCol + lhs);
             else 
       	     /* This is a serious error. Need to put in code */
	         	coeff = 0.0;
             /* end  if */
	     	/* To normalise the matrix coeff-1.0 is used so that I
			don't need to write extra code to multiply the row */
            arith_elem_row_op(Aug, lhs, lhs, (coeff-1.0));
      	 } /* end for */
      	 result = deaugment(Aug);
	  	 /* J1.14 now swap columns in reverse order to rows. J2.1 to understand this for
	  	    loop, see 2.1 comment above. This is a crap way of doing this */
	  	 for (i=swapstoreIndex; i != swapstore; i--)
	  		swap_elem_col_op(result, i->from, i->to); /* J2.1 */
     } else {
        fprintf(stderr, "Error in mat_invert\n");  /* J1.1 */
        result = (matrix_ptr) NULL;
     }
     free(swapstore);
     free_matrix(Aug);
     return result;
} /* end mat_invert */
/***********************************************************************/
/* Copy the matrix structure */
void mat_copy(matrix_ptr to, matrix_ptr mfrom)
{
    int i, no;
    int col = mfrom->col;
    int row = mfrom->row;
    MATTYPE* elem = to->element; /* j1.12 */
/* begin */
    to->row = row;
    to->col = col;
    no = row*col;
    free(to->element); /* J1.2 */
    to->element = (MATTYPE *) calloc(no, sizeof(MATTYPE)); /* J1.12 */
    for (i=0; i<no; i++)
        *(to->element+i) = *(mfrom->element+i);
    return;
} /* end mat_copy */
/*********************************************************************/
/* return a scalar value representing the matrix. 
 The infinity norm version J1.8 just finds the maximum absolute value of
 the element.
*/
MATTYPE matabs(matrix_ptr A) /* J1.12 */
{
    int i;
    int cols = A->col;
    int rows = A->row;
    MATTYPE elem; /* J1.12 */
    MATTYPE maxelem = 0.0;
/* begin */
    for (i=0; i<rows*cols; i++){
        elem = *(A->element+i);
        maxelem = myabs(elem)>maxelem ? myabs(elem) : maxelem;
    }
    return maxelem;
} /* end matabs */
/**********************************************************************/
/* J1.2 allocate memory and initialisr the element pointer */
matrix_ptr salloc(int n)
{
    matrix_ptr result;
    int i;
/* begin */
    switch (n) {
    case 0:
        result = (matrix_ptr) NULL;
        break; /* J2.2 found a bug. No break in case */
    case 1:
        result = (matrix_ptr) malloc(sizeof(struct matrix_type));
        result->element = (MATTYPE *) NULL; /* J1.12 */
        break;
    default:
        result = (matrix_ptr) calloc(n, sizeof(struct matrix_type));
        for (i=0; i<n; i++)
            (result+i)->element = (MATTYPE *) NULL; /* J1.12 */
        break;
    }
    return result;
} /* end salloc */
/***********************************************************************/
/* J1.4 free array of matricies. no_matricies is the same no as in salloc */
void free_matrix_array(matrix_ptr A, int no_matricies)
{
int i;
/* begin */
    for(i=(no_matricies-1); i>=0; i--)
        free((A+i)->element);
    free(A);
} /* end free_matrix_array */
/**************************************************************************/
/* J1.6 return absolute value of the MATTYPEing point number */
MATTYPE myabs(MATTYPE x) /* J1.12 */
{
    if (x<0)
        return -x;
    else
        return x;
}
/***************************************************************/
/* J1.12 This functions performs an elementry arithmetic row operation
   on a matrix. The sum is: A(lhs) = A(lhs) + coeff*A(rhs)
      where "A(i)" = row i of matrix A 
      It is important to note that the values of the passed
      matrix are changed in this routine 
      ALSO this does no checking to see if lhs or rhs are within 
      bounds */
/***************************************************************/
void arith_elem_row_op(matrix_ptr A, int lhs, int rhs, MATTYPE coeff)
{
int 	i, j; /* loop variables */
int	cols		= A->col;
/* begin */
    for (j=0; j<cols; j++)
	*(A->element + lhs*cols + j) += 
		     coeff * *(A->element + rhs*cols + j);
    return;
} /* end arith_elem_row_op */
/***************************************************************/
/* J1.12 This function returns the given matrix augmented with
   the identity matrix for finding the inverse matrix using
   row operations */
/***************************************************************/
matrix_ptr augment_with_I(matrix_ptr A)
{
matrix_ptr 	result 		= salloc(1);
int 		i, j, noElem;
int 		col 		= A->col;
int 		row 		= A->row;
/* begin */
    result->row 	= row;
    result->col 	= col*2; 
    /* augmented with I, so twice as many columns */
    noElem		= row*col*2;
    result->element 	= (MATTYPE *) calloc(noElem, sizeof(MATTYPE));
    
    for (i=0; i<row; i++)
	for (j=0; j<2*col; j++)
	    if (j<col) /* insert A element in augmented matrix */
		*(result->element + i*col*2 + j) = 
		*(A->element + i*col + j);
            else /* insert I matrix to augmented matrix */
		    *(result->element + i*col*2 + j) = 
		    (i == j - col) ? 1.0 : 0.0; /* 1 is diag, else 0 */
            /* end if */
        /* end for j */
    /* end for i */
    
    return result;
} /* end augment_with_I */
/***************************************************************/
/* J1.13 This function returns the given matrix deaugmented with
   the the last half. Used for finding the inverse matrix using
   row operations */
/***************************************************************/
matrix_ptr deaugment(matrix_ptr A)
{
matrix_ptr 	result 		= salloc(1);
int 		i, j, noElem;
int 		col 		= A->col;
int 		row 		= A->row;
int		ResCol		= col/2;
/* There are half as many columns in the resultant deaugmented matrix */
/* begin */
    result->row 	= row;
    result->col 	= ResCol; 
    /* augmented with I, so half as many columns */
    noElem		= row*ResCol;
    result->element 	= (MATTYPE *) calloc(noElem, sizeof(MATTYPE));
    
    for (i=0; i<row; i++)
	for (j=ResCol; j<col; j++)
            *(result->element + i*ResCol + (j-ResCol)) = 
	                        *(A->element + i*col + j);
        /* end for j */
    /* end for i */
    
    return result;
} /* end deaugment */
/****************************************************************
* J1.13 this function is the elementary row operation swap
* row x for row y
* it changes the matrix passed to it and returns nothing
****************************************************************/
void swap_elem_row_op(matrix_ptr A, int x, int y)
{
MATTYPE		temp;
int		col	= A->col;
int		row	= A->row;
int		i;
/* begin */
    for (i=0; i<row; i++)
    {
		temp		 		= *(A->element + x*col + i);
		*(A->element + x*col + i)	= *(A->element + y*col + i);
		*(A->element + y*col + i)	= temp;
    } /* end for i */
    return;
} /* end swap_elem_row_op */
/****************************************************************
* J2.1 this function is the elementary column operation swap
* col x for col y
* it changes the matrix passed to it and returns nothing
****************************************************************/
void swap_elem_col_op(matrix_ptr A, int x, int y)
{
MATTYPE		temp;
int		col	= A->col;
int		row	= A->row;
int		i;
/* begin */ 
    for (i=0; i<row; i++)
    {
		temp		 				= *(A->element + x + i*col);
		*(A->element + x + i*col)	= *(A->element + y + i*col);
		*(A->element + y + i*col)	= temp;
    } /* end for i */
    return;
} /* end swap_elem_row_op */
