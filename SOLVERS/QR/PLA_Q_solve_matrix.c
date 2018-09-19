/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Q_solve_matrix( int side, int trans, PLA_Obj A, PLA_Obj s, PLA_Obj B )
{
  if ( side == PLA_SIDE_LEFT ){
    if ( trans == PLA_NO_TRANSPOSE )
      PLA_QX_eq_B_solve_matrix( A, s, B );
    else if ( trans == PLA_TRANSPOSE )
      PLA_QtX_eq_B_solve_matrix( A, s, B );
    else
      PLA_Abort( "Case not yet implemented", __LINE__, __FILE__ );
  }
  else{
    PLA_Abort( "Case not yet implemented", __LINE__, __FILE__ );
  }
}

int PLA_QX_eq_B_solve_matrix( PLA_Obj A, PLA_Obj s, PLA_Obj B )
/*
  Purpose: Solve Q X = B where Q is stored as Householder vectors
           below the diagonal of A and in s.

  Input:  A       --   General mxn matrix A. Stores Householder vectors
                       below the diagonal
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  B       --   right-hand-side(s)
                       (MVECTOR or MATRIX of length m)

	  B       --   solution of Q X = B

  Return value: PLA_SUCCESS iff Q_solve is completed successfully
*/
{
  int       size, me,
            nb_alg1, nb_alg2, nb_alg;
  PLA_Template  templ;
  PLA_Obj   ABR         = NULL,
            A1          = NULL,     A2              = NULL,
            A1_mv       = NULL,     
            W_mv        = NULL,     
            Y_mv        = NULL,     
            sB          = NULL,     s1              = NULL,
            s1_dup      = NULL,     BB              = NULL;

  /* Determine algorithmic blocking size */
  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg1 ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg2 ); 
  nb_alg = ( nb_alg1 > nb_alg2 ? nb_alg1 : nb_alg2 );

  /* Initially partition A = / ATL |  *  \  where
                             \ ABL | ABR /  ATL is 0x0  */
  PLA_Obj_view_all( A, &ABR );

  /* Initially partition s = / sT \
                             \ sB / where sT is of length 0 */
  PLA_Obj_view_all( s, &sB );

  /* Initially partition B = / BT \  where
                             \ BB /  BT has 0 rows */
  PLA_Obj_view_all( B, &BB );

  while ( TRUE ) {
    /* Check if done */
    PLA_Obj_global_width( ABR , &size ); 
    if ( ( size = min( size, nb_alg ) ) == 0 ) break;

    /* Partition ABR = ( A1, A2 ) */
    PLA_Obj_vert_split_2( ABR,      size, &A1, &A2 );

    /* Partition sB  = ( s1, 
                         sB ) */
    PLA_Obj_horz_split_2( sB, size, &s1, 
                                    &sB );


    /* Create a duplicated copy of s1, in which to duplicate
       beta's */
    PLA_Mscalar_create_conf_to( s1, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s1_dup );
    PLA_Copy( s1, s1_dup );

    /* Copy the current panel to a multivector */
    PLA_Mvector_create_conf_to( A1, size, &A1_mv );

    PLA_Copy( A1, A1_mv );

    /* Create multivectors for W and Y and compute WY transform */
    PLA_Mvector_create_conf_to( A1, size, &W_mv );
    PLA_Mvector_create_conf_to( A1, size, &Y_mv );
    PLA_Compute_WY( A1_mv, s1_dup, W_mv, Y_mv );

    /* BB <- ( I + W Y^T ) BB */
    PLA_Apply_W_Y_transform ( PLA_SIDE_LEFT, PLA_TRANSPOSE,
                              W_mv, Y_mv, BB); 

    /* Update view ABR to view currently active part of A */
    PLA_Obj_horz_split_2( A2, size,          PLA_DUMMY,
                                             &ABR );

    /* Update view BB to view currently active part of B */
    PLA_Obj_horz_split_2( BB, size,          PLA_DUMMY,
                                             &BB );
  }

  /* Free temporary objects and views */
  PLA_Obj_free( &ABR ); 
  PLA_Obj_free( &A1 );            PLA_Obj_free( &A2 );
  PLA_Obj_free( &A1_mv );         
  PLA_Obj_free( &W_mv );          PLA_Obj_free( &Y_mv );          
  PLA_Obj_free( &sB );            PLA_Obj_free( &s1 );
  PLA_Obj_free( &s1_dup );
  PLA_Obj_free( &BB );

  return( PLA_SUCCESS );
}


int PLA_QtX_eq_B_solve_matrix( PLA_Obj A, PLA_Obj s, PLA_Obj B )
/*
  Purpose: Solve Q^T X = B where Q is stored as Householder vectors
           below the diagonal of A and in s.

  Input:  A       --   General mxn matrix A. Stores Householder vectors
                       below the diagonal
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  B       --   right-hand-side(s)
                       (MVECTOR or MATRIX of length m)

	  B       --   solution of Q^T X = B

  Return value: PLA_SUCCESS iff Q_solve is completed successfully
*/
{
  int       length, width, size, nb_alg1, nb_alg2, nb_alg;
  PLA_Template  templ;
  PLA_Obj   ATL         = NULL,     ABR            = NULL,
            BB          = NULL,
            AB1         = NULL,     AB1_mv         = NULL,     
            W_mv        = NULL,     Y_mv        = NULL,     
            sL          = NULL,
            s_cur       = NULL,     s_dup   = NULL;

  /* Determine algorithmic blocking size */
  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg1 ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg2 ); 
  nb_alg = ( nb_alg1 > nb_alg2 ? nb_alg1 : nb_alg2 );

  /* Apply the Householder vectors like 
                 ( H_1 ( ... ( H_n-1 ( H_n B ) ) ... ) ) */

  /* Initially partition A = / ATL |  *  \  where ATL is square and
                             \ ABL | ABR /  ABR has width 0  */

  PLA_Obj_global_length( A, &length );
  PLA_Obj_global_width ( A, &width );
  width = min( length, width );
  PLA_Obj_split_4( A, width, width,   &ATL,      PLA_DUMMY,
                                      PLA_DUMMY, &ABR );

  PLA_Obj_horz_split_2( B, width,     PLA_DUMMY, 
                                      &BB );

  PLA_Obj_horz_split_2( s, width,     &sL,       
                                      PLA_DUMMY );

  while ( TRUE ) {
    /* Check if done */

    PLA_Obj_global_width( ATL , &length ); 
    PLA_Obj_global_width( ATL , &width ); 
    if ( ( size = min( min( length, width), nb_alg ) ) == 0 ) break;

    /* Grow ABR by the panel from which to compute the next 
       WY transform.  Grow BB similarly */

    PLA_Obj_view_shift( ABR,    -size, 
                    -size,               0,
                                  0 );

    PLA_Obj_view_shift( BB,      -size,
                             0,           0,
                                    0 );

    /* Partition off AB1 from which to compute the next WY transform */
    PLA_Obj_vert_split_2( ABR, size, &AB1, PLA_DUMMY );

    /* Partition off the scaling factors and duplicate to all nodes */

    PLA_Obj_horz_split_2( sL,  -size,     &sL, 
                                          &s_cur );

    PLA_Mscalar_create_conf_to( s_cur, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s_dup );

    /* Redistribute AB1 as a multivector and compute W and Y */

    PLA_Mvector_create_conf_to( AB1, size, &AB1_mv );
    PLA_Mvector_create_conf_to( AB1, size, &W_mv );
    PLA_Mvector_create_conf_to( AB1, size, &Y_mv );

    PLA_Copy( AB1, AB1_mv );
    PLA_Copy( s_cur, s_dup );

    PLA_Compute_WY( AB1_mv, s_dup, W_mv, Y_mv );

    /* Update BB <- ( I + W Y^T ) BB */

    PLA_Apply_W_Y_transform ( PLA_SIDE_LEFT, PLA_NO_TRANSPOSE,
			      W_mv, Y_mv, BB ); 

    /* Update view of ATL */

    PLA_Obj_split_4( ATL, -size, -size,  &ATL,      PLA_DUMMY,
                                         PLA_DUMMY, PLA_DUMMY );
  }

  PLA_Obj_free( &ATL );           PLA_Obj_free( &ABR );            
  PLA_Obj_free( &BB ); 
  PLA_Obj_free( &AB1 );           PLA_Obj_free( &AB1_mv );
  PLA_Obj_free( &W_mv );          PLA_Obj_free( &Y_mv );          
  PLA_Obj_free( &s_cur );         PLA_Obj_free( &s_dup );
  PLA_Obj_free( &sL );

  return PLA_SUCCESS;
}
