/*
   PLAPACK Release 3.1
   
   Copyright (C) 2001
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trsm( int side, int uplo, int trans, int diag,
	       PLA_Obj alpha, PLA_Obj A, PLA_Obj B ) 
{
  int 
    value = PLA_SUCCESS,
    nb_alg;
  PLA_Template
    templ = NULL;
  PLA_Obj
    ONE = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Trsm_enter( side, uplo, trans, diag, alpha, A, B );

  PLA_Obj_template( A, &templ );

  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Create_constants_conf_to( A, NULL, NULL, &ONE );

  PLA_Scal( alpha, B );

  if ( side == PLA_SIDE_LEFT ) {
    if ( uplo == PLA_LOWER_TRIANGULAR ){
      if ( trans == PLA_NO_TRANSPOSE )
	PLA_Trsm_left_lower_notrans( diag, A, B, nb_alg );
      else if ( trans == PLA_TRANSPOSE )
	PLA_Trsm_left_lower_trans( diag, A, B, nb_alg );
      else
	PLA_Trsm_left_lower( trans, diag, ONE, A, B );
    }
    else{
      if ( trans == PLA_NO_TRANSPOSE )
	PLA_Trsm_left_upper_notrans( diag, A, B, nb_alg );
      else if ( trans == PLA_TRANSPOSE )
	PLA_Trsm_left_upper_trans( diag, A, B, nb_alg );
      else
	PLA_Trsm_left_upper( trans, diag, ONE, A, B );
    }
  }
  else { /* side == PLA_SIDE_RIGHT */
    if ( uplo == PLA_LOWER_TRIANGULAR ){
      if ( trans == PLA_NO_TRANSPOSE )
	PLA_Trsm_right_lower_notrans( diag, A, B, nb_alg );
      else if ( trans == PLA_TRANSPOSE )
	PLA_Trsm_right_lower_trans( diag, A, B, nb_alg );
      else
	PLA_Trsm_right_lower( trans, diag, ONE, A, B );
    }
    else{
      if ( trans == PLA_NO_TRANSPOSE )
	PLA_Trsm_right_upper_notrans( diag, A, B, nb_alg );
      else if ( trans == PLA_TRANSPOSE )
	PLA_Trsm_right_upper_trans( diag, A, B, nb_alg );
      else
	PLA_Trsm_right_upper( trans, diag, ONE, A, B );
    }
  }      

  PLA_Obj_free( &ONE );

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Trsm_exit( side, uplo, trans, diag, alpha, A, B );

  return value;
}
