/*
   PLAPACK Release 3.1
   
   Copyright (C) 2001
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symm( int side, int uplo, 
	      PLA_Obj alpha, PLA_Obj A, PLA_Obj B, 
	      PLA_Obj beta, PLA_Obj C ) 
{
  int 
    value = PLA_SUCCESS,
    nb_alg;
  PLA_Template
    templ = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Symm_enter( side, uplo, alpha, A, B, beta, C );

  PLA_Obj_template( A, &templ );

  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  if ( !value ){
    PLA_Scal( beta, C );

    if ( side == PLA_SIDE_LEFT ) {
      if ( uplo == PLA_LOWER_TRIANGULAR )
	value = PLA_Symm_left_lower( alpha, A, B, C, nb_alg );
      else 
	value = PLA_Symm_left_upper( alpha, A, B, C, nb_alg );
    }
    else{
      if ( uplo == PLA_LOWER_TRIANGULAR )
	value = PLA_Symm_right_lower( alpha, A, B, C, nb_alg );
      else 
	value = PLA_Symm_right_upper( alpha, A, B, C, nb_alg );
    }
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Symm_exit( side, uplo, alpha, A, B, beta, C );

  return value;
}

