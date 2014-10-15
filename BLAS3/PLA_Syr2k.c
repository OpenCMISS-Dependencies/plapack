/*
   PLAPACK Release 3.0
   
   Copyright (C) 2002
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Syr2k ( int uplo, int trans, PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                                     PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS,
    nb_alg;
  
  PLA_Template
    templ;

  if ( PLA_ERROR_CHECKING )
    value = PLA_Syr2k_enter ( uplo, trans, alpha, A, B, beta, C ); 

  PLA_Obj_template( C, &templ ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Syr2k_panpan ( nb_alg, uplo, trans, alpha, A, B, beta, C );

  if ( PLA_ERROR_CHECKING )
    value = PLA_Syr2k_exit ( uplo, trans, alpha, A, B, beta, C ); 

  return PLA_SUCCESS;
}
