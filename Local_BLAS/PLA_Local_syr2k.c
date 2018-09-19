/*
   PLAPACK Release 3.2
   
   Copyright (C) 2002
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_syr2k( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		     PLA_Obj beta,  PLA_Obj C )
{
  int 
    local_n, local_k,
    lda, ldb, ldc;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_B, *buf_C,
    *alphabuf;
  char 
    Trans[1], Uplo[1];
  
  PLA_Local_scal( beta, C );

  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( B, &buf_B);
  PLA_Obj_local_buffer( C, &buf_C);
  
  PLA_Obj_local_buffer( alpha, &alphabuf);
  
  PLA_Obj_local_length( C, &local_n );
  
  if ( PLA_LOWER_TRIANGULAR == uplo )
    Uplo[0] = 'L';
  else if ( PLA_UPPER_TRIANGULAR == uplo )
    Uplo[0] = 'U';

  if( PLA_NO_TRANS == trans) {
    Trans[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
  }
  else if( PLA_TRANS == trans ){
    Trans[0] = 'T';
    PLA_Obj_local_length( A, &local_k );
  }
  else if(PLA_CONJ == trans) {
    Trans[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
    PLA_Conjugate(A);
    PLA_Conjugate(B);
  }
  else /* if( PLA_CONJ_TRANS == trans ) */
    {
      Trans[0] = 'C';
      PLA_Obj_local_length(  A, &local_k );
    }
  
  PLA_Obj_local_ldim ( A, &lda);
  PLA_Obj_local_ldim ( B, &ldb);
  PLA_Obj_local_ldim ( C, &ldc);
  
  PLA_Obj_datatype(C, &datatype);
  
  if ( 0 != local_n && 0 != local_k ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;
      
      PLA_dsyr2k( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&d_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_ssyr2k( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&f_one, buf_C, &ldc); 
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_COMPLEX c_one = {1.0,0.0};
      
      PLA_csyr2k( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&c_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_DOUBLE_COMPLEX z_one = {1.0,0.0};
      
      PLA_zsyr2k( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&z_one, buf_C, &ldc); 
    }
  }
  
  if(PLA_CONJ == trans) {
    PLA_Conjugate(A);
    PLA_Conjugate(B);
  }
  
  return PLA_SUCCESS;
}
