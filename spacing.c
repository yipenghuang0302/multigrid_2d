#include "2d.h"

double * spacing (
  int n,
  double a,
  double b
) {

  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 ) { x[0] = ( a + b ) / 2.0; }
  else {
    for ( i = 0; i < n; i++ ) {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return x;
}