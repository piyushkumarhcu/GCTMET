#include "cordic.h"

void cordic(THETA_TYPE theta, COS_SIN_TYPE &s, COS_SIN_TYPE &c){
#pragma HLS PIPELINE II=9

  ap_uint<4> QUAD1 = 1;
  ap_uint<4> QUAD2 = 2;
  ap_uint<4> QUAD3 = 3;
  ap_uint<4> QUAD4 = 4;
  ap_uint<4> quadrant;
  THETA_TYPE z;


  COS_SIN_TYPE current_cos = 0.6073;
  COS_SIN_TYPE current_sin = 0.0;
  COS_SIN_TYPE factor = 1.0;
  COS_SIN_TYPE pi = 3.1415926;



  if ((theta >= 0) && (theta <= pi/2)){
      quadrant = QUAD1;
      z =theta;
   }
  else if ((theta > pi/2) && (theta<=pi)){
	  quadrant = QUAD2;
      z = (pi - theta);
   }

   /*else if (theta <= 4.71239)
      {
      quadrant = QUAD3;
      z = (theta - 3.14159);
      }
    else if (theta > 4.71239) && (theta <= 6.28319))
      {
      quadrant = QUAD4;
      z = (6.28319 - theta);
      }*/


  for (ap_uint<5> j = 0; j < NUM_ITERATIONS; j++) {


       #pragma HLS unroll
       ap_int<4> sigma = (z < 0) ? -1 : 1;



      COS_SIN_TYPE cos_shift = current_cos * sigma * factor;
      COS_SIN_TYPE sin_shift = current_sin * sigma * factor;


      current_cos = current_cos - sin_shift;

      current_sin = current_sin + cos_shift;


      z = z - sigma * cordic_phase[j];

      factor = factor/2;

  }

  //c = (quadrant==QUAD1 || quadrant==QUAD4) ? current_cos : -current_cos;
  //s = (quadrant==QUAD1 || quadrant==QUAD2) ? current_sin : -current_sin;

  s = current_sin;
  if(QUAD1 == quadrant){
	  c = current_cos;
  }
  else{
	  c = -current_cos;
  }


}
