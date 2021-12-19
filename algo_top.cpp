#include "algo_top_parameters.h"
#include "algo_top.h"
#include "cordic.h"
#include <algorithm>
#include <utility>

#include "../../../../include/objects.h"
#include "ap_fixed.h"
using namespace std;

TowersInEta unpackInputLink(ap_uint<576> &link) {
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS INLINE

TowersInEta tEta_;

  tEta_.towers[0]  = Tower(link( 31,   0));
  tEta_.towers[1]  = Tower(link( 63,  32));
  tEta_.towers[2]  = Tower(link( 95,  64));
  tEta_.towers[3]  = Tower(link(127,  96));
  tEta_.towers[4]  = Tower(link(159, 128));
  tEta_.towers[5]  = Tower(link(191, 160));
  tEta_.towers[6]  = Tower(link(223, 192));
  tEta_.towers[7]  = Tower(link(255, 224));
  tEta_.towers[8]  = Tower(link(287, 256));
  tEta_.towers[9]  = Tower(link(319, 288));
  tEta_.towers[10] = Tower(link(351, 320));
  tEta_.towers[11] = Tower(link(383, 352));
  tEta_.towers[12] = Tower(link(415, 384));
  tEta_.towers[13] = Tower(link(447, 416));
  tEta_.towers[14] = Tower(link(479, 448));
  tEta_.towers[15] = Tower(link(511, 480));
  tEta_.towers[16] = Tower(link(543, 512));

  return tEta_;
}

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]) {
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS INTERFACE ap_ctrl_hs port=return

  // Step 1: Unpack links
  // Input is 64 links carrying 32phix34eta towers
	TowersInEta towersInPosEta[TOWERS_IN_PHI];
	TowersInEta towersInNegEta[TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=towersInPosEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=towersInNegEta complete dim=0
     
  for (size_t ilink = 0; ilink < N_INPUT_LINKS/2; ilink++) {
#pragma LOOP UNROLL
	  size_t iPosEta = ilink;
	  size_t iNegEta = ilink + (N_INPUT_LINKS/2);
	  towersInPosEta[ilink] = unpackInputLink(link_in[iPosEta]);
	  towersInNegEta[ilink] = unpackInputLink(link_in[iNegEta]);
  }

   // Step 2: MET Algo goes here
  ap_uint<16> Exs[24];
  ap_uint<16> Eys[24];
	#pragma HLS ARRAY_PARTITION variable=Exs complete dim=0
	#pragma HLS ARRAY_PARTITION variable=Eys complete dim=0

	COS_SIN_TYPE sinphi[TOWERS_IN_PHI];
	COS_SIN_TYPE cosphi[TOWERS_IN_PHI];
	#pragma HLS ARRAY_PARTITION variable=sinphi complete dim=0
	#pragma HLS ARRAY_PARTITION variable=cosphi complete dim=0

	/*---sinphi and cosphi calculation through CORDIC---*/

  for(ap_uint<10> degree=0; degree < TOWERS_IN_PHI; degree++){
	#pragma HLS unroll
	ap_fixed<32, 20> theta = (degree*5+2.5)*0.01745329;  /*---degree to radian conversion---*/
	cordic(theta, sinphi[degree], cosphi[degree]);
  }

  for (ap_uint<6> b = 4; b < 28; b++) {
		#pragma hls unroll
	  	  ap_uint<16> Ey;
	  	  ap_uint<16> Ex;
	  	  ap_uint<16> j;
	  	  ap_uint<16> k, l, m, n, o, p, q, r, s, t, u;
	  	  ap_uint<16> v = 0;

  	 j = towersInPosEta[b].towers[0].tower_et() + towersInPosEta[b].towers[1].tower_et() + towersInPosEta[b].towers[2].tower_et() + towersInPosEta[b].towers[3].tower_et();
  	 k = towersInPosEta[b].towers[4].tower_et() + towersInPosEta[b].towers[5].tower_et() + towersInPosEta[b].towers[6].tower_et() + towersInPosEta[b].towers[7].tower_et();
  	 l = towersInPosEta[b].towers[8].tower_et() + towersInPosEta[b].towers[9].tower_et() + towersInPosEta[b].towers[10].tower_et() + towersInPosEta[b].towers[11].tower_et();
		 m = towersInPosEta[b].towers[12].tower_et() + towersInPosEta[b].towers[13].tower_et() + towersInPosEta[b].towers[14].tower_et() + towersInPosEta[b].towers[15].tower_et();
	   n = towersInPosEta[b].towers[16].tower_et() + towersInNegEta[b].towers[0].tower_et() + towersInNegEta[b].towers[1].tower_et() + towersInNegEta[b].towers[2].tower_et();
		 o = towersInNegEta[b].towers[3].tower_et() + towersInNegEta[b].towers[4].tower_et() + towersInNegEta[b].towers[5].tower_et() + towersInNegEta[b].towers[6].tower_et();
		 p = towersInNegEta[b].towers[7].tower_et() + towersInNegEta[b].towers[8].tower_et() + towersInNegEta[b].towers[9].tower_et() + towersInNegEta[b].towers[10].tower_et();
		 q = towersInNegEta[b].towers[11].tower_et() + towersInNegEta[b].towers[12].tower_et() + towersInNegEta[b].towers[13].tower_et() + towersInNegEta[b].towers[14].tower_et();
		 r = towersInNegEta[b].towers[15].tower_et() + towersInNegEta[b].towers[16].tower_et();
		 s = j + k + l + m;
		 t = n + o + p + q;
		 u = s + t + r + v;

		 Ey = sinphi[b]*u;
		 Eys[b-4] = Ey;

		 Ex = cosphi[b]*u;
		 Exs[b-4] = Ex;
  	}

  // Step 3: Pack the outputs

link_out[0] = 0;
size_t start = 0;
link0OutputLoop: for (size_t tower = 0; tower < 24; tower++) {
#pragma HLS UNROLL
    size_t end = start + 15;
    link_out[0].range(end, start) = Exs[tower];
    start += 16;
  }

link_out[1] = 0;
start = 0;
link1OutputLoop: for (size_t tower = 0; tower < 24; tower++) {
#pragma HLS UNROLL
    size_t end = start + 15;
    link_out[1].range(end, start) = Eys[tower];
    start += 16;
  }

}

