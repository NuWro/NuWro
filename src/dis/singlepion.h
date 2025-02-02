#ifndef _singlepion_h_
#define _singlepion_h_

#include "params.h"
#include "resevent_hybrid.h"
const size_t NSPPbins = 400;
const double SPP_MIN = 1210;
const double SPP_MAX = Wmax_hybrid/MeV;
const double NSPPSize = (SPP_MAX - SPP_MIN) / (NSPPbins -1);
extern double COUNTER[2][2][2][3][NSPPbins];
extern double OVERALL[2][2][2][NSPPbins];	//in OVERALL we do not distiguish channels
extern double SPP[2][2][2][3][NSPPbins];
void singlepion (params & p);//produce SPP table
#endif

