#include <cmath>
#include <algorithm>
#include <fstream>
#include "myheader.h"

float dist3D_pbc(Atom a, Atom b, float Lx, float Ly)
{
	float dis0 = sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
	float dis_pLx = sqrt( (a.x-(b.x+Lx))*(a.x-(b.x+Lx)) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
	float dis_nLx = sqrt( (a.x-(b.x-Lx))*(a.x-(b.x-Lx)) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z) );
	float dis_pLy = sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-(b.y+Ly))*(a.y-(b.y+Ly)) + (a.z-b.z)*(a.z-b.z) );
	float dis_nLy = sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-(b.y-Ly))*(a.y-(b.y-Ly)) + (a.z-b.z)*(a.z-b.z) );
	float dis_pLxpLy = sqrt( (a.x-(b.x+Lx))*(a.x-(b.x+Lx)) + (a.y-(b.y+Ly))*(a.y-(b.y+Ly)) + (a.z-b.z)*(a.z-b.z) );
	float dis_nLxnLy = sqrt( (a.x-(b.x-Lx))*(a.x-(b.x-Lx)) + (a.y-(b.y-Ly))*(a.y-(b.y-Ly)) + (a.z-b.z)*(a.z-b.z) );
	float dis_nLxpLy = sqrt( (a.x-(b.x-Lx))*(a.x-(b.x-Lx)) + (a.y-(b.y+Ly))*(a.y-(b.y+Ly)) + (a.z-b.z)*(a.z-b.z) );
	float dis_pLxnLy = sqrt( (a.x-(b.x+Lx))*(a.x-(b.x+Lx)) + (a.y-(b.y-Ly))*(a.y-(b.y-Ly)) + (a.z-b.z)*(a.z-b.z) );

	
	float dis[]={dis0,dis_pLx,dis_nLx,dis_pLy,dis_nLy,dis_pLxpLy,dis_nLxnLy,dis_nLxpLy,dis_pLxnLy};
	return *min_element(dis,dis+9);	
}
