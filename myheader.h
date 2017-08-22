using namespace std;

class Atom{
public:
int sn;
string res;
string name;
int idx;
double x,y,z,vx;
};

class MolInfo{
public:
string name;
double charge;
double mass;
};

const double pi = 3.1415926535897;
float dist_pbc (Atom a, Atom b, float Lx, float Ly);
float dist3D_pbc(Atom a, Atom b, float Lx, float Ly);
Atom cartesian_pbc(Atom a, Atom b, float Lx, float Ly);
double cos (Atom a);
double cosOd_Ion(Atom od, Atom ion);
double cosOdOH_SOL(Atom ow, Atom hw1, Atom hw2);
