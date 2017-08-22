#include <iostream>	// std::cout, std::endl
#include <iomanip>      // std::setw
#include <cstdlib> /* malloc, free, rand */
#include <sstream>
#include <cstring>
#include <time.h>
#include <cmath>
#include <vector>
#include <fstream>      // std::ifstream
#include "xdrfile_trr.h"
#include "myheader.h"
using namespace std;

#ifdef _OPENMP
#define OMP(x)  omp x
#else
#define OMP(x)  message("")
#endif

float dist3D_pbc(Atom a, Atom b, float Lx, float Ly);// dis3D_pbc is declared in dis3D_pbc.cpp

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}


int main(int argc, char* argv[])
{
	/*Open input and output files, excepting for .trr files*/
	ifstream inputSOL, molinfo;
	ofstream output1, output2, log;
	inputSOL.open(argv[2]);
	molinfo.open(argv[3]);
	output1.open(argv[4]);
	output2.open(argv[5]);
	log.open(argv[6]);
	
	time_t nStart = time(NULL);//Start timer
	int dbg=0;
	int nbin=800; // # of slices in z direction
	int skip=0; // # of ps to be skipped (1ps might contain 10 frames, check mdp settings)
	double comx,comy,comz;
	double dV, com1z, com2z, com1x, com2x, com1y;
	
	/*Read molecules information: atom names, mass, charges*/
        string line,word;
        stringstream linestream;
	int a;
	double d;
        getline(molinfo,line); //1st line "[ molecule ]"
        getline(molinfo,line); // 2nd line "Name"
        getline(molinfo,line);//3rd line
        linestream.clear(); linestream.str(line);
        linestream >> word;
        string Sysname = word;
        getline(molinfo,line); //"[Number of molecules]"
        getline(molinfo,line); // 5th line
        linestream.clear(); linestream.str(line);
        linestream >> a;
        int nSOL = a;
        getline(molinfo,line); // "[ Number of atoms per molecule]"
        getline(molinfo,line);
        linestream.clear(); linestream.str(line);
        linestream >> a;
        int nAtom = a;
        getline(molinfo,line); // "[ Time steps (ps) ]"
        getline(molinfo,line);
        linestream.clear(); linestream.str(line);
        linestream >> d;
        double dt = d; //Sampling time steps in ps
        getline(molinfo,line); // "[ atoms ]"
        getline(molinfo,line); //" Name   charge     mass"
        MolInfo atomtype;
        vector<MolInfo> info_ary;
        while(!molinfo.eof())
        {
                getline(molinfo,line);
                if ( (line.empty()) || (line.find("[")!=string::npos) )
                {
                        break;
                }
                else
                {
                        linestream.clear(); linestream.str(line);
                        linestream >> word;
                        atomtype.name = word;
                        linestream >> d;
                        atomtype.charge = d;
                        linestream >> d;
                        atomtype.mass = d;
                        info_ary.push_back(atomtype);
                }
        }
	getline(molinfo,line);//rcut value
	linestream.clear(); linestream.str(line);
	linestream >> word;
	float rcut = atof(word.c_str());

	//Output system information to log
        log << "System information: " << endl;
        log << "Molecule:       " << Sysname << endl;
        log << "nSOL:   " << nSOL << endl;
        log << "nAtom:  " << info_ary.size() << endl;
        log << "atom types: " ;
        for (int i=0; i<info_ary.size(); i++)
        {
                if (i<info_ary.size()-1)
                {
                        log << info_ary[i].name << " ";
                }
                else
                {
                        log <<  info_ary[i].name << endl;
                }
        }
        log << "atomic mass (amu): ";
        double MW=0;
        for (int i=0; i<info_ary.size(); i++)
        {
                if (i<info_ary.size()-1)
                {
                        log << info_ary[i].mass << " ";
                }
                else
                {
                        log <<  info_ary[i].mass << endl;
                }
                MW += info_ary[i].mass;
        }
        log << "atomic charges (e): ";
        for (int i=0; i<info_ary.size(); i++)
        {
                if (i<info_ary.size()-1)
                {
                        log << info_ary[i].charge << " ";
                }
                else
                {
                        log <<  info_ary[i].charge << endl;
                }
        }
        log << "Molecular Weight (amu): " << MW << endl;
        log << "Sampling Time step (ps): "<< dt << endl;
        log << "nbin:   " << nbin << endl;
	log << "rcut:	" << rcut << endl;
        log << "Skipped steps:  " << skip << "\n"<< endl;

        //Declare arrays to store number of molecules in each bin over time frames
        double *nDensOd_bin, *nDensOH_bin, *nDensEx_bin, *VelOd_bin, *VelOH_bin, *VelEx_bin, *Vel2Od_bin, *Vel2OH_bin, *Vel2Ex_bin;
        nDensOd_bin = (double *)malloc(nbin * sizeof(double));
	nDensOH_bin = (double *)malloc(nbin * sizeof(double));
	nDensEx_bin = (double *)malloc(nbin * sizeof(double));
        VelOd_bin = (double *)malloc(nbin * sizeof(double));
	VelOH_bin = (double *)malloc(nbin * sizeof(double));
	VelEx_bin = (double *)malloc(nbin * sizeof(double));
        Vel2Od_bin = (double *)malloc(nbin * sizeof(double));
	Vel2OH_bin = (double *)malloc(nbin * sizeof(double));
	Vel2Ex_bin = (double *)malloc(nbin * sizeof(double));
        for (int i=0; i<nbin; i++)
        {
                nDensOd_bin[i]=0; nDensOH_bin[i]=0; nDensEx_bin[i]=0;
		VelOd_bin[i]=0; VelOH_bin[i]=0; VelEx_bin[i]=0;
		Vel2Od_bin[i]=0; Vel2OH_bin[i]=0; Vel2Ex_bin[i]=0;
        }


	/*Reading the inputSOL to get index of target molecules*/
	Atom atom;
	vector<int> idx, idxOd, idxOH;
	vector<Atom> AtomGRO, OdGRO, OHGRO;
	int count=0, countOd=0, countOH=0;
	getline(inputSOL,line);//skip the 1st and 2nd line
	getline(inputSOL,line);
	while(!inputSOL.eof())
	{
		getline(inputSOL,line);
		if (line.empty())
		{
			break;
		}
		else
		{
			count++;
			atom.name = line.substr(10,5);
			linestream.clear(); linestream.str(atom.name);
			linestream >> word;
			atom.name = word;
			//Scan index of molecules whose velocity will be calculated
			for (int i=0; i<info_ary.size(); i++)
			{
				if(atom.name.find(info_ary[i].name) != string::npos)
				{
					AtomGRO.push_back(atom);	
					idx.push_back(count-1);//idx start from 0 in C
				}
			}
			//Scan index and coordinates of silica surface atoms, Od, OSil, and OGem
			if(atom.name == "Od")
			{
				atom.x = atof(line.substr(20,8).c_str());
				atom.y = atof(line.substr(28,8).c_str());
				atom.z = atof(line.substr(36,8).c_str());
				countOd++;
				OdGRO.push_back(atom);
				idxOd.push_back(countOd-1);
			}
			else if ( (atom.name == "OSil") || (atom.name == "OGem") )
			{
				atom.x = atof(line.substr(20,8).c_str());
				atom.y = atof(line.substr(28,8).c_str());
				atom.z = atof(line.substr(36,8).c_str());
				countOH++;
				OHGRO.push_back(atom);
				idxOH.push_back(countOH-1);
			}
		}
	}
	
	/*Open TRR file*/
	XDRFILE *TRR;
	TRR=xdrfile_open(argv[1],"r");
	int natoms, step, step2=10, nframes=0;
	float timeps, lambda;
	matrix box;
	read_trr_natoms(argv[1],&natoms);
	rvec *x = (rvec *)malloc(natoms * sizeof(rvec));
	rvec *v = (rvec *)malloc(natoms * sizeof(rvec));
	rvec *f = (rvec *)malloc(natoms * sizeof(rvec));
	int dtmd=100; //how many fs for each md step
	int stride=float(dt)/(float(dtmd)/1000.); //sample from every stride step

	/*Read data from TRR file*/
	vector<Atom> Atomcolumn;
	vector<vector<Atom> > AtomMatrix;
	for (int i=0; i<nAtom; i++)
	{
		AtomMatrix.push_back(Atomcolumn);
	}
	//Read the 1st frame (Store atoms of the initial time frame in the first nSOL row of AtomMatrix)
	read_trr(TRR, natoms, &step, &timeps, &lambda, box, x, v, f);
	nframes++;
	for (int i=0; i<nSOL; i++)
	{
		for (int j=0; j<nAtom; j++)
		{
			atom.name=AtomGRO[nAtom*i+j].name;
			atom.x = x[idx[nAtom*i+j]][0];
			atom.y = x[idx[nAtom*i+j]][1];
			atom.z = x[idx[nAtom*i+j]][2];
			atom.vx = v[idx[nAtom*i+j]][0];
			AtomMatrix.at(j).push_back(atom);
		}
		/*glue broken molecules*/
		comx=0; comy=0; comz=0;
		for (int j=0; j<nAtom; j++)
		{
			comx += AtomMatrix[j][i].x * info_ary[j].mass;
			comy += AtomMatrix[j][i].y * info_ary[j].mass;
			comz += AtomMatrix[j][i].z * info_ary[j].mass;
		}
		comx = comx/MW; comy=comy/MW; comz=comz/MW;
		for (int j=0; j<nAtom; j++)
		{
			AtomMatrix[j][i].x = AtomMatrix[j][i].x - round((AtomMatrix[j][i].x - comx)/box[0][0])*box[0][0];
			AtomMatrix[j][i].y = AtomMatrix[j][i].y - round((AtomMatrix[j][i].y - comy)/box[1][1])*box[1][1];
			AtomMatrix[j][i].z = AtomMatrix[j][i].z - round((AtomMatrix[j][i].z - comz)/box[2][2])*box[2][2];
		}
	}

	//Start reading from the 2nd frame. Delete the first nSOL row of AtomMatrix before the start of the next frame.
	while(1)
	{
		comx=0; comy=0; comz=0;
		read_trr(TRR, natoms, &step, &timeps, &lambda, box, x, v, f);
		if(step==step2)
		{
			cout << "\n";
			break;
		}
		else if (step/dtmd % int(stride) == 0)
		{
			step2=step;
			nframes++;
			if (int(timeps*10) %100 == 0)
			{
				cout << "\rReading step: " << step << " (" << timeps << " ps)" << flush;
				log << "Reading step: " << step << " (" << timeps << " ps)" << endl;
			}

			#pragma omp parallel num_threads(1)
			{	
				#pragma omp for private(atom,comx,comy,comz)
				for (int i=0; i<nSOL; i++)
				{
					for (int j=0; j<nAtom; j++)
					{
						atom.name=AtomGRO[nAtom*i+j].name;
						atom.x = x[idx[nAtom*i+j]][0];
						atom.y = x[idx[nAtom*i+j]][1];
						atom.z = x[idx[nAtom*i+j]][2];
						atom.vx = v[idx[nAtom*i+j]][0];
						#pragma omp critical
						{
							AtomMatrix.at(j).push_back(atom);
						}
						if (AtomMatrix[j][AtomMatrix[j].size()-1].name != info_ary[j].name)
						{
							cout << "The atomtype in the gro file:" << AtomMatrix[j][AtomMatrix[j].size()-1].name << " didn't match the atomtype in the mol file: " << info_ary[j].name << ". Job Abortted!" << endl;
							exit (EXIT_FAILURE);
						}
					}
					/*glue broken molecules*/
					for (int j=0; j<nAtom; j++)
					{
						comx += AtomMatrix[j][i+nSOL].x * info_ary[j].mass;
						comy += AtomMatrix[j][i+nSOL].y * info_ary[j].mass;
						comz += AtomMatrix[j][i+nSOL].z * info_ary[j].mass;
					}
					comx = comx/MW; comy=comy/MW; comz=comz/MW;
					for (int j=0; j<nAtom; j++)
					{
						AtomMatrix[j][i+nSOL].x = AtomMatrix[j][i+nSOL].x - round((AtomMatrix[j][i+nSOL].x - comx)/box[0][0])*box[0][0];
						AtomMatrix[j][i+nSOL].y = AtomMatrix[j][i+nSOL].y - round((AtomMatrix[j][i+nSOL].y - comy)/box[1][1])*box[1][1];
						AtomMatrix[j][i+nSOL].z = AtomMatrix[j][i+nSOL].z - round((AtomMatrix[j][i+nSOL].z - comz)/box[2][2])*box[2][2];
					}
				}

				//Check the number of atoms
				for (int n=0; n<nAtom; n++)
				{
					if (AtomMatrix[n].size()/2 != nSOL)
					{
						cout << "Number of " << info_ary[n].name << ":  " << AtomMatrix[n].size()/2 << " didn't match the number of nSOL:" << nSOL << endl;
						exit (EXIT_FAILURE);
					}
				}
			}	

			int N, binidx, binidx2, nsl, j, k;
			int Odid,OHid; //1 means the molecule is not classified as neighbors of Od/OH groups
			double dist_Od, dist_OH, r1, r2, d1, d2, d12, sl=box[2][2]/nbin, Vx; //define bin width
			//Velocity are calculated by displacement/dt, where the displacement  of each particle between the initial and final time frame
			#pragma omp parallel num_threads(4)
			{
				#pragma omp for private(j,k,N,nsl,binidx,binidx2,Odid,OHid,dist_Od,dist_OH,com1z,com1x,com2z,com2x,com1y,r1, r2, d1, d2, d12, Vx)
				for (int i=0; i<nSOL; i++)
				{
					com1z = 0; com1x=0; com1y=0;//initial time frame
					com2z = 0; com2x=0;//final time frame
					Vx=0;
					for (int n=0; n<nAtom; n++)
					{
						com1z += AtomMatrix[n][i].z * info_ary[n].mass;
						com2z += AtomMatrix[n][i+nSOL].z * info_ary[n].mass;
						com1x += AtomMatrix[n][i].x * info_ary[n].mass;
						com2x += AtomMatrix[n][i+nSOL].x * info_ary[n].mass;
						com1y += AtomMatrix[n][i].y * info_ary[n].mass;
						Vx += AtomMatrix[n][i].vx * info_ary[n].mass;
					}
					com1z = com1z/MW; com1x = com1x/MW; com1y = com1y/MW;
					//Calculate distOd and distOH
					Atom COM;
					COM.x = com1x; COM.y = com1y; COM.z = com1z;
					Odid=1, OHid=1;
					for (int iod=0; iod<idxOd.size(); iod++)
					{
						dist_Od = dist3D_pbc(OdGRO[iod], COM, box[0][0], box[1][1]);
						if (dist_Od < 0.3)
						{
							Odid=0;
							break;
						}
					}
					if (Odid ==1)
					{
						for (int iOH=0; iOH<idxOH.size(); iOH++)
						{
							dist_OH = dist3D_pbc(OHGRO[iOH], COM, box[0][0], box[1][1]);
							if (dist_OH < rcut)
							{
								OHid=0;
								break;
							}
						}
					}

					com2z = com2z/MW; com2x = com2x/MW;
					com2x = com2x + round((com1x-com2x)/box[0][0])*box[0][0]; //compute rectified trajectory in x direction
					Vx = Vx/MW;
					

					//Move all the z component of com into the box
					N = floor(abs(com1z)/box[2][2]);
					if (com1z > 0)
					{
						com1z += -N*box[2][2];
					}
					else if (com1z < 0)
					{
						com1z += (N+1)*box[2][2];
					}

					N = floor(abs(com2z)/box[2][2]);
					if (com2z > 0)
					{
						com2z += -N*box[2][2];
					}
					else if (com2z < 0)
					{
						com2z += (N+1)*box[2][2];
					}

					//Bin com of the ith molecule of the initial and final time frames
					for (j=0; j<nbin; j++)
					{
						if ( ( com1z >= sl*j ) && ( com1z < sl*(j+1) ) )
						{
							binidx = j;
							#pragma omp critical
							{
								if (Odid == 0)
								{
									nDensOd_bin[j] ++;
									Vel2Od_bin[j] += Vx;
								}
								else if (OHid == 0)
								{
									nDensOH_bin[j] ++;
									Vel2OH_bin[j] += Vx;
								}
								else
								{
									nDensEx_bin[j]++;
									Vel2Ex_bin[j] += Vx;
								}
							}
							break;
						}
					}
					for (j=0; j<nbin; j++)
					{
						if ( ( com2z >= sl*j ) && ( com2z < sl*(j+1) ) )
						{
							binidx2 = j;
							break;
						}
					}

					if (binidx == binidx2)
					{
						#pragma omp critical
						{
							if (Odid == 0)
							{
								VelOd_bin[binidx] += (com2x - com1x);
							}
							else if (OHid == 0)
							{
								VelOH_bin[binidx] += (com2x - com1x);
							}
							else
							{
								VelEx_bin[binidx] += (com2x - com1x);
							}
						}
					}
					else if (binidx2 > binidx)
					{
						d1 = sl - (com1z - sl*binidx);
						d2 = com2z - sl*binidx2;
						d12 = com2z - com1z;
						r1 = d1/d12;
						r2 = d2/d12;
						nsl = round ((d12-d1-d2)/sl);
						for (k=0; k<nsl; k++)
						{
							#pragma omp critical
							{
								if (Odid==0)
								{
									VelOd_bin[binidx+k+1] += (com2x - com1x)*sl/d12;
								}
								else if (OHid==0)
								{
									VelOH_bin[binidx+k+1] += (com2x - com1x)*sl/d12;
								}
								else
								{
									VelEx_bin[binidx+k+1] += (com2x - com1x)*sl/d12;
								}
							}
						}
						#pragma omp critical
						{
							if (Odid==0)
							{
								VelOd_bin[binidx] += (com2x - com1x)*r1;
								VelOd_bin[binidx2] += (com2x - com1x)*r2;
							}
							else if (OHid==0)
							{
								VelOH_bin[binidx] += (com2x - com1x)*r1;
								VelOH_bin[binidx2] += (com2x - com1x)*r2;
							}
							else
							{
								VelEx_bin[binidx] += (com2x - com1x)*r1;
								VelEx_bin[binidx2] += (com2x - com1x)*r2;
							}
						}
					}
					else
					{
						d1 = com1z - sl*binidx;
						d2 = sl - (com2z - sl*binidx2);
						d12 = com1z - com2z;
						r1 = d1/d12;
						r2 = d2/d12;
						nsl = round ((d12-d1-d2)/sl);
						for (k=0; k<nsl; k++)
						{
							#pragma omp critical
							{
								if(Odid==0)
								{
									VelOd_bin[binidx2+k+1] += (com2x - com1x)*sl/d12;
								}
								else if(OHid==0)
								{
									VelOH_bin[binidx2+k+1] += (com2x - com1x)*sl/d12;
								}
								else
								{
									VelEx_bin[binidx2+k+1] += (com2x - com1x)*sl/d12;
								}
							}
						}
						#pragma omp critical
						{
							if(Odid==0)
							{
								VelOd_bin[binidx] += (com2x - com1x)*r1;
								VelOd_bin[binidx2] += (com2x - com1x)*r2;
							}
							else if(OHid==0)
							{
								VelOH_bin[binidx] += (com2x - com1x)*r1;
								VelOH_bin[binidx2] += (com2x - com1x)*r2;
							}
							else
							{
								VelEx_bin[binidx2] += (com2x - com1x)*r2;
								VelEx_bin[binidx2] += (com2x - com1x)*r2;
							}
						}
					}

				}

			}// <--- End of loop over molecules

			//Empty AtomMatrix
			for (int n=0; n<nAtom; n++)
			{
				AtomMatrix[n].erase(AtomMatrix[n].begin(),AtomMatrix[n].begin()+nSOL);
			}
		}
	}/*<--- End of loop over time frames*/

        //Wrting outputs
        output1 << "# This file was created " << currentDateTime() << endl;
        output1 << "# Command line:" << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
        output1 << "@    title \"Number density of " << Sysname << "\"" << endl; //output header for OW number density
        output1 << "@    xaxis  label \"box (nm)\"" << endl;
        output1 << "@    yaxis  label \"Density (nm^-3)\"" << endl;
	output1 << "@TYPE xy \n@ view 0.15, 0.15, 0.75, 0.85" << endl;
	output1 << "@ legend on \n@ legend box on \n@ legend loctype view \n@ legend 0.78, 0.8 \n@ legend length 2" << endl;
	output1 << "@ s0 legend \"Od-type " << Sysname << " \"" << endl;
	output1 << "@ s1 legend \"OH-type " << Sysname << " \"" << endl;
	output1 << "@ s2 legend \"Other " << Sysname << " \"" << endl;
        output2 << "# This file was created " << currentDateTime() << endl;
        output2 << "# Command line:" << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
        output2 << "@    title \"Velocity profile of " << Sysname << "\"" << endl; //output header for OW number density
        output2 << "@    xaxis  label \"box (nm)\"" << endl;
        output2 << "@    yaxis  label \"V_x (nm/ps)\"" << endl;
        output2 << "@TYPE xy \n@ view 0.15, 0.15, 0.75, 0.85" << endl;
        output2 << "@ legend on \n@ legend box on \n@ legend loctype view \n@ legend 0.78, 0.8 \n@ legend length 2" << endl;
        output2 << "@ s0 legend \"Od-type " << Sysname << " \"" << endl;
	output2 << "@ s1 legend \"Od-type " << Sysname << " \"" << endl;
        output2 << "@ s2 legend \"OH-type " << Sysname << " \""<< endl;
	output2 << "@ s3 legend \"OH-type " << Sysname << " \"" << endl;
        output2 << "@ s4 legend \"Other " << Sysname << " \"" << endl;
	output2 << "@ s5 legend \"Other " << Sysname << " \"" << endl;

        double sum=0, sumVel=0;
	dV = box[0][0]*box[1][1]*box[2][2]/nbin;
        for (int i=0; i<nbin; i++)
        {
                output1 << setw(18) << setprecision(6) << fixed << box[2][2]/nbin*i;
                output1 << setw(18) << setprecision(6) << fixed << nDensOd_bin[i]/(nframes-1)/dV;
		output1 << setw(18) << setprecision(6) << fixed << nDensOH_bin[i]/(nframes-1)/dV; 
		output1 << setw(18) << setprecision(6) << fixed << nDensEx_bin[i]/(nframes-1)/dV << endl;
                output2 << setw(18) << setprecision(6) << fixed << box[2][2]/nbin*i;
                if (nDensOd_bin[i] != 0)
                {
                        output2 << setw(18) << setprecision(6) << fixed << VelOd_bin[i]/nDensOd_bin[i]/dt;
                        output2 << setw(18) << setprecision(6) << fixed << Vel2Od_bin[i]/nDensOd_bin[i];
                }
                else
                {
                        output2 << setw(18) << setprecision(6) << fixed << 0.0;
                        output2 << setw(18) << setprecision(6) << fixed << 0.0;
                }
                if (nDensOH_bin[i] != 0)
                {   
                        output2 << setw(18) << setprecision(6) << fixed << VelOH_bin[i]/nDensOH_bin[i]/dt;
                        output2 << setw(18) << setprecision(6) << fixed << Vel2OH_bin[i]/nDensOH_bin[i];
                }
                else
                {
                        output2 << setw(18) << setprecision(6) << fixed << 0.0; 
                        output2 << setw(18) << setprecision(6) << fixed << 0.0;
                }
                if (nDensEx_bin[i] != 0)
                {   
                        output2 << setw(18) << setprecision(6) << fixed << VelEx_bin[i]/nDensEx_bin[i]/dt;
                        output2 << setw(18) << setprecision(6) << fixed << Vel2Ex_bin[i]/nDensEx_bin[i] << endl;
                }
                else
                {
                        output2 << setw(18) << setprecision(6) << fixed << 0.0; 
                        output2 << setw(18) << setprecision(6) << fixed << 0.0 << endl;
                }



                //debug for mass density and polarization
                if (dbg==0)
                {
                        sum += (nDensOd_bin[i]+nDensOH_bin[i]+nDensEx_bin[i]);
                        sumVel += (VelOd_bin[i]+VelOH_bin[i]+VelEx_bin[i]);
                }
        }


	//Estimate cpu time for the calculations	
	time_t nEnd = time(NULL);//End timer
        cout<<"Elapsed time is :  "<< nEnd-nStart << " seconds " << endl;
        log << "Elapsed time is :  "<< nEnd-nStart << " seconds " << endl;

        //Output dbg inform
        if (dbg==0)
        {
                log << "\nDebug information: " << endl;
                log << "Number of frames: " << nframes << endl;
                log << "Total number of molecules: ";
                log << setw(12) << setprecision(6) << fixed << sum/(nframes-1) << endl;
                log << "Average velocity (nm/ps): ";
                log << setw(12) << setprecision(6) << fixed << sumVel/sum/dt << endl;
        }

	//Relase mem ry for dynamic arrays
        free(x); free(v); free(f); 
	free(nDensOd_bin); free(VelOd_bin); free(Vel2Od_bin);
	free(nDensOH_bin); free(VelOH_bin); free(Vel2OH_bin);
	free(nDensEx_bin); free(VelEx_bin); free(Vel2Ex_bin);

	/*Close Files*/
	xdrfile_close(TRR);
	inputSOL.close(); molinfo.close(); output1.close(); output2.close(); log.close();
	return 0;
}
