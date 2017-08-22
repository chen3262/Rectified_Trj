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
	double dV, com1z, com2z, com1x, com2x;
	
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
                if (line.empty())
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

	//Output system information to log
        log << "System information: " << endl;
        log << "Molecule:       " << Sysname << endl;
        log << "nSOL:   " << nSOL << endl;
        log << "nAtom:  " << info_ary.size() << endl;
        log << "atom tupes: " ;
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
        log << "Skipped steps:  " << skip << "\n"<< endl;

        //Declare arrays to store number of molecules in each bin over time frames
        double *nDens_bin, *Vel_bin, *Vel2_bin;
        nDens_bin = (double *)malloc(nbin * sizeof(double));
        Vel_bin = (double *)malloc(nbin * sizeof(double));
        Vel2_bin = (double *)malloc(nbin * sizeof(double));
        for (int i=0; i<nbin; i++)
        {
                nDens_bin[i]=0; Vel_bin[i]=0; Vel2_bin[i]=0;
        }


	/*Reading the inputSOL to get index of target molecules*/
	Atom atom;
	vector<int> idx;
	vector<Atom> AtomGRO;
	int count=0;
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
			for (int i=0; i<info_ary.size(); i++)
			{
				if(atom.name.find(info_ary[i].name) != string::npos)
				{
					AtomGRO.push_back(atom);	
					idx.push_back(count-1);//idx start from 0 in C
				}
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
			comz += AtomMatrix[j][i].z * info_ary[j].mass; //reuired in all directions for isotropic systems
		}
		comx = comx/MW; comy=comy/MW; comz=comz/MW;
		for (int j=0; j<nAtom; j++)
		{
			AtomMatrix[j][i].x = AtomMatrix[j][i].x - round((AtomMatrix[j][i].x - comx)/box[0][0])*box[0][0];
			AtomMatrix[j][i].y = AtomMatrix[j][i].y - round((AtomMatrix[j][i].y - comy)/box[1][1])*box[1][1];
			AtomMatrix[j][i].z = AtomMatrix[j][i].z - round((AtomMatrix[j][i].z - comz)/box[2][2])*box[2][2]; //reuired in all directions for isotropic systems
		}
	}

	//Start reading from the 2nd frame. Delete the first nSOL row of AtomMatrix before the start of the next frame.
	while(1)
	{
		read_trr(TRR, natoms, &step, &timeps, &lambda, box, x, v, f);
		if(step==step2)
		{
			cout << "\n";
			//cout << "Total steps: " << step << " (" << timeps << " ps)" << endl;
			//cout << "Frames used: " << nframes << endl;
			//cout << "natoms:        " << natoms << endl;
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
					comz += AtomMatrix[j][i+nSOL].z * info_ary[j].mass; //reuired in all directions for isotropic systems
				}
				comx = comx/MW; comy=comy/MW; comz=comz/MW;
				for (int j=0; j<nAtom; j++)
				{
					AtomMatrix[j][i+nSOL].x = AtomMatrix[j][i+nSOL].x - round((AtomMatrix[j][i+nSOL].x - comx)/box[0][0])*box[0][0];
					AtomMatrix[j][i+nSOL].y = AtomMatrix[j][i+nSOL].y - round((AtomMatrix[j][i+nSOL].y - comy)/box[1][1])*box[1][1];
					AtomMatrix[j][i+nSOL].z = AtomMatrix[j][i+nSOL].z - round((AtomMatrix[j][i+nSOL].z - comz)/box[2][2])*box[2][2]; ////reuired in all directions for isotropic systems
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

			int N, binidx, binidx2, nsl, j, k;
			double r1, r2, d1, d2, d12, sl=box[2][2]/nbin, Vx; //define bin width
			//Velocity are calculated by displacement/dt, where the displacement  of each particle between the initial and final time frame
			#pragma omp parallel num_threads(4)
			{
				#pragma omp for private(j,k,N,nsl,binidx,binidx2,com1z,com1x,com2z,com2x,r1, r2, d1, d2, d12, Vx)
				for (int i=0; i<nSOL; i++)
				{
					com1z = 0; com1x=0;//initial time frame
					com2z = 0; com2x=0;//final time frame
					Vx=0;
					for (int n=0; n<nAtom; n++)
					{
						com1z += AtomMatrix[n][i].z * info_ary[n].mass;
						com2z += AtomMatrix[n][i+nSOL].z * info_ary[n].mass;
						com1x += AtomMatrix[n][i].x * info_ary[n].mass;
						com2x += AtomMatrix[n][i+nSOL].x * info_ary[n].mass;
						Vx += AtomMatrix[n][i].vx * info_ary[n].mass;
					}
					com1z = com1z/MW; com1x = com1x/MW;
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
								nDens_bin[j] ++;
								Vel2_bin[j] += Vx;
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
							Vel_bin[binidx] += (com2x - com1x);
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
								Vel_bin[binidx+k+1] += (com2x - com1x)*sl/d12;
							}
						}
						#pragma omp critical
						{
							Vel_bin[binidx] += (com2x - com1x)*r1;
							Vel_bin[binidx2] += (com2x - com1x)*r2;
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
								Vel_bin[binidx2+k+1] += (com2x - com1x)*sl/d12;
							}
						}
						#pragma omp critical
						{
							Vel_bin[binidx] += (com2x - com1x)*r1;
							Vel_bin[binidx2] += (com2x - com1x)*r2;
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
        output2 << "# This file was created " << currentDateTime() << endl;
        output2 << "# Command line:" << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
        output2 << "@    title \"Velocity profile of " << Sysname << "\"" << endl; //output header for OW number density
        output2 << "@    xaxis  label \"box (nm)\"" << endl;
        output2 << "@    yaxis  label \"V_x (nm/ps)\"" << endl;

        double sum=0, sumVel=0;
	dV = box[0][0]*box[1][1]*box[2][2]/nbin;
        for (int i=0; i<nbin; i++)
        {
                output1 << setw(18) << setprecision(6) << fixed << box[2][2]/nbin*i;
                output1 << setw(18) << setprecision(6) << fixed << nDens_bin[i]/(nframes-1)/dV << endl;
                output2 << setw(18) << setprecision(6) << fixed << box[2][2]/nbin*i;
                if (nDens_bin[i] != 0)
                {
                        output2 << setw(18) << setprecision(6) << fixed << Vel_bin[i]/nDens_bin[i]/dt;// << endl;
                        output2 << setw(18) << setprecision(6) << fixed << Vel2_bin[i]/nDens_bin[i] << endl;
                }
                else
                {
                        output2 << setw(18) << setprecision(6) << fixed << 0.0; // << endl;
                        output2 << setw(18) << setprecision(6) << fixed << 0.0 << endl;
                }

                //debug for mass density and polarization
                if (dbg==0)
                {
                        sum += nDens_bin[i];
                        sumVel += Vel_bin[i];
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
        free(x); free(v); free(f); free(nDens_bin); free(Vel_bin); free(Vel2_bin);

	/*Close Files*/
	xdrfile_close(TRR);
	inputSOL.close(); molinfo.close(); output1.close(); output2.close(); log.close();
	return 0;
}
