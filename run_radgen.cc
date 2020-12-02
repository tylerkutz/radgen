#include <iostream>
#include <fstream>

#include "TMath.h"
#include "TVector3.h"

#include "RadGen.hh"

using namespace std;

int main(int argc, char ** argv){

	int numargs = 3;

	if (argc < numargs) {
		exit(-1);
	}

	RadGen* fRadGen = new RadGen();

	cout << "Made generator..." << endl;


	double Eb = 10.2; // GeV
	double Me = 0.000511; // GeV 
	double Mp = 0.938272;
	double deg = TMath::DegToRad();
	TVector3 k0 = TVector3(0.,0.,sqrt(Eb*Eb - Me*Me));

//	double minQ2 = 1.5;	// GeV^2
//	double maxQ2 = 10.0;

	double minQ2 = atof(argv[1]);
	double maxQ2 = atof(argv[2]);
	double dQ2 = 0.01;

	ofstream radOut("radgen.dat");

	for(double thisQ2 = minQ2; thisQ2 <= maxQ2; thisQ2+=dQ2) {

		double minW = 1.5;  	
		double maxW = sqrt(Mp*Mp + 2.*Mp*(Eb - 1.) - thisQ2);	
		double dW = 0.01;

		for(double thisW = minW; thisW <= maxW; thisW+=dW) {

		
			cout << "Q2 = " << thisQ2 << " W = " << thisW << endl;

			double nu = (thisW*thisW + thisQ2 - Mp*Mp)/(2.*Mp);
			double Ep = Eb - nu;
			double pe = sqrt(Ep*Ep - Me*Me);

			double cosTheta = 1. - thisQ2/(2.*Eb*Ep);
			double theta = acos(cosTheta);

			TVector3 electron;
			electron.SetMagThetaPhi(pe, theta, 10.*deg);

			TVector3 q = electron - k0;

			float vpgen[4] = {(float)q.x(), (float)q.y(), (float)q.z(), (float)nu};
			float vprad[4];
			float rprad[4];
			float Q2true;
			float Utrue;
			float weight;

			float ebf = (float)Eb;	
			fRadGen->Radiate(&ebf, vpgen, vprad, rprad, &Q2true, &Utrue, &weight); 

			radOut << thisQ2 << "\t" << thisW << "\t" << weight << "\n";


		}
	}

	radOut.close();

	cout << "Done!" << endl;
	return 0;

}
