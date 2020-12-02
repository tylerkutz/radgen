#include <iostream>

#include "RadGen.hh"

using namespace std;

extern "C" {

	void radgen_init_(int* CTARGET, float* ebeam, int* LST40, int* ixytst);
	void radgen_(float* e1, float* vpgen, float* vprad, float* phrad, float* q2tr, float*anutr, float* weight);

}

RadGen::RadGen() {

	cout << "Inside RadGen constructor..." << endl;
	
	target = 2;	// 1 = proton, 2 = deuterium, 3 = helium
	eBeam = 10.6; 	// GeV
	LST40 = 0;
	ixytst = -1;

	radgen_init_(&target, &eBeam, &LST40, &ixytst); 

}

RadGen::~RadGen() {

}

void RadGen::Radiate(float* e1, float* vpgen, float* vprad, float* phrad, float* q2tr, float* anutr, float* weight) {

	radgen_(e1, vpgen, vprad, phrad, q2tr, anutr, weight);

}


