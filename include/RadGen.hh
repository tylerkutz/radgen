#ifndef RADGEN_HH
#define RADGEN_HH

class RadGen{

public:
	RadGen();
	~RadGen();
	
	void Radiate(float*, float*, float*, float*, float*, float*, float*);

private:

	int target;
	float eBeam;
	int LST40;
	int ixytst;

};

#endif

