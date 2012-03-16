
class AWGN
{
public:
//	AWGN();
//	~AWGN();

	//
	double AWGN_sigma;
	double AWGN_mean;
	long int MBIG;//=0x0fffffff;
	long int MSEED;//=0x03c04ffd;
	int MZ;//=0;
	double FAC;
	long * seed_pointer;
	long seed_1;
	void Gauss_initial();
	double Crand(long *seed);
	double Gauss(long *seed);
	double AWGN_noise(double sigma, double mean);
	
protected:
private:
};