#include "Real.H"
#include "LargeMatrix.H"
#include "CH_Timer.H"



int main(){
	CH_TIMERS("MAIN");
	CH_TIMER("many_mult", t2);

#ifdef USE_BINARY_FILE
	int rows = 25 * 25;
	std::string mtx_name = "../parsing/matx_out_test.dat";
#endif

#ifndef USE_BINARY_FILE
	int rows = 1000;
	std::string mtx_name = "newmtx.txt";
	Mtx::LargeMatrix::GenerateRandom(mtx_name, rows);
#endif


	Mtx::LargeMatrix x(mtx_name, rows);
	Mtx::colVect a(rows, 0);
	std::cout << "A size: " << a.size() << std::endl;
	for (int i = 0 ; i < a.size() ; i++)
		{	
			a[i] = ((i + 7) * ((i+4) * 123)) % 37  * 111 % 10;
		}
	for (int i = 0 ; i < a.size() ; i++)
		{	
			std::cout << a[i] << " ";
		}
	std::cout << std::endl;

	Mtx::colVect b(a);
	std::cout << "B size: " << b.size() << std::endl;
	CH_START(t2);
	for (int i = 0 ; i < 1 ; i++)
		{
			b = x*b;
		}
	CH_STOP(t2);
	for (int i = 0 ; i < b.size() ; i++)
		{	
			std::cout << b[i] << " ";
		}
	std::cout << std::endl;

}
