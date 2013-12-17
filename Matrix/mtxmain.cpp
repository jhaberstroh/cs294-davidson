#include "Real.H"
#include "LargeMatrix.H"
#include "CH_Timer.H"


int main(){
	CH_TIMERS("MAIN");
	CH_TIMER("many_mult", t2);

	Mtx::Format fmt = Mtx::FMT_BIN;
	
	int rows;
	std::string mtx_name;

	if (fmt == Mtx::FMT_TXT)
		{
			rows = 1000;
			mtx_name = "newmtx.txt";
			Mtx::LargeMatrix::GenerateRandom(mtx_name, rows);
		}

	if (fmt == Mtx::FMT_BIN)
		{
			rows = 25 * 25;
			mtx_name = "../parsing/matx_out_test.dat";
		}
	

	Mtx::LargeMatrix x(mtx_name, rows);
	Mtx::colVect a(rows, 0);
	std::cout << "A size: " << a.size() << std::endl;
	for (int i = 0 ; i < a.size() ; i++)
		{	
			a[i] = ((i + 7) * ((i+4) * 123)) % 37  * 111 % 10;
		}
	Mtx::print(a);


	Mtx::colVect b(a);
	std::cout << "B size: " << b.size() << std::endl;
	CH_START(t2);
	for (int i = 0 ; i < 100 ; i++)
		{
			b = x*b;
		}
	CH_STOP(t2);
	Mtx::print(b);

}
