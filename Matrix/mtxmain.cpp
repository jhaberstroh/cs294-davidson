#include "Real.H"
#include "LargeMatrix.H"
#include "CH_Timer.H"
#include <fstream>
#include <iostream>

void mtxTest(Mtx::Format fmt)
	{
	
		int rows;
		std::string mtx_name;
	
		if (fmt == Mtx::FMT_TXT)
			{
				rows = 25*25;
				mtx_name = "newmtx.txt";
				Mtx::LargeMatrix::GenerateRandom(mtx_name, rows);
			}
	
		if (fmt == Mtx::FMT_BIN)
			{
				rows = 25 * 25;
				mtx_name = "../parsing/matx_out_test.dat";
			}
		
		if (fmt == Mtx::FMT_OPT_TEI)
			{
				rows = 25 * 25;
				mtx_name = "../parsing/matx_out_test.mot";
			}
		
	
		Mtx::LargeMatrix x(mtx_name, rows, fmt);
		Mtx::colVect a(rows, 0);
		std::cout << "A size: " << a.size() << std::endl;
		for (int i = 0 ; i < a.size() ; i++)
			{	
				a[i] = ((i + 7) * ((i+4) * 123)) % 37  * 111 % 10;
			}
		Mtx::print(a);
	
	
		Mtx::colVect b(a);
		std::cout << "B size: " << b.size() << std::endl;
		for (int i = 0 ; i < 5 ; i++)
			{
				b = x*b;
			}
		Mtx::print(b);
	}


int main()
	{
		CH_TIMERS("MAIN");
		CH_TIMER("many_mult_bin", t1);
		CH_TIMER("many_mult_opt", t2);
		CH_START(t1);
		mtxTest(Mtx::FMT_BIN);
		CH_STOP(t1);
		CH_START(t2);
		mtxTest(Mtx::FMT_OPT_TEI);
		CH_STOP(t2);


		std::ifstream test_file;
		test_file.open("../parsing/test_bin.bin", std::ifstream::binary);
		char buffer[8];
		test_file.read(buffer,8);
		std::cout << "The value of pi is "<<*(Real*)buffer <<std::endl;
		test_file.close();



	}
