#include "JDTemplate.H"
#include "LargeMatrix.H"
#include "lapacke.h"

#define N 5
#define LDA N

int main()
	{
		Mtx::Format fmt = Mtx::FMT_OPT_TEI;
		int rows = 25 * 25;
		std::string mtx_name = "../parsing/matx_out_test.mot";
		Mtx::LargeMatrix mt(mtx_name, rows, fmt, Mtx::FRead|Mtx::FWrite, 100000);


		JDRoutine(mt);
	}
