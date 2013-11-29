#include "Real.H"
#include "LargeMatrix.H"



int main(){
	int rows = 10;
	std::string mtx_name = "newmtx.txt";

	Mtx::LargeMatrix::GenerateRandom(mtx_name, rows);
	Mtx::LargeMatrix x(mtx_name, rows);
	Mtx::colVect a(rows, 0);
	for (int i = 0 ; i < a.size() ; i++)
		{	
			a[i] = ((i + 7) * ((i+4) * 123)) % 37  * 111 % 10;
		}
	for (int i = 0 ; i < a.size() ; i++)
		{	
			std::cout << a[i] << "\t";
		}
	std::cout << std::endl;

	Mtx::colVect b(x*a);
	for (int i = 0 ; i < b.size() ; i++)
		{	
			std::cout << b[i] << "\t";
		}
	std::cout << std::endl;

}
