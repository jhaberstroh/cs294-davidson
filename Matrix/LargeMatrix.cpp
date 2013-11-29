#include "LargeMatrix.H"
#include "Real.H"
#include <iostream>


Mtx::LargeMatrix::LargeMatrix(std::string a_filename, int a_iostate, int a_nums_in_memory)
	{
		m_filename = a_filename;
		m_iostate = a_iostate;
		assert(m_iostate | FRead);
		assert(m_iostate | FWrite);
		std::cout << "Succeeded assertions!"<<std::endl;
	}

int Mtx::LargeMatrix::EstimateMaxSize(int a_num_rows)
	{
		return a_num_rows * a_num_rows * 32;
	}


void Mtx::LargeMatrix::GenerateRandom(std::string a_filename, int a_num_rows)
	{
		std::cout << "Nothing yet implemented for GenerateRandom" << std::endl;
	}
