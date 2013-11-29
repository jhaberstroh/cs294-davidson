#ifndef __LARGEMATRIX_H_DEFINED__
#define __LARGEMATRIX_H_DEFINED__
#include <cassert>
#include <string>
#include <iostream>

namespace Mtx
{

enum RWFlags
	{
		FRead = 1,
		FWrite = 2
	};

inline RWFlags operator|(RWFlags a, RWFlags b)
	{
		return static_cast<RWFlags>(static_cast<int>(a) | static_cast<int>(b));
	}

class LargeMatrix
{
	private:
		std::string m_filename;
		int m_iostate;
		int m_nums_in_memory;
	public: 
		explicit LargeMatrix(std::string a_filename, int a_iostate = FRead|FWrite, int a_nums_in_memory = 1000000);

		static int EstimateMaxSize(int a_num_rows);
		static void GenerateRandom(std::string a_filename, int a_num_rows);

};

}


#endif//__LARGEMATRIX_H_DEFINED__
