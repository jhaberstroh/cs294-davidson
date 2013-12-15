#ifndef __LARGEMATRIX_H_DEFINED__
#define __LARGEMATRIX_H_DEFINED__
#include "Real.H"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#define USE_BINARY_FILE

namespace Mtx
{

typedef std::vector<Real> colVect;
typedef std::vector<Real> rowVect;

inline Real operator*(rowVect a, colVect b)
	{
		Real prod = 0;
		assert(a.size() == b.size());
		for (int i = 0 ; i < a.size() ; i++)
			{
				prod += a[i] * b[i];
			}
		return prod;
	}

enum RWFlags
	{
		FRead = 1,
		FWrite = 2,
		Overwrite = 4
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
		int m_num_rows;
	public: 
		explicit LargeMatrix(std::string a_filename, int a_num_rows, int a_iostate = FRead|FWrite, int a_nums_in_memory = 1000000);

		static int EstimateMaxSize(int a_num_rows);
		static void GenerateRandom(std::string a_filename, int a_num_rows, int a_ioflags = FWrite|Overwrite);

		colVect operator*(colVect vec);

};

}


#endif//__LARGEMATRIX_H_DEFINED__
