#ifndef __LARGEMATRIX_H_DEFINED__
#define __LARGEMATRIX_H_DEFINED__
#include "Real.H"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

namespace Mtx
{

typedef std::vector<Real> colVect;
typedef std::vector<Real> rowVect;

inline void print(colVect a)
	{
		for (int i = 0 ; i < a.size() ; i++)
			{
				if (i!=0)
					{
						std::cout << ",";
					}
				std::cout << a[i];
			}
		std::cout << std::endl;
	}


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
enum Format
	{
		FMT_TXT = 1, //.txt
		FMT_BIN = 2, //.dat
		FMT_OPT_OEI = 3,  //.moo
		FMT_OPT_TEI = 4   //.mot
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
		Format m_fmt;

		void TxtMult(colVect& output, const colVect& in) const;
		void BinMult(colVect& output, const colVect& in) const;
		void OptMult(colVect& output, const colVect& in) const;


	public: 
		explicit LargeMatrix(std::string a_filename, int a_num_rows, Format a_fmt = FMT_BIN, int a_iostate = FRead|FWrite, int a_nums_in_memory = 1000000);

		static int EstimateMaxSize(int a_num_rows);
		static void GenerateRandom(std::string a_filename, int a_num_rows, int a_ioflags = FWrite|Overwrite);

		colVect operator*(colVect vec) const;
		int size() const{ return m_num_rows; }

};

}


#endif//__LARGEMATRIX_H_DEFINED__
