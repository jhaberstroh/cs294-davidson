#include "LargeMatrix.H"
#include "Real.H"
#include "CH_Timer.H"
#include <boost/filesystem.hpp>
#include <gsl/gsl_rng.h>
#include <cblas.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

namespace Mtx
{

LargeMatrix::LargeMatrix(std::string a_filename, int a_num_rows, Format a_fmt, int a_iostate, int a_nums_in_memory)
	{
		m_filename = a_filename;
		m_fmt = a_fmt;
		m_iostate = a_iostate;
		m_num_rows = a_num_rows;
		m_nums_in_memory = a_nums_in_memory;
		assert(m_iostate | FRead);
		assert(m_iostate | FWrite);
		std::cout << "Succeeded assertions!"<<std::endl;
	}

int LargeMatrix::EstimateMaxSize(int a_num_rows)
	{
		return a_num_rows * a_num_rows * 32;
	}


void LargeMatrix::GenerateRandom(std::string a_filename, int a_num_rows, int a_ioflags)
	{
		//std::cout << "Nothing yet implemented for GenerateRandom" << std::endl;
		std::ofstream mtx_file;

		gsl_rng * r;
		const gsl_rng_type * T;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, 90210);
		std::vector<Real> rand_row(a_num_rows, 0);
		std::cout << "Generator type: " << gsl_rng_name(r) << std::endl;
		
		mtx_file.open(a_filename.c_str());
		for (int i = 0 ; i < a_num_rows ; i++)
			{	
				for (int j = 0 ; j < a_num_rows ; j++)
					{
						rand_row[j] = gsl_rng_uniform(r);
					}
				for (int j = 0 ; j < a_num_rows ; j++)
					{
						mtx_file << rand_row[j] << ",";
					}
				mtx_file << std::endl;
			}
		mtx_file.close();

		gsl_rng_free(r);
	}



void LargeMatrix::TxtMult(colVect& output, const colVect& arg)
	{

		CH_TIMERS("Multiplication");
		CH_TIMER("mult_file", t1);

		std::ifstream mtx_file;
		CH_START(t1);
		mtx_file.open(m_filename.c_str());
		//std::cout << "MATRIX LOADING: " << std::endl;
		std::stringstream ss;
		for (int i = 0 ; i < m_num_rows ; i++)
			{
				rowVect row(m_num_rows,0);
				std::string line;

				if (mtx_file.is_open())
					{
						getline(mtx_file,line);
					}
				ss.str(line);

				std::string remain;
				for (int j = 0 ; j < m_num_rows ; j++)
					{
						ss >> row[j];
						//std::cout << row[j] << "  \t";
						ss.ignore(256,',');
					}
				//std::cout << std::endl;
				
				output[i] = row*arg;
			}
		CH_STOP(t1);
		mtx_file.close();
	}

void LargeMatrix::BinMult(colVect& output, const colVect& arg)
	{
		CH_TIMERS("Multiplication");
		CH_TIMER("file_read", t1);
		CH_TIMER("line_splice", t2);
		CH_TIMER("vec_multiply", t3);
		//std::cout << "PERFORMING MULTIPLICATION!" << std::endl;
	
		std::ifstream mtx_file;
		mtx_file.open(m_filename.c_str(), std::ifstream::binary);
		//std::cout << "MATRIX LOADING: " << std::endl;
	
		int Nchar_double = 8;
		assert(m_num_rows = 25*25);
		int line_size = m_num_rows * Nchar_double;
		char buffer[line_size + 1];
		for (int i = 0 ; i < m_num_rows ; i++)
			{
				assert(mtx_file.is_open());
				CH_START(t1);
	
				rowVect row(m_num_rows,0);
				mtx_file.read(buffer, line_size + 1);
				if (!mtx_file)
					{
						mtx_file.close();
						std::cout << "FATAL ERROR: could not read line from file" << std::endl;
						abort();
					}
				assert(buffer[line_size] == '\n');
				CH_STOP(t1);
	
				CH_START(t2);
				for (int j = 0 ; j < m_num_rows ; j++)
					{
						int displace = Nchar_double * j;
	
						//assert(displace + Nchar_double  < line_size);
						
						row[j] = *(Real*) (buffer + Nchar_double * j);
					}
				CH_STOP(t2);
				
				CH_START(t3);
				output[i] = row*arg;
				CH_STOP(t3);
			}
		mtx_file.close();
	}

void LargeMatrix::OptMult(colVect& output, const colVect& arg)
	{
		CH_TIMERS("Multiplication");
		CH_TIMER("trivial", t0);
		CH_TIMER("file_read", t1);
		CH_TIMER("line_splice", t2);
		CH_TIMER("vec_multiply", t3);
		//std::cout << "PERFORMING MULTIPLICATION!" << std::endl;


		std::ifstream mtx_file;
		mtx_file.open(m_filename.c_str(), std::ifstream::binary);
		//std::cout << "MATRIX LOADING: " << std::endl;
	
		char header_buffer[19];
		mtx_file.read(header_buffer, 19);
		if (!mtx_file)
			{
				mtx_file.close();
				std::cout << "FATAL ERROR: could not read line from file" << std::endl;
				abort();
			}
		std::cout << "HEADER:: " << header_buffer << std::endl;
		unsigned int dat_size = *(int*) (header_buffer + 9);
		unsigned int num_row = *(int*) (header_buffer + 13);

		std::cout << "Dat Size:: " << dat_size << std::endl;
		std::cout << "Num rows:: " << num_row << std::endl;




		int max_lines_read = m_nums_in_memory / m_num_rows;
		//std::cout << "Max nums in memory: " << m_nums_in_memory << std::endl;
		//std::cout << "Max lines to read: " << max_lines_read << std::endl;

		int Nchar_double = 8;
		assert(m_num_rows = 25*25);
		int line_size = m_num_rows * Nchar_double;


		char buffer[ (line_size + 1) * max_lines_read];
		for (int i = 0 ; i < m_num_rows ; i+=max_lines_read)
			{
				CH_START(t0);
				assert(mtx_file.is_open());
				CH_STOP(t0);
				
				CH_START(t1);
				int num_lines_to_read = std::min(max_lines_read, m_num_rows - i);
				int num_bytes_to_read = (line_size + 1) * num_lines_to_read;

				//std::cout << "Reading " << num_lines_to_read << " lines.\n Buffer call to read " << num_bytes_to_read << " Bytes." << std::endl;
	
				rowVect row(m_num_rows,0);
				mtx_file.read(buffer, (line_size + 1) * num_lines_to_read);
				if (!mtx_file)
					{
						mtx_file.close();
						std::cout << "FATAL ERROR: could not read line from file" << std::endl;
						abort();
					}
				CH_STOP(t1);
	
				for (int lin = 0 ; lin < num_lines_to_read ; lin++)
					{
						CH_START(t2);
						int line_displace = lin * (line_size + 1);
						assert(buffer[line_displace + line_size] == '\n');
						//if (buffer[line_displace + line_size] != '\n')
						//	{
						//		std::cout << "Line "<< lin << " failed to pass newline assertion" << std::endl;
						//		abort();
						//	}
						Real* row_vect = (Real*) (buffer + line_displace);
						for (int j = 0 ; j < m_num_rows ; j++)
							{
								int float_displace = (j * Nchar_double) + line_displace;
	
								//assert(displace + Nchar_double  < line_size);
								row[j] = *(Real*) (buffer + float_displace);
								//if (i == 0)
								//	{
								//		if (j != 0)
								//			{
								//				std::cout << ",";
								//			}
								//		std::cout << row[j] ;
								//	}
							}
						if (i == 0)
							{
								//std::cout << std::endl;
							}
						CH_STOP(t2);
						
						CH_START(t3);
						output[i+lin] = cblas_ddot(m_num_rows, row_vect, 0, &*arg.begin(), 0);
						output[i+lin] = row*arg;
						CH_STOP(t3);
					}
			}
		mtx_file.close();
	}

colVect LargeMatrix::operator*(colVect arg)
	{
		assert(arg.size() == m_num_rows);
		colVect dot(arg.size(), 0);

		if (m_fmt == FMT_TXT)
			{
				TxtMult(dot, arg);
			}
		if (m_fmt == FMT_BIN)
			{
				BinMult(dot, arg);
			}
		if (m_fmt == FMT_OPT_TEI)
			{
				OptMult(dot, arg);
			}
		return dot;
	}


} //namespace Mtx 
