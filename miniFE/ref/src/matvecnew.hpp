#ifndef _matvecnew_hpp_
#define _matvecnew_hpp_

#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <SparseMatrix_functions.hpp>

namespace miniFE {

template<typename MatrixType, typename VectorType>
struct matvecnew {
	
	typedef typename MatrixType::ScalarType Scalar;
	typedef typename MatrixType::LocalOrdinalType LO;

	int off[27];
	Scalar cf[27];
	std::vector<std::pair<LO,LO>> runs;
	std::vector<LO> rest;

	explicit matvecnew(const MatrixType &A)
	{
		const LO n = A.rows.size();
		const LO *ro = &A.row_offsets[0];
		const auto *ci = &A.packed_cols[0];
		const Scalar *v = &A.packed_coefs[0];
		
		for (int i = 0; i < 27; ++i)
		{
			off[i] = 0;
			cf[i] = 0;
		}

		LO ref = -1;
		for (LO r = n/2; r < n && ref < 0; ++r)
		{
			if (ro[r+1] - ro[r] != 27)
				continue;

			bool ok = true;
			for (int i = 0; i < 27; ++i)
				ok = ok && (ci[ro[r]+i] < n);

			if (ok)
				ref = r;
		}

		Scalar vmax = 0;
		if (ref >= 0) {
			for (int i = 0; i < 27; ++i)
			{
				off[i] = (int)ci[ro[ref]+i] - ref;
				cf[i] = v[ro[ref]+i];
				vmax = std::max(vmax, std::abs(cf[i]));
			}
		}

		LO nst = 0;
		LO run_begin = -1;
		for (LO r = 0; r < n; ++r)
		{
			bool match = ref >= 0 && (ro[r+1]-ro[r] == 27);
			for (int i = 0; match && i < 27; ++i)
				match = ((int)ci[ro[r]+i] - r == off[i]) && (std::abs(v[ro[r]+i] - cf[i]) <= 1e-10 * vmax);

			if (match)
			{
				++nst;
				if (run_begin < 0)
					run_begin = r;
			}
			else
			{
				if (run_begin >= 0)
				{
					runs.push_back(std::make_pair(run_begin, r));
					run_begin = -1;
				}
				rest.push_back(r);
			}
		}
		if (run_begin >= 0)
			runs.push_back({run_begin, n});
		
		double frac = n > 0 ? double(nst) / double(n) : 0.0, minfrac = frac;
		int myproc = 0;

		#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
		MPI_Allreduce(&frac, &minfrac, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		#endif
		if (myproc == 0)
		{
			std::streamsize old = std::cout.precision(6);
			std::cout << "matvecnew: coverage rank0 = " << frac << ", min over ranks = " << minfrac << '\n';
			
			for (int i = 0; i < 27; ++i)
				std::cout << "  off " << std::setw(8) << off[i] << "  coef " << cf[i] << '\n';

			std::cout.precision(old);
		}
	}

	void operator()(MatrixType &A, VectorType &x, VectorType &y)
	{
		exchange_externals(A, x);
		const Scalar *__restrict xc = &x.coefs[0];
		Scalar *__restrict yc = &y.coefs[0];
		int oo[27];
		Scalar cc[27];
		for (int i = 0; i < 27; ++i)
		{
			oo[i] = off[i];
			cc[i] = cf[i];
		}

		for (size_t i = 0; i < runs.size(); ++i)
		{
			for (LO j = runs[i].first; j < runs[i].second; ++j)
			{
				Scalar s = 0;
				
				#pragma GCC unroll 27
				for (int k = 0; k < 27; ++k)
					s += cc[k] * xc[j + oo[k]];

				yc[j] = s;
			}
		}

		const LO *ro = &A.row_offsets[0];
		const auto *ci = &A.packed_cols[0];
		const Scalar *v = &A.packed_coefs[0];

		for (size_t i = 0; i < rest.size(); ++i)
		{
			LO r = rest[i];
			Scalar s = 0;

			for (LO j = ro[r]; j < ro[r+1]; ++j)
				s += v[j] * xc[ci[j]];

			yc[r] = s;
		}
	}

};

}

#endif
