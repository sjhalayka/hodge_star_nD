#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <sstream>
#include <algorithm>
#include <array>
#include <chrono>
#include <set>

using namespace std;


// https://claude.ai/chat/3caf4077-28b5-497f-b704-1b0c336a104d
// https://codereview.stackexchange.com/questions/295232/calculating-the-determinant-of-a-matrix-using-a-purely-analytical-method-that-in


template<class T, size_t N>
class Vector_nD
{
public:
	array<T, N> components;

	T operator[](size_t index) const
	{
		return components[index];
	}

	Vector_nD<T, N>(const array<T, N>& comps) : components(comps)
	{

	}

	Vector_nD<T, N>(void)
	{
		components.fill(0.0);
	}

	// Levi-Civita symbol, where i != j, 
	// and so it returns either 1 or -1
	static signed char permutation_sign(const array<int, (N - 1)> perm)
	{
		bool sign = true;

		for (size_t i = 0; i < (N - 2); i++)
		{
			for (size_t j = i + 1; j < (N - 1); j++)
			{
				if (perm[i] > perm[j])
				{
					// a swap would be needed here if sorting
					sign = !sign;
				}
			}
		}

		if (sign)
			return 1;
		else
			return -1;
	}

	// Hodge star operator, where k = n - 1
	static Vector_nD star(const vector<Vector_nD<T, N>>& vectors)
	{
		if (vectors.size() != (N - 1))
		{
			cout << "nD operation requires (n - 1) input vectors" << endl;
			return Vector_nD<T, N>();
		}

		array<T, N> result;

		for (size_t i = 0; i < N; i++)
			result[i] = 0.0;

		// These are the indices we'll use for each component calculation
		array<int, (N - 1)> base_indices;

		for (int i = 0; i < (N - 1); i++)
			base_indices[i] = i;

		// Fill cache
		set< array<int, (N - 1)>> cache_set_positive;

		do
		{
			signed char sign = Vector_nD<T, N>::permutation_sign(base_indices);
			
			if (sign == 1)
				cache_set_positive.insert(base_indices);

		} while (next_permutation(
			base_indices.begin(),
			base_indices.end()));


		// Skip k in our calculations - this is equivalent to removing the k-th column
		// For each permutation of the remaining (N - 1) indices
		for (int k = 0; k < N; k++)
		{
			do
			{
				// Do not use cache
				//signed char sign = Vector_nD<T, N>::permutation_sign(base_indices);

				// Use cache
				signed char sign = 0;

				auto ci = cache_set_positive.find(base_indices);

				if (ci != cache_set_positive.end())
					sign = 1;
				else
					sign = -1;

				// Calculate the product for this permutation
				T product = 1.0;
				ostringstream product_oss;

				for (int i = 0; i < (N - 1); i++)
				{
					const int col = base_indices[i];

					// Adjust column index if it's past k
					int actual_col = 0;

					if (col < k)
						actual_col = col;
					else
						actual_col = col + 1;

					product_oss << "v_{" << i << actual_col << "} ";

					product *= vectors[i][actual_col];
				}

				//if (sign == 1)
				//	cout << "x_{" << k << "} += " << product_oss.str() << endl;
				//else
				//	cout << "x_{" << k << "} -= " << product_oss.str() << endl;

				result[k] += sign * product;

			} while(next_permutation(
					base_indices.begin(), 
					base_indices.end()));
		}

		// Flip handedness
		for (size_t k = 0; k < N; k++)
			if (k % 2 == 1)
				result[k] = -result[k];

		cout << endl;

		for (int k = 0; k < N; k++)
			cout << "result[" << k << "] = " << result[k] << endl;

		cout << endl;

		if (N == 3)
		{
			// Demonstrate the traditional cross product too
			double x = vectors[0][0];
			double y = vectors[0][1];
			double z = vectors[0][2];

			double rhs_x = vectors[1][0];
			double rhs_y = vectors[1][1];
			double rhs_z = vectors[1][2];

			double cross_x = y * rhs_z - rhs_y * z;
			double cross_y = z * rhs_x - rhs_z * x;
			double cross_z = x * rhs_y - rhs_x * y;

			cout << cross_x << " " << cross_y << " " << cross_z << endl << endl;
		}

		return Vector_nD(result);
	}

	static T dot_product(const Vector_nD<T, N>& a, const Vector_nD<T, N>& b)
	{
		return inner_product(
			a.components.begin(), 
			a.components.end(), 
			b.components.begin(), 0.0);
	}
};


template <class T, typename size_t N>
T determinant_nxn(const MatrixX<T>& m)
{
	if (m.cols() != m.rows())
	{
		cout << "Matrix must be square" << endl;
		return 0;
	}

	// We will use this N-vector later, in the dot product operation
	Vector_nD<T, N> a_vector;

	for (size_t i = 0; i < N; i++)
		a_vector.components[i] = m(0, i);

	// We will use these (N - 1) N-vectors later, in the Hodge star operation
	vector<Vector_nD<T, N>> input_vectors;

	for (size_t i = 1; i < N; i++)
	{
		Vector_nD<T, N> b_vector;

		for (size_t j = 0; j < N; j++)
			b_vector.components[j] = m(i, j);

		input_vectors.push_back(b_vector);
	}

	// Compute the Hodge star operation using (N - 1) N-vectors
	Vector_nD<T, N> result = Vector_nD<T, N>::star(input_vectors);

	// Compute the dot product
	T det = Vector_nD<T, N>::dot_product(a_vector, result);

	// These numbers should match
	cout << "Determinant:       " << det << endl;
	cout << "Eigen Determinant: " << m.determinant() << endl << endl;

	return det;
}

int main(int argc, char** argv)
{
	srand(static_cast<unsigned int>(time(0)));

	const size_t N = 8; // Anything larger than 12 takes eons to solve for

	MatrixX<double> m(N, N);

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			m(i, j) = rand() / static_cast<double>(RAND_MAX);

			if (rand() % 2 == 0)
				m(i, j) = -m(i, j);
		}
	}

	std::chrono::high_resolution_clock::time_point start_time, end_time;
	start_time = std::chrono::high_resolution_clock::now();

	determinant_nxn<double, N>(m);

	end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float, std::milli> elapsed = end_time - start_time;

	cout << "Duration: " << elapsed.count() / 1000.0f << " seconds";

	return 0;
}
