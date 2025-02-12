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

	// Compute the Levi-Civita symbol using a pattern
	static signed char permutation_sign_pattern(const size_t term_index, const bool parity)
	{
		signed char sign = 0;

		if (term_index == 0)
		{
			sign = 1;
		}
		else
		{
			size_t local_term_index = term_index - 1;
			local_term_index /= 2;

			if (local_term_index % 2 != 0)
				sign = 1;
			else
				sign = -1;
		}

		if (false == parity)
			sign = -sign;

		return sign;
	}

	// Compute the Levi-Civita symbol using the determinant method
	static signed char permutation_sign_det(const array<int, (N - 1) >& indices)
	{
		size_t n = indices.size();

		MatrixX<double> m(n, n);
		m = MatrixX<double>::Zero(n, n);

		for (int i = 0; i < n; i++)
			m(i, indices[i]) = 1;

		// Compute the determinant of the Kronecker delta matrix
		return static_cast<signed char>(m.determinant());
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



		// Skip k in our calculations - this is equivalent to removing the k-th column
		// For each permutation of the remaining (N - 1) indices
		for (int k = 0; k < N; k++)
		{
			size_t term_index = 0;
			bool parity = false;

			if (N < 6)
				parity = true;

			//long signed int prev_coeff_index = -1;

			string prev_string = "";

			size_t swap_count = 0;
			size_t no_swap_count = 0;

			vector<string> prev_tokens;

			do
			{
				// Use pattern
				signed char sign = permutation_sign_pattern(term_index, parity);
				const signed char sign_perm = Vector_nD<T, N>::permutation_sign(base_indices);

				// Calculate the product for this permutation
				T product = 1.0;
				ostringstream product_oss;

				vector<string> tokens;

				for (int i = 0; i < (N - 1); i++)
				{
					const int col = base_indices[i];

					// Adjust column index if it's past k
					int actual_col = 0;

					if (col < k)
						actual_col = col;
					else
						actual_col = col + 1;

					if (i < N - 6)
						tokens.push_back(to_string(actual_col));

					product_oss << "v_{" << i << actual_col << "} ";

					product *= vectors[i][actual_col];

					if (N >= 6 && i == N - 6)
					{
						// Cheat if N == 10 or greater
						if (N >= 10)
						{
							if (sign != sign_perm)
							{
								parity = !parity;
								sign = -sign;
								term_index = 0;
							}
						}
						else
						{
							string str = to_string(actual_col);// product_oss.str();

							//							if( tokens.end() == find( tokens.begin(),  tokens.end(), str))
							if (prev_string != str)
							{
								// Calculate manually

								//cout << (int)sign << " " << (int)sign_perm << endl;
								//cout << prev_string << "     " << str << endl;

								//cout << "token sizes: " << tokens.size() << " " << prev_tokens.size() << endl;

								size_t different_count = 0;

								for (size_t j = 0; j < prev_tokens.size(); j++)
								{
									if (prev_tokens[j] != tokens[j])
										different_count++;
								}

								long signed int x = 1;

								if (different_count <= x)
								{
									parity = !parity;
									sign = -sign;
									term_index = 0;
								}

								prev_string = str;

							}
						}
					}
				}

				prev_tokens = tokens;

				term_index++;

				result[k] += sign * product;

				//if (sign == 1)
				//	cout << "x_{" << k << "} += " << product_oss.str() << endl;
				//else
				//	cout << "x_{" << k << "} -= " << product_oss.str() << endl;

			} while (next_permutation(
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

	const size_t N = 9;

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
