#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
#include <vector>
#include <numeric>
#include <mutex>
#include <thread>

using namespace std;


// https://claude.ai/chat/3caf4077-28b5-497f-b704-1b0c336a104d
// https://codereview.stackexchange.com/questions/295232/calculating-the-determinant-of-a-matrix-using-a-purely-analytical-method-that-in












template<class T, size_t N>
class Vector_nD
{
public:
	std::array<T, N> components;

	static void thread_func(const size_t k, double &result, const vector<Vector_nD<T, N>>& vectors)
	{
		result = 0;

		// These are the indices we'll use for each component calculation
		std::array<int, (N - 1)> base_indices;// = { 0, 1, 2, 3, 4, 5 };

		for (int i = 0; i < (N - 1); i++)
			base_indices[i] = i;

		// Skip k in our calculations - this is equivalent to removing the k-th column
		// For each permutation of the remaining (N - 1) indices
		do
		{
			// Calculate sign of this term
			const signed char sign = permutation_sign(base_indices);

			// Calculate the product for this permutation
			T product = 1.0;

			for (int i = 0; i < (N - 1); i++)
			{
				const int col = base_indices[i];

				// Adjust column index if it's past k
				int actual_col = 0;

				if (col < k)
					actual_col = col;
				else
					actual_col = col + 1;

				product *= vectors[i][actual_col];
			}

			result += sign * product;
		} 
		while (std::next_permutation(base_indices.begin(), base_indices.end()));
	}

	// Helper function to get the sign of permutation
	static signed char permutation_sign(const std::array<int, (N - 1)>& perm)
	{
		bool sign = true;

		for (int i = 0; i < (N - 2); i++)
			for (int j = i + 1; j < (N - 1); j++)
				if (perm[i] > perm[j])
					sign = !sign;

		if (sign)
			return 1;
		else
			return -1;
	}

	Vector_nD(const std::array<T, N>& comps) : components(comps) 
	{
	}

	Vector_nD(void)
	{
		components.fill(0.0);
	}

	T operator[](size_t index) const 
	{
		return components[index];
	}

	// Hodge star operator?
	static Vector_nD cross_product(const std::vector<Vector_nD<T, N>>& vectors)
	{
		if (vectors.size() != (N - 1))
		{
			cout << "nD cross product requires exactly (n - 1) input vectors" << endl;
			return Vector_nD<T, N>();
		}

		std::array<T, N> result;

		for (size_t i = 0; i < N; i++)
			result[i] = 0.0;

		vector<thread> threads;

		// For each component of the result vector
		for (int k = 0; k < N; k++)
			threads.push_back(thread(thread_func, k, ref(result[k]), ref(vectors)));

		for (size_t i = 0; i < N; i++)
			threads[i].join();

		return Vector_nD(result);
	}

	static T dot_product(const Vector_nD<T, N>& a, const Vector_nD<T, N>& b)
	{
		return inner_product(a.components.begin(), a.components.end(), b.components.begin(), 0.0);
	}
};


template <class T, typename size_t N>
T determinant_nxn(const MatrixXd& m)
{
	if (m.cols() != m.rows())
	{
		cout << "Matrix must be square" << endl;
		return 0;
	}

	// We will use this vector later, in the dot product operation
	Vector_nD<T, N> a_vector;

	for (size_t i = 0; i < N; i++)
		a_vector.components[i] = m(0, i);

	// We will use these (N - 1) vectors later, in the cross product operation
	std::vector<Vector_nD<T, N>> input_vectors;

	for (size_t i = 1; i < N; i++)
	{
		Vector_nD<T, N> b_vector;

		for (size_t j = 0; j < N; j++)
			b_vector.components[j] = m(i, j);

		input_vectors.push_back(b_vector);
	}

	// Compute the cross product using (N - 1) vectors
	Vector_nD<T, N> result = Vector_nD<T, N>::cross_product(input_vectors);

	// Flip handedness
	for (size_t i = 0; i < result.components.size(); i++)
		if (i % 2 == 1)
			result.components[i] = -result.components[i];

	T det = Vector_nD<T, N>::dot_product(a_vector, result);

	// These numbers should match
	cout << det << endl;
	cout << m.determinant() << endl << endl;

	return det;
}


int main(int argc, char** argv)
{
	srand(static_cast<unsigned int>(time(0)));

	const size_t N = 11; // Anything larger than 11 takes eons to solve for

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

	determinant_nxn<double, N>(m);

	return 0;
}
