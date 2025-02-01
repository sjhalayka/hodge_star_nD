#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
#include <vector>
using namespace std;


// https://claude.ai/chat/3caf4077-28b5-497f-b704-1b0c336a104d


template<size_t N>
class Vector_nD
{
public:
	std::array<double, N> components;

	// Helper function to get the sign of permutation
	static int permutationSign(const std::array<int, (N - 1)>& perm)
	{
		signed int sign = 1;

		for (int i = 0; i < (N - 2); i++)
			for (int j = i + 1; j < (N - 1); j++)
				if (perm[i] > perm[j])
					sign = -sign;

		return sign;
	}

	Vector_nD(const std::array<double, N>& comps) : components(comps) 
	{
	}

	Vector_nD(void)
	{
		for (size_t i = 0; i < N; i++)
			components[i] = 0;
	}

	double operator[](size_t index) const {
		if (index >= N) throw std::out_of_range("Index out of bounds");
		return components[index];
	}

	// Hodge star operator?
	static Vector_nD cross_product(const std::vector<Vector_nD<N>>& vectors)
	{
		if (vectors.size() != (N - 1))
			throw std::invalid_argument("nD cross product requires exactly (n - 1) vectors");

		std::array<double, N> result;// = { 0, 0, 0, 0, 0, 0, 0 };

		for (size_t i = 0; i < N; i++)
			result[i] = 0.0;

		// These are the indices we'll use for each component calculation
		std::array<int, (N - 1)> baseIndices;// = { 0, 1, 2, 3, 4, 5 };

		for (int i = 0; i < (N - 1); i++)
			baseIndices[i] = i;

		// For each component of the result vector
		for (int k = 0; k < N; k++)
		{
			// Skip k in our calculations - this is equivalent to removing the k-th column
			// For each permutation of the remaining (N - 1) indices
			do
			{
				// Calculate sign of this term
				const int sign = permutationSign(baseIndices);

				// Calculate the product for this permutation
				double product = 1.0;

				for (int i = 0; i < (N - 1); i++)
				{
					const int col = baseIndices[i];
					// Adjust column index if it's past k
					const int actualCol = (col >= k) ? col + 1 : col;
					product *= vectors[i][actualCol];
				}

				result[k] += sign * product;

			} while (std::next_permutation(baseIndices.begin(), baseIndices.end()));

			// Reset indices for next component
			for (int i = 0; i < (N - 1); i++)
				baseIndices[i] = i;
		}

		return Vector_nD(result);
	}

	static double dot_product(const Vector_nD<N>& a, const Vector_nD<N>& b)
	{
		double dot_prod = 0;

		for (size_t i = 0; i < N; i++)
			dot_prod += a[i] * b[i];

		return dot_prod;
	}
};



template <typename size_t N>
double determinant_nxn(const MatrixXd& m)
{
	if (m.cols() != m.rows())
	{
		cout << "Matrix must be square" << endl;
		return 0;
	}

	// We will use this vector later, in the dot product operation
	Vector_nD<N> a_vector;

	for (size_t i = 0; i < N; i++)
		a_vector.components[i] = m(0, i);

	// We will use these (N - 1) vectors later, in the cross product operation
	std::vector<Vector_nD<N>> input_vectors;

	for (size_t i = 1; i < N; i++)
	{
		Vector_nD<N> b_vector;

		for (size_t j = 0; j < N; j++)
			b_vector.components[j] = m(i, j);

		input_vectors.push_back(b_vector);
	}

	// Compute the cross product using (N - 1) vectors
	Vector_nD<N> result = Vector_nD<N>::cross_product(input_vectors);

	// Flip handedness
	for (size_t i = 0; i < result.components.size(); i++)
		if (i % 2 == 1)
			result.components[i] = -result.components[i];

	double det = Vector_nD<N>::dot_product(a_vector, result);

	// These numbers should match
	cout << det << endl;
	cout << m.determinant() << endl << endl;

	return det;
}



int main(int argc, char** argv)
{
	srand(static_cast<unsigned int>(time(0)));

	const size_t N = 11; // Anything larger than 11 takes eons to solve for

	MatrixXd m(N, N);

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			m(i, j) = rand() / static_cast<double>(RAND_MAX);

			if (rand() % 2 == 0)
				m(i, j) = -m(i, j);
		}
	}

	determinant_nxn<N>(m);

	return 0;
}