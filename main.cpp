#include <array>
#include <iostream>
#include <limits>

//#define DEBUG

template <typename T, size_t n, size_t m>
using Matrix = std::array<std::array<T, m>, n>;

template <typename T, size_t n>
using Vector = std::array<T, n>;

template <typename T, size_t n, size_t m>
void print_matrix(Matrix<T, n, m> matrix) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <typename T, size_t n>
void print_matrix(Vector<T, n> matrix) {
    for (size_t i = 0; i < n; ++i) {
        std::cout << matrix[i] << std::endl;
    }
    std::cout << std::endl;
}

// n is the number of xs, m is the number of constraints
template <size_t n, size_t m>
struct LinearProgram {
    // Objective function coefficients x_1 ... x_n ... x_m
    // x_n+1 ... x_m are the new slack variables
    Vector<float, n + m + 1> c;
    // A[i][j] is the coefficient of x_j in constraint i. A[i][m + n] is the constant
    Matrix<float, m, n + m + 1> A;

    Vector<size_t, m> loose;
    Vector<size_t, n> tight;

    LinearProgram() = delete;
    LinearProgram(Vector<float, n> c, Matrix<float, m, n> A, Vector<float, m> b) {
        // c
        for (size_t i = 0; i < n; ++i) {
            this->c[i] = c[i];
        }
        for (size_t i = n; i < n + m + 1; ++i) {
            this->c[i] = 0;
        }

        // A
        // actual coeff values
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                this->A[i][j] = -A[i][j];
            }
            
            for (size_t j = n; j < n + m + 1; ++j) {
                this->A[i][j] = 0;
            }

        }
        // constant values
        for (size_t i = 0; i < m; ++i) {
            this->A[i][n + m] = b[i];
        }

        // Tight and loose
        for (size_t i = 0; i < n; ++i) {
            this->tight[i] = i;
        }

        for (size_t i = n; i < n + m; ++i) {
            this->loose[i - n] = i;
        }
    }


    void solve() {

#ifdef DEBUG
        print_matrix(c);
        print_matrix(loose);
        print_matrix(tight);
        print_matrix(A);
#endif
        
        while (true) {
            // 1. Choose a variable to loosen (Dantzig's rule chooses the highest non-zero coefficient)
            size_t loosen = 0;
            float highest_coeff = -1;

            for (size_t i = 0; i < n; ++i) {
                size_t var = tight[i];
                float coeff = c[var];

                if (coeff > highest_coeff) {
                    loosen = var;
                    highest_coeff = coeff;
                }
            }

            if (highest_coeff <= 0) {
                return;
            }

            // 2. Tighten the variable with the largest non-positive ratio of constant and coefficient
            size_t tighten = 0;
            size_t constraint = 0;
            float highest_ratio = std::numeric_limits<float>::lowest(); // highest non-positive ratio

            for (size_t i = 0; i < m; ++i) {
                size_t var = loose[i];
                float constant = A[i][n + m];
                float coeff = A[i][loosen];
                if (coeff == 0) continue;
                float ratio = constant / coeff;

                if (ratio <= 0 && ratio > highest_ratio) {
                    tighten = var;
                    constraint = i;
                    highest_ratio = ratio;
                }
            }

            if (highest_ratio == std::numeric_limits<float>::lowest()) {
                std::cerr << "Something went wrong" << std::endl;
                return;
            }

        
            // 3. Fix the inequality so that all the basic (loose) variables are on the right

            float loosen_coeff = A[constraint][loosen];
#ifdef DEBUG
            std::cout << "Loosen: " << loosen << std::endl;
            std::cout << "Constraint: " << constraint << std::endl;
            std::cout << "LoosenCoeff: " << loosen_coeff << std::endl;
            std::cout << "Tighten: " << tighten << std::endl;
#endif
        
            // Move `loosen` to the left hand side
            loose[constraint] = loosen;
            A[constraint][loosen] = 0;

            // Move `tighten` to the right hand side
            A[constraint][tighten] = -1;

            // Divide by `loosen_coeff`
            for (size_t i = 0; i < n + m + 1; ++i) {
                A[constraint][i] = A[constraint][i] / -loosen_coeff;
            }

            // Remove `loosen` from the right side of everything else
            for (size_t i = 0; i < m; ++i) {
                float loosen_coeff = A[i][loosen];
                for (size_t j = 0; j < n + m + 1; ++j) {
                    A[i][j] += loosen_coeff * A[constraint][j];
                }
                A[i][loosen] = 0;
            }

            // Remove `loosen` from the objective function
            loosen_coeff = c[loosen];
            for (size_t i = 0; i < n + m + 1; ++i) {
                c[i] += loosen_coeff * A[constraint][i];
            }
            c[loosen] = 0;

#ifdef DEBUG
            print_matrix(c);
            print_matrix(loose);
            print_matrix(tight);
            print_matrix(A);
#endif
        }
    }
    
    float objective() const {
        return c[n + m];
    }
};

void example_1() {
    Matrix<float, 3, 2> A{{
            {1, 0},
            {0, 1},
            {1, 1},
        }};
    Vector<float, 3> b{{3000, 4000, 5000}};
    Vector<float, 2> c{{1.2, 1.7}};

    LinearProgram<2, 3> lp{c, A, b};
    lp.solve();
    std::cout << lp.objective() << std::endl;
}

void example_2() {
    Matrix<float, 3, 3> A{
        {
            {2, 3, 1},
            {4, 1, 2},
            {3, 4, 2},
        }
    };
    Vector<float, 3> b{{5, 11, 8}};
    Vector<float, 3> c{{5, 4, 3}};

    LinearProgram<3, 3> lp{c, A, b};
    lp.solve();
    std::cout << lp.objective() << std::endl;
};


int main() {
    example_1();
    example_2();
    return 0;
}
