#include "fraction.h"
#include <string>
#include <iostream>
#include <cassert>
using namespace std;

/***
 **  A simpleau implementation of the Simplex method algorithm for Standard LPs.
**  --- Made by Oh
**/

void printArrays(int n, int m, Fraction* b, Fraction** A, Fraction* c) {
    for (int i=0; i<n; i++)
        cout << c[i] << " ";
    cout << endl << endl;
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    for (int i=0; i<m; i++)
        cout << b[i] << " ";
    cout << endl;
}

void printArray(int* arr, int size) {
    for (int i=0; i<size; ++i) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

void printArray(Fraction* arr, int size) {
    for (int i=0; i<size; ++i) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

void printArray(int n, int m, Fraction** BA) {
    cout << "___" << endl;
    for (int y=0; y<m; y++) {
        for (int x=0; x<n; x++)
            cout << BA[y][x] << " ";
        cout << endl;
    }
    cout << "^^^" << endl;
}

void deleteArray(int size, Fraction* array) {
    for (int i=size-1; i>=0; --i) {
        array[i].~Fraction();
    }
    delete(array);
}

void deleteArray(int n, int m, Fraction** array) {
    for (int y=m-1; y>=0; --y) {
        for (int x=n-1; x>=0; --x) {
            array[y][x].~Fraction();
        }
        delete(array[y]);
    }
    delete(array);
}

Fraction* createFractionArray(int size) {
    return static_cast<Fraction*> (::operator new (sizeof(Fraction)* size));
}

void multiplyRow(int n, int m, Fraction** BA, int rw, Fraction f) {
    for (int i=0; i<n; i++) {
        BA[rw][i] *= f;
    }
}

void multiplyCol(int n, int m, Fraction** BA, int cl, Fraction f) {
    for (int i=0; i<m; i++) {
        BA[i][cl] *= f;
    }
}


// Returns the row index of the only positive (& non-zero) item in the column.
// If there is more than one non-zero item in the column, returns -1.
int findOnlyPositiveItemInColumn(int n, int m, Fraction** BA, int cl) {
    int nonZeroIndex = -1;
    for (int i=0; i<m; i++) {
        if (BA[i][cl].isPositive()) {
            if (nonZeroIndex == -1) {
                nonZeroIndex = i;
            } else {
                // More than one non-zero item.
                return -1;
            }
        }
    }
    return nonZeroIndex;
}

Fraction** extendTableau(int n, int m, int n_new, Fraction** BA) {
    // BA_2 is an extension of BA.
    Fraction** BA_2 = new Fraction*[m];
    for (int i=0; i<m; ++i) {
        BA_2[i] = createFractionArray(n_new);
        for (int j=0; j<n; ++j) {
            BA_2[i][j] = BA[i][j];
        }
    }
    return BA_2;
}

Fraction simplexMain(int n_curr, int n_orig, int m, Fraction** BA, Fraction* cBar,
        Fraction* basicSolution, int* basicVariables, Fraction neg_solution) {
    Fraction NEG_ONE = Fraction(-1);
    Fraction ONE = Fraction(1);

    while(true) {
        // >>> STEP 5A - Find a nonbasic variable with negative reduced cost.
        int selected_col = -1;
        int selected_row = -1;
        for (int x=0; x<n_curr; ++x) {
            if (cBar[x].isNegative()) {
                Fraction theta_bar = NEG_ONE;

                // >>> STEP 5B - Compute Theta_bar. Try to find a positive one.
                for (int y=0; y<m; ++y) {
                    if (BA[y][x].isPositive()) {
                        Fraction theta_y = basicSolution[y] / BA[y][x];
                        if (theta_bar.isNegative() || theta_y < theta_bar) {
                            theta_bar = theta_y;
                            selected_row = y;
                        }
                    }
                }
                if (theta_bar.isPositive()) {
                    selected_col = x;
                    x = n_curr; // break;
                } else if (theta_bar.isNegative()) {
                    // Solution is unbounded.
                    return Fraction(); // Invalid fraction -> unbounded
                }
            }
        }
        //cout << "SELECT " << selected_col << " ROW " << selected_row << endl;
        
        // selected_col is -1 when no column has been selected. Terminate...?
        // No variable of negative reduced cost.
        if (selected_col == -1) {

            // HOLD IT! DON'T BREAK YET! MAYBE THERE'S STILL SOME ZERO BASIC ARTIFICIAL VARIABLES!
            // original_n == -1 means it's the original LP, not the auxiliary LP.
            if (n_orig != -1 && neg_solution.isZero()) {
                selected_row = -1;
                for (int y=0; y<m; ++y) {
                    if (basicVariables[y] >= n_orig) {
                        selected_row = y;
                    }
                }
                if (selected_row != -1) {
                    // selected_row refers to a Zero Basic Aritificial Variable.
                    bool* isBasic = new bool[n_orig];
                    for (int x=0; x<n_orig; ++x) {
                        isBasic[x] = false;
                    }
                    for (int y=0; y<m; ++y) {
                        if (basicVariables[y] < n_orig) {
                            isBasic[basicVariables[y]] = true;
                        }
                    }
                    for (int x=0; x<n_orig; ++x) {
                        if (!isBasic[x]) {
                            selected_col = x;
                            x = n_orig; // break out of for loop.
                        }
                    }
                    // Don't break here! the simplex must continue!
                } else {
                    break;
                }
            } else {
                break; // End simplex method.
            }
        }

        // >>> STEP 5C - move in the selected direction (enter that variable)
        // We know that BA[selected_row][selected_col] is positive.
        // Replace the basic variable.
        basicVariables[selected_row] = selected_col;

        // Normalise the row selected_row.
        Fraction multiply = ONE / BA[selected_row][selected_col];
        multiplyRow(n_curr,m,BA, selected_row, multiply);
        basicSolution[selected_row] *= multiply;

        // zero out all the other rows in column selected_column.
        for (int y=0; y<m; ++y) {
            if (y != selected_row) {
                multiply = BA[y][selected_col];
                for (int x=0; x<n_curr; ++x) {
                    BA[y][x] -= BA[selected_row][x]*multiply;
                }
                basicSolution[y] -= basicSolution[selected_row]*multiply;
            }
        }
        // zero out the reduced costs row too.
        multiply = cBar[selected_col];
        for (int x=0; x<n_curr; ++x) {
            cBar[x] -= BA[selected_row][x]*multiply;
        }
        neg_solution -= basicSolution[selected_row]*multiply;


        //printArray(basicVariables, m);
        //printArray(basicSolution, m);
        //printArray(n_curr, m, BA);
        //printArray(cBar, n_curr);
        //cout << "----------------- " << neg_solution << endl;
    }
    return neg_solution;
}




// BA: B^-1*A
Fraction* runSimplex(int n, int m, Fraction* b, Fraction** BA, Fraction* c) {
    Fraction neg_solution = Fraction(0);
    int* basicVariables = new int[m]; // stores indices of basic variables.
    Fraction* basicSolution = createFractionArray(m);

    const Fraction ZERO = Fraction(0);
    const Fraction ONE = Fraction(1);
    const Fraction NEG_ONE = Fraction(-1);

    for (int i=0; i<m; i++) {
        basicVariables[i] = -1;
    }
    /**
    * ----- PART ZERO --------------------------------- +++++++++++++++++++++++
    * ------------------------- MAKING b POSITIVE ----- +++++++++++++++++++++++
    **/
    for (int i=0; i<m; ++i) {
        if (b[i].isNegative()) {
            for (int j=0; j<n; ++j) {
                BA[i][j] *= NEG_ONE;
            }
            b[i] *= NEG_ONE;
        }
    }



    /**
    * ----- PART ONE ---------------------------------- +++++++++++++++++++++++
    * --------------- GENERATING THE AUXILIARY LP ----- +++++++++++++++++++++++
    **/

    // >>> STEP 1 - Locate the rows without a corresponding elementary basis vector column.
    // If the column looks like this: 0 5 0 0 0, then divide the second row throughout by 5.
    for (int i=0; i<n; ++i) {
        int row = findOnlyPositiveItemInColumn(n,m,BA,i);
        if (row != -1) {
            if (basicVariables[row] == -1) {
                if (BA[row][i] != ONE) {
                    Fraction multiply = ONE / BA[row][i];
                    multiplyRow(n,m,BA,row,multiply);
                    b[row] *= multiply;
                }
                basicVariables[row] = i;
            }
        }
    }

    // STEP 1 - END

    //printArray(basicVariables, m);



    // >>> STEP 2 - Create auxiliary variables.
    int n_aux = n; // number of variables in auxiliary LP.
    for (int i=0; i<m; i++) {
        if (basicVariables[i] == -1) {
            n_aux++;
        }
    }

    // Resize the array BA.
    BA = extendTableau(n, m, n_aux, BA);
    {
        int col = n;
        for (int i=0; i<m; i++) {
            if (basicVariables[i] == -1) {
                basicVariables[i] = col;
                for (int j=0; j<m; j++) {
                    if (j == i) {
                        BA[j][col] = ONE;
                    } else {
                        BA[j][col] = ZERO;
                    }
                }
                col++;
            }
        }
        assert(col == n_aux);
    }

    // Use new cost vector.
    Fraction* c_aux = createFractionArray(n_aux);
    for (int i=0; i<n; ++i) {
        c_aux[i] = ZERO;
    }
    for (int i=n; i<n_aux; ++i) {
        c_aux[i] = ONE;
    }

    // STEP 2 - END

    //printArray(basicVariables, m);
    //printArrays(n_aux, m, b, BA, c_aux);



    /**
    * ----- PART TWO ---------------------------------- +++++++++++++++++++++++
    * ------------------ SOLVING THE AUXILIARY LP ----- +++++++++++++++++++++++
    **/

    // >>> STEP 3 - Get the first BFS.
    // Note: basic vectors are already chosen.
    // Set the current basicSolution.
    for (int i=0; i<m; i++) {
        basicSolution[i] = b[i];
    }
    // Set the current neg_solution value.
    for (int i=0; i<m; i++) {
        neg_solution -= basicSolution[i] * c_aux[basicVariables[i]];
    }

    // STEP 3 - END

    //cout << "SOLUTION = " << neg_solution << endl;


    // >>> STEP 4 - Compute reduced cost vector for aux LP. (for aux LP it's faster)
    // cBar is zero for basic variables.
    // for nonbasic variables, c[i] = -1 * dot(c_B,A_i)  where A_i is the ith column of A. 
    Fraction* cBar = createFractionArray(n_aux);
    {
        bool* isBasic = new bool[n_aux];
        for (int i=0; i<n_aux; ++i) {
            isBasic[i] = false;
        }
        for (int i=0; i<m; ++i) {
            isBasic[basicVariables[i]] = true;
        }
        for (int i=0; i<n_aux; ++i) {
            cBar[i] = ZERO;
            if (!isBasic[i]) {
                for (int j=0; j<m; j++) {
                    cBar[i] -= c_aux[basicVariables[j]] * BA[j][i];
                }
            }
        }
        delete(isBasic);
    }

    /// STEP 4 - END

    //printArray(cBar, n_aux);

        

    // >>> STEP 5 - Run the Simplex Algorithm!!!
    neg_solution = simplexMain(n_aux, n, m, BA, cBar, basicSolution,
                                basicVariables, neg_solution);
    
    deleteArray(n_aux, cBar);
    deleteArray(n_aux, c_aux);

    if (neg_solution.isNegative()) {
        // FAIL! Cleanup and return.
        deleteArray(m, basicSolution);
        delete(basicVariables);

        Fraction* fail = createFractionArray(1);
        fail[0] = NEG_ONE;
        return fail;
    }

    // STEP 5 - END


    //printArray(basicVariables, m);
    //printArrays(n_aux, m, b, BA, c_aux);
    //printArray(cBar, n_aux);

    

    /**
    * ----- PART THREE -------------------------------- +++++++++++++++++++++++
    * ------------------- SOLVING THE ORIGINAL LP ----- +++++++++++++++++++++++
    **/

    // >>> STEP 6 - Recompute the reduced costs.
    cBar = createFractionArray(n);
    {
        bool* isBasic = new bool[n];
        for (int x=0; x<n; ++x) {
            isBasic[x] = false;
            cBar[x] = c[x];
        }
        for (int y=0; y<m; ++y) {
            isBasic[basicVariables[y]] = true;
            cBar[basicVariables[y]] = ZERO;
        }
        for (int x=0; x<n; ++x) {
            if (!isBasic[x]) {
                for (int y=0; y<m; y++) {
                    cBar[x] -= c[basicVariables[y]] * BA[y][x];
                }
            }
        }
        delete(isBasic);
    }
    // Recompute the solution.
    neg_solution = ZERO;
    for (int i=0; i<m; i++) {
        neg_solution -= basicSolution[i] * c[basicVariables[i]];
    }

    // STEP 6 - END



    // >>> STEP 7 - Run the Simplex Algorithm on the original LP.
    neg_solution = simplexMain(n, -1, m, BA, cBar, basicSolution,
                                basicVariables, neg_solution);

    // Check for unbounded.
    if (neg_solution.isInvalid()) {
        // Unbounded Solution
        // Cleanup
        deleteArray(n, cBar);
        deleteArray(m, basicSolution);
        delete(basicVariables);

        Fraction* unbounded = createFractionArray(1);
        unbounded[0] = Fraction();
        return unbounded;
    }

    // STEP 7 - END
    

    //cout << "SOLUTION = " << neg_solution << endl;
    //printArray(basicVariables, m);
    //printArrays(n, m, b, BA, c);
    //printArray(cBar, n);
    //printArray(basicSolution, m);



    /**
    * ----- PART THREE -------------------------------- +++++++++++++++++++++++
    * ------------------- FORMATTING THE SOLUTION ----- +++++++++++++++++++++++
    **/

    // STEP 8 - FORMAT SOLUTION
    Fraction* optimalSolution = createFractionArray(n);
    for (int x=0; x<n; ++x) {
        optimalSolution[x] = ZERO;
    }
    for (int y=0; y<m; ++y) {
        optimalSolution[basicVariables[y]] = basicSolution[y];
    }

    // STEP 8 - END


    // Cleanup
    deleteArray(n, cBar);
    deleteArray(m, basicSolution);
    delete(basicVariables);

    return optimalSolution;
}

Fraction scanFraction() {
    string s;
    cin >> s;
    return parseFraction(s);
}


void interpretAndPrint(int n, Fraction* c, Fraction* optimalSolution) {

    if (optimalSolution[0].isInvalid()) {
        cout << "UNBOUNDED" << endl;
        cout << "Optimal Objective Value = -infinity";

    } else if (optimalSolution[0].isNegative()) {
        cout << "INFEASIBLE" << endl;

    } else {
        cout << "Optimal Solution: ";
        printArray(optimalSolution, n);
        Fraction value = Fraction(0);
        for (int x=0; x<n; ++x) {
            value += c[x]*optimalSolution[x];
        }
        cout << "Optimal Objective Value = " << value;

    }
}

int main(int argc, char** argv) {
    /**
    * Conditions for Standard LP:
    * min cTx
    * Ax = b
    * x >= 0
    *
    * Note: no condition that b >= 0. We need to convert it ourselves.
    */

    // Parsing - START
    int n, m;
    cin >> n;
    cin >> m;

    Fraction* b = createFractionArray(m);
    Fraction** A = new Fraction*[m];
    for (int i=0; i<m; i++) {
        A[i] = createFractionArray(n);
    }
    Fraction* c = createFractionArray(n);

    for (int i=0; i<n; ++i) {
        c[i] = scanFraction();
    }
    for (int i=0; i<m; ++i) {
        for (int j=0; j<n; ++j) {
            A[i][j] = scanFraction();
        }
    }
    for (int i=0; i<m; ++i) {
        b[i] = scanFraction();
    }
    // Parsing - END

    // Solving - START
    Fraction* optimalSolution = runSimplex(n, m, b, A, c);
    interpretAndPrint(n, c, optimalSolution);
    // Solving - END

    deleteArray(n, c);
    deleteArray(n, optimalSolution);
    deleteArray(n, m, A);
    deleteArray(m, b);
}