using namespace std;

using V = vector<double>;
using M = vector<V>;
using VI = vector<int>;


// Encontrar el mínimo de un vector (posición o valor).
int minVectorPos(V v);
double minVectorValue(V v);

// Matriz identidad de tamaño nxn.
M matrixIdentity(int n);

// Operaciones básicas con matrices, vectores y escalares.
M matrixMultiplication(M A, M B);
V vectorMatrixMultiplication(V v, M A);
V matrixVectorMultiplication(M A, V v);
V scalarVectorMultiplication(double M, V v);
V vectorSum(V v, V w);
V vectorDifference(V v, V w);

// Suma M veces la fila b-ésima a la fila a-ésima.
void combineMatrixRows(M &A, int a, int b, double M);
// Multiplica la fila i de la matriz por el escalar M.
void multiplyMatrixRow(M &A, int i, double M);

// Escoger diversos elementos de un vector o matriz.
V vectorChooseElements(V c, VI v);
M matrixChooseColumns(M A, VI v);
V matrixChooseColumn(M A, int j);

// Imprimir vector por pantalla.
void printVector(V v);
void printVector(VI v);
