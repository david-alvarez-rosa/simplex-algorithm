#include <iostream>
#include <vector>

using namespace std;


using V = vector<double>;
using M = vector<V>;
using VI = vector<int>;


// Devuelve la posición del mínimo de un vector.
int minVectorPos(V v) {
   int n = v.size();

   double min = v[0];
   int posMin = 0;
   for (int i = 1; i < n; ++i)
      if (v[i] < min) {
	 min = v[i];
	 posMin = i;
      }

   return posMin;
}


// Devuelve la posición del mínimo de un vector.
double minVectorValue(V v) {
   int n = v.size();

   double min = v[0];
   for (int i = 1; i < n; ++i)
      if (v[i] < min)
	 min = v[i];

   return min;
}


// Devuelve la matriz identidad (nxn).
M matrixIdentity(int n) {
   M I(n, V(n, 0));

   for (int k = 0; k < n; ++k)
      I[k][k] = 1;

   return I;
}


// Multiplica matrices (deben tener tamaños adecuados).
M matrixMultiplication(M A, M B) {
   int m = A.size();
   int n = A[0].size();
   int p = B[0].size();

   M C(m, V(p, 0));

   for (int i = 0; i < m; ++i)
      for (int j = 0; j < p; ++j)
	 for (int k = 0; k < n; ++k)
	    C[i][j] += A[i][k]*B[k][j];

   return C;
}


// Multiplica vector por matriz.
V vectorMatrixMultiplication(V v, M A) {
   int n = v.size();
   int p = A[0].size();
   
   V c(p, 0);

   for (int j = 0; j < p; ++j)
      for (int k = 0; k < n; ++k)
	 c[j] += v[k]*A[k][j];

   return c;
}


// Multiplica matriz por vector.
V matrixVectorMultiplication(M A, V v) {
   int m = A.size();
   int n = A[0].size();
   
   V c(m, 0);

   for (int i = 0; i < m; ++i)
      for (int k = 0; k < n; ++k)
	 c[i] += A[i][k]*v[k];

   return c;
}


// Multiplica vector por escalar.
V scalarVectorMultiplication(double M, V v) {
   int n = v.size();
   for (int i = 0; i < n; ++i) {
      v[i] *= M;
   }
   
   return v;
}


// Suma dos vectores.
V vectorSum(V v, V w) {
   int n = v.size();
   V t(n);

   for (int i = 0; i < n; ++i)
      t[i] = v[i] + w[i];

   return t;
}


// Resta dos vectores.
V vectorDifference(V v, V w) {
   int n = v.size();
   V t(n);

   for (int i = 0; i < n; ++i)
      t[i] = v[i] - w[i];

   return t;
}


// Suma M veces la fila b-ésima a la fila a-ésima.
void combineMatrixRows(M &A, int a, int b, double M) {
   int m = A[0].size();

   for (int j = 0; j < m; ++j)
      A[a][j] += M*A[b][j];
}


// Multiplica la fila i de la matriz por el escalar M.
void multiplyMatrixRow(M &A, int i, double M) {
   int n = A[0].size();
   for (int j = 0; j < n; ++j)
      A[i][j] *= M;
}


// Dados dos vectores, se queda con los elementos del primero indicados por el segundo.
V vectorChooseElements(V c, VI v) {
   int n = v.size();

   V cv(n);

   for (int i = 0; i < n; ++i)
      cv[i] = c[v[i] - 1];

   return cv;
}


// Dada matriz y vector se queda con las columnas indicadas por el vector.
M matrixChooseColumns(M A, VI v) {
   int m = A.size();
   int n = v.size();

   M C(m, V(n));
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
	 C[i][j] = A[i][v[j] - 1];

   return C;
}


// Devuelve la columna j-ésima.
V matrixChooseColumn(M A, int j) {
   int m = A.size();
   V v(m);

   for (int i = 0; i < m; ++i)
      v[i] = A[i][j - 1];

   return v;
}


// Imprime por pantalla un vector.
void printVector(V v) {
   int n = v.size();

   for (int i = 0; i < n; ++i)
	 cout << v[i] << '\t';
   cout << endl;
}


// Imprime por pantalla un vector.
void printVector(VI v) {
   int n = v.size();

   for (int i = 0; i < n; ++i)
	 cout << v[i] << '\t';
   cout << endl;
}
