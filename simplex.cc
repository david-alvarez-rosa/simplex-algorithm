#include <iostream>
#include <vector>
#include "functions.h"

using namespace std;


using V = vector<double>;
using M = vector<V>;
using VI = vector<int>;


const double TOLERANCE = 1e-8;


// Estructura de datos de número real, matriz, matriz.
struct dmm {
   double z;
   M B;
   M invB;
};


// Imprimir por pantalla la solución del Símplex (el óptimo si existiera).
void printSolution(VI vb, V xb, VI vn, V r, double z) {
   cout.setf(ios::fixed);
   cout.precision(4);
   
   cout << "VB*=" << endl << '\t';
   printVector(vb);
   cout << "xb*=" << endl << '\t';
   printVector(xb);
   cout << "VNB*=" << endl << '\t';
   printVector(vn);
   cout << "r*=" << endl << '\t';
   printVector(r);
   cout << "z*=" << endl;
   cout << '\t' << z << endl;
}


// Imprimir por pantalla información de la iteracción del símplex.
void printDataSimplexIter(int q, double rq, int Bp, double theta, double z) {
   cout << "q = " << q << ",\t";
   cout << "rq = " << rq << ",\t\t";
   cout << "B(p) = " << Bp << ",\t";
   cout << "theta* = " << theta << ", \t";
   cout << "z = " << z << endl;
}


// Calcula el coste reducido.
V reducedCost(V c, M A, M invB, VI vb, VI vn) {
   V cn = vectorChooseElements(c, vn);
   V cb = vectorChooseElements(c, vb);
   M An = matrixChooseColumns(A, vn);

   return
      vectorDifference(cn, vectorMatrixMultiplication(cb, matrixMultiplication(invB, An)));
}


// Escoge variable no básica para entrar a la base.
// Devuelve la posición de la no básica en el vector de no básicas.
int chooseNonBasicToEnter(VI vn, V r, int rule) {
   // Regla del coste más reducido.   
   if (rule == 1)
      return minVectorPos(r);

   // Regla de Bland.
   int qPos = -1;
   for (int i = 0; i < (int) r.size(); ++i)
      if (r[i] < 0 and (qPos == -1 or vn[i] < vn[qPos]))
         qPos = i;
   return qPos;
}


// Devuelve dirección básica de descenso asociada a la variable q.
V descentDirection(M invB, M A, int q) {
   V Aq = matrixChooseColumn(A, q);
   return
      scalarVectorMultiplication(-1, matrixVectorMultiplication(invB, Aq));
}


// Devuelve la variable básica para salir.
// Escojo según la regla de Bland.
int chooseBasicToExit(VI vb, V xb, V db) {
   int n = xb.size();
   V v, db_negative;

   int posMin = -1;
   double min;
   for (int i = 0; i < n; ++i)
      if (db[i] < 0) {
         if (posMin == -1) {
            posMin= i;
            min = -xb[i]/db[i];
         }
         else {
            double value = -xb[i]/db[i];
            if (value < min or (abs(value - min) < TOLERANCE and vb[i] <= vb[posMin])) {
               posMin = i;
               min = value;
            }
         }
      }

   return posMin + 1;
}


// Diversas actualizaciones tras la iteración del símplex.
void update(M A, M &B, M &invB, VI &vb, VI &vn, V &xb, double &z, double theta,
            V db, int p, int qPos, int q, V r) {
   // Actualizar la matriz B (Cambiar la columna B(p) de B por la columna q de A).
   for (int i = 0; i < (int) B.size(); ++i)
      B[i][p - 1] = A[i][q - 1];

   // Actualizar la inversa de B.
   multiplyMatrixRow(invB, p - 1, -1/db[p - 1]);
   for (int i = 0; i < (int) invB.size(); ++i)
      if (i != p - 1)
         combineMatrixRows(invB, i, p - 1, db[i]);

   // Actualizar variables básicas y no básicas.
   xb = vectorSum(xb, scalarVectorMultiplication(theta, db));
   xb[p - 1] = theta;
   vn[qPos] = vb[p - 1];
   vb[p - 1] = q;
   
   // Actualizar el valor de la función.
   double rq = r[qPos];
   z += theta*rq;
}


// Forzar que el vector b sea positivo (multiplicando por -1 cuando sea necesario).
void ensurePositivity(M &A, V &b) {
   int m = b.size();
   int n = A[0].size();
   for (int i = 0; i < m; ++i)
      if (b[i] < 0) {
         b[i] = -b[i];
         for (int j = 0; j < n; ++j)
            A[i][j] = -A[i][j];
      }
}


// Hace una iteración del símplex.
// Output: 1 si encontrada la óptima, 2 si problema ilimitado, 3 si es degenerado.
// Retorna 0 si hay que continuar iterando.
int simplexIter(V c, M A, VI &vb, VI &vn, V &xb, double &z, M &B, M &invB, int rule) {
   V r = reducedCost(c, A, invB, vb, vn);
   if (minVectorValue(r) >= 0) {
      printDataSimplexIter(0, 0, 0, 0, z);
      return 1;
   }

   int qPos = chooseNonBasicToEnter(vn, r, rule);
   int q = vn[qPos];
   
   V db = descentDirection(invB, A, q);
   if (minVectorValue(db) >= 0) {
      printDataSimplexIter(q, r[qPos], 0, 0, z);
      return 2;
   }
   
   int p = chooseBasicToExit(vb, xb, db);

   double theta = -xb[p - 1]/db[p - 1];

   int vbp = vb[p - 1];
   update(A, B, invB, vb, vn, xb, z, theta, db, p, qPos, q, r);
   printDataSimplexIter(q, r[qPos], vbp, theta, z);

   return 0;
}


// Fase I del símplex.
dmm simplexPhaseI(M A, V b, VI &vb, VI &vn, V &xb, int rule, int &iterCounter) {
   int m = A.size();
   int n = A[0].size();

   // Crear la matriz A prima del problema auxiliar.
   M Ap(m, V(n + m, 0));
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j)
         Ap[i][j] = A[i][j];
      Ap[i][i + n] = 1;
   }
   // Crear el vector c prima del problema auxiliar.
   V cp(n + m , 0);
   for (int i = n; i < n + m; ++i)
      cp[i] = 1;
   // Crear los vectores de variables básicas y no básicas para el problema auxiliar.
   VI vbp(m);
   for (int i = 0; i < m; ++i)
      vbp[i] = i + n + 1;
   VI vnp(n);
   for (int i = 0; i < n; ++i)
      vnp[i] = i + 1;
   V xbp = b;
   // Calcular el valor actual de la función auxiliar z prima.
   double zp = 0;
   for (int i = 0; i < m; ++i)
      zp += b[i];

   // Matriz B auxiliar y su inversa (ambas son la identidad).
   M Bp = matrixIdentity(m); M invBp = matrixIdentity(m);
   
   // Iterar el símplex para el problema auxiliar.
   int iout = 0;
   while (iout == 0) {
      cout << "\tIteración " << iterCounter << ": \t";
      iout = simplexIter(cp, Ap, vbp, vnp, xbp, zp, Bp, invBp, rule);
      ++iterCounter;
   }

   // Sacar variables artificiales de la base si fuera necesario.
   for (int i = 0; i < m; ++i)
      if (vbp[i] > n) {
         M auxM = matrixMultiplication(invBp, A);

         int cont = 0;
         bool finished = false;
         for (int j = 0; j < m and not finished; ++j)
            if (abs(auxM[i][j]) < TOLERANCE) // Si es igual a 0 (precisión de máquina).
               ++cont;
            else {
               V dbp = descentDirection(invBp, Ap, j + 1);
               update(Ap, Bp, invBp, vbp, vnp, xbp, zp, 0, dbp, i + 1, 0, j + 1, {0});
               finished = true;
            }
         if (cont == m)
            cout << "La restricción "  << i + 1 << "es redundante." << endl;
      }
      
   // Guardar la solución básica factible.
   vb = vbp;
   for (int i = 0; i < n; ++i)
      if (vnp[i] <= n)
         vn.push_back(vnp[i]);
   xb = xbp;

   // Devolver el valor de la función auxiliar y la mariz B auxiliar y su inversa.
   return {zp, Bp, invBp};
}


void simplexPhaseII(V c, M A, M B, M invB, V b, VI vb, VI vn, V xb, int rule,
                    int &iterCounter) {
   // Calcular el valor de la función z con la SBF inicial.
   double z = 0;
   for (int i = 0; i < (int) vb.size(); ++i)
      z += c[vb[i] - 1]*xb[i];


   int iout = 0;
   while (iout == 0) {
      cout << "\tIteración " << iterCounter << ": \t";
      iout = simplexIter(c, A, vb, vn, xb, z, B, invB, rule);
      ++iterCounter;
   }
   
   if (iout == 1) {
      cout.setf(ios::fixed);
      cout.precision(6);
      cout << endl << "\tSolución óptima encontrada, iteración "
           << iterCounter - 1 << ", z = " << z << "." << endl;
      cout << endl << endl << "Fin símplex primal." << endl << endl;

      // Mostrar solución por pantalla.
      V r = reducedCost(c, A, invB, vb, vn);
      printSolution(vb, xb, vn, r, z);
   }
   else if (iout == 2) {
      cout << endl << "\tEl problema es ILIMITADO." << endl;
      cout << endl << endl << "Fin símplex primal." << endl;
   }
}


// Función principal del símplex.
// Rule es el procedimiento de elección de variables: 1 cr más negativo o 2 Bland.
void simplex(V c, M A, V b, int rule) {
   if (rule == 1)
      cout << "Inicio símplex primal con la regla del coste reducido más negativo."
           << endl << endl;
   else 
      cout << "Inicio símplex primal con regla de Bland." << endl << endl;

   // Forzar que el vector b sea positivo (multiplicando por -1 cuando sea necesario).
   ensurePositivity(A, b);

   int iterCounter = 1;
   int m = A.size();
   
   // Fase I del símplex.
   // Tomar la solución básica factible encontrada.
   cout << endl << "    Fase I" << endl;
   
   VI vb(m), vn;
   V xb(m);
   dmm simplexPhaseIdata = simplexPhaseI(A, b, vb, vn, xb, rule, iterCounter);
   double zp = simplexPhaseIdata.z;
   M B = simplexPhaseIdata.B;
   M invB = simplexPhaseIdata.invB;
   
   if (zp > TOLERANCE) {
      cout << endl << "\tEl problema es INFACTIBLE." << endl;
      return;
   }
   cout << endl << "\tSolución básica factible encontrada, iteración "
        << iterCounter - 1 << "." << endl << endl << endl;

   // Fase II del símplex.
   cout << "    Fase II" << endl;
   simplexPhaseII(c, A, B, invB, b, vb, vn, xb, rule, iterCounter);
}


// Main para probar el código.
int main() {
   // Mostrar los números con 3 decimales.
   cout.setf(ios::fixed);
   cout.precision(3);

   // Leer datos.
   string aux;
   int m, n, rule;
   cin >> aux;   
   cin >> rule;
   cin >> aux;   
   cin >> m;
   cin >> aux;
   cin >> n;
   
   cin >> aux;
   V c(n);
   for (int i = 0; i < n; ++i) {
      double aux;
      cin >> aux;
      c[i] = aux;
   }

   cin >> aux;   
   M A(m, V(n));
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) {
         double aux;
         cin >> aux;
         A[i][j] = aux;
      }

   cin >> aux;
   V b(m);
   for (int i = 0; i < m; ++i) {
      double aux;
      cin >> aux;
      b[i] = aux;
   }

   // Ejecutar función símplex.
   simplex(c, A, b, rule);
}
