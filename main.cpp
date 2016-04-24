
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

double** newArray(const int size) {
  double** array = new double*[size];
  for (int i = 0; i < size; i++) {
    array[i] = new double[size];
    for (int j = 0; j < size; j++) {
      array[i][j] = 0.;
    }
  }
  return array;
}

double** copyArray(double** source, int size) {
  double** array = newArray(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      array[i][j] = source[i][j];
    }
  }
  return array;
}

void deleteArray(double** source, int size) {
  for (int i = 0; i < size; i++) {
    delete[] source[i];
  }
  delete[] source;
}

void printArray(double** array, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      cout << array[i][j] << " ";
    }
    cout << "\n";
  }
}

double Xi(double x, double delta) { return x * delta; }

double Yi(double y, double delta) { return y * delta; }

double ro(double x, double y) {
  double delta = 0.01;
  double R = 0.051;
  double x1 = 0.4;
  double y1 = 0.64;
  double x2 = 0.88;
  double y2 = 0.64;

  double Q = 80;

  double warunek1 = pow((Xi(x, delta) - x1), 2) + pow((Yi(y, delta) - y1), 2);
  double warunek2 = pow((Xi(x, delta) - x2), 2) + pow((Yi(y, delta) - y2), 2);

  if (warunek1 <= pow(R, 2)) {
    return Q;
  } else {
    if (warunek2 <= pow(R, 2)) {
      return -Q;
    } else {
      return 0;
    }
  }
}

void zadanie2(double** array, int size, double DELTA, double w) {
  stringstream ss;
  ss << "Zad2_" << w << ".dat";
  string nazwaPliku = ss.str();

  int k = 1;
  double sum_previous = 1.;
  double sum_current = 0.;
  double TOL = 10e-6;
  int counter = 0;
  std::ofstream file;
  file.open(nazwaPliku.c_str());
  file.precision(8);
  if (file.good()) {
    cout << fabs((sum_current - sum_previous) / sum_previous) << "\t" << TOL;
    while (fabs((sum_current - sum_previous) / sum_previous) > TOL) {
      // cout<<"wha";
      counter++;
      for (int i = 1; i < size - 1; i++) {
        for (int j = 1; j < size - 1; j++) {
          array[i][j] =
              ((1. - w) * array[i][j]) +
              (w * ((array[i + k][j] + array[i - k][j] + array[i][j + k] +
                     array[i][j - k] + pow((k * DELTA), 2) * ro(i, j)) /
                    4.));
        }
      }
      sum_previous = sum_current;
      sum_current = 0;
      for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size -1; ++j) {
          sum_current +=
              pow((k * DELTA), 2) *
              (0.5 * (pow(((array[i + k][j] - array[i][j]) / (k * DELTA)), 2) +
                      pow(((array[i][j + k] - array[i][j]) / (k * DELTA)), 2)) -
               (ro(i, j) * array[i][j]));
        }
      }
      file << counter << '\t' << sum_current << endl;
    }
  }
  // printArray(array,size);
  file.close();
}

void zadanie1(double** array, int size, double DELTA, double w) {
stringstream ss;
  ss << "Zad1_" << w << ".dat";
  string nazwaPliku = ss.str();

  int k = 1;
  double sum_previous = 1.;
  double sum_current = 0.;
  double TOL = 10e-6;
  int counter = 0;
  std::ofstream file;
  file.open(nazwaPliku.c_str());
  file.precision(8);
  if (file.good()) {
    cout << fabs((sum_current - sum_previous) / sum_previous) << "\t" << TOL;
    while (fabs((sum_current - sum_previous) / sum_previous) > TOL) {
      // cout<<"wha";
      counter++;
      double** copy = copyArray(array,size);

      for (int i = 1; i < size - 1; i++) {
        for (int j = 1; j < size - 1; j++) {
          array[i][j] =
              ((1. - w) * copy[i][j]) +
              (w * ((copy[i + k][j] + copy[i - k][j] + copy[i][j + k] +
                     copy[i][j - k] + pow((k * DELTA), 2) * ro(i, j)) /
                    4.));
        }
      }
      sum_previous = sum_current;
      sum_current = 0;
      for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size -1; ++j) {
          sum_current +=
              pow((k * DELTA), 2) *
              (0.5 * (pow(((array[i + k][j] - array[i][j]) / (k * DELTA)), 2) +
                      pow(((array[i][j + k] - array[i][j]) / (k * DELTA)), 2)) -
               (ro(i, j) * array[i][j]));
        }
      }
      file << counter << '\t' << sum_current << endl;
    }
  }
  // printArray(array,size);
  file.close();}

void zadanie3() {}
int main(int argc, char const* argv[]) {
  int SIZE = 64;
  double DELTA = 0.01;
  double** array = newArray(SIZE);
  // DIRICHLET
  for (int i = 1; i < (SIZE - 1); i++) {
    array[i][0] = 0.1;
  }

  for (int i = 1; i < SIZE - 1; i++) {
    array[i][SIZE - 1] = -0.1;
  }

  for (double i = 1.1; i <= 2.0; i += 0.1) {
    double** copy = copyArray(array, SIZE);
    zadanie2(copy, SIZE, DELTA, i);
    deleteArray(copy,SIZE);
  }
  // zadanie1(array,SIZE,DELTA);

  for (double i = 0.1; i <= 1.0; i += 0.1) {
    double** copy = copyArray(array, SIZE);
    zadanie1(copy, SIZE, DELTA, i);
    deleteArray(copy,SIZE);
  }
  zadanie3();
  // printArray(array,SIZE);
  deleteArray(array, SIZE);

  return 0;
}
