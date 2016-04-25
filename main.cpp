
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

double** newArray(const int x, const int y) {
  double** array = new double*[y];
  for (int i = 0; i < y; i++) {
    array[i] = new double[x];
    for (int j = 0; j < x; j++) {
      array[i][j] = 0.;
    }
  }
  return array;
}

double** copyArray(double** source, int x, int y) {
  double** array = newArray(x,y);
  for (int i = 0; i < y; i++) {
    for (int j = 0; j < x; j++) {
      array[i][j] = source[i][j];
    }
  }
  return array;
}

void deleteArray(double** source, int x, int y) {
  for (int i = 0; i < y; i++) {
    delete[] source[i];
  }
  delete[] source;
}

void printArray(double** array, int x, int y) {
  for (int i = 0; i < y; i++) {
    for (int j = 0; j < x; j++) {
      cout << array[i][j] << " ";
    }
    cout << "\n";
  }
}

double Xi(double x, int size) { return (x/size) - 0.5; }

double Yi(double y, int size) { return (y/size) - 0.5; }

double ro(double x, double y) {
  double delta = 0.01;
  int X = 100;
  int Y = 80;
  // double y1 = 0.64;
  // double x2 = 0.88;
  // double y2 = 0.64;

  // double Q = 80;

  // double warunek1 = pow((Xi(x, delta) - x1), 2) + pow((Yi(y, delta) - y1), 2);
  // double warunek2 = pow((Xi(x, delta) - x2), 2) + pow((Yijy,idelta) - y2), 2);

  // if (warunek1 <= pow(R, 2)) {
  //   return Q;
  // } else {
  //   if (warunek2 <= pow(R, 2)) {
  //     return -Q;
  //   } else {
  //     return 0;
  //   }
  // }
	return 100.*Xi(x,X)*Yi(y,Y)*exp(-50.*(pow(Xi(x,X),2)+pow(Xi(y,Y),2)));
}

double zadanie2(double** array, int x, int y, double DELTA, double w) {
  stringstream ss;
  ss << "Zad2_" << w << ".dat";
  string nazwaPliku = ss.str();

 cout<<"zadanie2"<<endl;

  int k = 1;
  double sum_previous = 1.;
  double sum_current = 0.;
  double TOL = 10e-8;
  int counter = 0;
  std::ofstream file;
  file.open(nazwaPliku.c_str());
  file.precision(8);
  if (file.good()) {
    while (fabs((sum_current - sum_previous) / sum_previous) > TOL) {
      // cout<<"wha";
      counter++;
      for (int i = 1; i < y - 1; i++) {
        for (int j = 1; j < x - 1; j++) {
          array[i][j] =
              ((1. - w) * array[i][j]) +
              (w * ((array[i + k][j] + array[i - k][j] + array[i][j + k] +
                     array[i][j - k] + pow((k * DELTA), 2) * ro(j, i)) /
                    4.));
        }
      }
      sum_previous = sum_current;
      sum_current = 0;
      for (int i = 0; i < y - 1; ++i) {
        for (int j = 0; j < x -1; ++j) {
          sum_current +=
              pow((k * DELTA), 2) *
              (0.5 * (pow(((array[i + k][j] - array[i][j]) / (k * DELTA)), 2) +
                      pow(((array[i][j + k] - array[i][j]) / (k * DELTA)), 2)) -
               (ro(j, i) * array[i][j]));
        }
      }
      file << counter << '\t' << sum_current << endl;
    }
  }
  std::ofstream file1;
  file1.open("zad2_gestosc.dat");
  file1.precision(8);
  for (int i = 0; i < y - 1; ++i) {
        for (int j = 0; j < x -1; ++j) {
        	file1 << j <<"\t"<< i <<"\t" << ro(j,i)<<endl;
        }
        file1 <<endl;
    }
  // printArray(array,size);
  file.close();
  file1.close();

    std::ofstream file2;
  file2.open("zad2_C.dat");
  file2.precision(8);
  for (int i = 1; i < y - 1; ++i) {
        for (int j = 1; j < x -1; ++j) {
        	double dx = (array[i+1][j] + array[i-1][j] - (2*array[i][j]))/pow(DELTA,2);
        	double dy = (array[i][j+1] + array[i][j-1] - (2*array[i][j]))/pow(DELTA,2);
        	file2 << j <<"\t"<< i <<"\t" << dx+dy+ro(j,i)<<endl;
        }
        file2 << endl;
    }
   file2.close();
    return counter;
}

double zadanie1(double** array, int x, int y, double DELTA, double w) {
stringstream ss;
  ss << "Zad1_" << w << ".dat";
  string nazwaPliku = ss.str();

  cout<<"zadanie1"<<endl;
  int counts[9];
  int k = 1;
  int temp = 0;
  double sum_previous = 1.;
  double sum_current = 0.;
  double TOL = 10e-8;
  int counter = 0;
  std::ofstream file;
  file.open(nazwaPliku.c_str());
  file.precision(8);
  if (file.good()) {
    // cout << fabs((sum_current - sum_previous) / sum_previous) << "\t" << TOL;
    while (fabs((sum_current - sum_previous) / sum_previous) > TOL) {
      counter++;
      double** copy = copyArray(array,x,y);

      for (int i = 1; i < y - 1; i++) {
        for (int j = 1; j < x - 1; j++) {
          array[i][j] =
              ((1. - w) * copy[i][j]) +
              (w * ((copy[i + k][j] + copy[i - k][j] + copy[i][j + k] +
                     copy[i][j - k] + pow((k * DELTA), 2) * ro(j, i)) /
                    4.));
        }
      }
      sum_previous = sum_current;
      sum_current = 0;
      for (int i = 0; i < y - 1; ++i) {
        for (int j = 0; j < x -1; ++j) {
          sum_current +=
              pow((k * DELTA), 2) *
              (0.5 * (pow(((array[i + k][j] - array[i][j]) / (k * DELTA)), 2) +
                      pow(((array[i][j + k] - array[i][j]) / (k * DELTA)), 2)) -
               (ro(j, i) * array[i][j]));
        }
      }
      file << counter << '\t' << sum_current << endl;
    }

  }
  // printArray(array,size);
  file.close();

  std::ofstream file1;
  file1.open("zad1_gestosc.dat");
  file1.precision(8);
  for (int i = 0; i < y - 1; ++i) {
        for (int j = 0; j < x -1; ++j) {
        	file1 << j <<"\t"<< i <<"\t" << ro(j,i)<<endl;
        }
        file1 << endl;
    }
   file1.close();

  std::ofstream file2;
  file2.open("zad1_C.dat");
  file2.precision(8);
  for (int i = 1; i < y - 1; ++i) {
        for (int j = 1; j < x -1; ++j) {
        	double dx = (array[i+1][j] + array[i-1][j] - (2*array[i][j]))/pow(DELTA,2);
        	double dy = (array[i][j+1] + array[i][j-1] - (2*array[i][j]))/pow(DELTA,2);
        	file2 << j <<"\t"<< i <<"\t" << dx+dy+ro(j,i)<<endl;
        }
        file2 << endl;
    }
   file2.close();

      return counter;
}

int counter = 0 ;

void zadanie3(
std::ofstream& file, double** array, int x,int y, double DELTA, int k) {
  cout<<"zadanie3"<<endl;

  int w = 1.9;
  double sum_previous = 1.;
  double sum_current = 0.;
  double TOL = 10e-8;
    while (fabs((sum_current - sum_previous) / sum_previous) > TOL) {
      //cout<<"wha";
      counter++;
      for (int i = k; i < y - k; i+=k) {
        for (int j = k; j < x - k; j+=k) {
          array[i][j] =
              ((1. - w) * array[i][j]) +
              (w * ((array[i + k][j] + array[i - k][j] + array[i][j + k] +
                     array[i][j - k] + pow((k * DELTA), 2) * ro(j, i)) /
                    4.));
        }
      }
      sum_previous = sum_current;
      sum_current = 0;
      for (int i = 0; i < y - k; i+=k) {
        for (int j = 0; j < x -k; j+=k) {
          sum_current +=
              pow((k * DELTA), 2) *
              (0.5 * (pow(((array[i + k][j] - array[i][j]) / (k * DELTA)), 2) +
                      pow(((array[i][j + k] - array[i][j]) / (k * DELTA)), 2)) -
               (ro(j, i) * array[i][j]));
        }
      }
      file << counter << '\t' << sum_current << endl;
    }
         for (int i = 0; i < y -k; i++) {
        for (int j = 0; j < x -k; j++) {
          array[(i+k)/2][(j+k)/2] = (array[i][j]+array[i+k][j]+array[i][j+k]+array[i+k][j+k])/4.;
          array[(i+k)/2][j] = (array[i][j]+array[i+k][j])/2.;
        }
      }
  // printArray(array,size);
  }


int main(int argc, char const* argv[]) {
  int X = 100;
  int Y = 80;
  double DELTA = 0.01;
  double** array = newArray(X,Y);
  // // DIRICHLET
  // for (int i = 1; i < (Y - 1); i++) {
  //   array[i][0] = 0.1;
  // }

  // for (int i = 1; i < Y - 1; i++) {
  //   array[i][X - 1] = -0.1;
  // }
/*
    std::ofstream file3;
  	file3.open("zad2_count.dat");
  	file3.precision(8);
  for (double i = 1.1; i <= 2.0; i += 0.1) {
    double** copy = copyArray(array, X,Y);
    int count = zadanie2(copy, X,Y, DELTA, i);
    file3<< i << '\t' << count << endl;
    deleteArray(copy,X,Y);
  }
  file3.close();


  // zadanie1(array,X,Y,DELTA);
    std::ofstream file4;
  	file4.open("zad1_count.dat");
  	file4.precision(8);
  for (double i = 0.1; i <= 1.0; i += 0.1) {
    double** copy = copyArray(array, X,Y);
    int count = zadanie1(copy, X,Y, DELTA, i);
    	file4<< i << '\t' << count << endl;
    deleteArray(copy,X,Y);
  }
  	file4.close();

*/
std::ofstream file;
  file.open("Zad3.dat");
  file.precision(8);
  if (file.good()) {
 
  	double** copy = copyArray(array, X,Y);
  for (int k = 16; k >= 1; k /= 2) {
    zadanie3(file,copy, X,Y, DELTA, k);
   

  }
  deleteArray(copy,X,Y);
}
file.close();
  
  // printArray(array,X,Y);
  deleteArray(array, X,Y);

  return 0;
}
