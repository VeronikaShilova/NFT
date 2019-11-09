#include "ProgramParams.h"

const string path = "Output/";

int darkSolitonPotential(const unsigned int N) {
  const string Mode = "TIB";
  const double T1 = 0.;
  const double T2 = 3.;
  const double alpha = 1.;
  const double nu = 0.8;
  darkSoliton signal(alpha, nu, T1, T2);
  double h = (T2 - T1) / N; 
  ofstream file;
  vector<complex<double>> a = signal.get_analytical(N);
  file.open(path + "Analytical_solution_alpha=" + to_string(alpha) + "_nu=" 
            + to_string(nu) + "_" + Mode + ".txt");
  for (int j = 0; j < a.size(); j++) {
    file << T1 + j*h << " " << a[j].real() << endl;
  }
  file.close();
  vector<complex<double>> b = InnerBordering(N, signal);
  file.open(path + "Numerical_solution_alpha=" + to_string(alpha) + "_nu=" 
            + to_string(nu) + "_" + Mode + ".txt");
  for (int j = 0; j < b.size(); j++) {
    file << T1 + j*h << " " << b[j].real() << endl;
  }
  file.close();
  double c = ErrorL2(a, b, h);
  file.open(path + "Error_L2_alpha=" + to_string(alpha) + "_nu=" + to_string(nu) + "_" 
            + Mode + ".txt", ios::app);
  file << N << " " << c << endl;
  file.close();
  return 0;
}

int SatsumaPotential(const unsigned int N) {
  const double T1 = -10.;
  const double T2 = 10.;
  const double A = 5.25;
  Satsuma signal(A, T1, T2);
  double timeStep = (T2 - T1) / N;
  const int contPoints = N;
  const int xiPoints = N;
  double xiStep = (signal.Xi2_sp() - signal.Xi1_sp()) / xiPoints;
  ofstream file;
  vector<complex<double>> analytical_solution = signal.get_analytical(N);
  file.open(path + "Analytical_solution_SY.txt");
  for (int j = 0; j < analytical_solution.size(); j++) {
    file << T1 + j*timeStep << " " << analytical_solution[j].real() << endl;
  }
  file.close();
 vector<complex<double>> cont_sp(xiPoints);
  for (int j = 0; j < xiPoints; j++) {
    vector<complex<double>> a = BOsborne(signal, signal.Xi1_sp() + j*xiStep, N);
    cont_sp[j] = conj(a[1]) / a[0];
  }
  signal.set_cont_sp(move(cont_sp));
  cout << "Continuous spectrum is done!" << endl;
  signal.set_disc_sp(aDL(signal, RectangularContour(-1., 1., 0.1, 5., contPoints), contPoints, N));
  cout << "Discrete spectrum is done!" << endl;
  for (const auto item : signal.disc_sp()) {
    cout << "Eigenvalue: " << item << endl;
  }
  vector<complex<double>> numerical_solution = InnerBordering(N, signal);
  file.open(path + "Numerical_solution_SY.txt");
  for (int j = 0; j < numerical_solution.size(); j++) {
    file << T1 + j*timeStep << " " << numerical_solution[j].real() << endl;
  }
  file.close();
  double c = ErrorL2(analytical_solution, numerical_solution, timeStep);
  file.open(path + "Error_L2_SY_A="+to_string(A)+".txt", ios::app);
  file << N << " " << c << endl;
  file.close();
  return 0;
}

int SolitonPotential(const unsigned int N, const double numberOfSolitons) {
  const string Mode = "TIB";
  nSoliton signal(numberOfSolitons, 15, 1.0);
  const double T1 = signal.T1();
  const double T2 = signal.T2();
  double h = (T2 - T1) / N;
  const int contPoints = N;
  const int xiPoints = N;
  double xiStep = (signal.Xi2_sp() - signal.Xi1_sp()) / xiPoints;
  ofstream file;
  vector<complex<double>> a = signal.get_analytical(N);
  file.open(path + "Analytical_solution_" + to_string(int(numberOfSolitons)) + "soliton.txt");
  for (int j = 0; j < a.size(); j++) {
    file << T1 + j * h << " " << abs(a[j]) << endl;
  }
  file.close();
  vector<complex<double>> cont_sp(xiPoints);
  for (int j = 0; j < xiPoints; j++) {
    vector<complex<double>> a = BOsborne(signal, signal.Xi1_sp() + j * xiStep, N);
    cont_sp[j] = conj(a[1]) / a[0];
  }
  signal.set_cont_sp(move(cont_sp));
  cout << "Continuous spectrum is done!" << endl;
  signal.set_disc_sp(aDL(signal, RectangularContour(-1., 1., 0.1, 5., contPoints), contPoints, N));
  cout << "Discrete spectrum is done!" << endl;
  for (const auto item : signal.disc_sp()) {
    cout << "Eigenvalue: " << item << endl;
  }
  vector<complex<double>> b = InnerBordering(N, signal);
  file.open(path + "Numerical_solution_" + to_string(int(numberOfSolitons)) + "soliton.txt");
  for (int j = 0; j < b.size(); j++) {
    file << T1 + j * h << " " << abs(b[j]) << endl;
  }
  file.close();
  double c = ErrorL2(a, b, h);
  file.open(path + "Error_L2_" + to_string(int(numberOfSolitons)) + "soliton.txt", ios::app);
  file << N << " " << c << endl;
  file.close();
  return 0;
}


int main() {
  const unsigned int N = 2048;
  //darkSolitonPotential(N);
  //SatsumaPotential(N);

  SolitonPotential(N, 3.);
  
  /*const double numberOfSolitons = 2.;
  nSoliton signal(numberOfSolitons, 30, 1.0);
  const double T1 = signal.T1();
  const double T2 = signal.T2();
  const int contPoints = N;
  const int xiPoints = N;
  signal.set_disc_sp(aDL(signal, RectangularContour(-1., 1., 0.1, 5., contPoints), contPoints, N));
  cout << "Discrete spectrum is done!" << endl;
  for (const auto item : signal.disc_sp()) {
    cout << "Eigenvalue: " << item << endl;
  }*/
  cout << "end"; 
  getchar();
  return 0;
}