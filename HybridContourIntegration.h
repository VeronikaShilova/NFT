#pragma once

#include "ProgramParams.h"

vector<complex<double>> aDL(Potential const & signal, RectangularContour const & contour, 
                            const int contourPoints, const int timePoints) {
  int numberOfZeros;
  complex<double> integral(0, 0);
  vector<complex<double>> F(contourPoints), z(contourPoints);
  
  //z = contour.getPoints(contourPoints);

  for (int m = 0; m < contourPoints; m++) {
    vector<complex<double>> a(2);
    //a = BOsborne(signal, z[m], timePoints);
    a = BOsborne(signal, contour.z(m), timePoints);
    F[m] = a[0];
  }
   
  for (int m = 1; m < contourPoints; m++) {
    integral += (1.0 - F[m - 1] / F[m]);
  }
  numberOfZeros = (int)round((integral / (2.0 * pi*i)).real());
  vector<complex<double>> roots(numberOfZeros);

  if (numberOfZeros != 0) {

    vector<complex<double>> s(numberOfZeros, 0. + i*0.), sigma(numberOfZeros, 0. + i*0.);

    for (int k = 0; k < numberOfZeros; k++) {
      for (int m = 1; m < contourPoints; m++) {
        //s[k] += pow(z[m], k + 1)*(1.0 - F[m - 1] / F[m]);
        s[k] += pow(contour.z(m), k + 1) * (1.0 - F[m - 1] / F[m]);
      }
      s[k] = s[k] / (2.0*pi*i);
    }
    sigma[0] = -s[0];

    for (int k = 1; k < numberOfZeros; k++) {
      for (int j = 1; j <= k; j++) {
        sigma[k] += sigma[j - 1] * s[k - j];
      }
      sigma[k] += s[k];
      sigma[k] = -sigma[k] / (double)(k + 1);
    }

    Polynom poly(sigma);
    roots = findRoots(numberOfZeros, sigma);
    //roots = findRoots(poly);

    for (int j = 0; j < numberOfZeros; j++) {
      roots[j] = Secant(signal, roots[j], timePoints);
    }

    sort(begin(roots), end(roots),
      [](complex<double> item1, complex<double> item2) {return item1.imag() < item2.imag(); });

    for (int j = 0; j < numberOfZeros; j++) {
      vector<complex<double>> a = BOsborne(signal, roots[j], timePoints);
      complex<double> aPrime = BOsbornePrime(signal, roots[j], timePoints);
      roots.push_back(1. / (a[1] * aPrime));
    }
  }
  return roots;
}
