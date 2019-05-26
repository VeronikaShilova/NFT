#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <cassert>

using namespace std;

const complex<double> i(0, 1); //unit imaginary number
const double pi = 3.14159265358979;

#include "Potential.h"
#include "ContourRectangle.h"
#include "BOsborne.h"
#include "PolynomRootsFinder.h"
#include "IterativeMethodsForDiscreteSp.h"
#include "HybridContourIntegration.h"
#include "InnerBordering.h"
#include "DifferentialEquationsINFT.h"
#include "Error.h"