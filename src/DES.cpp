/*
        dy/dx = ax - by,
        y(0) = d;

        a = 0.3
        b = 0.2
        d = 1

        C = d + a/b^2 = 1 + (0.3)/0.2^2 = 8.5

        Solve on an interval of 0.01 with a step of 0.001
*/

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double getC(double a, double b, double d) { return d + a / (b * b); }

double f(double x, double y) { return (-2 * x - 0.2 * y); }

// Const step
double runge4(double t, double y0, double step) {
  double k1, k2, k3, k4;
  k1 = f(t, y0);
  k2 = f(t + step / 2, y0 + step * k1 / 2);
  k3 = f(t + step / 2, y0 + step * k2 / 2);
  k4 = f(t + step, y0 + step * k3);
  double y = y0 + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
  return y;
}

// Autostep
double runge45(double t, double y0, double *step, double tol) {
  double k1, k2, k3, k4, k5, k6;
  double h2 = *step / 2.0;
  double h6 = *step / 6.0;

  if (*step < 1e-6) {
    return y0;
  }

  k1 = *step * f(t, y0);
  k2 = *step * f(t + h2, y0 + k1 / 2.0);
  k3 = *step * f(t + h2, y0 + k2 / 2.0);
  k4 = *step * f(t + *step, y0 + k3);
  k5 = *step * f(t + *step / 2.0, y0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 8.0);
  k6 = *step * f(t + *step / 2.0, y0 + (k1 - 2.0 * k2 + 2.0 * k3 - k4) / 8.0);

  double y5 = y0 + (k1 + 4.0 * k2 + 6.0 * k3 + 4.0 * k4 + k5) * h6;
  double y4 = y0 + (k1 + 4.0 * k2 + 6.0 * k3 + 4.0 * k4 + k6) * h6;

  double delta =
      abs(y5 - y4) / (pow(tol, 4) + pow(tol * max(abs(y4), abs(y5)), 4));

  if (delta > 1.0) {
    // Reduce step size if error is too high
    *step *= 0.5 * pow(delta, -0.2);
  } else if (delta > 0.01) {
    // Reduce step size more gradually if error is moderate
    *step *= 0.9;
  } else if (delta > 0.0001) {
    // Increase step size if error is small
    *step *= 1.5;
  } else {
    // Increase step size more if error is very small
    *step *= 2.0;
  }

  y0 = y5;
  t += *step;

  //  cout << "step = " << *step << endl;
  return y0;
}

double calcMaxDelta(double C, std::vector<double> y, std::vector<double> time,
                    double step) {
  cout << endl
       << "Analytical        "
       << "Runge-Kutta            "
       << "Delta               "
       << "Step\n\n";

  double max_delta = 0;

  for (int i = time.size() - 10; i < time.size(); i++) {
    double u = 0.3 / 0.2 * (time[i] - 1 / 0.2) + C * exp(-0.2 * time[i]);
    double delta = abs(u - y[i]);
    max_delta = max(delta, max_delta);

    cout << "u(t) =" << setw(11) << u << " y(t) =" << setw(11) << y[i]
         << setw(15) << delta << setw(15) << step << endl;
  }
  return max_delta;
}

int main() {
  double initVal = 1;
  double step = 0.001; // starting step size
  double t = 0;        // initial time
  double T = 0.01;     // end time
  double tol = 1e-4;   // tolerance for error control
  double C = getC(0.3, 0.2, 1);

  std::vector<double> time, Time, y, Y;
  time.push_back(t);
  Time.push_back(t);
  y.push_back(initVal);
  Y.push_back(initVal);

  // Without step adjustment
  while (t < T) {
    double dy = runge4(Time.back(), Y.back(), step);
    Time.push_back(t);
    Y.push_back(dy);
    t += step;
  }

  double max_delta = calcMaxDelta(C, Y, Time, step);
  cout << "\nMaximum deviation without step adjustment: " << max_delta
       << "\nStep = " << step << "\n\n";

  // With step adjustment
  t = 0;
  int max_iter = 15;
  int iter = 0; // iteration counter
  while (t < T && iter < max_iter) {
    // calculation of the time step and the Runge-Kutta method
    double dy = runge45(t, y.back(), &step, tol);
    y.push_back(dy);
    t += step;
    time.push_back(t);

    // checking to reach time T
    if (t + step > T) {
      step = T - t;
    }

    iter++; // iteration counter++

    // message about reaching the maximum number of iterations
    if (iter == max_iter) {
      cout << "Maximum number of iterations reached!" << endl;
    }
  }
  max_delta = calcMaxDelta(C, y, time, step);
  cout << "\nMaximum deviation with step adjustment: " << max_delta << "\n"
       << "Autostep \n\n";

  return 0;
}
