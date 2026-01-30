#include "roots.hpp"
#include <cmath>
#include <algorithm>

const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 1000000;

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{

    // Check if a and b have opposite signs
    if (f(a) * f(b) >= 0.0)
    {
        return false;
    }

    double c;
    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        c = (a + b) / 2.0;
        double fc = f(c);

        // Check if we found the root
        if (std::abs(fc) < TOLERANCE || std::abs(b - a) < TOLERANCE)
        {
            *root = c;
            return true;
        }

        // Determine which half contains the root
        if (f(a) * fc < 0.0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }

    *root = (a + b) / 2.0;
    return true;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
  
    // Check if a and b have opposite signs
    if (f(a) * f(b) >= 0.0)
    {
        return false;
    }

    double fa = f(a);
    double fb = f(b);
    double c;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        // Compute the false position
        c = a - (fa * (b - a)) / (fb - fa);

        // Check if c is within bounds
        if (c < std::min(a, b) || c > std::max(a, b))
        {
            return false;
        }

        double fc = f(c);

        // Check if we found the root
        if (std::abs(fc) < TOLERANCE)
        {
            *root = c;
            return true;
        }

        // Determine which half contains the root
        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    *root = c;
    return true;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
 
    double x = c;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        double fx = f(x);
        double gx = g(x);

        // Check if derivative is zero
        if (std::abs(gx) < TOLERANCE)
        {
            return false;
        }

        // Newton-Raphson step
        double x_next = x - fx / gx;

        // Check if we left the interval
        if (x_next < std::min(a, b) || x_next > std::max(a, b))
        {
            return false;
        }

        // Check convergence
        if (std::abs(x_next - x) < TOLERANCE)
        {
            *root = x_next;
            return true;
        }

        x = x_next;
    }

    *root = x;
    return true;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    
    double x_prev = c;
    double x_curr = c + 0.001; // Small perturbation for initial second point

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);

        // Check if denominator is too small
        if (std::abs(f_curr - f_prev) < TOLERANCE)
        {
            return false;
        }

        // Secant method step
        double x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

        // Check if we left the interval
        if (x_next < std::min(a, b) || x_next > std::max(a, b))
        {
            return false;
        }

        // Check convergence
        if (std::abs(x_next - x_curr) < TOLERANCE)
        {
            *root = x_next;
            return true;
        }

        x_prev = x_curr;
        x_curr = x_next;
    }

    *root = x_curr;
    return true;
}
