import scipy.integrate as integrate
import numpy as np

# Integral a
def integrand_a(z):
    return (np.exp(z) + 1) / (np.exp(z) + z)

result_a, error_a = integrate.quad(integrand_a, 0, 1)
print(f"Integral a: {result_a}")

# Integral b
def integrand_b(x):
    return (x * np.sin(x)) / (1 + np.cos(x)**2)

result_b, error_b = integrate.quad(integrand_b, 0, np.pi)
print(f"Integral b: {result_b}")
import random

# Monte Carlo Simulation for Integral a
def monte_carlo_a(num_points):
    total = 0
    for _ in range(num_points):
        z = random.uniform(0, 1)
        total += (np.exp(z) + 1) / (np.exp(z) + z)
    return total / num_points

result_mc_a = monte_carlo_a(100000)
print(f"Monte Carlo Result for Integral a: {result_mc_a}")

# Monte Carlo Simulation for Integral b
def monte_carlo_b(num_points):
    total = 0
    for _ in range(num_points):
        x = random.uniform(0, np.pi)
        total += (x * np.sin(x)) / (1 + np.cos(x)**2)
    return total * np.pi / num_points

result_mc_b = monte_carlo_b(100000)
print(f"Monte Carlo Result for Integral b: {result_mc_b}")
import numpy as np
from scipy import integrate


a, b = 0, 2


def f_y(y):
    return y**2

volume_y = np.pi * integrate.quad(lambda y: f_y(y)**2, a, b)[0]
print(f"Volume about y-axis: {volume_y}")

def f_x(x):
    return x**3

volume_x = np.pi * integrate.quad(lambda x: f_x(x)**2, a, b)[0]
print(f"Volume about x-axis: {volume_x}")

import random

# Monte Carlo Simulation
def monte_carlo_lune(num_points):
    count = 0
    for _ in range(num_points):
        x = random.uniform(-5, 5)
        y = random.uniform(-5, 5)
        if x**2 + y**2 <= 25 and x**2 + y**2 >= 9:
            count += 1
    area = (count / num_points) * (10 * 10)  # Bounding box area
    return area

area_lune = monte_carlo_lune(100000)
print(f"Area of Lune: {area_lune}")
import sympy as sp

# Define symbols
x, a, b, L, lambda_, epsilon_0 = sp.symbols('x a b L lambda epsilon_0')

# Define the integrand
integrand = (lambda_ * b) / (4 * sp.pi * epsilon_0 * (x**2 + b**2)**(3/2))

# Evaluate the integral
E_P = sp.integrate(integrand, (x, -a, L - a))
print(f"Electric Field E(P): {E_P}")