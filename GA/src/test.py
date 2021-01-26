import numpy as np
import LinReg

lin_reg = LinReg.LinReg()
data = np.genfromtxt("data.csv", delimiter=',')
x = data[:, :-1]
y = data[:, -1:]
print(lin_reg.get_fitness(x, y))
x = data[:, 10:-1]
print(lin_reg.get_fitness(x, y))
x = data[:, 20:-1]
print(lin_reg.get_fitness(x, y))
x = data[:, 30:-1]
print(lin_reg.get_fitness(x, y))
