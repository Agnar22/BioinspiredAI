import LinReg
import math
import param
import numpy as np

reg = LinReg.LinReg()
pre_computed = dict()

def sin(value: int):
  return math.sin(value) + 1


def linear(value: int):
  return value

def lin_reg(genes: str, actual=False):
  if genes in pre_computed:
    return pre_computed[genes]
  x = reg.get_columns(reg.x, genes)
  fitness = reg.get_fitness(x, reg.y)
  pre_computed[genes] = np.exp(-fitness)
  if actual:
    return fitness
  return np.exp(-fitness)