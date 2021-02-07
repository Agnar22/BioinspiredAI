import LinReg
import math

reg = LinReg.LinReg()
pre_computed = dict()


def sin(value: int):
  return math.sin(value), math.sin(value) + 1


def linear(value: int):
  return value, value


def lin_reg(genes: str):
  if genes in pre_computed:
    return pre_computed[genes]
  x = reg.get_columns(reg.x, genes)
  fitness = reg.get_fitness(x, reg.y)
  pre_computed[genes] = (fitness, 1 / fitness)
  return fitness, 1 / fitness
