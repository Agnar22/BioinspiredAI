import math
import param


def sin(value: int):
  #return math.sin(int(genes, 2) / ((2 ** param.BITSTRING_SIZE) - 1) * 128) + 1
  return math.sin(value) + 1


def linear(value: int):
  #return int(genes, 2) / ((2 ** param.BITSTRING_SIZE) - 1)
  return value
