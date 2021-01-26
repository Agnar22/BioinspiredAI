import numpy as np
import matplotlib.pyplot as plt


def draw_sin(x: np.ndarray, y: np.ndarray):
  x_sin = np.arange(0, 128, 0.1)
  y_sin = np.sin(x_sin)
  plt.plot(x_sin, y_sin)
  plt.scatter(x, y-1, color='orange')
  plt.show()
