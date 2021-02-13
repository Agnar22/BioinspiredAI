import numpy as np
import matplotlib.pyplot as plt


def draw_sin(x: np.ndarray, y: np.ndarray, epoch):
  x_sin = np.arange(0, 128, 0.1)
  y_sin = np.sin(x_sin)
  plt.plot(x_sin, y_sin)
  plt.scatter(x, y, color='orange')
  plt.title(f'Epoch {epoch}')
  plt.show()

def draw(*args):
  for x, y, label in args:
    plt.plot(x, y, label=label)
  plt.legend()
  plt.show()