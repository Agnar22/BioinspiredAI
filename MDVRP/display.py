import matplotlib.pyplot as plt
import os
import random
from typing import List

prob_dir = 'Data files project 2/Testing Data/Data Files/'
sol_dir = 'Data files project 2/Testing Data/Solution files/'


def read_problem_file(filename: str) -> List[List[float]]:
  customers = []
  depots = []

  with open(filename, mode='r') as f:
    lines = f.readlines()
    vcl_per_dep, num_cust, num_dep = map(int, lines[0].split())

    for line in range(num_dep + 1, num_cust + num_dep + 1):
      _, x, y, *_ = map(int, lines[line].split())
      customers.append([x, y])
    for line in range(num_cust + num_dep + 1, len(lines)):
      _, x, y, *_ = map(int, lines[line].split())
      depots.append([x, y])
  return customers, depots


def read_solution_file(filename: str) -> List[List[float]]:
  routes = []
  with open(filename, mode='r') as f:
    lines = f.readlines()
    for line in lines[1:]:
      dep_no, vcl_no = map(int, line.split()[:2])
      route_dur = float(line.split()[2])
      load = int(line.split()[3])
      customers = list(map(int, line.split()[5:]))
      route = [dep_no, *customers, dep_no]
      routes.append(route)
  return routes


def display_problem(ax: plt.axes, customers: List[List[float]], depots: List[List[float]], display=True, title=None):
  if title is not None:
    plt.title(title)
  cust_x, cust_y = map(list, zip(*customers))
  ax.scatter(cust_x, cust_y, marker='o')

  dep_x, dep_y = map(list, zip(*depots))
  ax.scatter(dep_x, dep_y, marker='x')
  if display:
    plt.show()


def display_problem_and_solution(customers: List[List[float]], depots: List[List[float]], routes: List[List[int]],
                                 title=None):
  fig = plt.figure(1)
  ax = fig.add_subplot(111)
  if title is not None:
    plt.title(title)
  display_problem(ax, customers, depots, display=False)
  depot_color = {}
  for route in routes:
    x, y = [], []
    for num, customer in enumerate(route):
      customer -= 1
      if num == 0 or num == len(route) - 1:
        x.append(depots[customer][0])
        y.append(depots[customer][1])
      else:
        x.append(customers[customer][0])
        y.append(customers[customer][1])
    #if route[0] in depot_color:
    depot_color[route[0]] = [random.random(), random.random(), random.random()]
    ax.plot(x, y, color=depot_color[route[0]])
  plt.show()


def display_problem_from_file(filename: str):
  display_problem(plt.figure(1).add_subplot(111), *read_problem_file(prob_dir + filename), display=True, title=filename)


def display_problem_and_solution_from_file(filename: str):
  display_problem_and_solution(*read_problem_file(prob_dir + filename), read_solution_file(sol_dir + filename + '.res'),
                               title=filename)


if __name__ == '__main__':
  #for filename in os.listdir(prob_dir):
  #  display_problem_from_file(filename)
  #  display_problem_and_solution_from_file(filename)
  problem_name = "p07"
  sol_dir = './build/'
  display_problem_and_solution_from_file(problem_name)
  sol_dir = 'Data files project 2/Testing Data/Solution files/'
  display_problem_and_solution_from_file(problem_name)
