import param
import copy
import random
import fitness
import draw_graph
import math
import numpy as np
import matplotlib.pyplot as plt
import LinReg
from functools import reduce


class Individual:
  def __init__(self, gene_sequence: str = None):
    self.fitness = None
    self.actual_fitness = None
    if gene_sequence is not None:
      self.genes = gene_sequence
    else:
      self.genes = ''.join([str(random.randint(0, 1)) for _ in range(param.BITSTRING_SIZE)])


class SimpleGA:
  def __init__(self, population_size: int, fitness_function, genes_to_value_function, crowding):
    self.fitness_function = fitness_function
    self.genes_to_value_function = genes_to_value_function
    self.population_size = population_size
    self.population = []
    self.init_population()
    self.entropies = []
    self.all_time_best = math.inf
    self.crowding = crowding

  def init_population(self):
    for x in range(self.population_size):
      self.population.append(Individual())

  def simulate(self, mod_draw_fitness_graph, hamming):
    self.test_fitness(self.population)

    for x in range(param.EPOCHS):
      if mod_draw_fitness_graph is not None and x % mod_draw_fitness_graph == 0:
        self.draw_fitness_graph(x)
      self.add_entropy()

      selected_pairs = self.select()
      children = self.crossover(selected_pairs)
      children = self.mutate(children)
      if self.crowding:
        children = self.survive(self.population, selected_pairs, children, hamming)
        self.population = self.kill_weakest(children)
      else:
        self.population = self.kill_weakest(children)
      self.test_fitness(self.population, verbose=True)

      # if x % 5 == 0 and self.crowding:
      #   self.test_fitness(self.population, verbose=True)
        # print(f'Population size: {len(self.population)}')

  def print_population(self):
    for indiv in self.population:
      print(f'{indiv.genes} {self.fitness_function(self.genes_to_value_function(indiv.genes))}')

  def test_fitness(self, population, verbose=False):
    tot_fitness = 0
    best_fitness = math.inf if param.MINIMIZE_FITNESS else -math.inf
    for indiv in population:
      if indiv.fitness is None:
        indiv.actual_fitness, indiv.fitness = self.fitness_function(self.genes_to_value_function(indiv.genes))
        actual_fitness = indiv.actual_fitness
      else:
        actual_fitness = indiv.actual_fitness
      tot_fitness += actual_fitness
      best_fitness = min(actual_fitness, best_fitness) if param.MINIMIZE_FITNESS else max(actual_fitness, best_fitness)
    if verbose:
      print(f'Max fitness: {best_fitness:.5f}, average fitness: {tot_fitness / len(population):.5f}')
    self.all_time_best = min(self.all_time_best, best_fitness) if param.MINIMIZE_FITNESS else max(self.all_time_best,
                                                                                                  best_fitness)

  def draw_fitness_graph(self, num):
    x = [self.genes_to_value_function(indiv.genes) for indiv in self.population]
    y = [self.fitness_function(val)[0] for val in x]
    draw_graph.draw_sin(np.array(x), np.array(y), num)

  def select(self):
    tot_fitness = 0
    for indiv in self.population:
      tot_fitness += indiv.fitness

    selected_pairs = []
    for x in range(param.CHILDREN):
      parent_one = self.get_pos_from_fitness(random.random() * tot_fitness)
      parent_two = self.get_pos_from_fitness(random.random() * tot_fitness)
      while parent_one == parent_two:
        parent_one = self.get_pos_from_fitness(random.random() * tot_fitness)
        parent_two = self.get_pos_from_fitness(random.random() * tot_fitness)
      selected_pairs.append((parent_one, parent_two))
    return selected_pairs

  def get_pos_from_fitness(self, fitness):
    cumulative_sum = 0
    for x in range(len(self.population)):
      cumulative_sum += self.population[x].fitness
      if cumulative_sum >= fitness:
        return x
    raise ("Fitness not within cumulative sum range.")

  def crossover(self, selected_pairs):
    children = []
    for parents in selected_pairs:
      if random.random() >= param.CROSSOVER_PROB:
        children.append(copy.deepcopy(self.population[parents[0]]))
        children.append(copy.deepcopy(self.population[parents[1]]))
      else:
        idx = random.randint(1, param.BITSTRING_SIZE - 1)
        children.append(Individual(self.population[parents[0]].genes[0:idx] + self.population[parents[1]].genes[idx:]))
        children.append(Individual(self.population[parents[1]].genes[0:idx] + self.population[parents[0]].genes[idx:]))
    return children

  @staticmethod
  def mutate(children):
    for indiv in children:
      for pos in range(param.BITSTRING_SIZE):
        if random.random() < param.MUTATION_PROB:
          genes = list(indiv.genes)
          genes[pos] = '1' if indiv.genes[pos] == '0' else '0'
          indiv.genes = ''.join(genes)
    return children

  def survive(self, parents, selected_individuals, children, hamming):
    next_gen = []
    for pos, child in enumerate(children):
      child.actual_fitness, child.fitness = self.fitness_function(self.genes_to_value_function(child.genes))
      parent1 = selected_individuals[pos // 2][0]
      parent2 = selected_individuals[pos // 2][1]
      most_similar = reduce(
        lambda x, y: x if SimpleGA.distance(parents[x], child, hamming) < SimpleGA.distance(parents[y], child, hamming) else y,
        [parent1, parent2]
      )
      if parents[most_similar].fitness < child.fitness:
        #parents[most_similar] = child
        next_gen.append(copy.deepcopy(child))
      else:
        next_gen.append(copy.deepcopy(parents[most_similar]))
    return next_gen

  def add_entropy(self):
    on = [0 for _ in range(param.BITSTRING_SIZE)]
    for indiv in self.population:
      for pos, bit in enumerate(indiv.genes):
        if bit == "1":
          on[pos] += 1
    prob = np.array(on) / len(self.population)
    epsilon = 1e-8
    self.entropies.append(-np.sum((prob + epsilon) * np.log(prob + epsilon)))

  def get_entropies(self):
    return self.entropies

  @staticmethod
  def distance(indiv_1: Individual, indiv_2: Individual, hamming=False):
    if not hamming:
      return abs(int(indiv_1.genes) - int(indiv_2.genes))
    hamming_distance = 0
    for gene_1, gene_2 in zip(indiv_1.genes, indiv_2.genes):
      if gene_1 == gene_2:
        hamming_distance += 1
    return hamming_distance

  def kill_weakest(self, population):
    for indiv in population:
      indiv.actual_fitness, indiv.fitness = self.fitness_function(self.genes_to_value_function(indiv.genes))
    population.sort(key=lambda x: x.fitness, reverse=True)
    # self.population = self.population[:int((1 - param.KILL_LOWEST) * self.population_size)]
    return population[:param.POPULATION_SIZE]


def genes_to_value(start, stop, return_bitstring=False):
  def value_scaled(genes):
    if return_bitstring:
      return genes
    interval = stop - start
    return int(genes, 2) / ((2 ** param.BITSTRING_SIZE) - 1) * interval + start

  return value_scaled


def hyperparam_search(search_params):
  for kill_lowest, mut_prob, cross_prob, pop_size in search_params:
    param.CROWDING_FACTOR = kill_lowest
    param.MUTATION_PROB = mut_prob
    param.CROSSOVER_PROB = cross_prob
    param.CHILDREN = 1

    param.BITSTRING_SIZE = 101

    ga_cw_lin = SimpleGA(param.POPULATION_SIZE, fitness.lin_reg, genes_to_value(0, 128, return_bitstring=True), True)
    ga_cw_lin.simulate()
    ga_cw_lin.print_population()

    with open("hyperparameters.txt", "a+") as file:
      file.write(
        f'{param.EPOCHS},\t{kill_lowest},\t{mut_prob:.3f},\t{cross_prob:.2f},\t{pop_size}:\t{ga_cw_lin.all_time_best:.5f}\n')


if __name__ == '__main__':
  if param.HYPERPARAMETER_SEARCH:
    hyperparam_search(
      [
        [random.randint(1, 12), random.random() * 0.03, random.random() * 0.3 + 0.7, random.randint(5, 60)] for _ in
        range(200)
      ]
    )

  # Task 1A - 1E
  print("Task 1A - 1E")
  param.BITSTRING_SIZE = 32
  param.MINIMIZE_FITNESS = False
  param.MOD_DRAW_FITNESS = 5
  param.EPOCHS, param.POPULATION_SIZE, param.CHILDREN, param.CROSSOVER_PROB, param.MUTATION_PROB = 30, 400, 1200, 0.7, 0.02
  ga_sin = SimpleGA(param.POPULATION_SIZE, fitness.sin, genes_to_value(0, 128, return_bitstring=False), crowding=False)
  ga_sin.simulate(param.MOD_DRAW_FITNESS, False)
  ga_sin.print_population()

  # Task 1F
  print("Task 1F")
  param.BITSTRING_SIZE = 101
  param.MINIMIZE_FITNESS = True
  param.MOD_DRAW_FITNESS = None
  param.EPOCHS, param.POPULATION_SIZE, param.CHILDREN, param.CROSSOVER_PROB, param.MUTATION_PROB = 15, 150, 450, 0.7, 0.02
  lin_reg = LinReg.LinReg()
  print(f'Loss with all columns: {lin_reg.get_fitness(lin_reg.x, lin_reg.y):.5f}')
  ga_reg = SimpleGA(param.POPULATION_SIZE, fitness.lin_reg, genes_to_value(0, 128, return_bitstring=True), crowding=False)
  ga_reg.simulate(param.MOD_DRAW_FITNESS, True)
  print(f'All time best: {ga_reg.all_time_best}')

  # Task 1G
  print("Task 1G - Sin")
  param.BITSTRING_SIZE = 32
  param.MINIMIZE_FITNESS = False
  param.MOD_DRAW_FITNESS = 5
  param.EPOCHS, param.POPULATION_SIZE, param.CHILDREN, param.CROSSOVER_PROB, param.MUTATION_PROB = 30, 400, 200, 0.7, 0.02
  # Sin
  ga_cw_sin = SimpleGA(param.POPULATION_SIZE, fitness.sin, genes_to_value(0, 128, return_bitstring=False), crowding=True)
  ga_cw_sin.simulate(param.MOD_DRAW_FITNESS, False)
  ga_cw_sin.print_population()

  draw_graph.draw(
    [[x for x in range(len(ga_sin.get_entropies()))], ga_sin.get_entropies(), "sga_sin"],
    [[x for x in range(len(ga_cw_sin.get_entropies()))], ga_cw_sin.get_entropies(), "sga_cw_sin"]
  )


  # LinReg
  print("Task 1G - linreg")
  param.BITSTRING_SIZE = 101
  param.MINIMIZE_FITNESS = True
  param.MOD_DRAW_FITNESS = None
  param.EPOCHS, param.POPULATION_SIZE, param.CHILDREN, param.CROSSOVER_PROB, param.MUTATION_PROB = 15, 400, 300, 0.99, 0.02
  ga_cw_reg = SimpleGA(param.POPULATION_SIZE, fitness.lin_reg, genes_to_value(0, 128, return_bitstring=True),
                       crowding=True)
  ga_cw_reg.simulate(param.MOD_DRAW_FITNESS, True)
  ga_cw_reg.print_population()

  draw_graph.draw(
    [[x for x in range(len(ga_reg.get_entropies()))], ga_reg.get_entropies(), "sga_reg"],
    [[x for x in range(len(ga_cw_reg.get_entropies()))], ga_cw_reg.get_entropies(), "sga_cw_reg"]
  )
  print(f'All time best: {ga_cw_reg.all_time_best}')
