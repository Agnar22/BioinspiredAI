import param
import copy
import random
import fitness
import draw_graph
import math
import numpy as np


class Individual:
  def __init__(self, gene_sequence: str = None):
    self.fitness = 0
    if gene_sequence is not None:
      self.genes = gene_sequence
    else:
      self.genes = ''.join([str(random.randint(0, 1)) for _ in range(param.BITSTRING_SIZE)])


class SimpleGA:
  def __init__(self, population_size: int, fitness_function, genes_to_value_function):
    self.fitness_function = fitness_function
    self.genes_to_value_function = genes_to_value_function
    self.population_size = population_size
    self.population = []
    self.all_time_best = -math.inf
    for x in range(self.population_size):
      self.population.append(Individual())

  def simulate(self):
    self.test_fitness()
    for x in range(param.EPOCHS):
      selected_pairs = self.select()
      children = self.crossover(selected_pairs)
      self.population = self.mutate(children)
      self.test_fitness()
      self.kill_weakest()

  def print_population(self):
    print(map(lambda x: x.genes, self.population))

  def test_fitness(self):
    tot_fitness = 0
    max_fitness = -math.inf
    for indiv in self.population:
      indiv.fitness = self.fitness_function(self.genes_to_value_function(indiv.genes))
      tot_fitness += -np.log(indiv.fitness)
      max_fitness = max(np.log(indiv.fitness), max_fitness)
    print(f'Max fitness: {max_fitness}, average fitness: {tot_fitness / len(self.population)}')
    self.all_time_best = max(max_fitness, max_fitness)

  def draw_fitness_graph(self):
    x = [self.genes_to_value_function(indiv.genes) for indiv in self.population]
    y = [self.fitness_function(val) for val in x]
    draw_graph.draw_sin(np.array(x), np.array(y))

  def select(self):
    tot_fitness = 0
    for indiv in self.population:
      tot_fitness += indiv.fitness

    selected_pairs = []
    for x in range(self.population_size // 2):
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

  def mutate(self, children):
    for indiv in children:
      for pos in range(param.BITSTRING_SIZE):
        if random.random() < param.MUTATION_PROB:
          genes = list(indiv.genes)
          genes[pos] = '1' if indiv.genes[pos] == '0' else '0'
          indiv.genes = ''.join(genes)
    return children

  def kill_weakest(self):
    for indiv in self.population:
      indiv.fitness = indiv.fitness
    self.population.sort(key=lambda x: x.fitness, reverse=True)
    self.population = self.population[:int((1-param.KILL_LOWEST) * self.population_size)]


def genes_to_value(start, stop, return_bitstring=False):
  def value_scaled(genes):
    if return_bitstring:
      return genes
    interval = stop - start
    return int(genes, 2) / ((2 ** param.BITSTRING_SIZE) - 1) * interval + start
  return value_scaled

def hyperparam_search(search_params):
  for kill_lowest, mut_prob, cross_prob, pop_size in search_params:
    param.KILL_LOWEST = kill_lowest
    param.MUTATION_PROB = mut_prob
    param.CROSSOVER_PROB = cross_prob
    param.POPULATION_SIZE = pop_size

    ga = SimpleGA(param.POPULATION_SIZE, fitness.lin_reg, genes_to_value(0, 128, return_bitstring=True))
    ga.simulate()

    with open("hyperparameters.txt", "a+") as file:
      file.write(f'{param.EPOCHS},\t{kill_lowest},\t{mut_prob},\t{cross_prob},\t{pop_size}:\t{ga.all_time_best:.5f}\n')


if __name__ == '__main__':
  # hyperparam_search(
  #   [
  #     [0.9, 0.04, 0.6, 150],
  #     [0.75, 0.02, 0.7, 150],
  #     [0.7, 0.02, 0.75, 150],
  #     [0.8, 0.02, 0.7, 150],
  #     [0.75, 0.02, 0.65, 150],
  #     [0.75, 0.03, 0.7, 150],
  #     [0.75, 0.01, 0.7, 150],
  #   ]
  # )
  ga = SimpleGA(param.POPULATION_SIZE, fitness.lin_reg, genes_to_value(0, 128, return_bitstring=True))
  ga.simulate()
