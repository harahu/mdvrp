import random, math, point, numpy, os
import matplotlib.pyplot as plt 

script_dir = os.path.dirname(os.path.realpath('__file__'))
p23 = os.path.join(script_dir, 'Data', 'Data Files', 'pr02') 

class Chromosome(object):
    """docstring for Chromosome"""
    def __init__(self, genes, C, D):
        super(Chromosome, self).__init__()
        self.genes = genes
        self.C = C
        self.D = D
        self.solution = None

    def get_solution(self):
        if self.solution == None:
            self.solution = MDVRPSolution(self)
        return self.solution

    def get_copy(self):
        genes_copy = []
        for gene in self.genes:
            genes_copy.append(gene[:])
        return Chromosome(genes_copy, self.C, self.D)

class MDVRPSolution(object):
    """Docstring"""
    def __init__(self, chromosome):
        super(MDVRPSolution, self).__init__()
        self.tours = [[] for i in range(len(chromosome.genes))]
        self.C = chromosome.C
        self.D = chromosome.D
        for i in range(len(chromosome.genes)):
            self.tours[i] = self._schedule_tours(chromosome.genes[i], self.D[i])
        self.total_distance = 0
        for depot in self.tours:
            for tour in depot:
                self.total_distance += self._tour_dist(tour)

    def reschedule_depot(self, depot_cluster, i):
        for tour in self.tours[i]:
            self.total_distance -= self._tour_dist(tour)
        self.tours[i] = self._schedule_tours(depot_cluster, self.D[i])
        for tour in self.tours[i]:
            self.total_distance += self._tour_dist(tour)

    def _schedule_tours(self, depot_cluster, d):
        routes = []
        length = 0
        load = 0
        tour = [d]
        for c in depot_cluster:
            total_duration = length + c.distance(tour[-1]) + c.duration + d.distance(c)
            if load + c.demand <= d.max_load and d.duration_check(total_duration):
                tour.append(c)
                length += c.distance(tour[-1]) + c.duration
                load += c.demand
            else:
                tour.append(d)
                routes.append(tour)
                length = 0
                load = 0
                tour = [d]
                total_duration = length + c.distance(tour[-1]) + c.duration + d.distance(c)
                assert load + c.demand <= d.max_load and d.duration_check(total_duration)
                tour.append(c)
                length += c.distance(tour[-1]) + c.duration
                load += c.demand
        tour.append(d)
        routes.append(tour)
        #build second alternative
        routes2 = [route for route in routes]
        for i in range(len(routes2)):
            routes2[i] = routes2[i][1:-1]
        for i in range(1, len(routes2)):
            routes2[i] = [routes2[i-1][-1]] + routes2[i]
            routes2[i-1] = routes2[i-1][:-1]
        for i in range(len(routes2)):
            routes2[i] = [d] + routes2[i] + [d]
        valid = True
        for route in routes2:
            if not (d.duration_check(self._tour_dist(route)) and self._tour_load(route) <= d.max_load):
                valid = False
                break
        if valid and self._routes_dist(routes) >= self._routes_dist(routes2):
            return routes2
        return routes

    def _tour_dist(self, tour):
        distance = 0
        for i in range(1, len(tour)):
            distance += tour[i].distance(tour[i-1])
        return distance

    def _routes_dist(self, tours):
        distance = 0
        for tour in tours:
            distance += self._tour_dist(tour)
        return(distance)

    def _tour_load(self, tour):
        load = 0
        for customer in tour[1:-1]:
            load += customer.demand
        return load

    def write_to_file(self):
        with open('p01_s', 'w') as f:
            f.write("{:.2f}".format(self.total_distance))
            f.write('\n')
            for i in range(len(self.tours)):
                depot = self.tours[i]
                for j in range(len(depot)):
                    tour = depot[j]
                    f.write(str(i+1)+' ')
                    f.write(str(j+1)+' ')
                    f.write("{:.2f}".format(self._tour_dist(tour))+' ')
                    f.write('0 ')
                    for pnt in tour[1:-1]:
                        f.write(str(pnt.id)+' ')
                    f.write('0')
                    f.write('\n')

    def plot(self):
        customer_x = [c.x for c in self.C]
        customer_y = [c.y for c in self.C]
        depot_x = [d.x for d in self.D]
        depot_y = [d.y for d in self.D]
        fig, ax = plt.subplots()
        ax.scatter(customer_x, customer_y, marker='x')
        for depot in self.tours:
            for tour in depot:
                xs = [point.x for point in tour]
                ys = [point.y for point in tour]
                ax.plot(xs, ys, c=numpy.random.rand(3,1))
        ax.scatter(depot_x, depot_y, marker='o', c = 'r')
        plt.show()

def get_problem_set(filename):
    with open(filename, 'r') as f:
        dataset = [line.split() for line in f]
    m, n, t = tuple([int(i) for i in dataset[0]])
    #depot list
    D = []
    for i in range(t):
        d = [i] + [float(j) for j in dataset[1+t+n+i][1:3]] + [float(j) for j in dataset[1+i]] + [m]
        d = point.Depot(*tuple(d))
        D.append(d)
    #customer list
    C = []
    for line in dataset[1+t:1+t+n]:
        c = [int(line[0])] + [float(i) for i in line[1:5]]
        c = point.Customer(*tuple(c))
        C.append(c)
    return C, D

def min_dist_randomized_chromosome(C, D):
    genes = [[] for _ in range(len(D))]
    for c in C:
        min_dist_d = 0
        min_dist = c.distance(D[0])
        for d in D[1:]:
            curr_dist = c.distance(d)
            if curr_dist < min_dist:
                min_dist_d = d.id
                min_dist = curr_dist
        genes[min_dist_d].append(c)
    for gene in genes:
        random.shuffle(gene)
    chromosome = Chromosome(genes, C, D)
    return chromosome

def initial_population(C, D, size):
    return [min_dist_randomized_chromosome(C, D) for _ in range(size)]

def tournament_selection(population, tourney_size):
    tourney = random.sample(population, 2)
    if random.random() < 0.9:
        return min(tourney, key=lambda chromosome: chromosome.get_solution().total_distance)
    return max(tourney, key=lambda chromosome: chromosome.get_solution().total_distance)

def best_cost_route_crossover(genes, depot, tour, C, D):
    for d in genes:
        for c in tour:
            if c in d:
                d.remove(c)
    stripped_solution = Chromosome(genes, C, D).get_solution()
    for c in tour:
        stripped_cost = stripped_solution.total_distance
        insertion_costs = []
        for i in range(len(genes[depot])+1):
            genes[depot].insert(i, c)
            stripped_solution.reschedule_depot(genes[depot], depot)
            insertion_costs.append(stripped_solution.total_distance - stripped_cost)
            del genes[depot][i]
        genes[depot].insert(insertion_costs.index(min(insertion_costs)), c)
        stripped_solution.reschedule_depot(genes[depot], depot)
    return genes

def reversal_mutation(gene):
    cutpoints = random.sample(range(len(gene)), 2)
    cutpoints.sort()
    gene[cutpoints[0]:cutpoints[1]] = gene[cutpoints[0]:cutpoints[1]][::-1]

def swap_mutation(gene):
    swap_points = random.sample(range(len(gene)), 2)
    gene[swap_points[0]], gene[swap_points[1]] = gene[swap_points[1]], gene[swap_points[0]]

def mutate_genes(genes):
    gene = random.choice(genes)
    if random.random() <= 0.5:
        reversal_mutation(gene)
    else:
        swap_mutation(gene)

def recombination(p1, p2, mutate):
    c1_genes = p1.get_copy().genes
    c2_genes = p2.get_copy().genes
    if random.random() <= 0.6:
        depot = random.randrange(0, len(p1.genes))
        p1_tour = random.choice(p1.get_solution().tours[depot])[1:-1]
        p2_tour = random.choice(p2.get_solution().tours[depot])[1:-1]
        c2_genes = best_cost_route_crossover(c2_genes, depot, p1_tour, p1.C, p1.D)
        c1_genes = best_cost_route_crossover(c1_genes, depot, p2_tour, p1.C, p1.D)
    if mutate:
        mutate_genes(c1_genes)
        mutate_genes(c2_genes)
    return Chromosome(c1_genes, p1.C, p1.D), Chromosome(c2_genes, p1.C, p1.D)

def genetic_algorithm(population, next_gen_size):
    new_pop = []
    #inserting elites
    for i in range(next_gen_size//100):
        new_pop.append(population[i])
    #filling in with children
    mutate = False
    if random.random() <= 0.2:
        mutate = True
    while len(new_pop) < next_gen_size:
        p1 = tournament_selection(population, 2)
        p2 = tournament_selection(population, 2)
        children = recombination(p1, p2, mutate)
        new_pop.extend(children)
    sort_pop(new_pop)
    return new_pop

def pop_info(population):
    min_fitness = population[0].get_solution().total_distance
    avg_fitness = 0
    for i in range(len(population)):
        avg_fitness += population[i].get_solution().total_distance
    avg_fitness /= len(population)
    print('Avg Fit:', avg_fitness)
    print('Min Fit:', min_fitness)

def sort_pop(population):
    population.sort(key = lambda chromosome: chromosome.get_solution().total_distance)

def main():
    C, D = get_problem_set(p23)
    population = initial_population(C, D, 100)
    sort_pop(population)
    for i in range(10):
        population = genetic_algorithm(population, 100)
        print('Gen ', i+1)
        pop_info(population)
        print('-----------------')
    population[0].get_solution().write_to_file()    
    population[0].get_solution().plot()

if __name__ == '__main__':
    main()
