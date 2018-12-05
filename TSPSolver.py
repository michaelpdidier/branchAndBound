#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
from math import inf, isinf
from statistics import mean

class TSPSolver:
    def __init__( self, gui_view ):
        self._scenario = None

    def setupWithScenario( self, scenario ):
        self._scenario = scenario


    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution, 
        time spent to find solution, number of permutations tried during search, the 
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def defaultRandomTour( self, time_allowance=60.0 ):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time()-start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation( ncities )
            route = []
            # Now build the route using the random permutation
            for i in range( ncities ):
                route.append( cities[ perm[i] ] )
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results


    ''' <summary>
        This is the entry point for the greedy solver, which you must implement for 
        the group project (but it is probably a good idea to just do it for the branch-and
        bound project as a way to get your feet wet).  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found, the best
        solution found, and three null values for fields not used for this 
        algorithm</returns> 
    '''

    def greedy( self,time_allowance=60.0 ):
        pass


    ''' 
        <summary>This is the function used to reduce rows and columns in a state matrix.</summary>
        <returns>Total cost of reduction.</returns> 
    '''

    def reduceMatrix(self, state, citiesLen):
        # Save reduction cost
        totalCost = 0

        # Perform row reduction
        for departingCity in range(citiesLen):
            minCost = min(state[departingCity, :])
            if isinf(minCost):
                minCost = 0
            totalCost += minCost
            state[departingCity, :] -= minCost

        # Perform column reduction
        for arrivingCity in range(citiesLen):
            minCost = min(state[:, arrivingCity])
            if isinf(minCost):
                minCost = 0
            totalCost += minCost
            state[:, arrivingCity] -= minCost

        return totalCost


    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints: 
        max queue size, total number of states created, and number of pruned states.</returns> 
    '''

    def branchAndBound( self, time_allowance=60.0 ):

        # Find initial BSSF by calling the random algorithm and setting that as BSSF
        bssf = self.defaultRandomTour()['soln']

        # Create variables to track max size of priority queue, # of bssf updates, total states created, total states pruned, all initialized to 0
        maxHeapSize = bssfUpdateCount = childrenCreated = childrenKilled = 0

        # Create set of cities to be used in determining potential children later and list of cities to create state matrix
        citiesList = self._scenario.getCities()
        citiesSet = set(self._scenario.getCities())
        citiesLen = len(citiesList)
        maxDepth = citiesLen

        # Start timer
        startTime = time.time()

        # Create first state matrix
        state = np.zeros((citiesLen, citiesLen))

        for departingCity in range(citiesLen):
            for arrivingCity in range(citiesLen):
                state[departingCity][arrivingCity] = citiesList[departingCity].costTo(citiesList[arrivingCity])

        # Reduce initial matrix
        nodeCost = self.reduceMatrix(state, citiesLen)

        totalCost = totalCount = 0

        # Find average cost
        for row in state:
            for column in row:
                if not isinf(column):
                    totalCost += column
                    totalCount += 1

        averageCost = totalCount/totalCount

        # Choose first city to be root node and put matrix into node along with cost. Set parents as empty list. Set current node = this node
        currentNode = stateNode(0, nodeCost, citiesList[0], state, 0, [], set())

        # Create priority queue to hold child nodes for potential exploration
        childHeap = []

        # While loop terminates if 60 seconds have passed or if the BSSF has been found
        while time.time()-startTime < time_allowance and currentNode is not None:
            state = currentNode.state
            currentNodeIndex = currentNode.city._index

            # Get set of unvisited cities
            childSet = citiesSet - currentNode.ancestorSet

            for child in childSet:
                childIndex = child._index  # ['city']._index
                interCityCost = state[currentNodeIndex][childIndex]

                # check if route is possible
                if not isinf(interCityCost):
                    # Create cost variable, aggregate the following costs into it: the current node's cost, the cost of travelling from current city to child city
                    nodeCost = currentNode.nodeCost + interCityCost

                    # Create parent list and set for child node from current node and add current node into it
                    childAncestorList = currentNode.ancestorList.copy()
                    childAncestorList.append(currentNode.city)
                    childAncestorSet = currentNode.ancestorSet.copy()
                    childAncestorSet.add(currentNode.city)

                    # Record depth of child node for later use in childHeap
                    childDepth = currentNode.depth + 1

                    # Create the state matrix by copying the current node's state matrix
                    childState = state.copy()

                    # Iterate over parent set of current node and set entries travelling from child to those cities to inf
                    for ancestor in childAncestorList:
                        childState[childIndex][ancestor._index] = inf

                    # Set row corresponding to current city and column corresponding to child city to inf
                    childState[currentNodeIndex, :] = inf
                    childState[:, childIndex] = inf

                    # Row and column reduce matrix again, tracking additional cost and add to cost variable
                    nodeCost += self.reduceMatrix(childState, citiesLen)

                    if nodeCost <= bssf.cost:
                        if len(childSet) == 1:
                            # Replace BSSF
                            bssf = TSPSolution(childAncestorList.append(child))
                            bssfUpdateCount += 1

                        else:
                            # Add child to priority queue using node cost as key in tuple
                            childNode = stateNode(nodeCost+(averageCost*(maxDepth-childDepth)), nodeCost, child, childState, childDepth, childAncestorList, childAncestorSet)

                            #Create priority key by adding the average cost * number of cities remaining to current cost. This should shift the balance between breadth and depth the search
                            heapq.heappush(childHeap, childNode)

                    else:
                        childrenKilled += 1

                    # Increment total states created count
                    childrenCreated += 1

            # Update max heap size
            if len(childHeap) > maxHeapSize:
                maxHeapSize = len(childHeap)

            # Check for next possible path
            nodeFound = False
            currentNode = None

            while not nodeFound and len(childHeap) > 0:
                candidateNode = heapq.heappop(childHeap)

                # Check for nodes with worse paths than BSSF
                if candidateNode.nodeCost < bssf.cost:
                    currentNode = candidateNode
                    nodeFound = True
                # Prune node
                else:
                    childrenKilled += 1


        endTime = time.time()
        childrenKilled += len(childHeap)

        results = {'cost': bssf.cost}
        results['time'] = endTime - startTime
        results['count'] = bssfUpdateCount
        results['soln'] = bssf
        results['max'] = maxHeapSize
        results['total'] = childrenCreated
        results['pruned'] = childrenKilled
        return results




    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found during search, the 
        best solution found.  You may use the other three field however you like.
        algorithm</returns> 
    '''

    def fancy( self,time_allowance=60.0 ):
        pass




