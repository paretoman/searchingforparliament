#!/usr/bin/env python

from gurobipy import *

import math
import StringIO
import numpy
from pyvotecore.stv import STV
from pyvotecore.plurality_at_large import PluralityAtLarge
from pyvotecore.schulze_stv import SchulzeSTV
from openstv.ballots import Ballots
from openstv.plugins import getMethodPlugins, getReportPlugins
import json
import community
from collections import defaultdict
import networkx
import random
import seaborn
import pandas
from scipy.cluster import hierarchy
from sklearn.manifold import SpectralEmbedding,LocallyLinearEmbedding
from sklearn import manifold
import csv
from subprocess import Popen, PIPE, STDOUT
from sklearn.neighbors import NearestNeighbors
import os.path

# example of problem data
# voters = [[c1,c2] for c1 in range(10) for c2 in range(10)]
# reps = [[f1*3+1.5,f2*3+1.7] for f1 in range(3) for f2 in range(3)]
# options = {"keepsmultiplier":1,"normalizeBallots":0,"oneOverDistanceBallots":1,"exponentialBallots":0,"thresholdBallots":0,"seatsPlusOne":1,"cosineSimilarity":1,"l1Similarity":0,"jaccardSimilarity":0,"numberOfWinners":5}

def mycallback(model, where):
    if where == GRB.callback.MESSAGE:
        print >>model.__output, model.cbGet(GRB.callback.MSG_STRING),

def optimize(voters, reps, options, output=False):
    
    nWinners = options['numberOfWinners']
    
    options_str = '\n'.join(["%s - %s" % (options[i],i) for i in options])
    g_log = "" # for logging
    
    numReps = len(reps)
    numVoters = len(voters)

    if 0:
        # make a random sample of voters
        numVoters0 = numVoters
        voters0 = voters
        if numVoters > 140 and options['phragmen']:
            rv = numpy.array(random.sample(range(numVoters),140))
            rv.sort() # there is actually a good ordering to the data already. neat.
            numVoters = 140
            voters = numpy.array(voters)[rv,:].tolist()
            print(voters)
    
    # Add constants
    d = numpy.zeros((numVoters,numReps))
    b = numpy.zeros((numVoters,numReps))
    s = numpy.zeros((numReps,numReps))
    t = numpy.zeros(numReps)

    def distance(a,b):
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        return math.sqrt(dx*dx + dy*dy)
        
    for i in range(numVoters):
        for j in range(numReps):
            d[i,j] = distance(voters[i], reps[j])
    
    if options["exponentialBallots"]:
        b = numpy.exp(- d/10 )
    elif options["oneOverDistanceBallots"]:
        b = 1 /( d/10 + 1 )
    else: #if options["linearBallots"]:
        b = (numpy.max(d) - d) / 400
    
    if options['normalizeBallots']:
        for i in range(numVoters):
            normalizer = 1 / max(b[i,:])
            for j in range(numReps):
                b[i,j] *= normalizer
    
    
    #
    # We're all set up to do the computations.
    # This next section will get s and t, which are for computing just the BQP problem.  s and t are also needed for a nice election result output table.
    # 
    
    if options['seatsPlusOne']:
        keep = options['keepsmultiplier'] * sum(numpy.max(b,1)) / (1+nWinners)
        #keep = options['keepsmultiplier'] * sum(b) / (1 + 5) # might not work
    elif options['seatsPlusHalf']:
        keep = options['keepsmultiplier'] * sum(numpy.max(b,1)) / (.5+nWinners)
    else:
        keep = options['keepsmultiplier'] * sum(numpy.max(b,1)) / (nWinners)
    
    def cosine_similarity(a,b):
        return numpy.dot(a,b) / numpy.sqrt(numpy.sum(a**2) * numpy.sum(b**2))
    def l1_similarity(a,b): # well sorensen-dice coefficient
        return 2*numpy.sum(numpy.minimum(a,b)) / (numpy.sum(a) + numpy.sum(b))
    def multiply_support_not_min(a,b): # instead of finding the amount of support shared by two candidates, multiply their support ... maybe
        return numpy.dot(a,b) / (numpy.sum(a) * numpy.max(b))
    def jaccard_similarity(a,b):
        return numpy.sum(numpy.minimum(a,b)) / numpy.sum(numpy.maximum(a,b))
    def both_out_of_one(a,b):
        return numpy.sum(numpy.minimum(a,b)) / numpy.sum(a)
    def oneFromBoth(a,b):
        return numpy.sum(a) / numpy.sum(numpy.maximum(a,b))
    def simultaneous(a,b):
        return both_out_of_one(a,b)-.5*keep*numpy.sum(numpy.minimum(a,b)) / (numpy.sum(a) * numpy.sum(b))
    def integrateKeeps(a,b):
        av = numpy.sum(a)
        bv = numpy.sum(b)
        cv = numpy.sum(numpy.minimum(a,b))
        if 0:
            return cv * (1-numpy.exp(-(1/av+1/bv)*keep)) / (2 * keep)
        elif 0:
            if keep>=cv:
                return 1000000  # don't allow this to happen
            keep_star = -cv/2*numpy.log(1-keep/cv)
            return cv * (1-numpy.exp(-(1/av+1/bv)*keep_star)) / (2 * keep)
        elif 0:
            if (cv/av+cv/bv)*keep >cv:
                return 100000 # don't allow the pair if they don't have enough votes (if the venn diagram would go negative)
            else:
                return cv/av;
        elif 0:
            if (cv/av+cv/bv)*keep >cv:
                extra = (cv/av+cv/bv)*keep - cv
                return cv/av + extra/(2*keep) # double cost for borrowing to meet the keep
            else:
                return cv/av;
        elif 0:
            return .5 * cv * min(1/av+1/bv,1/keep) # no negatives in the venn diagram
        elif 0:
            return min(cv/(2*keep),1)  # separate voter groups
        elif 1:
            return cv/av * (numVoters-av-bv+cv)/numVoters  # reduce number of total voters
        else:
            jg = 1/2*cv/(2*(av+bv-cv)-numVoters) #hmm doesn't work
            kg = 1+jg
            return jg
            
    if options['l1Similarity']:
        similarity = l1_similarity
    elif options['cosineSimilarity']:
        similarity = cosine_similarity
    elif options["bothOutOfOne"]: 
        similarity = both_out_of_one
    elif options["oneFromBoth"]: 
        similarity = oneFromBoth
    elif options["multiplySupport"]: 
        similarity = multiply_support_not_min
    elif options["simultaneous"]: 
        similarity = simultaneous
    elif options["integrateKeeps"]: 
        similarity = integrateKeeps
    else: #options["jaccardSimilarity"]:
        similarity = jaccard_similarity
            
    for j in range(numReps):
        for k in range(numReps):
            s[j,k] = similarity(b[:,j],b[:,k])
        t[j] = sum(b[:,j])
    
    
    

    
    #
    # there are a few different computations, depending on the option
    #
    
    if options['phragmen'] or options['Phragmen bid'] or options['RRV max easy'] or options['RRV max'] or options["RRVloadbalance"] or options['RRVbid'] or options["RRVloadbalanceEasy"] or options['computeClustering'] or options['computeMaxRRV'] or options['computeBQP']: # uses gurobi
    
        m = Model()
        if not output:
            m.params.OutputFlag = 0
        m.setParam('TimeLimit', 100)
        
        # Add variables
        x = {}
        for j in range(numReps):
            x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)
            
        rrvfudgefactor = 1
        
        if options['phragmen']:
        
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            # huh, I didn't even need y<x because the optimization already took care of that.
            #m.addConstr(nWinners == quicksum(quicksum( y[(i,j)] for j in range(numReps) ) for i in range(numVoters) ) )
            #m.addConstr(nWinners == quicksum(quicksum( b[i,j] * y[(i,j)] for j in range(numReps) ) for i in range(numVoters) ) )
            for j in range(numReps):
                m.addConstr( x[j] == quicksum(b[i,j]*y[(i,j)] for i in range(numVoters)) ) 
            for i in range(numVoters):
                m.addConstr( Z >= quicksum( y[(i,j)] for j in range(numReps) ) ) 
            m.update()
            m.setObjective(Z,GRB.MINIMIZE)
        
        
        elif options["RRVloadbalance"]:
            # fast? way
            f = {}
            for i in range(numVoters): # should this have an upper bound of 1?
                f[i] = m.addVar(lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="f%d" % (i))
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            for j in range(numReps):
                m.addConstr( x[j] == x[j] * (quicksum((b[i,j])*f[i] for i in range(numVoters))  ) ) # equal seats
            if 0: # why did I think I needed this?
                for i in range(numVoters): # setting f
                    m.addConstr( f[i] * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )
            for i in range(numVoters):
                m.addConstr( Z >= nWinners * f[i] ) # minimize
            m.update()
            m.setObjective(Z,GRB.MINIMIZE)

        elif options["RRVloadbalanceEasy"]:
            # easy  way
            f = {}
            for i in range(numVoters): # should this have an upper bound of 1?
                f[i] = m.addVar(lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="f%d" % (i))
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, ub=1000, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            #m.addConstr(nWinners == quicksum(quicksum( y[(i,j)] for j in range(numReps) ) for i in range(numVoters) ) )
            #m.addConstr(nWinners == quicksum(quicksum( b[i,j] * y[(i,j)] for j in range(numReps) ) for i in range(numVoters) ) )
            for j in range(numReps):
                m.addConstr( x[j] == quicksum(b[i,j]*y[(i,j)] for i in range(numVoters)) ) 
            for i in range(numVoters):
                m.addConstr( Z >= quicksum( y[(i,j)] for j in range(numReps) ) ) 
            for i in range(numVoters):
                for j in range(numReps):
                    m.addConstr( y[(i,j)] == f[i] * x[j] )
            m.update()
            m.setObjective(Z,GRB.MINIMIZE)
            
        elif options['RRVbid']:
            f = {}
            for i in range(numVoters):
                f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % (i))
                # an option might be to put ub=0.
                # runs in about 20 seconds
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            for j in range(numReps):
                m.addConstr( Z  <= 10000* (1-x[j]) + quicksum(b[i,j]*f[i] for i in range(numVoters)) ) # equal seats
            for i in range(numVoters): # setting f
                m.addConstr( f[i] * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )
                #m.addConstr( f[i] + 1/nWinners * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )
            m.update()
            m.setObjective(Z,GRB.MAXIMIZE)
            
        elif options['Phragmen bid']:
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            for i in range(numVoters):
                for j in range(numReps):
                    m.addConstr( y[(i,j)] <= x[j] ) # maybe i don't need this, and the optimization will take care of it.
            for j in range(numReps):
                m.addConstr( Z  <= numVoters*10*(1-x[j]) + quicksum(b[i,j]*y[(i,j)] for i in range(numVoters)) ) # equal seats
            for i in range(numVoters): # setting f
                m.addConstr( quicksum( y[(i,j)] for j in range(numReps) ) == 1 ) # this gives the same phragmen choices.
                #m.addConstr( quicksum( b[i,j] * y[(i,j)] for j in range(numReps) ) == 1 ) # really weird choices
                #m.addConstr( f[i] + 1/nWinners * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )  # interesting TDON method
            m.update()
            m.setObjective(Z,GRB.MAXIMIZE)
            
        elif options['RRV max easy']:
            f = {}
            for i in range(numVoters):
                f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % (i))
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            for i in range(numVoters):
                for j in range(numReps):
                    m.addConstr( y[(i,j)] == f[i] * b[i,j] * x[j] ) # maybe i don't need this, and the optimization will take care of it.
            for j in range(numReps):
                m.addConstr( Z  <= numVoters*10*(1-x[j]) + quicksum(b[i,j]*y[(i,j)] for i in range(numVoters)) ) # equal seats
            for i in range(numVoters): # setting f
                m.addConstr( quicksum( y[(i,j)] for j in range(numReps) ) == 1 ) # this gives the same phragmen choices.
                #m.addConstr( quicksum( b[i,j] * y[(i,j)] for j in range(numReps) ) == 1 ) # really weird choices
                #m.addConstr( f[i] + 1/nWinners * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )  # interesting TDON method
            m.update()
            m.setObjective(Z,GRB.MAXIMIZE)
            
        elif options['RRV max']:
            f = {}
            for i in range(numVoters):
                #ubf = 1000
                ubf = 1/(numpy.min(b[i,:])*nWinners)
                f[i] = m.addVar(lb=0, ub=ubf,vtype=GRB.CONTINUOUS, name="f%d" % (i))
                # really f shouldn't be bound here because it is bound elsewhere but gurobi needs it
                # we could say f < 1/ min b
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            for j in range(numReps):
                m.addConstr( Z  <= numVoters*10*(1-x[j]) + quicksum(b[i,j]*b[i,j]*f[i] for i in range(numVoters)) ) # equal seats
            for i in range(numVoters): # setting f
                m.addConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) <= 1 ) # this gives the same phragmen choices.
                #m.addConstr( quicksum( b[i,j] * y[(i,j)] for j in range(numReps) ) == 1 ) # really weird choices
                #m.addConstr( f[i] + 1/nWinners * quicksum( b[i,j] * x[j] for j in range(numReps) ) == 1 )  # interesting TDON method
            m.update()
            m.setObjective(Z,GRB.MAXIMIZE)
            
            
        elif options['computeClustering']:
                
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            m.update()

            # Add constraints
            for i in range(numVoters):
                for j in range(numReps):
                    m.addConstr(y[(i,j)] <= x[j])

            if 1:
                for i in range(numVoters):
                    m.addConstr(quicksum(y[(i,j)] for j in range(numReps)) == 1)
            else:
                for i in range(numVoters):
                    m.addConstr(quicksum(y[(i,j)]*b[i,j] for j in range(numReps)) == 1)
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
                
            if "min total distance" == options['loadType']:
                m.setObjective( quicksum( quicksum(d[(i,j)]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ), GRB.MINIMIZE)
            elif "max total ballot" == options['loadType']:
                m.setObjective( quicksum( quicksum(b[(i,j)]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ), GRB.MAXIMIZE)
            else: 
                Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
                m.update()
                # Add constraints
                for j in range(numReps):
                
                    # oops had i and j in the sum below
                    if 0:
                        m.addConstr( Z >= quicksum( quicksum(d[(i,j)]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ) ) # the winning candidates's distances are bound from above by Z
                    elif 0:
                        m.addConstr( Z >= quicksum( quicksum(1/b[i,j]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ) ) 
                        # this actually worked.  London had 2 reps.
                    elif 0:
                        m.addConstr( Z >= quicksum( quicksum((1-b[i,j])*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ) ) 
                        # this worked too. London gets 2 reps.
                        
                    # fix
                    elif 0:
                        m.addConstr( Z >= quicksum(d[(i,j)]*y[(i,j)] for i in range(numVoters))  ) # eh, London only gets 2 when it has 3/5.
                    elif 0:
                        m.addConstr( Z >= quicksum(1/b[i,j]*y[(i,j)] for i in range(numVoters))  ) 
                        # eh, same
                    elif "minimax 1-b" == options['loadType']:
                        m.addConstr( Z >= quicksum((1-b[i,j])*y[(i,j)] for i in range(numVoters))  ) 
                        # oh cool, london actually gets 4, and this is unnormalized ballots.  This is really good.
                    
                    # fix again because ... forgot
                
                # both seem to work
                if "max min-total-ballot (jefferson)" == options['loadType']: # jefferson
                    m.setObjective( Z, GRB.MINIMIZE)
                else: # "min diff-total-ballot (webster)" == options['loadType']:
                    Y = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="y")
                    m.update()
                    for j in range(numReps):
                        m.addConstr(  Y>= quicksum((b[i,j])*y[(i,j)] for i in range(numVoters))  ) # 
                    m.setObjective( Z+Y, GRB.MINIMIZE)
                
                
        
        elif options['computeBQP']:
            
            m.update()
            
            # Add constraints
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
                
            d_obj = LinExpr()
            for j in range(numReps):
                d_obj += t[j]*x[j]
                for k in range(numReps):
                    if k != j:
                        d_obj += -.5*keep*s[j,k]*x[j]*x[k]

            m.setObjective( d_obj , GRB.MAXIMIZE)
            


        elif 0: # piggyback on computeMaxRRV when it is selected
        
        
            # this is actually an iterative process and we only stop when we reach a particular condition.  Let's just do one iteration
            
            # initialization
            
            xx = numpy.zeros(numVoters)
            xx[0] = nWinners
            
            
            #
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
                    m.addConstr( y[(i,j)] <= 1) # constraint 1
                    m.addConstr(y[(i,j)] <= x[j]) # y is zero if the rep is not selected. constraint 2
            
            m.addConstr(quicksum( quicksum(d[(i,j)]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ) == nWinners) # constraint 3
            
            for i in range(numVoters):
                m.addConstr(quicksum(y[(i,j)] for j in range(numReps)) == 1) # constraint 4 but with additional constraint that doesn't change things.
                
            # add variables
            e = {}
            s = {}
            t = {}
            for i in range(numVoters):
                for j in range(numVoters):
                    e[(i,j)] = m.addVar(vtype=GRB.BINARY, name="t%d,%d" % (i,j)) # constraint 7
                s[i] = m.addVar(vtype=GRB.BINARY, name="s%d" % i) # constraint 8
                t[i] = m.addVar(vtype=GRB.BINARY, name="t%d" % i) # constraint 9
            e_epsilon = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="e_epsilon")
            
            # constrain
            for i in range(numVoters):
                m.addConstr( s[i] + quicksum( e[(i,j)] for j in range(numVoters) ) == 1 ) # constraint 10
            for j in range(numVoters):
                m.addConstr( t[j] + quicksum( e[(i,j)] for i in range(numVoters) ) <= 1 ) # constraint 11
            m.addConstr( quicksum( t[j] for j in range(numVoters) ) == 1 ) # constraint 12
            
            for i in range(numVoters):
                for j in range(numVoters):
                    m.addConstr( quicksum( y[(i,cc)] for cc in range(numReps) ) - nWinners * ( 1 - e[(i,j)] ) <= xx[j] ) # constraint 13
                    m.addConstr( quicksum( y[(i,cc)] for cc in range(numReps) ) - nWinners * ( 2 - s[i] - t[j] ) <= xx[j] - e_epsilon ) #  constraint 14

            # ok now maximize e_epsilon.
            m.setObjective( e_epsilon, GRB.MAXIMIZE)   
                 

            m.update()
            
        elif 0:
            # apportionment method
            m = Model("qcp")
            if not output:
                m.params.OutputFlag = 0
            m.setParam('TimeLimit', 100)
            
            # Add variables
            x = {}
            for j in range(numReps):
                x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)
            f = {}
            for i in range(numVoters):
                f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % i)
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            # Add constraints
            for j in range(numReps):
                m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters) ) >= Z * x[j] ) # the winning candidates's scores are bound from below by Z
            for i in range(numVoters):
                m.addQConstr( quicksum( f[i] * (b[i,j] * x[j]) for j in range(numReps) ) <= 1 )
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            m.update()
            if 1: # jefferson
                m.setObjective( Z, GRB.MAXIMIZE)
            else: # webster
                Y = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="y")
                m.update()
                for j in range(numReps):
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters) ) <= Y * x[j] ) # the winning candidates's scores are bound from above by Y
                m.setObjective( Y-Z, GRB.MINIMIZE)
            
            m.update()
            # Gurobi can't do this because, as it says, "the q matrix is not positive semidefinite".
        
        elif 0:
            # 1/b attampt at apportionment
            m = Model("qcp")
            if not output:
                m.params.OutputFlag = 0
            m.setParam('TimeLimit', 100)
            
            # Add variables
            x = {}
            for j in range(numReps):
                x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)
            f = {}
            for i in range(numVoters):
                f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % i)
            Z = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="z")
            m.update()
            # Add constraints
            for j in range(numReps):
                if 0:
                    m.addQConstr( quicksum( f[i] * 1/b[i,j] * x[j] for i in range(numVoters) ) <= Z  ) # the winning candidates's scores are bound from below by Z
                elif 0:
                    m.addQConstr( quicksum( f[i] * (1-b[i,j]) * x[j] for i in range(numVoters) ) <= Z  ) # the winning candidates's scores are bound from below by Z
                else:
                    m.addQConstr( quicksum( f[i] * (d[i,j]) * x[j]+1 for i in range(numVoters) ) <= Z  ) # the winning candidates's scores are bound from below by Z
            for i in range(numVoters):
                if 0:
                    m.addQConstr( quicksum( f[i] * (1/b[i,j] * x[j]) for j in range(numReps) ) >= 1 )
                elif 0:
                    m.addQConstr( quicksum( f[i] * ((1-b[i,j]) * x[j]) for j in range(numReps) ) >= 1 )
                else:
                    m.addQConstr( quicksum( f[i] * ((d[i,j]) * x[j]-1) for j in range(numReps) ) >= 1 )
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            m.update()
            if 1: # jefferson
                m.setObjective( Z, GRB.MINIMIZE)
            else: # webster
                Y = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="y")
                m.update()
                for j in range(numReps):
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters) ) <= Y * x[j] ) # the winning candidates's scores are bound from above by Y
                m.setObjective( Y-Z, GRB.MINIMIZE)
            
            m.update()
            # Gurobi can't do this because, as it says, "the q matrix is not positive semidefinite".
        
        elif 0:  # working
            y = {}
            for i in range(numVoters):
                for j in range(numReps):
                    y[(i,j)] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            m.update()

            # Add constraints
            for i in range(numVoters):
                for j in range(numReps):
                    m.addConstr(y[(i,j)] <= x[j])

            for i in range(numVoters):
                m.addConstr(quicksum(y[(i,j)] for j in range(numReps)) == 1)
        
        elif options['computeMaxRRV']:
            
            # still to do
            
            m = Model("qcp")
            if not output:
                m.params.OutputFlag = 0
            m.setParam('TimeLimit', 100)
            

            # Add variables
            x = {}
            for j in range(numReps):
                x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)
                
            # Add variables
            f = {}
            for i in range(numVoters):
                f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % i)
                
            # Add constraints
                
            if 0:
                for i in range(numVoters):
                    m.addConstr(f[i] <= 1000) # maybe not needed ... turns out we needed it.
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) <= 1 ) 
                    # I know it should be == 1 but gurobi won't allow it.
                    # still gurobi doesn't work with <= 1 because it says the problem is not positive semidefinite.
                    # maybe the <=1 is okay because we are maximizing, and there is no reason to make f smaller than the optimal f, where the equality holds. 
            
            elif 0:  # try this
                won1 = {}
                for i in range(numVoters):
                    m.addConstr(f[i] <= 1000) # maybe not needed ... turns out we needed it.
                    won1[i] = m.addVar(vtype=GRB.BINARY, name="won1%d" % i)
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) == won1[i] )
                    # it looks like equality is allowed if there is a constraint on f.  f <= 1 or even f <= 1000.          
            
            elif 0:
                for i in range(numVoters):
                    m.addConstr(f[i] <= 1000)
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) == 1 )
                    
            elif 0: 
                for i in range(numVoters):
                    m.addConstr(f[i] <= 1)
                    m.addQConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) <= 1 )      
            elif 0: 
                for i in range(numVoters):
                    m.addConstr(f[i] <= 1)
                    m.addQConstr( quicksum( f[i] * (b[i,j]+1) * x[j] for j in range(numReps) ) <= 1 ) # just trying out this b+1
            elif 1:  # try this
                for i in range(numVoters):
                    # hey I don't need the extra f <= 1 constraints.
                    m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] + 1) for j in range(numReps) ) <= 1 ) # this is a better place for the + 1.  Maybe we should subtract the average score.
            else:  # try this
                for i in range(numVoters):
                    # hey I don't need the extra f <= 1 constraints.
                    nw0 = (nWinners-1) / nWinners # for some reason, computations go a lot faster with this factor.
                    m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] * nw0 + 1) for j in range(numReps) ) <= 1 ) # Let's also subtract the average score to try to match RRV better.
            m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            
            
            m.setObjective( quicksum( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters)) for j in range(numReps) ), GRB.MAXIMIZE)
            
        # optimize and wrap up output and x[j]
        output = StringIO.StringIO()
        m.__output = output
        m.optimize(mycallback) 
        g_log += output.getvalue()
        if (m.status != 2):
            return ["error",g_log]       
        xo = numpy.zeros(numReps)
        for j in range(numReps):
            if (x[j].X > .5):
                xo[j] = 1


    elif options['computeSTV'] or options["computePluralityMultiwinner"] or options["computeSchulzeSTV"] or options["MeeksSTV"] or options["openstv"]:  # uses openstv or pyvotecore
        
        # convert ballot format
        bSTV = []
        bPlur = []
        b_text = ""
        for i in range(numVoters):
            sb=b[i,:].argsort()[::-1]
            ssb = ["%d" % n_name for n_name in sb]
            ssb0 = ["%d" % sb[0]]
            dssb = {"count":1,"ballot":ssb}
            dssb0 = {"count":1,"ballot":ssb0}
            bSTV.append(dssb)
            bPlur.append(dssb0)
            a_line = " ".join(ssb)
            b_text += a_line + "\n"
        b_text = b_text[:-1] # remove last newline
        
        if options['computeSTV'] or options["computePluralityMultiwinner"] or options["computeSchulzeSTV"]:
            if options['computeSTV']:
                outputSTV = STV(bSTV, required_winners=nWinners).as_dict()
            elif options["computePluralityMultiwinner"]:
                outputSTV = PluralityAtLarge(bPlur, required_winners=nWinners).as_dict()
            elif options["computeSchulzeSTV"]:
                outputSTV = SchulzeSTV(bSTV, required_winners=nWinners, ballot_notation=SchulzeSTV.BALLOT_NOTATION_GROUPING).as_dict()
            winSet = outputSTV['winners'] # set of winners
            g_log += "PyVoteCore (not Gurobi) \n\n"
            
        elif options["MeeksSTV"] or options["openstv"]:
            if options["MeeksSTV"]:
                name="MeekSTV"
            elif options["openstv"]:
                name=options["stvtype"]
            methods = getMethodPlugins("byName", exclude0=False)
            dirtyBallots = Ballots()
            dirtyBallots.loadText(b_text)
            dirtyBallots.numSeats = nWinners
            cleanBallots = dirtyBallots.getCleanBallots()
            e = methods[name](cleanBallots)
            e.runElection()
            winSet = e.getWinnerList()
            g_log += "OpenSTV (not Gurobi) \n\n"
        
        
        xo = numpy.zeros(numReps)
        for j in range(numReps):
            if "%d" % j in winSet:
                xo[j] = 1
                

    elif options['computeRRV'] or options['computeRRV-TDON']:
        
        winBool = numpy.zeros(numReps)
        winList = []
        voterWeight = numpy.ones(numVoters)
        voterWinSum = numpy.zeros(numVoters)
        voterMax = numpy.max(b,1)
        runt = t # running tally is updated
        for i in range(nWinners):
            runt = numpy.matmul(voterWeight,b)
            runt[winList]=0
            winner = numpy.argmax(runt)
            winBool[winner]=1
            winList += [winner]
            voterWinSum += b[:,winner]
            if options['computeRRV']:
                voterWeight = 1/(1+voterWinSum/voterMax)
            else: #method from TDON
                voterWeight = 1+i-voterWinSum # or 1-voterWinSum/(i+1) or 1-voterWinSum/(maximum possible voterWinSum)
        
        g_log += "RRV (not Gurobi) \n\n"
        xo = numpy.zeros(numReps)
        for j in range(numReps):
            if j in winList:
                xo[j] = 1
                
    
    
    ### dealing with output
    
    solution1 = []
    solution2 = []
    
    solutionfid = numpy.zeros(numReps)
    i=0
    for j in range(numReps):
        if (xo[j]):
            solution1.append(j)
            solutionfid[j] = i
            i+=1
    
    votercolor = numpy.zeros(numVoters)
    yo = numpy.zeros([numVoters,numReps])
    for i in range(numVoters):
        maxj = 0
        maxb = 0
        for j in solution1:
            if options['phragmen'] or options['Phragmen bid'] or options['RRV max easy'] or options["RRVloadbalanceEasy"]:
                yo_here = y[(i,j)].X
                b_close = b[i,j] * y[(i,j)].X
                if 0: # old way
                    b_close = y[(i,j)].X
            elif options["RRVloadbalance"] or options['RRVbid'] or options['RRV max']:
                yo_here = b[i,j]*f[i].X * x[j].X
                b_close = b[i,j]*b[i,j]* f[i].X * x[j].X
                #y[(i,j)].X = b_close # might not need this
                if 0: # old way
                    b_close = f[i].X * x[j].X
            else:
                yo_here = b[i,j] # not sure about this, but I don't think any method that fell into this else would use the y term... except for the clustering method.
                b_close = b[i,j]
            yo[i,j] = yo_here
            if b_close > maxb:
                maxj = j
                maxb = b_close
                # maybe I should assign colors ot the biggest y*b?
        if 0:
            if options['phragmen']:
                if maxb > numReps/numVoters:
                    solution2.append((i,maxj))
                else:
                    solution2.append((i,-1))
            else:
                solution2.append((i,maxj))
        else:
            solution2.append((i,maxj))
        votercolor[i] = solutionfid[maxj]
        
    if 0:
        # do a similarity cluster map
        bp = pandas.DataFrame(b)
        if 0:  # this is basically the easy way to do what is below, but we want the orderings
            scm = seaborn.clustermap(bp)
            scm.savefig("scm.png")
        row_linkage = hierarchy.linkage(bp)
        col_linkage = hierarchy.linkage(bp.transpose())
        row_order = hierarchy.leaves_list(row_linkage)
        col_order = hierarchy.leaves_list(col_linkage)
        scm2 = seaborn.clustermap(bp, row_linkage=row_linkage, col_linkage=col_linkage)
        scm2.savefig("scm2.png")
    
    def getspecorder(voters):
        spec=SpectralEmbedding(n_components=1,affinity='rbf',gamma=30)
        ordd = spec.fit_transform(numpy.array(voters))
        specorder = numpy.argsort(ordd[:,0])
        return specorder
    def getspecorder2(aff):
        spec=SpectralEmbedding(n_components=1,affinity='precomputed')
        ordd = spec.fit_transform(aff)
        specorder = numpy.argsort(ordd[:,0])
        return specorder
    def getspecorder3(voters):
        spec=manifold.LocallyLinearEmbedding(int(min(numVoters2*.8,20*numVoters/numVoters2)),1)
        ordd = spec.fit_transform(numpy.array(voters))
        specorder = numpy.argsort(ordd[:,0])
        return specorder
    def getspecorder4(voters):
        spec=manifold.MDS(1)
        #spec=manifold.Isomap(int(numVoters2*.5),1)
        ordd = spec.fit_transform(numpy.array(voters))
        specorder = numpy.argsort(ordd[:,0])
        return specorder
        
    def getseriated(voters):
        csvdata = voters.tolist()
        a_command = "Rscript"
        #a_command = 'K:\\install\\R\\R-3.1.0\\bin\\Rscript.exe'
        path2script = 'seriate.R'

        
        stringarg = StringIO.StringIO()
        writer = csv.writer(stringarg)
        writer.writerows(csvdata)
        onearg = stringarg.getvalue()

        p = Popen([a_command ,path2script], stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
        grep_stdout = p.communicate(input=onearg)[0]
        #print(grep_stdout.decode())

        f = StringIO.StringIO(grep_stdout)
        reader = csv.reader(f, delimiter=',')
        s_strings =  [row for row in reader]
        for i_start in range(len(s_strings)):
            if s_strings[i_start][0] == "this is the start of the data":
                break
                
        print(i_start)
        print("R says:")
        print(s_strings[0:i_start])
        #g_data = [[float(x) for x in row] for row in s_strings]
        g_ord = [int(row[0]) for row in s_strings[i_start+1:]]

        return g_ord
        
    wholetable = -.5*keep*s
    for j in range(numReps):
        wholetable[j,j] = t[j]
    solutiontable = wholetable[solution1][:,solution1]
    
    tableslog = "\n\nProblem:\n"
    for i in wholetable:
        for j in i:
            tableslog += "%6i" % j
        tableslog += '\n'
    tableslog += '\nSolution:\n'
    for i in solutiontable:
        for j in i:
            tableslog += "%6i" % j
        tableslog += '\n'
    
    sudoku = numpy.zeros([numReps,numReps])
    runsum = 0
    for j in range(numReps):
        for k in range(numReps):
            if 1:
                if j != k:
                    a_entry = -.5*keep*s[j,k] + .5*1/(numReps-1)*(t[j]+t[k])
                    sudoku[j,k] = a_entry
                    runsum += a_entry 
            else:
                if j == k:
                    sudoku[j,k] = .5*1/numReps*(t[j]+t[k])
                else:
                    sudoku[j,k] = -.5*keep*s[j,k] + .5*1/numReps*(t[j]+t[k])
    a_avg = runsum / ( numReps * (numReps-1) )
    for j in range(numReps):
        sudoku[j,j] = a_avg
    def normalize_a(a):
        ma = numpy.max(numpy.max(a))
        mi = numpy.min(numpy.min(a))
        return (a-mi)/(ma-mi)
    nsudoku = normalize_a(sudoku)
    
    # find the communities within the adjacency matrix
    
    if 1:
        G=networkx.from_numpy_matrix(s)
    else:
        G=networkx.from_numpy_matrix(nsudoku)
    # Run louvain community finding algorithm
    louvain_community_dict = community.best_partition(G)
    # Convert community assignmet dict into list of communities
    louvain_comms = defaultdict(list)
    for node_index, comm_id in louvain_community_dict.iteritems():
        louvain_comms[comm_id].append(node_index)
    louvain_comms = louvain_comms.values()
    nodes_louvain_ordered = [node for comm in louvain_comms for node in comm]
    orderedtable = wholetable[nodes_louvain_ordered][:,nodes_louvain_ordered]
    
    nordsudoku = nsudoku[nodes_louvain_ordered][:,nodes_louvain_ordered]
    ords = s[nodes_louvain_ordered][:,nodes_louvain_ordered]
    ordt = t[nodes_louvain_ordered]
    orderedxo = xo[nodes_louvain_ordered]
    orderedsolutionfid = solutionfid[nodes_louvain_ordered]
    
    # add a table of numbers to solve, like a sudoku
    tableslog += "\nOrdered Table:\n"
    for i in range(numReps):
        for j in range(numReps):
            tableslog += "%6i" % orderedtable[i,j]
            if orderedxo[i] and orderedxo[j]:
                tableslog += "*"
            else:
                tableslog += " "
        tableslog += '\n'
    
    
    # see what a similarity measure would look like for voters - if there are communities
    do_setup = not os.path.isfile("vorder.txt") # see if we already have a starting file ready
    print("setup?")
    print(do_setup)
    
    if options["Calculate Voter Communities"]:
        if options["stateInitialMap"] and not do_setup:
            # don't need to do new calculations
            vorder = numpy.genfromtxt('vorder.txt', delimiter=',',dtype='int').tolist()
            rv = numpy.genfromtxt('rv.txt', delimiter=',',dtype='int').tolist()
            votercomcolor = votercolor[vorder]
            votercommunities = numpy.genfromtxt('voterCom.txt', delimiter=',')
            votercommunities = votercommunities.reshape((len(vorder),len(vorder)))
        else:
            # make a random sample of voters
            if numVoters > 140:
                if 0: #test
                    numVoters2 = numVoters
                else:    
                    numVoters2 = 140
                rv = numpy.array(random.sample(range(numVoters),numVoters2))
                rv.sort() # there is actually a good ordering to the data already. neat.
            else:
                numVoters2 = numVoters
                rv = numpy.array(range(numVoters))
            vd = numpy.zeros((numVoters2,numVoters2))
            vb = numpy.zeros((numVoters2,numVoters2))
            vs = numpy.zeros((numVoters2,numVoters2))
            for i in range(numVoters2):
                for j in range(numVoters2):
                    vd[i,j] = distance(voters[rv[i]], voters[rv[j]])
            if options["exponentialBallots"]:
                vb = numpy.exp(- vd/10 )
            elif options["oneOverDistanceBallots"]:
                vb = 1 /( vd/10 + 1 )
            else: #if options["linearBallots"]:
                vb = (numpy.max(vd) - vd) / 400
            for j in range(numVoters2):
                for k in range(numVoters2):
                    vs[j,k] = similarity(vb[:,j],vb[:,k])
            # find the communities within the adjacency matrix
            def getcommunityorder(a):
                G=networkx.from_numpy_matrix(a)
                # Run louvain community finding algorithm
                louvain_community_dict = community.best_partition(G)
                # Convert community assignmet dict into list of communities
                louvain_comms = defaultdict(list)
                for node_index, comm_id in louvain_community_dict.iteritems():
                    louvain_comms[comm_id].append(node_index)
                louvain_comms = louvain_comms.values()
                nodes_louvain_ordered = [node for comm in louvain_comms for node in comm]
                return nodes_louvain_ordered
            def gethierarchyorder(a):
                bp = pandas.DataFrame(a)
                if 0:  # this is basically the easy way to do what is below, but we want the orderings
                    scm = seaborn.clustermap(bp)
                    scm.savefig("scm.png")
                row_linkage = hierarchy.linkage(bp)
                col_linkage = hierarchy.linkage(bp.transpose())
                row_order = hierarchy.leaves_list(row_linkage)
                col_order = hierarchy.leaves_list(col_linkage)
                scm2 = seaborn.clustermap(bp, row_linkage=row_linkage, col_linkage=col_linkage)
                scm2.savefig("scm3.png")
                return row_order
            if 0:
                nvorder = getcommunityorder(vs)
            elif 0:
                nvorder = gethierarchyorder(vs)
            elif 0:
                nvorder = getspecorder(numpy.array(voters)[rv])
            elif 0:
                nvorder = getspecorder2(vs)
            elif 0:
                nvorder = getspecorder3(numpy.array(voters)[rv])
            elif 0:
                nvorder = getspecorder4(numpy.array(voters)[rv])
            else:
                nvorder = getseriated(numpy.array(voters)[rv])
            votercommunities = vs[nvorder][:,nvorder]
            vorder = rv[nvorder].tolist()
            votercomcolor = votercolor[vorder]
        
    else:
        votercommunities = numpy.zeros(2)
        votercomcolor = numpy.zeros(2)
        vorder = options['vorder']
        
    if options['findnearestneighborrorder']:
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(numpy.array(voters)[vorder])
        a_distances, a_indices = nbrs.kneighbors(reps)
        a_indices = [x[0] for x in a_indices]
        print(a_indices)
        ord_can = numpy.argsort(a_indices).tolist()
        print(ord_can)
        orb = b[vorder][:,ord_can] # b[vorder] same as b[rv][nvorder]
    else:
        orb = numpy.zeros(2)
        ord_can = [0,0]
    
    showy=0
    oryo = yo
    oryo = yo[vorder][:,ord_can]
    noryo = yo[vorder][:,ord_can] / numpy.max(yo)
    if options['phragmen'] or options['Phragmen bid'] or options['RRV max easy'] or options['RRV max'] or options["RRVloadbalance"] or options['RRVbid'] or options["RRVloadbalanceEasy"]:
        showy = 1
    
    # just once for set up.
    if do_setup:
        f = open('voterCom.txt','w')
        votercommunities.tofile(f,",","%1.2f")
        #f.write( [["%1.3f" % a for a in b ] for b in votercommunities.tolist()] )
        f.close()
        f = open('vorder.txt','w')
        st = ''
        for a in vorder:
            st += "%i," % a
        st = st[:-1]
        f.write(st)
        f.close()
        f = open('rv.txt','w')
        st = ''
        for a in vorder:
            st += "%i," % a
        st = st[:-1]
        f.write(st)
        f.close()
        # if doing setup, remember to delete the extra commas and then change the filename to erase the 1.
        
    # also, let's see if there are voters who were evenly split between two winners.
    # y > .99
    # well, this might be more complicated because we have to exchange pairs.  Also, there probably aren't many of these ties.  Unless we are going to do a threshold ballot.
        
    def norm1(x):
        return x* 1/numpy.max(x)
        
    return [solution1, 
    solution2, 
    g_log + tableslog + options_str,
    ords.tolist(),
    orderedxo.tolist(),
    ordt.tolist(),
    nordsudoku.tolist(),
    nodes_louvain_ordered,
    orderedsolutionfid.tolist(),
    keep,
    options["Calculate Voter Communities"],
    votercommunities.tolist(),
    votercomcolor.tolist(),
    votercolor.tolist(),
    vorder,
    orb.tolist(),
    ord_can,
    showy,
    noryo.tolist(),
    solutionfid.tolist(),
    xo.tolist(),
    norm1((orb*oryo)).tolist(),
    norm1(numpy.sum(orb*oryo,1)).tolist(),
    norm1(numpy.sum(oryo,1)).tolist(),
    oryo.tolist(),
    (orb*oryo).tolist(),
    d[vorder][:,ord_can].tolist()]

def handleoptimize(jsdict):
    if 'clients' in jsdict and 'facilities' in jsdict and 'charge' in jsdict:
        optionsValues = jsdict['charge']
        optionsNames =  ["numberOfWinners","keepsmultiplier","stvtype","loadType","Calculate Voter Communities","vorder","findnearestneighborrorder","stateInitialMap","normalizeBallots","oneOverDistanceBallots","linearBallots","exponentialBallots","thresholdBallots","phragmen",'Phragmen bid','RRV max',"computeBQP","computeSTV","MeeksSTV","computeRRV","computeRRV-TDON","openstv","computePluralityMultiwinner","computeSchulzeSTV","computeClustering","computeMaxRRV","RRVloadbalance","RRVloadbalanceEasy",'RRVbid','RRV max easy',"seatsPlusZero","seatsPlusHalf","seatsPlusOne","jaccardSimilarity","bothOutOfOne","oneFromBoth","simultaneous","integrateKeeps","cosineSimilarity","l1Similarity","multiplySupport"]
        options = dict(zip(optionsNames,optionsValues))
        solution = optimize(jsdict['clients'], jsdict['facilities'], options)
        return {'solution': solution }

if __name__ == '__main__':
    jsdict = json.load(sys.stdin)
    jsdict = handleoptimize(jsdict)
    print 'Content-Type: application/json\n\n'
    print json.dumps(jsdict)


