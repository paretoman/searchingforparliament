#!/usr/bin/env python

from gurobipy import *

import math
import StringIO
import numpy

# example of problem data
# voters = [[c1,c2] for c1 in range(10) for c2 in range(10)]
# reps = [[f1*3+1.5,f2*3+1.7] for f1 in range(3) for f2 in range(3)]
# options = {"keepsmultiplier":1,"normalizeBallots":0,"oneOverDistanceBallots":1,"exponentialBallots":0,"thresholdBallots":0,"seatsPlusOne":1,"cosineSimilarity":1,"l1Similarity":0,"jaccardSimilarity":0,"numberOfWinners":5}

def mycallback(model, where):
    if where == GRB.callback.MESSAGE:
        print >>model.__output, model.cbGet(GRB.callback.MSG_STRING),

def optimize(voters, reps, options, output=False):

    nWinners = options['numberOfWinners']
    m = Model()
    if not output:
        m.params.OutputFlag = 0
    m.setParam('TimeLimit', 100)
    
    numReps = len(reps)
    numVoters = len(voters)

    # Add variables
    x = {}
    for j in range(numReps):
        x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)

    
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
                b[i,j] = numpy.exp(- d[i,j]/10 )
            else:
                b[i,j] = 1 /( d[i,j]/10 + 1 )
    
    if options['normalizeBallots']:
        for i in range(numVoters):
            normalizer = 1 / max(b[i,:])
            for j in range(numReps):
                b[i,j] *= normalizer

    if options['computeClustering']:
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
    
        m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            
        m.setObjective( quicksum( quicksum(d[(i,j)]*y[(i,j)] for i in range(numVoters)) for j in range(numReps) ), GRB.MINIMIZE)


    elif options['computeSTV']:
        
        from pyvotecore.stv import STV

        bSTV = []
        for i in range(numVoters):
            sb=b[i,:].argsort()[::-1]
            ssb = ["%d" % n_name for n_name in sb]
            dssb = {"count":1,"ballot":ssb}
            bSTV.append(dssb)
        print(bSTV)
        outputSTV = STV(bSTV, required_winners=nWinners).as_dict()
        winSet = outputSTV['winners'] # set of winners
        
                
        solution1 = []
        solution2 = []
        
        for j in range(numReps):
            if "%d" % j in winSet:
                solution1.append(j)
        
        for i in range(numVoters):
            maxj = 0
            maxb = 0
            for j in solution1:
                if b[i,j] > maxb:
                    maxj = j
                    maxb = b[i,j]
            solution2.append((i,maxj))

        return [solution1, solution2, "STV"] 

    elif options['computeMaxRRV']:
        
        # still to do
        
        # Add variables
        f = {}
        for i in range(numVoters):
            f[i] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="f%d" % i)
            #m.addConstr(f[i] <= 1) # maybe not needed
        
        # Add constraints
        m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            
        for i in range(numVoters):
            m.addConstr( quicksum( f[i] * b[i,j] * x[j] for j in range(numReps) ) == 1 ) 
            # I know it should be == 1 but gurobi won't allow it.
            # still doesn't work with <= 1
        
        m.setObjective( quicksum( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters))
                                 for j in range(numReps) ), GRB.MAXIMIZE)
        

    
    else : # computeBQP
        def cosine_similarity(a,b):
            return numpy.dot(a,b) / numpy.sqrt(numpy.sum(a**2) * numpy.sum(b**2))
        def l1_similarity(a,b):
            return numpy.dot(a,b) / (numpy.sum(a) * numpy.sum(b))

        for j in range(numReps):
            for k in range(numReps):
                if options['l1Similarity']:
                    s[j,k] = l1_similarity(a,b)
                else:
                    s[j,k] = cosine_similarity(b[:,j],b[:,k])
            t[j] = sum(b[:,j])

        if options['seatsPlusOne']:
            keep = options['keepsmultiplier'] * sum(numpy.max(b,1)) / (1+nWinners)
            #keep = options['keepsmultiplier'] * sum(b) / (1 + 5) # might not work
        else:
            keep = options['keepsmultiplier'] * sum(numpy.max(b,1)) / (nWinners)
        m.update()
        
        # Add constraints
        m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
            
        d_obj = LinExpr()
        for j in range(numReps):
            d_obj += t[j]*x[j]
            for k in range(numReps):
                if k != j:
                    d_obj += -keep*s[j,k]*x[j]*x[k]

        m.setObjective( d_obj , GRB.MAXIMIZE)
        
    output = StringIO.StringIO()
    m.__output = output

    m.optimize(mycallback)

    if (m.status != 2):
        return ["error"]

    solution1 = []
    solution2 = []

    for j in range(numReps):
        if (x[j].X > .5):
            solution1.append(j)
    
    for i in range(numVoters):
        maxj = 0
        maxb = 0
        for j in solution1:
            if b[i,j] > maxb:
                maxj = j
                maxb = b[i,j]
        solution2.append((i,maxj))

    return [solution1, solution2, output.getvalue()]

def handleoptimize(jsdict):
    if 'clients' in jsdict and 'facilities' in jsdict and 'charge' in jsdict:
        optionsValues = jsdict['charge']
        optionsNames =  ["keepsmultiplier","normalizeBallots","oneOverDistanceBallots","exponentialBallots","thresholdBallots","seatsPlusOne","cosineSimilarity","l1Similarity","jaccardSimilarity","numberOfWinners","computeBQP","computeSTV","computeClustering","computeMaxRRV"]
        options = dict(zip(optionsNames,optionsValues))
        solution = optimize(jsdict['clients'], jsdict['facilities'], options)
        return {'solution': solution }

if __name__ == '__main__':
    import json
    jsdict = json.load(sys.stdin)
    jsdict = handleoptimize(jsdict)
    print 'Content-Type: application/json\n\n'
    print json.dumps(jsdict)


