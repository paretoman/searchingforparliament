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

                
    # there are a few different computations, depending on the option
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

    elif 0:
    
    
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
        elif 0:  # try this
            for i in range(numVoters):
                m.addConstr(f[i] <= 1)
                m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] + 1) for j in range(numReps) ) <= 1 ) # this is a better place for the + 1.  Maybe we should subtract the average score.
        else:  # try this
            for i in range(numVoters):
                m.addConstr(f[i] <= 1)
                nw0 = (nWinners-1) / nWinners
                m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] * nw0 + 1) for j in range(numReps) ) <= 1 ) # Let's also subtract the average score to try to match RRV better.
        m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
        
        
        m.setObjective( quicksum( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters)) for j in range(numReps) ), GRB.MAXIMIZE)
        

    
    else : # computeBQP
        def cosine_similarity(a,b):
            return numpy.dot(a,b) / numpy.sqrt(numpy.sum(a**2) * numpy.sum(b**2))
        def l1_similarity(a,b):
            return numpy.dot(a,b) / (numpy.sum(a) * numpy.sum(b))

        for j in range(numReps):
            for k in range(numReps):
                if options['l1Similarity']:
                    s[j,k] = l1_similarity(b[:,j],b[:,k])
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


