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
    options_str = '\n'.join(["%s - %s" % (options[i],i) for i in options])
    
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


    elif options['computeSTV'] or options["computePluralityMultiwinner"] or options["computeSchulzeSTV"] or options["MeeksSTV"] or options["openstv"]:
        
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
        
        if options['computeSTV']:
            outputSTV = STV(bSTV, required_winners=nWinners).as_dict()
            winSet = outputSTV['winners'] # set of winners
        elif options["computePluralityMultiwinner"]:
            outputSTV = PluralityAtLarge(bPlur, required_winners=nWinners).as_dict()
            winSet = outputSTV['winners'] # set of winners
        elif options["computeSchulzeSTV"]:
            outputSTV = SchulzeSTV(bSTV, required_winners=nWinners, ballot_notation=SchulzeSTV.BALLOT_NOTATION_GROUPING).as_dict()
            winSet = outputSTV['winners'] # set of winners
        elif options["MeeksSTV"]:
            methods = getMethodPlugins("byName", exclude0=False)
            name="MeekSTV"
            dirtyBallots = Ballots()
            dirtyBallots.loadText(b_text)
            dirtyBallots.numSeats = nWinners
            cleanBallots = dirtyBallots.getCleanBallots()
            e = methods[name](cleanBallots)
            e.runElection()
            winSet = e.getWinnerList()
        elif options["openstv"]:
            methods = getMethodPlugins("byName", exclude0=False)
            name=options["stvtype"]
            dirtyBallots = Ballots()
            dirtyBallots.loadText(b_text)
            dirtyBallots.numSeats = nWinners
            cleanBallots = dirtyBallots.getCleanBallots()
            e = methods[name](cleanBallots)
            e.runElection()
            winSet = e.getWinnerList()
            
        
                
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

        return [solution1, solution2, "STV" + options_str] 

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
        
    elif options['computeRRV']:
        
        winBool = numpy.zeros(numReps)
        winList = []
        voterWeight = numpy.ones(numVoters)
        voterWinSum = numpy.zeros(numVoters)
        voterMax = numpy.max(b,1)
        for i in range(nWinners):
            t = numpy.matmul(voterWeight,b)
            t[winList]=0
            winner = numpy.argmax(t)
            winBool[winner]=1
            winList += [winner]
            voterWinSum += b[:,winner]
            voterWeight = 1/(1+voterWinSum/voterMax)
        
                
        solution1 = []
        solution2 = []
        
        for j in range(numReps):
            if j in winList:
                solution1.append(j)
        
        for i in range(numVoters):
            maxj = 0
            maxb = 0
            for j in solution1:
                if b[i,j] > maxb:
                    maxj = j
                    maxb = b[i,j]
            solution2.append((i,maxj))

        return [solution1, solution2, "STV" + options_str] 
    
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
                # hey I don't need the extra f <= 1 constraints.
                m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] + 1) for j in range(numReps) ) <= 1 ) # this is a better place for the + 1.  Maybe we should subtract the average score.
        else:  # try this
            for i in range(numVoters):
                # hey I don't need the extra f <= 1 constraints.
                nw0 = (nWinners-1) / nWinners # for some reason, computations go a lot faster with this factor.
                m.addQConstr( quicksum( f[i] * (b[i,j] * x[j] * nw0 + 1) for j in range(numReps) ) <= 1 ) # Let's also subtract the average score to try to match RRV better.
        m.addConstr(quicksum(x[j] for j in range(numReps)) == nWinners)
        
        
        m.setObjective( quicksum( quicksum( f[i] * b[i,j] * x[j] for i in range(numVoters)) for j in range(numReps) ), GRB.MAXIMIZE)
        

    
    else : # computeBQP
        
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
        
    output = StringIO.StringIO()
    m.__output = output

    m.optimize(mycallback)

    if (m.status != 2):
        return ["error"]

        
    solution1 = []
    solution2 = []

    solutionb = numpy.zeros(numReps)
    solutionfid = numpy.zeros(numReps)
    i=0
    for j in range(numReps):
        if (x[j].X > .5):
            solution1.append(j)
            solutionb[j] = 1
            solutionfid[j] = i
            i+=1
    
    wholetable = -.5*keep*s
    for j in range(numReps):
        wholetable[j,j] = t[j]
    solutiontable = wholetable[solution1][:,solution1]
    
    out3 = "\n\nProblem:\n"
    for i in wholetable:
        for j in i:
            out3 += "%6i" % j
        out3 += '\n'
    out3 += '\nSolution:\n'
    for i in solutiontable:
        for j in i:
            out3 += "%6i" % j
        out3 += '\n'
    
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
    orderedsolutionb = solutionb[nodes_louvain_ordered]
    orderedsolutionfid = solutionfid[nodes_louvain_ordered]
    
    # add a table of numbers to solve, like a sudoku
    out3 += "\nOrdered Table:\n"
    for i in range(numReps):
        for j in range(numReps):
            out3 += "%6i" % orderedtable[i,j]
            if orderedsolutionb[i] and orderedsolutionb[j]:
                out3 += "*"
            else:
                out3 += " "
        out3 += '\n'
        
    options_str += out3
    
    votercolor = numpy.zeros(numVoters)
    for i in range(numVoters):
        maxj = 0
        maxb = 0
        for j in solution1:
            if b[i,j] > maxb:
                maxj = j
                maxb = b[i,j]
        solution2.append((i,maxj))
        votercolor[i] = solutionfid[maxj]
        
    # see what a similarity measure would look like for voters - if there are communities
    if options["Calculate Voter Communities"]:
        vd = numpy.zeros((numVoters,numVoters))
        vb = numpy.zeros((numVoters,numVoters))
        vs = numpy.zeros((numVoters,numVoters))
        for i in range(numVoters):
            for j in range(numVoters):
                vd[i,j] = distance(voters[i], voters[j])
        if options["exponentialBallots"]:
            vb = numpy.exp(- vd/10 )
        elif options["oneOverDistanceBallots"]:
            vb = 1 /( vd/10 + 1 )
        else: #if options["linearBallots"]:
            vb = (numpy.max(vd) - vd) / 400
        for j in range(numVoters):
            for k in range(numVoters):
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
        vorder = getcommunityorder(vs)
        votercommunities = vs[vorder][:,vorder]
        votercomcolor = votercolor[vorder]
    else:
        votercommunities = numpy.zeros(2)
        votercomcolor = numpy.zeros(2)
        vorder = [0,0]
    
    # just once for set up.
    if 0:
        f = open('voterCom.txt','w')
        votercommunities.tofile(f,",","%1.2f")
        #f.write( [["%1.3f" % a for a in b ] for b in votercommunities.tolist()] )
        f.close()
        f = open('vorder.txt','w')
        f.write("[")
        for a in vorder:
            f.write( "%i," % a )
        f.write("]")
        f.close()
        
    return [solution1, solution2, output.getvalue() + options_str,ords.tolist(),orderedsolutionb.tolist(),ordt.tolist(),nordsudoku.tolist(),nodes_louvain_ordered,orderedsolutionfid.tolist(),keep,options["Calculate Voter Communities"],votercommunities.tolist(),votercomcolor.tolist(),votercolor.tolist(),vorder]

def handleoptimize(jsdict):
    if 'clients' in jsdict and 'facilities' in jsdict and 'charge' in jsdict:
        optionsValues = jsdict['charge']
        optionsNames =  ["numberOfWinners","keepsmultiplier","stvtype","Calculate Voter Communities","seatsPlusZero","seatsPlusHalf","seatsPlusOne","normalizeBallots","oneOverDistanceBallots","linearBallots","exponentialBallots","thresholdBallots","jaccardSimilarity","bothOutOfOne","oneFromBoth","simultaneous","integrateKeeps","cosineSimilarity","l1Similarity","multiplySupport","computeBQP","computeSTV","MeeksSTV","computeRRV","computePluralityMultiwinner","computeSchulzeSTV","openstv","computeClustering","computeMaxRRV"]
        options = dict(zip(optionsNames,optionsValues))
        solution = optimize(jsdict['clients'], jsdict['facilities'], options)
        return {'solution': solution }

if __name__ == '__main__':
    jsdict = json.load(sys.stdin)
    jsdict = handleoptimize(jsdict)
    print 'Content-Type: application/json\n\n'
    print json.dumps(jsdict)


