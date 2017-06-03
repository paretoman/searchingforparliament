#!/usr/bin/env python

from gurobipy import *
import math
import StringIO

clients = [[100,200], [150,250], [650, 200], [50, 300]]

facilities = []; charge = []

for i in range(10):
    for j in range(10):
        facilities.append([i*70, j*50])
        charge.append(0)  # here I made an edit so that charge is 0 instead of 1
        # there is no additional cost for having a facility.  Though this doesn't really make any difference.

def mycallback(model, where):
    if where == GRB.callback.MESSAGE:
        print >>model.__output, model.cbGet(GRB.callback.MSG_STRING),

def distance(a,b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return math.sqrt(dx*dx + dy*dy)

def optimize(clients, facilities, charge, output=False):
    numFacilities = len(facilities)

    numClients = len(clients)

    m = Model()

    if not output:
        m.params.OutputFlag = 0

    #m.setParam('TimeLimit', 10000)

    # Add variables
    x = {}
    y = {}
    d = {} # Distance matrix (not a variable)

    for j in range(numFacilities):
        x[j] = m.addVar(vtype=GRB.BINARY, name="x%d" % j)
        # so there is a variable x1, x2, etc for each possible facility
        # I want to make these sum to a set number, like 5

    for i in range(numClients):
        for j in range(numFacilities):
            y[(i,j)] = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name="t%d,%d" % (i,j))
            d[(i,j)] = distance(clients[i], facilities[j])
        # there are also y1,1 y1,2 etc for each pair of client and facility
        # there is also distance between each pair
    m.update()

    # Add constraints
    for i in range(numClients):
        for j in range(numFacilities):
            m.addConstr(y[(i,j)] <= x[j])
            # oh, the two facilities can supply the same client.
            # the facility has to be on.
            # the amount of supply is a fraction if the facility is on.
            
    for i in range(numClients):
        m.addConstr(quicksum(y[(i,j)] for j in range(numFacilities)) == 1)
        # each client is served fully, ==1
        
    for j in range(numFacilities):
        m.addConstr(quicksum(y[(i,j)] for i in range(numClients)) == .2*numClients)
        # make each facility equally supported by the same number of clients.
    
    # I will add my own constraint
    m.addConstr(quicksum(x[j] for j in range(numFacilities)) == 5)
    # There are 5 facilities total.
        
    m.setObjective( quicksum( charge[j]*x[j] + quicksum(d[(i,j)]*y[(i,j)] for i in range(numClients))
                             for j in range(numFacilities) ), GRB.MINIMIZE)
    # all done, now we know how many facilities we're going to have and we're going to minimize the distance between the facilities and the customers.  And a facility can partially cover a customer.  Although, I guess you wouldn't need to.  Oh wait, does a facility have a capacity?
    
    output = StringIO.StringIO()
    m.__output = output

    m.optimize(mycallback)

    if (m.status != 2):
        return ["error"]

    solution1 = []
    solution2 = []

    for j in range(numFacilities):
        if (x[j].X > .5):
            solution1.append(j)

    for i in range(numClients):
        for j in range(numFacilities):
            if (y[(i,j)].X > .5):
                solution2.append([i,j])

    return [solution1, solution2, output.getvalue()]

def handleoptimize(jsdict):
    if 'clients' in jsdict and 'facilities' in jsdict and 'charge' in jsdict:
        solution = optimize(jsdict['clients'], jsdict['facilities'], jsdict['charge'])
        return {'solution': solution }

if __name__ == '__main__':
    import json
    jsdict = json.load(sys.stdin)
    jsdict = handleoptimize(jsdict)
    print 'Content-Type: application/json\n\n'
    print json.dumps(jsdict)

