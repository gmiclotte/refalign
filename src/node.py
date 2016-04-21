class Node:
        def __init__(self, id, seq):
                self.id = id
                self.inCount = 0
                self.inArcs = []
                self.outCount = 0
                self.outArcs = []
                self.seq = seq

        def __str__(self):
                s = ''
                s += 'id: ' + str(self.id) + '\n'
                s += 'in count: ' + str(self.inCount) + '\n'
                s += 'in arcs: ' + ' '.join(str(i) for i in self.inArcs) + '\n'
                s += 'out count: ' + str(self.outCount) + '\n'
                s += 'out arcs: ' + ' '.join(str(i) for i in self.outArcs) + '\n'
                s += 'sequence: ' + self.seq + '\n'
                return s
