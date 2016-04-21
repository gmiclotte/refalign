#!/usr/bin/python

from util import fasta_parse, sam_parse, cigar_parse
from settings import Settings
from node import Node

settings = Settings()
settings.graphName = 'data/DBGraph.fasta'
settings.samName = 'data/sorted.sam'
settings.k = 63

nodes = []

def get_out_arcs(nodeID):
        if nodeID < 0:
                return [-a for a in nodes[-nodeID - 1].inArcs]
        else:
                return nodes[nodeID - 1].outArcs

for meta, seq in fasta_parse(settings.graphName, allmeta = True):
        node = Node(int(meta[0]), seq)
        node.size = int(meta[1])
        node.inCount = int(meta[2])
        j = 3
        while j < node.inCount + 3:
                node.inArcs.append(int(meta[j]))
                j += 1
        node.outCount = int(meta[j])
        j += 1
        while j < node.outCount + node.inCount + 4:
                node.outArcs.append(int(meta[j]))
                j += 1
        print(node)
        nodes.append(node)

class Entry:
        def __init__(self, entry):
                if (int(entry[1]) & 0x10):
                        self.nodeID = -int(entry[0])
                else:
                        self.nodeID = int(entry[0])
                self.pos = int(entry[3])
                self.len = 0
                self.cigar = entry[5]
                for c in cigar_parse(self.cigar):
                        if c[1] == 'D':
                                self.len -= c[0]
                        if c[1] == 'S' or c[1] == 'H':
                                self.len += c[0]
                        else:
                                self.len += c[0]

        def __str__(self):
                s = ''
                s += 'id: ' + str(self.nodeID) + ' '
                s += '[' + str(self.pos) + ' - ' + str(self.pos + self.len) + '] '
                s += self.cigar + '\n'
                return s

        def __repr__(self):
                return str(self)

current_paths = []
finished_paths = []
for t, l in sam_parse(settings.samName):
        if t == 'meta':
                print('meta: ' + l)
        else:
                entry = Entry(l)
                next_paths = []
                appended = False
                for cp in current_paths:
                        if cp[-1].pos + cp[-1].len - (settings.k - 1) < entry.pos:
                                finished_paths.append(cp)
                        elif cp[-1].pos + cp[-1].len - (settings.k - 1) == entry.pos \
                             and entry.nodeID in get_out_arcs(cp[-1].nodeID):
                                cp.append(entry)
                                next_paths.append(cp)
                                appended = True
                if not appended:
                        next_paths.append([entry])
                current_paths = next_paths
for cp in current_paths:
        finished_paths.append(cp)

for p in finished_paths:
        if 1 or len(p) > 1:
                print(p)



