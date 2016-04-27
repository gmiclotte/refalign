#!/usr/bin/python3

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
        nodes.append(node)

class Entry(object):
        def __init__(self, entry):
                if (int(entry[1]) & 0x10):
                        self.nodeID = -int(entry[0])
                else:
                        self.nodeID = int(entry[0])
                self.pos = int(entry[3])
                self.len = 0
                self.cigar = entry[5]
                for c in cigar_parse(self.cigar):
                        if c[1] == 'I':
                                self.len += 0
                        elif c[1] == 'D':
                                self.len -= c[0]
                        elif c[1] == 'S' or c[1] == 'H':
                                self.len += 0
                        else:
                                self.len += c[0]
        def __str__(self):
                s = ''
                s += 'id: ' + str(self.nodeID) + ' '
                s += '[' + str(self.pos) + ' - ' + str(self.end()) + '] '
                s += self.cigar
                return s
        def __repr__(self):
                return str(self)
        def __len__(self):
                return self.len
        def end(self):
                return self.pos + self.len - 1
        def shifted_end(self):
                return self.end() - settings.k + 1

class SAMNode(object):
        def __init__(self, entry_, id_):
                self.entries = [entry_]
                self.id = id_
                self.next = []
                self.deleted = False
                self.next_matches = 0
                self.max_matches = 0
                self.update_matches()
        def count_matches(self):
                count = 0
                for idx in range(0, len(self.entries)):
                        entry = self.entries[idx]
                        if idx == 0:
                                ignore = 0
                        else:
                                ignore = settings.k - 1
                        for c, t in cigar_parse(entry.cigar):
                                if ignore > 0:
                                        m = min(ignore, c)
                                        ignore -= m
                                        c -= m
                                if t == 'M':
                                        count += c
                                else:
                                        count -= 0 * c
                return count
        def get_max_matches(self):
                self.update_matches()
                return self.max_matches
        def update_next_matches(self, next):
                self.next_matches = next
                self.update_matches()
        def update_matches(self):
                self.max_matches = self.count_matches() + self.next_matches
        def collapse(self, next_node):
                self.next = next_node.next
                for entry in next_node.entries:
                        self.entries.append(entry)
                self.update_next_matches(next_node.next_matches)
        def container(self, other):
                if self.entries[0].pos == other.entries[0].pos:
                        if other.entries[-1].end() < self.entries[-1].end():
                                return self.id
                        elif self.entries[-1].end() < other.entries[-1].end():
                                return other
                        else:
                                if other.get_max_matches() < self.get_max_matches():
                                        return self.id
                                elif self.get_max_matches() < other.get_max_matches():
                                        return other.id
                                else:
                                        return None
                elif self.entries[0].pos < other.entries[0].pos:
                        if other.entries[-1].end() <= self.entries[-1].end():
                                return self.id
                        else:
                                return None
                else:
                        if self.entries[-1].end() <= other.entries[-1].end():
                                return other
                        else:
                                return None

class SAMGraph(object):
        def __init__(self):
                self.nodes = []
                self.position = 0
        def __len__(self):
                return len(self.nodes)
        def add_node(self, entry, id_):
                if not id_ == len(self.nodes):
                        print('An error occured while building the SAM graph.')
                        exit()
                self.nodes.append(SAMNode(entry, id_))
        def first(self):
                return self.get(self.position)
        def get(self, idx):
                if idx < len(self):
                        return self.nodes[idx]
        def match_entry(self, entry):
                id_ = len(self.nodes)
                self.add_node(entry, id_)
                while self.position < len(self) and self.first().entries[-1].shifted_end() < entry.pos:
                        self.position += 1
                for idx in range(self.position, len(self)):
                        node = self.get(idx)
                        if node.entries[-1].shifted_end() == entry.pos and entry.nodeID in get_out_arcs(node.entries[-1].nodeID):
                                node.next.append(id_)
                                self.nodes[idx] = node
        def update_potentials(self):
                for node in reversed(self.nodes):
                        next = 0
                        for nxt in node.next:
                                self.nodes[nxt].update_matches()
                                count = self.nodes[nxt].get_max_matches()
                                next = count if count > next else next
                        node.update_next_matches(next)
        def collapse_linear(self):
                for node in reversed(self.nodes):
                        if not node.deleted:
                                while len(node.next) == 1:
                                        #print('collapsing ' + str(node.id) + ' ' + str(node.next[0]))
                                        next = node.next[0]
                                        node.collapse(self.nodes[next])
                                        self.nodes[next].deleted = True
        def collapse_greedy(self):
                self.update_potentials()
                for node in reversed(self.nodes):
                        node.next = [n.id for n in reversed(sorted([self.nodes[next] for next in node.next], key=SAMNode.get_max_matches))]
                        if len(node.next) > 0:
                                node.next = [node.next[0]]
        def filter(self):
                for node in self.nodes:
                        if not node.deleted:
                                for other in self.nodes:
                                        if not other.deleted:
                                                container = node.container(other)
                                                if container == node.id:
                                                        other.deleted = True
                                                elif container == other.id:
                                                        node.deleted = True
                                                        break
        def __str__(self):
                s = 'Contig results:\n'
                for n in self.nodes:
                        if not n.deleted:
                                s += str(n.id) + ' ' + str(n.get_max_matches()) + ' ' + str(n.entries[0]) + ' ' + str(n.entries[-1]) + ' ' + str(len(n.entries)) + '\n'
                return s

sg = []
current_reference = ''
for t, l in sam_parse(settings.samName):
        if t == 'meta':
                print('meta: ' + l)
        else:
                entry = Entry(l)
                if l[2] != current_reference:
                        current_reference = l[2]
                        sg.append(SAMGraph())
                sg[-1].match_entry(entry)

for g in sg:
        g.collapse_greedy()
        g.collapse_linear()
        g.filter()
        print(g)

exit()
