#!/usr/bin/python3
import sys
import getopt
from util import fasta_parse, sam_parse, cigar_parse, md_parse
from settings import Settings
from node import Node

def m2eqx(cigar, md):
        md = md.split(':')[-1]
        mds = [m for m in md_parse(md)]
        cigars = [c for c in cigar_parse(cigar)]
        c_idx = 0
        m_idx = 0
        cigar = ''
        for c_idx in range(len(cigars)):
                c = cigars[c_idx]
                if c[1] in ['I', 'N', 'S', 'H', 'P']:
                        #masked bases or insertions don't occur in the MD tag #TODO
                        cigar += str(c[0]) + c[1]
                elif c[1] in ['D', 'X', '=']:
                        #deletion or cigar already in X/= format
                        if c[1] == 'X':
                                m_idx += 2 * c[0] - 1
                        else:
                                m_idx += 1
                        cigar += str(c[0]) + c[1]
                elif c[1] == 'M':
                        m_count = 0
                        while m_count < c[0]:
                                m = mds[m_idx]
                                if m[0] in [str(i) for i in range(10)]:
                                        #match
                                        
                                        
                                        if c[0] < m_count + int(m):
                                                cigar += str(c[0] - m_count) + '='
                                                m = str(int(m) - c[0] + m_count)
                                                mds[m_idx] = m
                                                m_count += c[0] - m_count
                                        elif c[0] == m_count + int(m):
                                                cigar += str(c[0] - m_count) + '='
                                                m_count += c[0] - m_count
                                                m_idx += 1
                                        else:
                                                cigar += m + '='
                                                m_count += int(m)
                                                m_idx += 1
                                else:
                                        #mismatch
                                        size = 1
                                        while m_idx + 2 < len(mds) and mds[m_idx + 1] == '0':
                                                size += 1
                                                m_idx += 2
                                        cigar += str(size) + 'X'
                                        m_count += size
                                        m_idx += 1
                else:
                        print('Unknown cigar entry!', c)
                        exit()
        return cigar

class DBGraph:
        def __init__(self):
                self.nodes = []
        def get_out_arcs(self, nodeID):
                if nodeID < 0:
                        return [-a for a in self.nodes[-nodeID - 1].inArcs]
                else:
                        return self.nodes[nodeID - 1].outArcs

def fix_cigar(err):
        cigar = ''
        count = 0
        curr = None
        for c, t in cigar_parse(err):
                if t == curr:
                        count += c
                else:
                        if count > 0:
                                cigar += str(count) + curr
                        curr = t
                        count = c
        if count > 0:
                cigar += str(count) + curr
        return cigar

def merge_entry_cigars(settings, entries):
        cigar = ''
        count = 0
        curr = None
        for idx in range(0, len(entries)):
                if idx == 0:
                        ignore = 0
                else:
                        ignore = settings.k - 1
                for c, t in cigar_parse(entries[idx].cigar):
                        if ignore > 0:
                                m = min(ignore, c)
                                ignore -= m
                                c -= m
                        if t == curr:
                                count += c
                        else:
                                if count > 0:
                                        cigar += str(count) + curr
                                curr = t
                                count = c
        if count > 0:
                cigar += str(count) + curr
        return fix_cigar(cigar)

class Entry(object):
        def __init__(self, settings_, entry):
                self.settings = settings_
                if (int(entry[1]) & 0x10):
                        self.nodeID = -int(entry[0])
                else:
                        self.nodeID = int(entry[0])
                self.pos = int(entry[3])
                self.len = 0
                self.cigar = m2eqx(entry[5], entry[12])
                for c in cigar_parse(self.cigar):
                        if c[1] in ['M', '=', 'X', 'D', 'N']:
                                self.len += c[0]
                        elif c[1] in ['I']:
                                self.len -= c[0]
                        elif c[1] in ['S', 'H']:
                                self.len += 0
                        else:
                                print('Unknown cigar.')
                                exit()
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
                return self.end() - self.settings.k + 1

class SAMNode(object):
        def __init__(self, settings_, entry_, id_):
                self.settings = settings_
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
                                ignore = self.settings.k - 1
                        for c, t in cigar_parse(entry.cigar):
                                if ignore > 0:
                                        m = min(ignore, c)
                                        ignore -= m
                                        c -= m
                                if t == '=' or t == 'M':
                                        count += c
                                elif t in ['H']:
                                        count += 0
                                else:
                                        count -= 100 * c
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
        def __init__(self, settings_, dbg_):
                self.settings = settings_
                self.dbg = dbg_
                self.nodes = []
                self.position = 0
        def __len__(self):
                return len(self.nodes)
        def add_node(self, entry, id_):
                if not id_ == len(self.nodes):
                        print('An error occured while building the SAM graph.')
                        exit()
                self.nodes.append(SAMNode(self.settings, entry, id_))
        def first(self):
                return self.get(self.position)
        def get(self, idx):
                if idx < len(self):
                        return self.nodes[idx]
        def match_entry(self, entry):
                id_ = len(self.nodes)
                self.add_node(entry, id_)
                while self.position < len(self) and self.first().entries[-1].shifted_end() < entry.pos - 1:
                        self.position += 1
                for idx in range(self.position, len(self)):
                        node = self.get(idx)
                        if entry.nodeID in self.dbg.get_out_arcs(node.entries[-1].nodeID):
                                if node.entries[-1].shifted_end() == entry.pos - 1:
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
        def cutoff(self, c):
                for node in self.nodes:
                        node.deleted = node.deleted or (node.entries[-1].end() - node.entries[0].pos + 1) < c
        def __str__(self):
                s = ''
                for n in self.nodes:
                        if not n.deleted:
                                qname = 'path'
                                flag = '0'
                                rname = 'contig_name'
                                pos = str(n.entries[0].pos)
                                mapq = '0'
                                cigar = merge_entry_cigars(self.settings, n.entries)
                                rnext = '*'
                                pnext = '*'
                                tlen = '0'
                                seq = '*'
                                qual = ''
                                s += qname + ' ' + flag + ' ' + rname + ' ' + pos + ' ' + mapq + ' ' + cigar + ' ' + rnext + ' ' + pnext + ' ' + tlen + ' ' + seq + ' ' + qual + '\n'
                return s


def align_to_graph(settings):

        dbg = DBGraph()
        ofile = open('output.sam', 'w')
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
                dbg.nodes.append(node)

        sg = []
        current_reference = ''
        for t, l in sam_parse(settings.samName):
                if t == 'meta':
                        ofile.write(l)
                        print('meta: ' + l[:-1])
                else:
                        entry = Entry(settings, l)
                        if l[2] != current_reference:
                                current_reference = l[2]
                                sg.append(SAMGraph(settings, dbg))
                        sg[-1].match_entry(entry)

        print('Finished reading data.')

        for g in sg:
                g.collapse_greedy()
                g.collapse_linear()
                print('Deleted', len([n for n in g.nodes if n.deleted]), len(g.nodes))
                g.filter()
                print('Deleted', len([n for n in g.nodes if n.deleted]), len(g.nodes))
                g.cutoff(5000)
                print('Deleted', len([n for n in g.nodes if n.deleted]), len(g.nodes))
                print('\nContig results:')
                print(g)
                ofile.write(str(g))

def main(argv = None):
        if argv is None:
                argv = sys.argv
        try:
                opts, args = getopt.getopt(argv[1:], 'hg:s:k:', ['help', 'graph=', 'sam=', 'kmersize='])
        except getopt.error:
                print >>sys.stderr, 'For help use --help'
                return 2
        settings = Settings()
        for opt, val in opts:
                if opt == '-h' or opt == '--help':
                        print('Help message should come here.') #TODO
                elif opt == '-g' or opt == '--graph':
                        settings.graphName = val
                elif opt == '-s' or opt == '--sam':
                        settings.samName = val
                elif opt == '-k' or opt == '--kmersize':
                        settings.k = int(val)
        align_to_graph(settings)

if __name__ == "__main__":
        sys.exit(main())
