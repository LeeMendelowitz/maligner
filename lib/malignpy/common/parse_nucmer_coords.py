#!/bin/bash
import pandas as pd

class Record(object):
    num_fields = 13

    keys = ['s1', 'e1', 's2', 'e2', 'l1', 'l2','idy','lenr','lenq','covr','covq', 'r','q']
    index = keys + ['orientation']

    def __init__(self, s1, e1, s2, e2, l1, l2, idy, lenr, lenq, covr, covq, r, q):
        self.s1 = int(s1)
        self.e1 = int(e2)
        self.s2 = int(s2)
        self.e2 = int(e2)
        self.l1 = int(l1)
        self.l2 = int(l2)
        self.idy = float(idy)
        self.lenr = int(lenr)
        self.lenq = int(lenq)
        self.covr = float(covr)
        self.covq = float(covq)
        self.r = r
        self.q = q

    @classmethod
    def from_line(cls, line):
        fields = line.strip().split()
        fields = [f for f in fields if f != '|']
        if len(fields) != cls.num_fields:
            return None
        return cls(*fields)
        
    @property
    def is_forward(self):
        return self.e2 > self.s2
    
    @property
    def is_reverse(self):
        return not self.is_forward

    def to_dict(self):
        ret = {k:getattr(self,k,None) for k in self.keys}
        ret['orientation'] = 'F' if self.is_forward else 'R'
        return pd.Series(ret, index = self.index, dtype=object)


def parse_file(f):
    with open(f) as fin:
        # Skip 5 line header
        for i in range(5):
            fin.next()
        records = [Record.from_line(l).to_dict() for l in fin]
    d = pd.DataFrame(records)
    # Reorder columns
    d = d[Record.index]
    return d

def run(f, pfx):
    MIN_QUERY_LENGTH=1000
    MIN_COVQ=90
    d = parse_file(f)
    d.to_csv('%s.nucmer.coords.tsv'%pfx, sep='\t', index=False)

    d2 = d[(d.covq >= MIN_COVQ) & (d.lenq >= MIN_QUERY_LENGTH)]
    d2 = d2[['q', 's1', 'e1', 'orientation']]
    d2.to_csv('%s.nucmer.locations.tsv'%pfx, sep='\t', index=False)
    

if __name__ == '__main__':
    f = sys.argv[1]
    pfx = sys.argv[2]
    run(f, pfx)
 
