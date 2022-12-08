import sys
from ImageD11.indexing import indexer

def test():
    i = indexer()
    i.readgvfile('simAl.gve')
    i.score_all_pairs()
    assert len(i.ubis)==20