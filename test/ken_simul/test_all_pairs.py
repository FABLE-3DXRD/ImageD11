import sys, os
from ImageD11.indexing import indexer

def test():
    i = indexer()
    i.readgvfile(os.path.join(os.path.dirname(__file__),'simAl.gve'))
    i.score_all_pairs()
    assert len(i.ubis)==20