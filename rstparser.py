## main.py
## Author: Yangfeng Ji
## Date: 09-25-2015
## Time-stamp: <yangfeng 09/26/2015 00:10:59>

from code.evalparser import evalparser
from cPickle import load
import gzip, sys

def main(coref_path, save_path, path, draw=True, use_entities=True):
    with gzip.open("resources/bc3200.pickle.gz") as fin:
        print 'Load Brown clusters for creating features ...'
        bcvocab = load(fin)
    evalparser(coref_path, save_path, path=path, report=False, draw=draw, use_entities=use_entities, 
               bcvocab=bcvocab,
               withdp=False)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        coref_path = sys.argv[1]
        path = sys.argv[2]
        print 'Read files from: {}'.format(path)
        main(coref_path, path)
    elif len(sys.argv) == 4:
        coref_path = sys.argv[1]
        path = sys.argv[2]
        draw = eval(sys.argv[3])
        print 'Read files from {}'.format(path)
        main(coref_path, path, draw)
    elif len(sys.argv) == 6:
        coref_path = sys.argv[1]
        save_path = sys.argv[2]
        path = sys.argv[3]
        draw = eval(sys.argv[4])
        use_entities = eval(sys.argv[5])
        print 'Read files from {}'.format(path)
        main(coref_path, save_path, path, draw, use_entities)
    else:
        print "Usage: python rstparser.py coref_path file_path [draw_rst_tree] [use_entities]"
        print "\tcoref_path - path to the coref grids"
        print "\tfile_path - path to the segmented file"

