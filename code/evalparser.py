## evalparser.py
## Author: Yangfeng Ji
## Date: 11-05-2014
## Time-stamp: <yangfeng 09/25/2015 16:32:42>

from model import ParsingModel
from tree import RSTTree
from docreader import DocReader
from evaluation import Metrics
from os import listdir
from os.path import join as joinpath
from util import drawrst

try:
	import cPickle as pickle
except ImportError:
	import pickle
	
import pip
#pip.main(['install', '--user', 'pandas'])

import pandas as pd
from pandas.io.common import EmptyDataError

parse_tree = None
tokendict = None
edudict = None
coarse_to_fine_mapper = {}
doc_to_edu_relations = {}

def parse(pm, doc):
    """ Parse one document using the given parsing model

    :type pm: ParsingModel
    :param pm: an well-trained parsing model

    :type fedus: string
    :param fedus: file name of an document (with segmented EDUs) 
    """
    pred_rst = pm.sr_parse(doc)
    return pred_rst


def writebrackets(fname, brackets):
    """ Write the bracketing results into file
    """
    print 'Writing parsing results into file: {}'.format(fname)
    with open(fname, 'w') as fout:
        for item in brackets:
            fout.write(str(item) + '\n')


def evalparser(coref_path, save_path, path='./examples', report=False, 
               bcvocab=None, draw=True, use_entities=True,
               withdp=False, fdpvocab=None, fprojmat=None):
    """ Test the parsing performance

    :type path: string
    :param path: path to the evaluation data

    :type report: boolean
    :param report: whether to report (calculate) the f1 score
    """
    # ----------------------------------------
    # Load the parsing model
    print 'Load parsing model ...'
    pm = ParsingModel(withdp=withdp,
        fdpvocab=fdpvocab, fprojmat=fprojmat)
    pm.loadmodel("model/parsing-model.pickle.gz")
    # ----------------------------------------
    # Evaluation
    met = Metrics(levels=['span','nuclearity','relation'])
    # ----------------------------------------
    # Read all files from the given path
    doclist = [joinpath(path, fname) for fname in listdir(path) if fname.endswith('.merge')]
    doc_to_edu_relations = {}
    for fmerge in doclist:
        # ----------------------------------------
        # Read *.merge file
        dr = DocReader()
        doc = dr.read(fmerge)
        # ----------------------------------------
        # Parsing
        pred_rst = pm.sr_parse(doc, bcvocab)
        
        global parse_tree
        parse_tree = pred_rst
        
        if use_entities:
             global tokendict
             global edudict
             tokendict = doc.tokendict
             edudict = doc.edudict

             #Read coref file
             coref_fname = fmerge[fmerge.rfind("/")+1:fmerge.rfind(".")] + "_grid.tsv"
             try:
                  df_corefs = pd.read_csv(joinpath(coref_path,coref_fname), sep="\t", header=0, index_col=False, dtype=str)
                  df_grid = df_corefs.applymap(get_relations)#, args=(pred_rst,))
                  save_fname = fmerge[fmerge.rfind("/")+1:fmerge.rfind(".")] + "_rst.csv"
                  df_grid.to_csv(joinpath(coref_path,save_fname), encoding='utf-8')
                  print ('Wrote results to ', joinpath(coref_path,save_fname))
             except EmptyDataError:
                  print "Coref grid file is empty! No rst grid will be created.", coref_fname
        else:
             build_relation_mapper()
             leaves = get_leaves(parse_tree)
             assert (len(leaves) == len(doc.edudict)), "The number of leaf nodes in the parse tree is not equal to the number of EDUs!  leaf nodes: %r, EDUs: %r" % (len(leaves), len(doc.edudict))
             rels = get_relations_for_leaves(leaves)
             doc = fmerge[fmerge.rfind("/")+1:fmerge.rfind(".")]
             doc_to_edu_relations[doc] = rels
             
        if draw:
            strtree = pred_rst.parse()
            drawrst(strtree, fmerge.replace(".merge",".ps"))
        # Get brackets from parsing results
        pred_brackets = pred_rst.bracketing()
        fbrackets = fmerge.replace('.merge', '.brackets')
        # Write brackets into file
        writebrackets(fbrackets, pred_brackets)
        # ----------------------------------------
        # Evaluate with gold RST tree
        if report:
            fdis = fmerge.replace('.merge', '.dis')
            gold_rst = RSTTree(fdis, fmerge)
            gold_rst.build()
            gold_brackets = gold_rst.bracketing()
            met.eval(gold_rst, pred_rst)
    if report:
        met.report()
    
    if not use_entities:
        pickle.dump(doc_to_edu_relations, open(save_path, "wb" ))       

def get_relations(spans):
    relations = []
    spans = spans.strip('[]').split(",")
    for span in spans:
        span_relations = []
        span = span.strip().split()
        #check span is not empty
        if span:
            #print 'going to look for span: ', span
            matching_node = parse_tree.get(span)
            matching_num_sents = get_num_sents(matching_node)
            #print ('edu span: ', matching_node.eduspan)
            #print ('relation: ', matching_node.relation)
            #print ('nucedu: ', matching_node.nucedu)
            #print ('prop: ', matching_node.prop)
            
            if matching_node.prop == "Nucleus":
                #traverse up all the nucleus EDUs
                current_node = matching_node
                current_num_sents = matching_num_sents
                while (current_node.prop == "Nucleus" and current_num_sents <= 2):
                    relation = current_node.relation
                    #in case of span, get actual relation from sister node
                    if (current_node.relation == "span"):
                        if (current_node.pnode.rnode.prop == "Satellite"):
                            relation = current_node.pnode.rnode.relation
                        else:
                            relation = current_node.pnode.lnode.relation
                    span_relations.append(relation+".N")
                    current_node = current_node.pnode
                    current_num_sents = get_num_sents(current_node)
                    
                #now get the last satellite relation (unless we traversed all the way up to the root)
                if(current_node.relation and current_num_sents <= 2):
                    span_relations.append(current_node.relation+".S")
            elif matching_num_sents <= 2:
                span_relations.append(matching_node.relation+".S")
            #print ('final relations list for one span: ', span_relations)
        relations.append(span_relations)    
    return relations
    
def get_num_sents(node):
    #print ('edu span: ', node.eduspan)
    if node.eduspan[0] != node.eduspan[1]:
        begin_span = edudict[node.eduspan[0]]
        end_span = edudict[node.eduspan[1]]
        begin_sent = tokendict[begin_span[0]].sidx
        end_sent = tokendict[end_span[len(end_span)-1]].sidx
        #print('begin edu to token indx: ', begin_span)
        #print('end edu to token indx: ', end_span)
        #print ('first edu sidx: ', begin_sent)
        #print ('second edu sidx: ', end_sent)
        return abs(begin_sent - end_sent) + 1
    else:
        return 1

def get_leaves(root):
    leafs = []
    def _get_leaf_nodes( node):
        if node is not None:
            if node.lnode is None and node.rnode is None:
                leafs.append(node)
            for n in [node.lnode, node.rnode]:
                _get_leaf_nodes(n)
    _get_leaf_nodes(root.tree)
    return leafs

def get_relations_for_leaves(leaves):
    global coarse_to_fine_mapper
    relations = []
    for leaf in leaves:
        if (leaf.relation == "span"):
            if (leaf.pnode.rnode.prop == "Satellite"):
                relation = leaf.pnode.rnode.relation
            else:
                relation = leaf.pnode.lnode.relation
        else:
            relation  = leaf.relation
        if leaf.prop == "Nucleus":
            relations.append(coarse_to_fine_mapper[relation + ".N"])
        else:
            relations.append(coarse_to_fine_mapper[relation + ".S"])
    return relations
    
def build_relation_mapper():
        global coarse_to_fine_mapper
        coarse_to_fine_mapper['attribution.N'] = "Attribution.N"
        coarse_to_fine_mapper['attribution.S'] = "Attribution.S"

        coarse_to_fine_mapper['background.N'] = "Background.N"
        coarse_to_fine_mapper['circumstance.N'] = "Background.N"
        coarse_to_fine_mapper['background.S'] = "Background.S"
        coarse_to_fine_mapper['circumstance.S'] = "Background.S"

        coarse_to_fine_mapper['cause.N'] = "Cause.N"
        coarse_to_fine_mapper['result.N'] = "Cause.N"
        coarse_to_fine_mapper['consequence.N'] = "Cause.N"
        coarse_to_fine_mapper['cause.S'] = "Cause.S"
        coarse_to_fine_mapper['result.S'] = "Cause.S"
        coarse_to_fine_mapper['consequence.S'] = "Cause.S"

        coarse_to_fine_mapper['comparison.N'] = "Comparison.N"
        coarse_to_fine_mapper['preference.N'] = "Comparison.N"
        coarse_to_fine_mapper['analogy.N'] = "Comparison.N"
        coarse_to_fine_mapper['proportion.N'] = "Comparison.N"
        coarse_to_fine_mapper['comparison.S'] = "Comparison.S"
        coarse_to_fine_mapper['preference.S'] = "Comparison.S"
        coarse_to_fine_mapper['analogy.S'] = "Comparison.S"

        coarse_to_fine_mapper['condition.N'] = "Condition.N"
        coarse_to_fine_mapper['hypothetical.N'] = "Condition.N"
        coarse_to_fine_mapper['contingency.N'] = "Condition.N"
        coarse_to_fine_mapper['otherwise.N'] = "Condition.N"
        coarse_to_fine_mapper['condition.S'] = "Condition.S"
        coarse_to_fine_mapper['hypothetical.S'] = "Condition.S"
        coarse_to_fine_mapper['contingency.S'] = "Condition.S"
        coarse_to_fine_mapper['otherwise.S'] = "Condition.S"

        coarse_to_fine_mapper['contrast.N'] = "Contrast.N"
        coarse_to_fine_mapper['concession.N'] = "Contrast.N"
        coarse_to_fine_mapper['antithesis.N'] = "Contrast.N"
        coarse_to_fine_mapper['concession.S'] = "Contrast.S"
        coarse_to_fine_mapper['antithesis.S'] = "Contrast.S"

        coarse_to_fine_mapper['elaboration.N'] = "Elaboration.N"
        coarse_to_fine_mapper['example.N'] = "Elaboration.N"
        coarse_to_fine_mapper['definition.N'] = "Elaboration.N"
        coarse_to_fine_mapper['elaboration.S'] = "Elaboration.S"
        coarse_to_fine_mapper['example.S'] = "Elaboration.S"
        coarse_to_fine_mapper['definition.S'] = "Elaboration.S"
        
        coarse_to_fine_mapper['purpose.N'] = "Enablement.N"
        coarse_to_fine_mapper['enablement.N'] = "Enablement.N"
        coarse_to_fine_mapper['purpose.S'] = "Enablement.S"
        coarse_to_fine_mapper['enablement.S'] = "Enablement.S"

        coarse_to_fine_mapper['evaluation.N'] = "Evaluation.N"
        coarse_to_fine_mapper['interpretation.N'] = "Evaluation.N"
        coarse_to_fine_mapper['conclusion.N'] = "Evaluation.N"
        coarse_to_fine_mapper['evaluation.S'] = "Evaluation.S"
        coarse_to_fine_mapper['interpretation.S'] = "Evaluation.S"
        coarse_to_fine_mapper['conclusion.S'] = "Evaluation.S"

        coarse_to_fine_mapper['evidence.N'] = "Explanation.N"
        coarse_to_fine_mapper['explanation.N'] = "Explanation.N"
        coarse_to_fine_mapper['reason.N'] = "Explanation.N"
        coarse_to_fine_mapper['evidence.S'] = "Explanation.S"
        coarse_to_fine_mapper['explanation.S'] = "Explanation.S"
        coarse_to_fine_mapper['reason.S'] = "Explanation.S"

        coarse_to_fine_mapper['list.N'] = "Joint.N"
        coarse_to_fine_mapper['disjunction.N'] = "Joint.N"

        coarse_to_fine_mapper['manner.N'] = "Manner-Means.N"
        coarse_to_fine_mapper['means.N'] = "Manner-Means.N"
        coarse_to_fine_mapper['manner.S'] = "Manner-Means.S"
        coarse_to_fine_mapper['means.S'] = "Manner-Means.S"

        coarse_to_fine_mapper['problem.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['question.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['statement.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['topic.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['comment.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['rhetorical.N'] = "Topic-Comment.N"
        coarse_to_fine_mapper['problem.S'] = "Topic-Comment.S"
        coarse_to_fine_mapper['question.S'] = "Topic-Comment.S"
        coarse_to_fine_mapper['statement.S'] = "Topic-Comment.S"
        coarse_to_fine_mapper['topic.S'] = "Topic-Comment.S"
        coarse_to_fine_mapper['comment.S'] = "Topic-Comment.S"
        coarse_to_fine_mapper['rhetorical.S'] = "Topic-Comment.S"

        coarse_to_fine_mapper['summary.N'] = "Summary.N"
        coarse_to_fine_mapper['restatement.N'] = "Summary.N"
        coarse_to_fine_mapper['summary.S'] = "Summary.S"
        coarse_to_fine_mapper['restatement.S'] = "Summary.S"

        coarse_to_fine_mapper['temporal.N'] = "Temporal.N"
        coarse_to_fine_mapper['sequence.N'] = "Temporal.N"
        coarse_to_fine_mapper['inverted.N'] = "Temporal.N"
        coarse_to_fine_mapper['temporal.S'] = "Temporal.S"

        coarse_to_fine_mapper['same_unit.N'] = "Same-unit.N"

        coarse_to_fine_mapper['textualorganization.N'] = "Textual-organization.N"
