#!/usr/bin/env python
#-*-coding:utf-8-*-

import os,sys

def de_bruijn_ize(st, k):
    """ generate the edges and nodes of the k-mer de_bruijn graph of string """
    nodes = set()
    edges = []
    for i in range(len(st)-k+1):
        edges.append( (st[i:i+k-1],st[i+1:i+k]) )
        nodes.add(st[i:i+k-1])
        nodes.add(st[i:i+k-1])

    return nodes,edges

def visualize_de_bruijn(st, k):
    """ visualize a directed multigraph using graphviz. """
    nodes, edges = de_bruijn_ize(st, k)
    dot_str = 'digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str += ' %s [label="%s"] ;\n ' %(node,node)
    for src,dst in edges:
        dot_str += ' %s -> %s ;\n' %(src,dst)
    return dot_str + '}\n'





