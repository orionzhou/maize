#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
From synteny blocks, reconstruct ancestral order by interleaving the genes in
between the anchors. This is the bottom-up method used first in Bowers (2003),
and in Tang (2010), to reconstruct pre-alpha and pre-rho order, respectively.
"""
from __future__ import print_function

import sys
import logging

from math import sqrt
from six.moves import zip_longest

from jcvi.compara.synteny import AnchorFile, check_beds
from jcvi.formats.bed import Bed
from jcvi.utils.grouper import Grouper
from jcvi.apps.base import OptionParser, ActionDispatcher

def main():

    actions = (
        ('mergechrom', 'merge synteny blocks on the same chrom'),
        ('pairs', 'convert anchorsfile to pairsfile'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def pairs(args):
    """
    %prog pairs anchorsfile prefix

    Convert anchorsfile to pairsfile.
    """
    p = OptionParser(pairs.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    anchorfile, prefix = args
    outfile = prefix + ".pairs"
    fw = open(outfile, "w")

    af = AnchorFile(anchorfile)
    blocks = af.blocks
    pad = len(str(len(blocks)))
    npairs = 0
    for i, block in enumerate(blocks):
        block_id = "{0}{1:0{2}d}".format(prefix, i + 1, pad)
        lines = []
        for q, s, score in block:
            npairs += 1
            score = score.replace('L', '')
            lines.append("\t".join((q, s, score, block_id)))
        print("\n".join(sorted(lines)), file=fw)

    fw.close()
    logging.debug("A total of {0} pairs written to `{1}`.".
                  format(npairs, outfile))

def get_collinear(block):
    # block contains (gene a, gene b, score)
    asc_score, asc_chain = print_chain(block)
    desc_score, desc_chain = print_chain(block, ascending=False)
    return asc_chain if asc_score > desc_score else desc_chain

def print_chain(block, ascending=True):

    scope = 50  # reduce search complexity
    if not ascending:
        block = [(a, -b, c) for (a, b, c) in block]

    block.sort()
    bsize = len(block)
    fromm = [-1] * bsize
    scores = [score_convert(c) for (a, b, c) in block]

    for i, (a, b, c) in enumerate(block):
        for j in range(i + 1, i + scope):
            if j >= bsize:
                break

            d, e, f = block[j]

            # Ensure strictly collinear
            if d == a or b >= e:
                continue

            this_score = scores[i] + score_convert(f)
            if this_score > scores[j]:
                fromm[j] = i
                scores[j] = this_score

    scoresfromm = list(zip(scores, fromm))
    maxchain = max(scoresfromm)
    chainscore, chainend = maxchain
    solution = [scoresfromm.index(maxchain), chainend]
    last = chainend
    while True:
        _last = fromm[last]
        if _last == -1:
            break
        last = _last
        solution.append(last)

    solution.reverse()
    solution = [block[x] for x in solution]
    if not ascending:
        solution = [(a, -b, c) for (a, b, c) in solution]
    return chainscore, solution

def mergechrom(args):
    """
    %prog mergechrom a.b.anchors

    merge synteny blocks on the same chromosome
    """
    p = OptionParser(mergechrom.__doc__)
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    af = AnchorFile(anchorfile)
    newanchorfile = anchorfile.rsplit(".", 1)[0] + ".mergechrom.anchors"
    fw = open(newanchorfile, "w")

    qchrom_dic = dict((b.accn,b.seqid) for b in qbed)
    schrom_dic = dict((b.accn,b.seqid) for b in sbed)
    block_dic = dict()
    blocks = af.blocks
    for (i,block) in enumerate(blocks):
        q, s, score = block[0]
        qchrom, schrom = qchrom_dic[q], schrom_dic[s]
        k = "%s_%s" % (qchrom, schrom)
        if k not in block_dic: block_dic[k] = []
        block_dic[k].append(i)

    for (k, idxs) in block_dic.items():
        print("#" * 3, file=fw)
        for i in idxs:
            for q, s, score in blocks[i]:
                print("\t".join((q, s, str(score))), file=fw)

    fw.close()
    print("%d blocks merged to %d" % (len(blocks), len(block_dic.keys())))

if __name__ == '__main__':
    main()
