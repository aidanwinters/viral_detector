from numpy.random import choice


def getProportions(seq):

    props = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    for s in seq:
        props[s] += 1

    for k in props.keys():
        props[k] = props[k]/len(seq)

    return list(props.values()), len(seq)

def genSequence(proportions, length):
    draw = choice(['A', 'T', 'C', 'G'], length, p=proportions)
    return (''.join(draw))
