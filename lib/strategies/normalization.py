def default(difference, reference):
    return difference / reference

def nonzero(difference, reference):
    nonzero = reference > 0
    return difference[nonzero] / reference[nonzero]

def plusone(difference, reference):
    return difference / (reference + 1)