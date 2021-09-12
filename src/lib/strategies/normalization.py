class base(object):
    def __init__(self, matrix_type):
        self.matrix_type = matrix_type


class standard(base):
    def normalize(self, difference, reference):
        return difference / reference

    def get_label(self):
        return "$\Delta$ %s/%s" % (self.matrix_type, self.matrix_type)

    def get_prefix(self):
        return "standard"


class nonzero(base):
    def normalize(self, difference, reference):
        nonzero = reference > 0
        return difference[nonzero] / reference[nonzero]

    def get_label(self):
        return "$\Delta$ %s/%s (%s>0)" % (self.matrix_type, self.matrix_type, self.matrix_type)

    def get_prefix(self):
        return "nonzero"


class plusone(base):
    def normalize(self, difference, reference):
        return difference / (reference + 1)

    def get_label(self):
        return "$\Delta$ %s/(%s+1)" % (self.matrix_type, self.matrix_type)

    def get_prefix(self):
        return "plusone"


class none(base):
    def normalize(self, difference, reference):
        return difference

    def get_label(self):
        return "$\Delta$ %s" % self.matrix_type

    def get_prefix(self):
        return ""
