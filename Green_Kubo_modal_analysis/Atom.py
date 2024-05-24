
class Atom:
    def __init__(self, id, type, x, y, z, vx, vy, vz):
        assert isinstance(x, float)
        self.id = id
        self.type = type
        #self.masss = mass
        self.x = [x, y, z]
        self.v = [vx, vy, vz]

