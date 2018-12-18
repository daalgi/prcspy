

class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return "Point({}, {})".format(self.x, self.y)

    def move(self, dx, dy):
        self.x += dx
        self.y += dy

class Polygon(object):
    def __init__(self, *args):
        self.points = args

    def __repr__(self):
        return 'Polygon(' + ', '.join(map(lambda p: str(p), self.points)) + ')'

    def move(self, dx, dy):
        for p in self.points:
            p.move(dx, dy)