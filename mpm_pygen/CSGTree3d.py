#!/usr/bin/env python
import Primitives3d

class Node:
    def __init__(self, primitive=None, left=None, right=None, operation=None):
        self.primitive = primitive
        self.left = left
        self.right = right
        self.operation = operation

    def is_leaf(self):
        return (self.left == None and self.right == None)

    def evaluate(self, point):
        if (self.is_leaf()):
            return self.primitive.encompasses(point)
        else:
            # lv = self.left.evaluate(point)
            # rv = self.right.evaluate(point)
            # return self.operation(lv, rv)
            return self.operation(self.left.evaluate(point), self.right.evaluate(point))

if __name__ == '__main__':
    lcenter = Primitives.Point(0,1)
    rcenter = Primitives.Point(1,0)
    nl = Node(Primitives.Circle(lcenter, 1.2))
    nr = Node(Primitives.Circle(rcenter, 0.5))

    nl2 = Node(None, nl, nr, lambda x,y: x and y)

    a = Primitives.Point(0, 0)
    b = Primitives.Point(1, 0)
    c = Primitives.Point(0, 1)
    tri = Primitives.ConvexPolygon([a,b,c])

    nr2 = Node(tri)

    tree = Node(None, nl2, nr2, lambda x,y: x or y)

    N = 1000
    for i in range(0,N):
        for j in range(0,N):
            x = 3 * float(j) / N - 1.0
            y = 3 * float(i) / N - 1.0
            p = Primitives.Point(x, y)
            if tree.evaluate(p):
               print str(p.x)+', '+str(p.y) 
    '''
            if tree.evaluate(p):
                print 'X',
            else:
                print '.',
        print ''
    '''
