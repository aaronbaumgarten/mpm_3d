#!/usr/bin/env python
import sys
import Primitives3d
import CSGTree3d

class CartesianPointGrid:
    def __init__(self, Lx, Ly, Lz, num_linear_elements, linear_element_material_point_density):
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.Lz = float(Lz)
        self.Ne = num_linear_elements
        self.Nn = self.Ne + 1
        self.lmpdensity = linear_element_material_point_density
        self.primitive = None

    def generate_point_array(self,primitive=None):
        self.primitive = primitive
        self.body = CSGTree3d.Node(primitive)

        # cell spacing for material point
        cs = 1.0 / self.lmpdensity
        s = map(lambda k: (0.5 + k) * cs, range(0, self.lmpdensity))

        # material point positions inside cell
        xyz_s = [(x,y,z) for x in s for y in s for z in s]

        # element spacing
        dx = self.Lx / self.Ne
        dy = self.Ly / self.Ne
        dz = self.Lz / self.Ne

        ijk_array = [(i, j, k) for i in range(0, self.Ne) for j in range(0, self.Ne) for k in range(0, self.Ne)]

        xyz_nodes = map(lambda ijk: (self.Lx*float(ijk[0])/self.Ne, self.Ly*float(ijk[1])/self.Ne, self.Lz*float(ijk[2])/self.Ne), ijk_array)

        xyz_material_point_candidates = []
        i=0
        iEnd=len(xyz_nodes)*len(xyz_s)
        for n in xyz_nodes:
            for sp in xyz_s:
                xyz_material_point_candidate = Primitives3d.Point(n[0]+dx*sp[0], n[1]+dy*sp[1], n[2]+dz*sp[2])
                if (self.primitive==None or self.body.evaluate(xyz_material_point_candidate)):
                    xyz_material_point_candidates.append(xyz_material_point_candidate)
                i+=1
                sys.stdout.write("\r" + str(i) + " / " + str(iEnd))
        sys.stdout.write("\n")

        self.point_array = xyz_material_point_candidates
        self.material_point_volume = cs*cs*cs*dx*dy*dz
        self.element_volume = dx*dy*dz

    def write_grid_file(self, filename):
        with open(filename, 'w') as f:
            f.write('{}\n{}\n'.format(self.Nn, self.Lx))

if __name__ == '__main__':
    g = CartesianPointGrid(1, 1, 1, 100, 1)
    print g.point_array
