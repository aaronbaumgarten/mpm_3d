#!/usr/bin/env python
import sys
import Primitives3d
import CSGTree3d

class CartesianPointGrid:
    def __init__(self, Lx, Ly, Lz, num_linear_elements_x, num_linear_elements_y, num_linear_elements_z, linear_element_material_point_density):
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.Lz = float(Lz)
        self.Nex = num_linear_elements_x
        self.Ney = num_linear_elements_y
        self.Nez = num_linear_elements_z
        self.Nnx = self.Nex + 1
        self.Nny = self.Ney + 1
        self.Nnz = self.Nez + 1
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
        dx = self.Lx / self.Nex
        dy = self.Ly / self.Ney
        dz = self.Lz / self.Nez

        ijk_array = [(i, j, k) for i in range(0, self.Nex) for j in range(0, self.Ney) for k in range(0, self.Nez)]

        xyz_nodes = map(lambda ijk: (self.Lx*float(ijk[0])/self.Nex, self.Ly*float(ijk[1])/self.Ney, self.Lz*float(ijk[2])/self.Nez), ijk_array)

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
            f.write('{} {} {}\n{} {} {}\n'.format(self.Nnx,self.Nny,self.Nnz,self.Lx,self.Ly,self.Lz))

if __name__ == '__main__':
    g = CartesianPointGrid(1, 1, 1, 100, 1)
    print g.point_array
