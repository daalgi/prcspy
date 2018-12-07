from shapely.geometry import LineString, MultiLineString, Polygon, Point
from shapely import affinity
import math


class Section(object):
    distance_precision = 1e-6
    infinite_distance = 1e8
    width = 0.0

    def __init__(self, polygon):
        self.polygon = polygon
        self.angle = 0.0
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return self.polygon.area

    def extract_poly_coords(self, polygon):
        if polygon.type == 'Polygon':
            exterior_coords = polygon.exterior.coords[:]
            interior_coords = []
            for i in polygon.interiors:
                interior_coords += i.coords[:]
        elif polygon.type == 'MultiPolygon':
            exterior_coords = []
            interior_coords = []
            for part in polygon:
                epc = self.extract_poly_coords(part)  # Recursive call
                exterior_coords += epc['exterior_coords']
                interior_coords += epc['interior_coords']
        else:
            raise ValueError('Unhandled geometry type: ' + repr(polygon.type))
        return {'exterior_coords': exterior_coords,
                'interior_coords': interior_coords}

    def eval_static_moment_x(self, polygon):
        sum = 0
        vertices = self.extract_poly_coords(polygon)
        x0, y0 = vertices['exterior_coords'][-2]
        for x1, y1 in vertices['exterior_coords'][:-1]:
            temp = y0 * y0 + y0 * y1 + y1 * y1
            sum += temp * (x1 - x0)
            x0 = x1; y0 = y1
        return -sum / 6.0

    def eval_static_moment_y(self, polygon):
        sum = 0
        vertices = self.extract_poly_coords(polygon)
        x0, y0 = vertices['exterior_coords'][-2]
        for x1, y1 in vertices['exterior_coords'][:-1]:
            temp = x0 * x0 + x0 * x1 + x1 * x1
            sum += temp * (y1 - y0)
            x0 = x1; y0 = y1
        return -sum / 6.0

    def eval_inertia_xx(self, polygon):
        sum = 0
        vertices = self.extract_poly_coords(polygon)
        x0, y0 = vertices['exterior_coords'][-2]
        for x1, y1 in vertices['exterior_coords'][:-1]:
            temp = y0 * y0 * y0 + y0 * y0 * y1 + y0 * y1 * y1 + y1 * y1 * y1
            sum += temp * (x1 - x0)
            x0 = x1; y0 = y1
        return -sum / 12.0

    def eval_inertia_xy(self, polygon):
        sum = 0
        vertices = self.extract_poly_coords(polygon)
        x0, y0 = vertices['exterior_coords'][-2]
        for x1, y1 in vertices['exterior_coords'][:-1]:
            sum += y0 * y0 * (x1 * x1 - x0 * x0)
            x0 = x1; y0 = y1
        return -sum / 4.0

    def eval_inertia_yy(self, polygon):
        sum = 0
        vertices = self.extract_poly_coords(polygon)
        x0, y0 = vertices['exterior_coords'][-2]
        for x1, y1 in vertices['exterior_coords'][:-1]:
            temp = x0 * x0 * x0 + x0 * x0 * x1 + x0 * x1 * x1 + x1 * x1 * x1
            sum += temp * (y1 - y0)
            x0 = x1; y0 = y1
        return +sum / 12.0

    def eval_centroid(self):
        return self.polygon.centroid

    def update_parameters(self):
        self.minx = self.polygon.bounds[0]
        self.miny = self.polygon.bounds[1]
        self.maxx = self.polygon.bounds[2]
        self.maxy = self.polygon.bounds[3]
        self.height = self.maxy - self.miny
        self.max_width = self.eval_max_width()
        self.centroid = self.eval_centroid()
        self.x0 = self.centroid.coords.xy[0][0]
        self.y0 = self.centroid.coords.xy[1][0]
        self.area = self.eval_area()
        self.static_moment_x = self.eval_static_moment_x(self.polygon)
        self.static_moment_y = self.eval_static_moment_y(self.polygon)
        self.inertia_xx = self.eval_inertia_xx(self.polygon)
        self.inertia_xy = self.eval_inertia_xy(self.polygon)
        self.inertia_yy = self.eval_inertia_yy(self.polygon)
        self.section_modulus_x = self.inertia_xx / max(self.maxy, abs(self.miny))
        self.section_modulus_y = self.inertia_yy / max(self.maxx, abs(self.minx))

    def width_at(self, vertical_coord):
        horizontal_line = LineString([(self.minx, vertical_coord), (self.maxx, vertical_coord)])
        intersection = self.polygon.intersection(horizontal_line)

        if isinstance(intersection, MultiLineString):
            width = 0
            for line in intersection:
                width += line.length
            return width

        elif isinstance(intersection, LineString):
            return intersection.length

        else:
            return 0

    def eval_max_width(self):
        return 0

    def max_x_at(self, vertical_coord):
        horizontal_line = LineString([(self.minx, vertical_coord), (self.maxx, vertical_coord)])
        intersection = self.polygon.intersection(horizontal_line)

        if isinstance(intersection, MultiLineString) or isinstance(intersection, LineString):
            return intersection.bounds[2]
        else:
            return 0

    def min_x_at(self, vertical_coord):
        horizontal_line = LineString([(self.minx, vertical_coord), (self.maxx, vertical_coord)])
        intersection = self.polygon.intersection(horizontal_line)

        if isinstance(intersection, MultiLineString) or isinstance(intersection, LineString):
            return intersection.bounds[0]
        else:
            return 0

    def center_section_at_centroid(self):
        x = (self.polygon.centroid.coords.xy[0])[0]
        y = (self.polygon.centroid.coords.xy[1])[0]
        self.displace_origin_by(x, y)

    def displace_origin_by(self, x, y):
        self.polygon = affinity.translate(self.polygon, -x, -y)
        self.update_parameters()

    def rotate(self, angle):
        self.angle = angle
        self.polygon = affinity.rotate(self.polygon, angle=angle,     # angle in degrees
                                       origin=self.polygon.centroid)
        self.update_parameters()

    def is_centroid_at_origin(self):
        circle_buffer = Point(0, 0).buffer(self.distance_precision)
        return self.polygon.centroid.within(circle_buffer)


class RectangularSection(Section):

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.polygon = Polygon([(0, 0), (width, 0), (width, height), (0, height)])
        self.angle = 0.0
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return self.width * self.height

    def width_at(self, vertical_coord):
        if self.angle == 0.:
            if self.miny <= vertical_coord <= self.maxy:
                return self.width
            else:
                return 0.
        else:
            return Section.width_at(self, vertical_coord)

    def eval_inertia_xx(self, polygon):
        """if self.angle == 0. and abs(self.maxy) == abs(self.miny):
        ixx = self.width * self.height ** 3. / 12.
        return ixx"""
        if self.angle == 0.:
            ixx = self.width * self.height ** 3. / 12.
            if abs(self.maxy) == abs(self.miny):
                return ixx
            else:
                return ixx + self.area * self.y0 ** 2.
        else:
            return Section.eval_inertia_yy(self, polygon)

    def eval_inertia_yy(self, polygon):
        """if self.angle == 0. and abs(self.maxy) == abs(self.miny):
        iyy = self.width ** 3. * self.height / 12.
        return iyy"""
        if self.angle == 0.:
            iyy = self.width ** 3. * self.height / 12.
            if abs(self.maxy) == abs(self.miny):
                return iyy
            else:
                return iyy + self.area * self.x0 ** 2.
        else:
            return Section.eval_inertia_yy(self, polygon)


class TshapedSection(Section):

    def __init__(self, flange_width, flange_thickness, web_height, web_thickness):
        self.flange_width = flange_width
        self.flange_thickness = flange_thickness
        self.web_height = web_height
        self.web_thickness = web_thickness
        self.polygon = Polygon([
            (-web_thickness / 2, -web_height), (+web_thickness / 2, -web_height),
            (+web_thickness / 2, 0.0), (+flange_width / 2, 0.0),
            (+flange_width / 2, +flange_thickness), (-flange_width / 2, +flange_thickness),
            (-flange_width / 2, 0.0), (-web_thickness / 2, 0.0),
            ])
        self.angle = 0.0
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return self.flange_thickness * self.flange_width \
               + self.web_thickness * self.web_height


class IshapedSection(Section):
    def __init__(self, upper_flange_width, upper_flange_thickness,
                 web_height, web_thickness,
                 lower_flange_width, lower_flange_thickness):
        self.upper_flange_width = upper_flange_width
        self.upper_flange_thickness = upper_flange_thickness
        self.web_height = web_height
        self.web_thickness = web_thickness
        self.lower_flange_width = lower_flange_width
        self.lower_flange_thickness = lower_flange_thickness
        self.polygon = Polygon([
            (-web_thickness / 2, -web_height), (-lower_flange_width / 2, -web_height),
            (-lower_flange_width / 2, -web_height - lower_flange_thickness),
            (+lower_flange_width / 2, -web_height - lower_flange_thickness),
            (+lower_flange_width / 2, -web_height),
            (+web_thickness / 2, -web_height),
            (+web_thickness / 2, 0.0), (+upper_flange_width / 2, 0.0),
            (+upper_flange_width / 2, +upper_flange_thickness),
            (-upper_flange_width / 2, +upper_flange_thickness),
            (-upper_flange_width / 2, 0.0),
            (-web_thickness / 2, 0.0)
        ])
        self.angle = 0.0
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return self.upper_flange_thickness * self.upper_flange_width \
               + self.web_thickness * self.web_height \
               + self.lower_flange_thickness * self.lower_flange_width


class RectangularHollowSection(Section):
    def __init__(self, width, height, web_thickness, flange_thickness):
        self.width = width
        self.height = height
        self.web_thickness = web_thickness
        self.flange_thickness = flange_thickness
        ext = [(0, 0), (width, 0), (width, height), (0, height)]
        int = [(+web_thickness, +flange_thickness),
               (width - web_thickness, +flange_thickness),
               (width - web_thickness, height - flange_thickness),
               (+web_thickness, height - flange_thickness)]
        self.angle = 0.0
        self.polygon = Polygon(ext, [int])
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return self.width * self.height - \
               (self.width - 2 * self.web_thickness) * \
               (self.height - 2 * self.flange_thickness)


class CircularSection(Section):
    def __init__(self, diameter):
        self.diameter = diameter
        self.width = diameter
        self.radius = diameter / 2.0
        self.polygon = Point(0, 0).buffer(self.radius)
        self.angle = 0.0
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return math.pi * self.radius ** 2

    def width_at(self, vertical_coord):
        if self.miny < vertical_coord < self.maxy:
            return 2 * math.sqrt(self.radius ** 2
                                 - (vertical_coord - self.y0) ** 2)
        else:
            return 0

    def eval_max_width(self):
        return self.diameter

    def eval_inertia_xx(self, polygon):
        ixx = math.pi * self.diameter ** 4. / 64.
        if abs(self.maxy) == abs(self.miny):
            return ixx
        else:
            return ixx + self.area * self.y0 ** 2.

    def eval_inertia_xy(self, polygon):
        if abs(self.maxy) == abs(self.miny) or abs(self.maxx) == abs(self.minx):
            return 0.
        else:
            return Section.eval_inertia_xy(self, polygon)

    def eval_inertia_yy(self, polygon):
        iyy = math.pi * self.diameter ** 4. / 64.
        if abs(self.maxx) == abs(self.minx):
            return iyy
        else:
            return iyy + self.area * self.x0 ** 2.


class CircularHollowSection(Section):
    def __init__(self, outer_diameter, inner_diameter):
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        self.thickness = (outer_diameter - inner_diameter) / 2.
        self.outer_radius = outer_diameter / 2.
        self.inner_radius = inner_diameter / 2.
        ext = Point(0, 0).buffer(self.outer_radius)
        int = Point(0, 0).buffer(self.inner_diameter)
        self.polygon = Polygon(ext, [int])
        self.angle = 0.
        self.center_section_at_centroid()
        self.update_parameters()

    def eval_area(self):
        return math.pi * (self.outer_radius ** 2. - self.inner_radius ** 2.)

    def eval_inertia_xx(self, polygon):
        ixx = self.eval_inertia()
        if abs(self.maxy) == abs(self.miny):
            return ixx
        else:
            return ixx + self.area * self.y0 ** 2.

    def eval_inertia_xy(self, polygon):
        if abs(self.maxy) == abs(self.miny) or abs(self.maxx) == abs(self.minx):
            return 0.
        else:
            return Section.eval_inertia_xy(self, polygon)

    def eval_inertia_yy(self, polygon):
        iyy = self.eval_inertia()
        if abs(self.maxx) == abs(self.minx):
            return iyy
        else:
            return iyy + self.area * self.x0 ** 2.

    def eval_inertia(self):
        return math.pi / 64. * \
               (self.outer_diameter ** 4. - self.inner_diameter ** 4.)

    def width_at(self, vertical_coord):
        if self.miny < vertical_coord < self.maxy:
            ext = 2. * math.sqrt(self.outer_radius ** 2. - (vertical_coord - self.y0) ** 2.)
            int = 0.
            if -self.inner_radius <= vertical_coord - self.y0 <= self.inner_radius:
                int = 2. * math.sqrt(self.inner_radius ** 2. - (vertical_coord - self.y0) ** 2 )
            return ext - int
        else:
            return 0.


class CircularElement(object):
    """
    Circular element with origin of coordinates at the center of the circle of its base
    """

    def thickness_at(self, radius):
        return NotImplemented

    def width_at(self, radius):
        return NotImplemented

    def volume(self):
        pass

    def weight(self):
        return self.volume() * self.material.specific_weight


class CircularSlab(CircularElement):

    def __init__(self, material, bottom_diameter, edge_thickness, top_diameter=None, slope_thickness=None,
                 tag=''):
        self.material = material
        self.bottom_diameter = bottom_diameter
        self.bottom_radius = bottom_diameter / 2.
        self.edge_thickness = edge_thickness
        self.top_diameter = top_diameter if top_diameter is not None else bottom_diameter
        self.top_radius = self.top_diameter / 2.
        self.slope_thickness = slope_thickness if slope_thickness is not None else 0
        self.thickness = self.edge_thickness + self.slope_thickness
        self.tag = tag

    def thickness_at(self, radius):

        if 0 < radius <= self.top_radius:

            return self.thickness

        elif self.top_radius < radius <= self.bottom_radius:

            h1 = self.edge_thickness
            h2 = self.slope_thickness
            r_top = self.top_radius
            r_bot = self.bottom_radius
            return (h1 + h2) + (abs(radius) - r_top) * ((h1) - (h1 + h2)) / (r_bot - r_top)

        else:
            return 0

    def width_at(self, radius):

        if -self.bottom_radius < radius < +self.bottom_radius:

            return 2 * math.sqrt(self.bottom_radius**2 - radius**2)

        else:

            return 0

    def volume(self):
        # Cylinder: pi * R ^ 2 * H
        vol = math.pi * self.bottom_diameter ** 2 / 4. * self.edge_thickness
        # Truncated cone: pi * (R1 ^ 2 + R1 * R2 + R2 ^ 2) * H / 3
        vol += math.pi / 4. * (self.bottom_diameter ** 2 +
                               self.top_diameter ** 2 +
                               self.bottom_diameter * self.top_diameter) / 3. * self.slope_thickness
        return vol