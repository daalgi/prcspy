import abc
from material import ReinforcingSteel, PrestressingSteel
from shapely.geometry import Point
import math


diameters = [0.008, 0.010, 0.012, 0.016, 0.020, 0.025, 0.032, 0.040]
default_nominal_cover = 0.05
preferred_diameter = 0.016

def design_longitudinal_bars_horizontal_row(
        area, width_available, effective_depth, preferred_diameter=0.016,
        steel=ReinforcingSteel()):
    """
    :param area: needed reinforcement area
    :param width_available: width where the bars can be placed in the 
    concrete section: width_available = span - diameter
    :return: list of LongitudinalBar objects with the longitudinal 
    bar row configuration
    """
    if area > 0:
        span = width_available - preferred_diameter
        num_bars = num_bars_needed(area, preferred_diameter)
        # Reduce diameter if only one bar is needed
        while num_bars == 1:
            preferred_diameter = previous_diameter(preferred_diameter)
            span = width_available - preferred_diameter
            num_bars = num_bars_needed(area, preferred_diameter)

        clear_distance = bar_clear_distance(
            num_bars, span, preferred_diameter)

        # Greater diameter if the clear distance between bar is too small
        while clear_distance < ec2_bar_min_clear_distance(preferred_diameter):
            preferred_diameter = next_diameter(preferred_diameter)
            span = width_available - preferred_diameter
            num_bars = num_bars_needed(area, preferred_diameter)
            clear_distance = bar_clear_distance(
                num_bars, span, preferred_diameter)
        # print('d:', preferred_diameter)
        # print('n_bars:', num_bars)
        # print('b_dist:', bar_clear_distance)
        r = LongitudinalHorizontalBarsRow(
            num_bars=num_bars, diameter=preferred_diameter,
            effective_depth=effective_depth, steel=steel)
        r.update_row_configuration(num_bars, preferred_diameter, span)
        return r


def design_longitudinal_bars_circular_row(
        needed_reinforcement_area, section_diameter, effective_depth,
        steel=ReinforcingSteel()):
    pass


def design_circular_footing_radial_reinforcement(
        area_per_m, mfooting, preferred_diameter=0.032):
    if area_per_m > 0:
        area = area_per_m * math.pi * mfooting.pedestal_geometry.diameter
        num_bars = num_bars_needed(area, preferred_diameter)
        span = math.pi * mfooting.pedestal_geometry.diameter * \
               (num_bars - 1) / num_bars
        # Reduce diameter if only one bar is needed
        while num_bars == 1:
            preferred_diameter = previous_diameter(preferred_diameter)
            span = math.pi * mfooting.pedestal_geometry.diameter * \
                   (num_bars - 1) / num_bars
            num_bars = num_bars_needed(area, preferred_diameter)

        clear_distance = bar_clear_distance(
            num_bars, span, preferred_diameter)

        # Greater diameter if the clear distance between bar is too small
        while clear_distance < ec2_bar_min_clear_distance(preferred_diameter):
            preferred_diameter = next_diameter(preferred_diameter)
            span = math.pi * mfooting.pedestal_geometry.diameter * \
                   (num_bars - 1) / num_bars
            num_bars = num_bars_needed(area, preferred_diameter)
            clear_distance = bar_clear_distance(
                num_bars, span, preferred_diameter)

        effective_depth = sum(mfooting.buried_heights[0:2]) - \
                          mfooting.concrete.nominal_cover - preferred_diameter / 2.
        r = LongitudinalHorizontalBarsRow(
            num_bars=num_bars, diameter=preferred_diameter,
            effective_depth=effective_depth, steel=mfooting.reinforcing_steel)
        r.update_row_configuration(num_bars, preferred_diameter, span)
        return r


def num_bars_needed(area, diameter):
    return math.ceil(area / (math.pi * diameter ** 2 / 4))


def bar_clear_distance(num_bars, span, diameter):
    if num_bars > 1.0 and span > 0:
        return span / (num_bars - 1) - diameter
    else:
        return 0


def ec2_bar_min_clear_distance(diameter, max_aggregate_size=0.020):
    # EN 1992-1-1:2004, 8.2. Spacing of bars
    return max(diameter, max_aggregate_size + 0.005, 0.020)


def next_diameter(diameter):
    if diameter < max(diameters):
        return [x for x in diameters if x > diameter][0]
    else:
        # TODO: design a configuration of two layers
        return 0


def previous_diameter(diameter):
    if diameter > min(diameters):
        return [x for x in diameters if x < diameter][-1]
    else:
        # TODO: handle this case
        return 0


class LongitudinalBarsRow(object):
    """diameters = {       # m
        'EC2': [0.008, 0.010, 0.012, 0.016, 0.020, 0.025, 0.032, 0.040]
    }"""
    diameters = [0.008, 0.010, 0.012, 0.016, 0.020, 0.025, 0.032, 0.040]
    default_nominal_cover = 0.05
    preferred_diameter = 0.016
    bar = []

    @abc.abstractmethod
    def align_row_in_cross_section(self):
        pass

    def eval_num_bars(self, area, diameter):
        return math.ceil(area / (math.pi * diameter ** 2 / 4))

    def eval_bar_clear_distance(self, num_bars, span, diameter):
        if num_bars > 1.0 and span > 0:
            return span / (num_bars - 1) - diameter
        else:
            return 0

    def eval_bar_effective_distance(self, num_bars, span):
        if num_bars > 1.0 and span > 0:
            return span / (num_bars - 1)
        else:
            return 0

    def ec2_bar_min_clear_distance(self, diameter, max_aggregate_size=0.020):
        # EN 1992-1-1:2004, 8.2. Spacing of bars
        return max(diameter, max_aggregate_size + 0.005, 0.020)

    def next_diameter(self, diameter):
        if diameter < max(self.diameters):
            return [x for x in self.diameters if x > diameter][0]
        else:
            # TODO: design a configuration of two layers
            return 0

    def previous_diameter(self, diameter):
        if diameter > min(self.diameters):
            return [x for x in self.diameters if x < diameter][-1]
        else:
            # TODO: handle this case
            return 0


class LongitudinalHorizontalBarsRow(LongitudinalBarsRow):
    """::param bars_separation: distance between the centroid of consecutive bars in the row
    """
    def __init__(self, num_bars, diameter, effective_depth, span=0,
                 steel=ReinforcingSteel(), bar=[], area=None):
        """
        :param num_bars: number of bars 
        :param diameter (m): bar diameter
        :param effective_depth (m): effective depth of the bars of the row
        :param span (m): distance between the centroid of the extreme bars 
        in the row 
        :param steel: material
        :param bar: list of ReinforcementElement objects. 
        :param area (m2): row area. If defined, the diameter will be adjusted so
        for the given number of bars the area will be equal to 'area' 
        By default the vertical coordinates of the row are Y = -effective_depth"""
        self.effective_depth = effective_depth
        self.steel = steel
        if isinstance(area, (int, float)):
            diameter = math.sqrt(4. * area / (num_bars * math.pi))
        self.update_row_configuration(num_bars, diameter, span, bar)

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

    def update_row_configuration(self, num_bars, diameter, span, bar=[]):
        self.num_bars = num_bars
        self.diameter = diameter
        self.area = num_bars * math.pi * diameter ** 2. / 4.
        self.bar_effective_distance = self.eval_bar_effective_distance(num_bars, span)
        self.bar_clear_distance = self.eval_bar_clear_distance(num_bars, span, diameter)

        if span == 0:
            self.bar = bar
        elif span > 0:
            self.bar = self.update_bar_coords(
                cross_section_max_y=0,
                cross_section_min_x_at_effective_depth=-span / 2 - diameter / 2,
                cross_section_max_x_at_effective_depth=+span / 2 + diameter / 2,
                nominal_cover=0)
        else:
            self.span = 0
            self.bar = []

    def update_bar_coords(self, cross_section_max_y,
                          cross_section_min_x_at_effective_depth,
                          cross_section_max_x_at_effective_depth,
                          nominal_cover):
        """:param cross_section_max_y: 
        :param cross_section_min_x_at_effective_depth: 
        :param cross_section_max_x_at_effective_depth: 
        :param nominal_cover: distance between the section edge and the closest bar edge
        (may include the concrete cover and the transverse reinforcement)
        :return: list with LongitudinalBar objects"""
        self.span = cross_section_max_x_at_effective_depth \
                    - cross_section_min_x_at_effective_depth \
                    - 2 * nominal_cover - self.diameter
        self.bar_effective_distance = self.eval_bar_effective_distance(
            self.num_bars, self.span)
        self.bar_clear_distance = self.eval_bar_clear_distance(
            self.num_bars, self.span, self.diameter)
        return [LongitudinalBar(diameter=self.diameter,
                                steel=self.steel,
                                point=Point(cross_section_min_x_at_effective_depth
                                            + nominal_cover
                                            + self.diameter / 2
                                            + i * self.bar_effective_distance,
                                            cross_section_max_y - self.effective_depth)
                                ) for i in range(0, self.num_bars)]

    def align_row_in_cross_section(self, cross_section_max_y,
                                   cross_section_min_x_at_effective_depth,
                                   cross_section_max_x_at_effective_depth,
                                   nominal_cover):
        #span = cross_section_max_x_at_effective_depth \
        #       - cross_section_min_x_at_effective_depth \
        #       - 2 * nominal_cover - self.diameter
        #self.bar_clear_distance = span / (self.num_bars - 1)
        self.bar = self.update_bar_coords(
            cross_section_max_y=cross_section_max_y,
            cross_section_min_x_at_effective_depth=cross_section_min_x_at_effective_depth,
            cross_section_max_x_at_effective_depth=cross_section_max_x_at_effective_depth,
            nominal_cover=nominal_cover)

    def design_row(self, area, width_available, preferred_diameter=0.016):
        """
        :param area: needed reinforcement area
        :param width_available: width where the bars can be placed in the concrete section
         width_available = span - diameter
        :return: list of LongitudinalBar objects with the longitudinal bar row configuration
        """
        if area > 0:
            span = width_available - preferred_diameter
            num_bars = self.eval_num_bars(area, preferred_diameter)
            # Reduce diameter if only one bar is needed
            while num_bars == 1:
                preferred_diameter = self.previous_diameter(preferred_diameter)
                span = width_available - preferred_diameter
                num_bars = self.eval_num_bars(area, preferred_diameter)

            bar_clear_distance = self.eval_bar_clear_distance(
                num_bars, span, preferred_diameter)

            # Greater diameter if the clear distance between bar is too small
            while bar_clear_distance < self.ec2_bar_min_clear_distance(preferred_diameter):
                preferred_diameter = self.next_diameter(preferred_diameter)
                span = width_available - preferred_diameter
                num_bars = self.eval_num_bars(area, preferred_diameter)
                bar_clear_distance = self.eval_bar_clear_distance(
                    num_bars, span, preferred_diameter)
            #print('d:', preferred_diameter)
            #print('n_bars:', num_bars)
            #print('b_dist:', bar_clear_distance)
            self.update_row_configuration(num_bars, preferred_diameter, span)


class LongitudinalCircularBarsRow(LongitudinalBarsRow):
    """::param bars_separation: distance between the centroid of consecutive bars in the row"""
    def __init__(self, num_bars, diameter, effective_depth, csec_diameter,
                 steel=ReinforcingSteel()):
        """
        :param num_bars: number of bars 
        :param diameter: bar diameter
        :param effective_depth: defines the location of the row in the cross-section 
        :param csec_diameter: diameter of the concrete cross-section 
        :param steel: material
        :param bar: list of LongitudinalBar objects.
         By default the circular row it's centered at the origin (0,0)"""
        self.effective_depth = effective_depth
        self.steel = steel
        self.update_row_configuration(num_bars, diameter, csec_diameter)

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

    def update_row_configuration(self, num_bars, diameter, csec_diameter, csec_centroid=Point(0,0)):
        self.num_bars = num_bars
        self.diameter = diameter
        self.area = num_bars * math.pi * diameter ** 2 / 4.0
        self.csec_diameter = csec_diameter
        self.bar = self.update_bar_coords(csec_centroid, csec_diameter=csec_diameter)

    def update_bar_coords(self, csec_centroid, csec_diameter):
        """
        :param csec_centroid: Point object of the centroid of the concrete cross-section  
        :param csec_diameter: diameter of the concrete section
        :return: list with LongitudinalBar objects"""
        phi = [2 * math.pi / self.num_bars * i for i in range(0, self.num_bars)]
        r = self.effective_depth - csec_diameter / 2.0
        x = csec_centroid.coords.xy[0][0]
        y = csec_centroid.coords.xy[1][0]
        return [LongitudinalBar(diameter=self.diameter,
                                steel=self.steel,
                                point=Point(x + r * math.cos(p), y + r * math.sin(p)))
                                    for p in phi]

    def align_row_in_cross_section(self, csec_centroid=Point(0,0)):
        self.bar = self.update_bar_coords(csec_centroid, self.csec_diameter)


class LongitudinalTendonRow(object):
    diameters = [0.008, 0.010, 0.012, 0.016, 0.020, 0.025, 0.032, 0.040]
    default_nominal_cover = 0.05
    preferred_diameter = 0.016
    tendon = []

    @abc.abstractmethod
    def align_row_in_cross_section(self):
        pass

    def eval_num_tendons(self, area, diameter):
        return math.ceil(area / (math.pi * diameter ** 2 / 4))

    def eval_tendon_clear_distance(self, num_bars, span, diameter):
        if num_bars > 1.0 and span > 0:
            return span / (num_bars - 1) - diameter
        else:
            return 0

    def eval_tendon_effective_distance(self, num_bars, span):
        if num_bars > 1.0 and span > 0:
            return span / (num_bars - 1)
        else:
            return 0

    def ec2_tendon_min_clear_distance(self, diameter, max_aggregate_size=0.020):
        # EN 1992-1-1:2004, 8.2. Spacing of bars ----- CHEEEEEECK
        return max(diameter, max_aggregate_size + 0.005, 0.020)

    def next_diameter(self, diameter):
        if diameter < max(self.diameters):
            return [x for x in self.diameters if x > diameter][0]
        else:
            # TODO: design a configuration of two layers
            return 0

    def previous_diameter(self, diameter):
        if diameter > min(self.diameters):
            return [x for x in self.diameters if x < diameter][-1]
        else:
            # TODO: handle this case
            return 0


class LongitudinalHorizontalTendonsRow(LongitudinalTendonRow):
    def __init__(self, num_tendons, effective_depth,
                 diameter, hole_diameter, prestressing_force=0, initial_strain=0,
                 span=0, steel=PrestressingSteel(), tendon=[]):
        """
        :param num_tendons: number of tendons
        :param effective_depth: effective depth of the tendons of the row
        :param prestressing_force: in kN
        :param span: distance between the centroid of the extreme bars in the row 
        :param steel: material
        :param tendon: list of ReinforcementElement objects. 
        By default the vertical coordinates of the row are Y = -effective_depth"""
        self.effective_depth = effective_depth
        self.steel = steel
        self.update_row_configuration(num_tendons, diameter, hole_diameter,
                                      prestressing_force, initial_strain, span, tendon)

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

    def update_row_configuration(self, num_tendons, diameter, hole_diameter,
                                 prestressing_force, initial_strain, span, tendon=[]):
        self.num_tendons = num_tendons
        self.diameter = diameter
        self.area = num_tendons * math.pi * diameter ** 2 / 4.0
        self.hole_diameter = hole_diameter
        self.hole_area = num_tendons * math.pi * hole_diameter ** 2 / 4.0
        self.prestressing_force = prestressing_force
        self.initial_strain = initial_strain
        self.span = span

        if num_tendons == 1:
            self.tendon = self.update_tendon_coords(
                cross_section_max_y=0,
                row_min_x=0,
                row_max_x=0)
        elif span == 0:
            self.tendon = tendon
        elif span > 0:
            self.tendon = self.update_tendon_coords(
                cross_section_max_y=0,
                row_min_x=-span / 2,
                row_max_x=+span / 2)
        else:
            self.span = 0
            self.tendon = []

    def update_tendon_coords(self, cross_section_max_y, row_min_x, row_max_x):
        """
        :param cross_section_max_y: 
        :param row_min_x: Minimum X coord. of the centroid of the tendons in the row
        :param row_max_x: Maximum X coord. of the centroid of the tendons in the row
        :return: list with LongitudinalBar objects"""
        if self.num_tendons == 1:
            return [LongitudinalTendon(
                diameter=self.diameter, hole_diameter=self.hole_diameter,
                prestressing_force=self.prestressing_force,
                initial_strain=self.initial_strain,
                steel=self.steel,
                point=Point(row_min_x, cross_section_max_y - self.effective_depth))]

        elif self.num_tendons > 1:
            self.span = row_max_x - row_min_x
            self.tendon_effective_distance = self.eval_tendon_effective_distance(
                self.num_tendons, self.span)
            return [LongitudinalTendon(
                diameter=self.diameter, hole_diameter=self.hole_diameter,
                prestressing_force=self.prestressing_force,
                initial_strain=self.initial_strain,
                steel=self.steel,
                point=Point(row_min_x + i * self.tendon_effective_distance,
                            cross_section_max_y - self.effective_depth)
            ) for i in range(0, self.num_tendons)]
        else:
            return []

    def align_row_in_cross_section(self, cross_section_max_y, row_min_x, row_max_x):
        self.tendon = self.update_tendon_coords(
            cross_section_max_y=cross_section_max_y,
            row_min_x=row_min_x, row_max_x=row_max_x)


class LongitudinalCircularTendonsRow(LongitudinalTendonRow):
    """::param bars_separation: distance between the centroid of consecutive bars in the row"""
    def __init__(self, num_tendons, diameter, hole_diameter,
                 effective_depth, csec_diameter,
                 prestressing_force=0, initial_strain=0,
                 steel=PrestressingSteel()):
        """
        :param num_tendons: number of tendons 
        :param diameter (m): tendon diameter
        :param hole_diameter (m)
        :param effective_depth: defines the location of the row in the cross-section 
        :param csec_diameter: diameter of the concrete cross-section
        :param prestressing_force (kN)
        :param initial_strain (m/m)
        :param steel: material
         By default the circular row it's centered at the origin (0,0)"""
        self.effective_depth = effective_depth
        self.steel = steel
        self.update_row_configuration(num_tendons, diameter, hole_diameter,
                                      csec_diameter, Point(0,0),
                                      prestressing_force, initial_strain,
                                      )

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

    def update_row_configuration(self, num_tendons, diameter, hole_diameter,
                                 csec_diameter, csec_centroid=Point(0, 0),
                                 prestressing_force=0, initial_strain=0):
        self.num_tendons = num_tendons
        self.diameter = diameter
        self.area = num_tendons * math.pi * diameter ** 2 / 4.0
        self.hole_diameter = hole_diameter
        self.hole_area = num_tendons * math.pi * hole_diameter ** 2 / 4.0
        self.prestressing_force = prestressing_force
        self.initial_strain = initial_strain
        self.csec_diameter = csec_diameter
        self.tendon = self.update_bar_coords(csec_centroid, csec_diameter=csec_diameter)

    def update_bar_coords(self, csec_centroid, csec_diameter):
        """
        :param csec_centroid: Point object of the centroid of the concrete cross-section  
        :param csec_diameter: diameter of the concrete section
        :return: list with LongitudinalTendon objects"""
        phi = [2 * math.pi / self.num_tendons * i for i in range(0, self.num_tendons)]
        r = self.effective_depth - csec_diameter / 2.0
        x = csec_centroid.coords.xy[0][0]
        y = csec_centroid.coords.xy[1][0]
        return [LongitudinalTendon(
            diameter=self.diameter, hole_diameter=self.hole_diameter,
            prestressing_force=self.prestressing_force, initial_strain=self.initial_strain,
            steel=self.steel, point=Point(x + r * math.cos(p), y + r * math.sin(p)))
            for p in phi]

    def align_row_in_cross_section(self, csec_centroid=Point(0,0)):
        self.tendon = self.update_bar_coords(csec_centroid, self.csec_diameter)


class ReinforcementElement(object):
    diameter = 0.0
    area = 0.0
    x = 0.0
    y = 0.0

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__

    def update_coordinates(self, point):
        self.point = point  # Point object from the shapely library
        self.coords = point.coords
        self.x = point.coords.xy[0][0]
        self.y = point.coords.xy[1][0]

    def rotate(self, angle):
        angle = angle * math.pi / 180.0  # radians
        self.update_coordinates(Point(
            self.x * math.cos(angle) - self.y * math.sin(angle),
            self.x * math.sin(angle) + self.y * math.cos(angle)
        ))

    def effective_depth(self, csec_y_max):
        return csec_y_max - self.y


class LongitudinalBar(ReinforcementElement):
    def __init__(self, diameter, steel=ReinforcingSteel(), point=Point()):
        self.diameter = diameter
        self.area = math.pi * diameter ** 2 / 4.0
        self.steel = steel
        self.update_coordinates(point)


class LongitudinalTendon(ReinforcementElement):
    def __init__(self, diameter, hole_diameter=0,
                 prestressing_force=0, initial_strain=0,
                 steel=PrestressingSteel(), point=Point()):
        """
        :param diameter: 
        :param hole_diameter: 
        :param prestressing_force: in kN 
        :param initial_strain: positive value, algorithm converts the sign
        :param steel: 
        :param point: 
        """
        self.diameter = diameter
        self.area = math.pi * diameter ** 2 / 4.0
        self.hole_diameter = hole_diameter
        self.hole_area = math.pi * hole_diameter ** 2 / 4.0
        if prestressing_force == 0:
            self.initial_strain = initial_strain
            self.prestressing_force = initial_strain * steel.young_modulus * self.area
            self.prestressing_force *= 1e3  # kN
        elif initial_strain == 0:
            self.prestressing_force = prestressing_force
            self.initial_strain = prestressing_force / (steel.young_modulus * self.area)
            self.initial_strain *= 1e-3     # m/m
        self.steel = steel
        self.update_coordinates(point)