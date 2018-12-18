import material
import math

class Connection(object):
    pass


class EmbeddedSteelRing(Connection):
    pass


class AnchorCage(Connection):

    def __init__(self, bolts_arrangement, flange, upper_anchor_plate, grout, lower_anchor_plate, weight=0):
        self.bolts_arrangement = bolts_arrangement
        self.bolt = bolts_arrangement.bolt
        self.num_bolts = bolts_arrangement.num_bolts
        self.flange = flange
        self.upper_anchor_plate = upper_anchor_plate
        self.grout = grout
        self.lower_anchor_plate = lower_anchor_plate
        self.weight = weight
        self.eval_properties()

    def eval_properties(self):
        self.flange.eval_properties(
            num_bolts=self.bolts_arrangement.num_bolts,
            num_rows=self.bolts_arrangement.num_rows,
            hole_diameter=self.bolt.hole_diameter)
        if self.upper_anchor_plate is not None:
            self.upper_anchor_plate.eval_properties(
                num_bolts=self.bolts_arrangement.num_bolts,
                num_rows=self.bolts_arrangement.num_rows,
                hole_diameter=self.bolt.hole_diameter)
        self.grout.eval_properties(
            num_bolts=self.bolts_arrangement.num_bolts,
            num_rows=self.bolts_arrangement.num_rows,
            hole_diameter=self.bolt.hole_diameter)
        self.lower_anchor_plate.eval_properties(
            num_bolts=self.bolts_arrangement.num_bolts,
            num_rows=self.bolts_arrangement.num_rows,
            hole_diameter=self.bolt.hole_diameter)
        self.average_diameter = self.flange.average_diameter
        self.num_divisions = self.flange.num_divisions
        self.num_rows = self.bolts_arrangement.num_rows
        self.depth_per_division = math.pi * self.average_diameter / self.num_divisions

    def bolt_prestress_force_per_division(self, load_partial_safety_factor=1.1):
        return load_partial_safety_factor * self.bolt.prestress_force * self.num_rows


class BoltsCircularArrangement(object):

    def __init__(self, bolt, num_bolts, outer_bolt_circle_diameter=0, inner_bolt_circle_diameter=0):
        """
        :param bolt: Bolt object or standarized bolt name, i.e., M36-8.8
        :param num_bolts: total number of bolts
        :param outer_bolt_circle_diameter: diameter of the outer circle of bolts, if any
        :param inner_bolt_circle_diameter: diameter of the inner circle of bolts, if any
        """
        self.bolt = bolt
        self.num_bolts = num_bolts
        self.outer_bolt_circle_diameter = outer_bolt_circle_diameter
        self.inner_bolt_circle_diameter = inner_bolt_circle_diameter
        self.bolts_distance_in_division = (outer_bolt_circle_diameter - inner_bolt_circle_diameter) / 2.
        self.num_rows = self.check_num_rows()

    def check_num_rows(self):
        if all([self.outer_bolt_circle_diameter, self.inner_bolt_circle_diameter]) > 0:
            return 2
        elif any([self.outer_bolt_circle_diameter, self.inner_bolt_circle_diameter]) > 0:
            return 1
        else:
            return 0


class Bolt(object):

    # TODO library of bolts M36 --> diamete, shaft area, stress area

    def __init__(self, diameter, shaft_area=None, stress_area=None, effective_free_length=None,
                 hole_diameter=None, steel=material.Steel(), prestress_force=None, washer=None):
        """

        :param diameter: in m
        :param shaft_area: shaft area
        :param stress_area: stress area at the thread of the bolt
        :param effective_free_length: free length, i.e., between bottom of grout and top of embedded flange ring
        :param hole_diameter: if not specified, equal to diameter + 6 mm
        :param steel:
        :param prestress_force:
        :param washer: Washer object
        """
        self.diameter = diameter
        self.shaft_area = shaft_area
        self.stress_area = stress_area
        self.effective_free_length = effective_free_length
        self.hole_diameter = diameter + 0.006 if hole_diameter is None else hole_diameter
        self.hole_area = math.pi * self.hole_diameter ** 2. / 4.
        self.steel = steel
        self.prestress_force = prestress_force
        self.ultimate_prestress_force = steel.fuk * 1e3 * stress_area / 1.1
        self.washer = washer
        if washer is None:
            self.washer = Washer(diameter=diameter+0.05, hole_diameter=diameter+0.01, thickness=0.01, steel=steel)

    def __str__(self):
        return "M%.0f - %s" % (self.diameter * 1e3, self.steel.strength_grade)

    def prestress_ratio(self):
        return self.prestress_force / self.ultimate_prestress_force \
            if self.prestress_force is not None else None


class Washer(object):
    def __init__(self, diameter, hole_diameter, thickness, steel=material.Steel()):
        self.diameter = diameter
        self.hole_diameter = hole_diameter
        self.thickness = thickness
        self.steel = steel


class CircularSteelPlate(object):

    def __init__(self, outer_diameter, inner_diameter, thickness, steel=material.Steel()):
        self.set_plate_parameters(outer_diameter, inner_diameter, thickness, steel)

    def set_plate_parameters(self, outer_diameter, inner_diameter, thickness, steel=material.Steel()):
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        self.thickness = thickness
        self.steel = steel

    def eval_properties(self, num_bolts=0, num_rows=0, hole_diameter=0):
        self.gross_area = math.pi * (self.outer_diameter ** 2. - self.inner_diameter ** 2.) / 4.
        self.net_area = self.gross_area
        self.gross_volume = self.gross_area * self.thickness
        self.net_volume = self.net_area * self.thickness
        self.average_diameter = (self.outer_diameter + self.inner_diameter) / 2.
        self.num_divisions = 0
        self.net_area_per_division = 0.
        self.width_per_division = 0
        self.average_length_per_division = 0
        if num_rows > 0:
            self.num_divisions = (num_bolts / num_rows)
            self.net_area -= num_bolts * math.pi * hole_diameter ** 2. / 4
            self.net_area_per_division = self.net_area / self.num_divisions
            self.width_per_division = (self.outer_diameter - self.inner_diameter) / 2.
            self.average_length_per_division = math.pi * self.average_diameter / self.num_divisions
            self.net_volume = self.net_area * self.thickness
        self.weight = self.net_area * self.steel.specific_weight


class Flange(CircularSteelPlate):

    def __init__(self, outer_diameter, inner_diameter, thickness, steel=material.Steel(), flange_type='T'):
        super(Flange, self).__init__(outer_diameter, inner_diameter, thickness, steel)
        self.type = flange_type


class GroutLayer(object):
    pass


class GroutCircularLayer(GroutLayer):

    def __init__(self, outer_diameter, inner_diameter, thickness, material=material.Concrete(fck=80.)):
        self.outer_diameter = outer_diameter
        self.inner_diameter = inner_diameter
        self.thickness = thickness
        self.material = material

    def eval_properties(self, num_bolts=0, num_rows=0, hole_diameter=0):
        self.gross_area = math.pi * (self.outer_diameter ** 2. - self.inner_diameter ** 2.) / 4.
        self.net_area = self.gross_area
        self.gross_volume = self.gross_area * self.thickness
        self.net_volume = self.net_area * self.thickness
        self.average_diameter = (self.outer_diameter + self.inner_diameter) / 2.
        self.width_per_division = (self.outer_diameter - self.inner_diameter) / 2.
        self.net_area_per_division = 0.
        self.num_divisions = 0
        if num_rows > 0:
            self.num_divisions = (num_bolts / num_rows)
            self.net_area -= num_bolts * math.pi * hole_diameter ** 2. / 4.
            self.net_area_per_division = self.net_area / self.num_divisions
            self.average_length_per_division = math.pi * self.average_diameter / self.num_divisions
            self.net_volume = self.net_area * self.thickness
        self.weight = self.net_area * self.material.specific_weight