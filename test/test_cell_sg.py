
import logging
import unittest
from ImageD11 import unitcell, parameters

RHOMB = [5., 5., 5., 70., 70., 70.]
HEX = [4.9, 4.9, 13., 90., 90., 120.]
CUBIC = [4.05, 4.05, 4.05, 90., 90., 90.]


class Catch(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self, logging.WARNING)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


def warnings_from(*args, **kwds):
    """ unitcell(*args, **kwds) and the warnings it logged """
    catch = Catch()
    logging.getLogger().addHandler(catch)
    try:
        cell = unitcell.unitcell(*args, **kwds)
    finally:
        logging.getLogger().removeHandler(catch)
    return cell, catch.messages


def hkls(cell, dsmax=1.0):
    return set(h for d, h in cell.gethkls(dsmax))


class test_cell_sg(unittest.TestCase):

    def test_rhombohedral_axes_from_number(self):
        # 166 on rhombohedral axes : no centring absences
        cell, w = warnings_from(RHOMB, 166)
        self.assertEqual(cell.cell_choice, "rhombohedral")
        self.assertEqual(w, [])
        self.assertEqual(hkls(cell), hkls(unitcell.unitcell(RHOMB, "P", sgname="R-3mr")))

    def test_hexagonal_axes_from_number(self):
        cell, w = warnings_from(HEX, 166)
        self.assertEqual(cell.cell_choice, "standard")
        self.assertEqual(w, [])
        self.assertEqual(hkls(cell), hkls(unitcell.unitcell(HEX, "R")))

    def test_more_specific_name(self):
        cell, w = warnings_from(CUBIC, "F", sgname="Fm-3m")
        self.assertEqual(w, [])
        self.assertEqual(cell.spacegroup, 225)
        self.assertEqual(cell.symmetry, "F")
        self.assertEqual(hkls(cell), hkls(unitcell.unitcell(CUBIC, 225)))

    def test_incompatible_name(self):
        for symmetry in ("I", 229):
            cell, w = warnings_from(CUBIC, symmetry, sgname="Fm-3m")
            self.assertEqual(len(w), 1, symmetry)
            self.assertEqual(cell.spacegroup, 225)  # the name wins
            self.assertEqual(hkls(cell), hkls(unitcell.unitcell(CUBIC, 225)))

    def test_cell_without_the_symmetry(self):
        for cell, sg in [([3., 4., 5., 90., 90., 100.], 14),     # unique axis c
                         ([4.05, 4.06, 4.05, 90., 90., 90.], 225),
                         ([4.9, 4.9, 13., 90., 90., 90.], 166)]:
            self.assertEqual(len(warnings_from(cell, sg)[1]), 1, sg)
        for cell, sg in [([3., 4., 5., 90., 100., 90.], 14), (CUBIC, 225), (HEX, 194)]:
            self.assertEqual(warnings_from(cell, sg)[1], [], sg)
        self.assertEqual(warnings_from([3., 4., 5., 90., 90., 100.], "P")[1], [])

    def test_par_round_trip(self):
        pars = unitcell.unitcell(RHOMB, "P", sgname="R-3mr").to_par_dict()
        self.assertEqual(pars["cell_sg"], "R-3mr")
        for p in (pars, parameters.parameters(**pars)):
            cell = unitcell.unitcell_from_parameters(p)
            self.assertEqual(cell.sgname, "R-3mr")
            self.assertEqual(cell.cell_choice, "rhombohedral")
        # old files without cell_sg
        pars = unitcell.unitcell(CUBIC, 225).to_par_dict()
        self.assertNotIn("cell_sg", pars)
        for p in (pars, parameters.parameters(**pars)):
            self.assertEqual(unitcell.unitcell_from_parameters(p).sgname, None)


if __name__ == "__main__":
    unittest.main()
