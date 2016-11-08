import unittest
from .. import pls


class TestInput(unittest.TestCase):
    def testInput(self):
        y_col = "lreJB_kpc"
        x1_col = "lsigma_cor"
        x2_col = "lIeJB_cor"
        fits_table, total_galaxies = pls.read_clusters("comafit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        fits_table1, total_galaxies1 = pls.read_clusters("comafit.fits", False, "coma", "group", y_col, x1_col, x2_col)
        self.assertNotEquals(fits_table, fits_table1)
        self.assertEquals(total_galaxies, total_galaxies1)
        self.assertEqual(fits_table["group"], fits_table1["group"])


class TestProgram(unittest.TestCase):
    def setUp(self):
        y_col = "lreJB_kpc"
        x1_col = "lsigma_cor"
        x2_col = "lIeJB_cor"
        coma_table_three, coma_galaxies_three = pls.read_clusters("comafit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        coma_table_two, coma_galaxies_two = pls.read_clusters("comafit.fits", False, "coma", "group", y_col, x1_col, x2_col)
        y_col = "lreJB_kpc_DEV"
        x1_col = "lsigma_cor"
        x2_col = "lIeJB_DEV"
        rxj0152_table, rxj0152_galaxies = pls.read_clusters("rxj0152allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        rxj1226_table, rxj1226_galaxies = pls.read_clusters("rxj1226allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        y_col = "lML_JB_DEV"
        x1_col = "lMass_DEV"
        rxj1226_table, rxj1226_galaxies = pls.read_clusters("rxj1226allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)

    def testZeropoint(self):
        return 0

    def testResiduals(self):
        return 0

    def testBootstrap(self):
        return 0

    def testChangeCoefficients(self):
        return 0

    def testUncertainty(self):
        return 0

    def testMinimization(self):
        return 0

    def testTfitlin(self):
        return 0