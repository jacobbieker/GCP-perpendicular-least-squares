import unittest
import pls


class TestInput(unittest.TestCase):
    def testInput(self):
        y_col = "lreJB_kpc".upper()
        x1_col = "lsigma_cor".upper()
        x2_col = "lIeJB_cor".upper()
        fits_table, total_galaxies = pls.read_clusters(["comafit.fits"], True, "GALAXY", "GROUP", y_col, x1_col, x2_col)
        fits_table1, total_galaxies1 = pls.read_clusters(["comafit.fits"], False, "GALAXY", "GROUP", y_col, x1_col, x2_col)
        #self.assertNotEquals(fits_table["GALAXY"], fits_table1["GALAXY"])
        self.assertEquals(total_galaxies, total_galaxies1)
        #self.assertEqual(fits_table["GROUP"], fits_table1["GROUP"])


class TestProgram(unittest.TestCase):
    def setUp(self):
        '''
        y_col = "lreJB_kpc_DEV"
        x1_col = "lsigma_cor"
        x2_col = "lIeJB_DEV"
        rxj0152_table, rxj0152_galaxies = pls.read_clusters("rxj0152allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        rxj1226_table, rxj1226_galaxies = pls.read_clusters("rxj1226allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        y_col = "lML_JB_DEV"
        x1_col = "lMass_DEV"
        rxj1226_table, rxj1226_galaxies = pls.read_clusters("rxj1226allfit.fits", True, "coma", "group", y_col, x1_col, x2_col)
        '''
    def testZeropoint(self):
        cluster = 1
        y_col = "lreJB_kpc".upper()
        x1_col = "lsigma_cor".upper()
        x2_col = "lIeJB_cor".upper()
        type_solution = ["median", "mean"]
        res_choice = ["per", "x1", "x2"]
        a_factor = 0.05
        b_factor = 0.02
        coma_table_three, coma_galaxies_three = pls.read_clusters(["comafit.fits"], True, "GALAXY", "GROUP", y_col, x1_col, x2_col)
        #coma_table_two, coma_galaxies_two = pls.read_clusters(["comafit.fits"], False, "GALAXY", "GROUP", y_col, x1_col, x2_col)

        results_1 = pls.zeropoint(coma_table_three, cluster, type_solution[0], res_choice[0], y_col, x1_col, x2_col, a_factor, b_factor, True)
        #results_2 = pls.zeropoint(coma_table_two, cluster, type_solution[1], res_choice[1], y_col, x1_col, x2_col, a_factor, b_factor)

        #self.assertNotEqual(results_1, results_2)
        #self.assertEqual(coma_table_three['res'], coma_table_two['res'])

    def testResiduals(self):
        cluster = 1
        y_col = "lreJB_kpc".upper()
        x1_col = "lsigma_cor".upper()
        x2_col = "lIeJB_cor".upper()
        type_solution = ["median", "mean"]
        res_choice = ["per", "x1", "x2"]
        a_factor = 0.05
        b_factor = 0.02
        coma_table_three, coma_galaxies_three = pls.read_clusters(["comafit.fits"], True, "GALAXY", "GROUP", y_col, x1_col, x2_col)
        #coma_table_two, coma_galaxies_two = pls.read_clusters(["comafit.fits"], False, "GALAXY", "GROUP", y_col, x1_col, x2_col)

        results_1 = pls.zeropoint(coma_table_three, cluster, type_solution[0], res_choice[0], y_col, x1_col, x2_col, a_factor, b_factor, True)
        #results_2 = pls.zeropoint(coma_table_two, cluster, type_solution[1], res_choice[1], y_col, x1_col, x2_col, a_factor, b_factor)

        #residuals_1 = pls.residuals(results_1[4], results_1[0], results_1[1], results_1[2], results_1[3])
        #residuals_2 = pls.residuals(results_2[4], results_2[0], results_2[1], results_2[2], results_2[3])

        #self.assertNotEqual(residuals_1['res'], None)
        #self.assertNotEqual(residuals_2['res'], None)

        #self.assertNotEqual(residuals_1['r' + str(results_1[1])], not None)
        #self.assertNotEqual(residuals_2['r' + str(results_2[1])], not None)

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

    def testSetup(self):
        #pls.startup(filename="comafit.fits", y_col="lreJB_kpc", x1_col="lsigma_cor", x2_col="lIeJB_cor",
        #            factor_change_a=1.1860, factor_change_b=-0.8020, zeropoint_choice="median")
        pls.startup(filename="comafit.fits", y_col="lML_JB", x1_col="lMass")