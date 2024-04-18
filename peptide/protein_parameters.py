import numpy as np
import matplotlib.pyplot as plt

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    compositionAnalyzed = False

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        self.originalInput = protein
        self.protein, self.gappedProtein = self.clean()
        self.aaCounts = {aa:0 for aa in self.aa2mw.keys()}
        self.aaComposition()

    def clean (self):
        """Cleans the initially input protein sequence.

        Args:
            None

        Returns:
            cleanedProtein: The cleaned protein sequence
            gappedProtein: The gapped protein sequence

        Raises:
            None
        """
        # Gapped protein for alignment purposes
        gappedProtein = "".join([residue if residue in self.aa2mw else "-"
                                 for residue in self.originalInput.upper()])
        # Actual cleaned protein sequence
        cleanedProtein = "".join([residue
                                  for residue in gappedProtein if residue != "-"])
        # print(self.originalInput)
        # print(gappedProtein)
        # print(cleanedProtein)
        return cleanedProtein, gappedProtein

    def aaCount (self):
        """Find the number of amino acids in the protein.

        Args:
            None

        Returns:
            The length of the protein sequence.

        Raises:
            None
        """

        return len(self.protein)

    def pI (self, algorithm = "naive"):
        """Computes the theoretical isolelectric point of the protein.

        Args:
            algorithm: Which algorithm to use. Can be
                "naive" (brute force - default)
                "binary" (binary search)
                "gd" (gradient decent) - not working properly
                "metropolis" (Metropolis-Hastings) - not implemented
        """
        # Used for selection of which algorithm to use
        pIAlogirthms = {
            "naive": self.naivepI(),
            "binary": self.binarySearchpI(),
            "gd": self._gradientDecentpI_(), # FIXME - recalculate gradient
            "metropolis": self._metropolisHastingspI_() # FIXME - needs to be implemented
            }
        pH, charges, pHValues = pIAlogirthms[algorithm]
        # self.plotpHandCharge(pH, charges, pHValues)
        return pH

    def plotpHandCharge(self, minChargepH, charges, pHValues):
        """Plots pH versus charge and marks the location of the isoelectric point.

        Args:
            minChargepH: The pH at which the lowest charge was found.
            charges: A list of charges corresponding to each queried pH.
            pHValues: A list of pHs corresponding to each resulting charge.

        Returns:
            None

        Raises:
            None

        """
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 6))  # 1 row, 2 columns, and optional figure size

        # True Charge
        ax0.plot(pHValues, charges)
        ax0.set_xlabel("pH")
        ax0.set_ylabel("Charge")
        ax0.axhline(y=0, color='orange', linestyle='dotted')
        ax0.axvline(x=minChargepH, color='orange', linestyle='dotted')
        ax0.text(minChargepH, 0, f' pH={minChargepH:.2f}',
                 verticalalignment='bottom',
                 horizontalalignment='left',
                 color='black')

        # Absolute Charges
        absCharges = [abs(charge) for charge in charges]
        ax1.plot(pHValues, absCharges)
        ax1.set_xlabel("pH")
        ax1.set_ylabel("Absolute Charge")
        ax1.axhline(y=0, color='orange', linestyle='dotted')
        ax1.axvline(x=minChargepH, color='orange', linestyle='dotted')
        ax1.text(minChargepH, 0, f' pH={minChargepH:.2f}',
                 verticalalignment='bottom',
                 horizontalalignment='left',
                 color='black')

        plt.show()


    def naivepI(self, interval=0.01):
        """Loops through a range of pH values and finds the pH at which
        the charge is closest to zero.

        Args:
            interval: The interal with which to query the pH.

        Returns:
            The pH at which the charge is closest to zero.

        Raises:
            None
        """

        charges = []
        pHValues = []
        for i in list(np.arange(0.0,14.0,interval)):
            charge = self._charge_(i)
            charges.append(charge)
            pHValues.append(i)

        minCharge = min(charges, key=abs) # Find the minimum charge closest to zero
        minChargeIndex = charges.index(minCharge) # Get the index of this charge
        pH = pHValues[minChargeIndex] # Retrieve the minimum pH

        return pH, charges, pHValues

    def _gradientDecentpI_(self, initialpH=7.0, learningRate=0.01,
                     tolerance=0.001, maxIterations=1000):
        """Uses gradient decent to find the pH at which the charge is
        closest to zero.

        Args:
            initialpH: The initial pH to start at.
            learningRate: The rate at which to update the pH.
            tolerance: The tolerance with which to stop the algorithm.
            maxIterations: The maximum number of iterations to run.

        Returns:
            The pH at which the charge is closest to zero.

        Raises:
            None
        """
        pH = initialpH
        charges = []
        pHValues = []

        for i in range(maxIterations):
            current_charge = self._charge_(pH)
            gradient = self._chargeGradient(pH)
            nextpH = pH - learningRate * gradient

            if abs(current_charge - self._charge_(nextpH)) < tolerance:
                break

            pH = nextpH
            charges.append(current_charge)
            pHValues.append(pH)

        return pH, charges, pHValues

    def _metropolisHastingspI_(self):
        """Metropolis Hastings algorithm to find the pH at which the charge is
        closest to zero.

        NOT IMPLEMENTED!

        Args:
            None

        Returns:
            The pH at which the charge is closest to zero.

        Raises:
            None

        """
        pass

    def _chargeGradient(self, pH): # FIXME: Not sure this is correctly calculated.
        """Computes the gradient of the charge with respect to pH.

        Args:
            pH: The pH at which to calculate the gradient.

        Returns:
            The gradient of the charge with respect to pH.

        Raises:
            None
        """

        logBase10 = np.log(10)

        # The first derivative of the first sum of the
        positiveGradient = sum([
            -self.aaCounts[aa] * (10 ** self.aa2chargePos[aa])
            * logBase10 * (10 ** pH) /
            ((10 ** self.aa2chargePos[aa] + 10 ** pH) ** 2)
            for aa in ['K', 'R', 'H']
        ])

        # Adding N-terminus contribution
        positiveGradient += -logBase10 * (10 ** pH) * (10 ** self.aaNterm) / \
                            ((10 ** self.aaNterm + 10 ** pH) ** 2)

        # Gradient for the negatively charged amino acids
        negativeGradient = sum([
            self.aaCounts[aa] * (10 ** pH) * logBase10
            / (10 ** self.aa2chargeNeg[aa] + 10 ** pH) -
            ((10 ** pH) ** 2) * logBase10
            / ((10 ** self.aa2chargeNeg[aa] + 10 ** pH) ** 2)
            for aa in ['D', 'E', 'C', 'Y']
        ])

        # Adding C-terminus contribution
        negativeGradient += (10 ** pH) * logBase10 / (10 ** self.aaCterm + 10 ** pH) - \
                            ((10 ** pH) ** 2) * logBase10 / ((10 ** self.aaCterm + 10 ** pH) ** 2)

        net_gradient = positiveGradient + negativeGradient

        return net_gradient

    def binarySearchpI(self, initialpH=7, iterations = 1000):
        """Uses binary search to find the pH at which the charge is
        closest to zero.

        Args:
            initialpH: The initial pH to start at.
            iterations: The number of iterations to run.

        Returns:
            The pH at which the charge is closest to zero.

        Raises:
            None
        """

        pH = initialpH
        charges = []
        pHValues = []
        for i in range(iterations):
            charge = self._charge_(pH)
            charges.append(charge)
            pHValues.append(pH)
            if charge > 0:
                pH += pH / 2
            elif charge < 0:
                pH -= pH / 2
            else:
                return pH, charges, pHValues
        return pH, charges, pHValues

    def aaComposition (self) :
        """Computes the composition of amino acids in the protein.

        Args:
            None

        Returns:
            A dictionary containing the composition of amino acids
            in the protein.

        Raises:
            None
        """
        for aa in self.protein:
            self.aaCounts.update({aa:self.aaCounts[aa]+1})

        self.compositionAnalyzed = True
        return self.aaCounts


    def _charge_(self, pH=7.0):
        """Computes the net charge of the protein.

        Args:
            pH: The pH at which to calculate the charge (default is 7.0).


        Returns:
            The net charge of the protein.

        Raises:
            None
        """

        # Positive charge contributions from basic amino acids (K, R, H)
        netPositiveCharge = sum([self.aaCounts[aa] * (10 ** self.aa2chargePos[aa]) /
                                (10 ** self.aa2chargePos[aa] + 10 ** pH)
                                for aa in ['K', 'R', 'H']])

        # Adding the N-terminus charge contribution
        netPositiveCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)

        # Negative charge contributions from acidic amino acids (D, E, C, Y)
        netNegativeCharge = sum([self.aaCounts[aa] * 10 ** pH /
                                (10 ** self.aa2chargeNeg[aa] + 10 ** pH)
                                for aa in ['D', 'E', 'C', 'Y']])

        # Adding the C-terminus charge contribution
        netNegativeCharge += 10 ** pH / (10 ** self.aaCterm + 10 ** pH)

        # Net charge is the difference between positive and negative charges
        return netPositiveCharge - netNegativeCharge

    def molarExtinction (self, cystine=True):
        """Computes the molar extinction coefficient of the protein.

        Args:
            cystine: Whether to include cystine in the calculation (default is True).

        Returns:
            The molar extinction coefficient of the protein.

        Raises:
            None
        """

        if cystine:
            return sum([self.aaCounts[aa]
                        * self.aa2abs280[aa] for aa in ['Y', 'W', 'C']])
        else:
            return sum([self.aaCounts[aa]
                        * self.aa2abs280[aa] for aa in ['Y', 'W']])

    def massExtinction (self, cystine=True):
        """Computes the mass extinction coefficient of the protein.

        Args:
            cystine: Whether to include cystine in the calculation (default is True).

        Returns:
            The mass extinction coefficient of the protein.

        Raises:
            None
        """
        if cystine:
            myMW =  self.molecularWeight()
            return self.molarExtinction() / myMW if myMW else 0.0
        else:
            myMW =  self.molecularWeight(True)
            return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self, customSequence = False):
        """Computes the molecular weight of the protein.

        Args:
            customSequence: A custom sequence to use for the calculation (default is False).

        Returns:
            The molecular weight of the protein.

        Raises:
            None
        """
        if customSequence != None:
            return sum([self.aa2mw[residue] for residue in self.protein]) \
                    - self.mwH2O * (self.aaCount() - 1)
        else:
            return sum([self.aa2mw[residue] for residue in self.protein.replace('C','')]) \
                    - self.mwH2O * (self.aaCount() - 1)

    def _abs280_ (self):
        pass



# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        # print ("net charge: {:.2f}".format(myParamMaker._charge_()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI("naive")))
        print ("Amino acid composition:")

        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present

        for aa,n in sorted(myParamMaker.aaComposition().items(),
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))

        inString = input('protein sequence?')

if __name__ == "__main__":
     main()
