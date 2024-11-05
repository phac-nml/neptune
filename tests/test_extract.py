#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2024

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""

import os
import sys
import io

from tests.TestingUtility import *

from neptune.ExtractSignatures import *

import unittest

"""
# =============================================================================

CALCULATE PROBABILITY HOMOLOGOUS BASES MUTATE AND MATCH

# =============================================================================
"""
class TestCalculateProbHBMM(unittest.TestCase):

    """ 
    # =============================================================================

    test_zero

    PURPOSE:
        Tests a GC-content of zero.

    INPUT:
        0: gc = float(0.0)
        1: gc = int(0)

    EXPECTED:
        0: 1.0
        1: 1.0

    # =============================================================================
    """
    def test_zero(self):

        # 0: gc = 0.0 (float)
        result = calculateProbHBMM(float(0.0))
        self.assertAlmostEqual(result, 1.0, 5)

        # 1: gc = 0 (integer)
        result = calculateProbHBMM(int(0))
        self.assertAlmostEqual(result, 1.0, 5)

    """ 
    # =============================================================================

    test_one

    PURPOSE:
        Tests a GC-content of 1.0. 

    INPUT:
        0: gc = float(1.0)
        1: gc = int(1)

    EXPECTED:
        0: 1.0
        1: 1.0

    # =============================================================================
    """
    def test_one(self):
            
        # 0: gc = 1.0 (float)
        result = calculateProbHBMM(float(1.0))
        self.assertAlmostEqual(result, 1.0, 5)

        # 1: gc = 1 (integer)
        result = calculateProbHBMM(int(1))
        self.assertAlmostEqual(result, 1.0, 5)

    """ 
    # =============================================================================

    test_50

    PURPOSE:
        Tests a GC-content of 0.50.

    INPUT:
        0: gc = 0.50

    EXPECTED:
        0: 1/3

    # =============================================================================
    """
    def test_50(self):

        # 0: gc = 0.50
	    result = calculateProbHBMM(0.50)
	    self.assertAlmostEqual(result, 1.0/3.0, 5)

    """ 
    # =============================================================================

    test_25

    PURPOSE:
        Tests a GC-content of 0.25.

    INPUT:
        0: gc = 0.25
        1: gc = 0.75

    EXPECTED:
        0: 0.426939
        1: 0.426939

    # =============================================================================
    """
    def test_25(self):

        # 0: gc = 0.25
        result25 = calculateProbHBMM(0.25)

        # 1: gc = 0.75
        result75 = calculateProbHBMM(0.75)

        self.assertAlmostEqual(result25, 0.426939, 5)
        self.assertAlmostEqual(result75, 0.426939, 5)
        self.assertAlmostEqual(result25, result75, 5)

    """ 
    # =============================================================================

    test_10

    PURPOSE:
        Tests a GC-content of 0.10.

    INPUT:
        0: gc = 0.10
        1: gc = 0.90

    EXPECTED:
        0: 0.662508
        1: 0.662508

    # =============================================================================
    """
    def test_10(self):

        # 0: gc = 0.10
        result10 = calculateProbHBMM(0.10)

        # 1: gc = 0.90
        result90 = calculateProbHBMM(0.90)

        self.assertAlmostEqual(result10, 0.662508, 5)
        self.assertAlmostEqual(result90, 0.662508, 5)
        self.assertAlmostEqual(result10, result90, 5)

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the bounds of the function.

    INPUT:
        0: gc = 1.1
        1: gc = -0.1

    EXPECTED:
        0: RuntimeError
        1: RuntimeError

    # =============================================================================
    """
    def test_bounds(self):

	    # 0: gc = 1.1 -- high gc
	    with self.assertRaises(RuntimeError):
		    result = calculateProbHBMM(1.1)

	    # 1: gc = -0.1 -- low gc
	    with self.assertRaises(RuntimeError):
		    result = calculateProbHBMM(-0.1)

"""
# =============================================================================

CALCULATE PROBABILITY HOMOLOGOUS BASES MATCH

# =============================================================================
"""
class TestCalculateProbHBM(unittest.TestCase):

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the bounds of the function.

    INPUT:
        0: mutationRate = 0.01, gc = 1.1
        1: mutationRate = 0.01, gc = -0.1
        2: mutationRate = 1.1, gc = 0.5
        3: mutationRate = -0.1, gc = 0.5

    EXPECTED:
        0: RunetimeError
        1: RunetimeError
        2: RunetimeError
        3: RunetimeError

    # =============================================================================
    """
    def test_bounds(self):
	
        # 0: mutationRate = 0.01, gc = 1.1 -- gc high
        with self.assertRaises(RuntimeError):
            result = calculateProbHBM(0.01, 1.1)

        # 1: mutationRate = 0.01, gc = -0.1 -- gc low
        with self.assertRaises(RuntimeError):
            result = calculateProbHBM(0.01, -0.1)

        # 2: mutationRate = 1.1, gc = 0.5 -- mutation high
        with self.assertRaises(RuntimeError):
            result = calculateProbHBM(1.1, 0.5)

        # 3: mutationRate = -0.1, gc = 0.5 -- mutation low
        with self.assertRaises(RuntimeError):
            result = calculateProbHBM(-0.1, 0.5)

    """ 
    # =============================================================================

    test_standard

    PURPOSE:
        Tests a standard example.

    INPUT:
        0: mutationRate = 0.01, gc = 0.5

    EXPECTED:
        0: 0.980133

    # =============================================================================
    """
    def test_standard(self):

	    result = calculateProbHBM(0.01, 0.5)
	    self.assertAlmostEqual(result, 0.980133, 5)

    """ 
    # =============================================================================

    test_high_gc

    PURPOSE:
        Tests a high, but acceptable, GC-content.

    INPUT:
        0: mutation = 0.01, gc = 0.9

    EXPECTED:
        0: 0.980166

    # =============================================================================
    """
    def test_high_gc(self):

        result = calculateProbHBM(0.01, 0.9)
        self.assertAlmostEqual(result, 0.980166, 5)

    """ 
    # =============================================================================

    test_high_mutation

    PURPOSE:
        Tests a high, but acceptable, mutation rate.

    INPUT:
        0: mutation = 0.5, gc = 0.5

    EXPECTED:
        0: 0.333333

    # =============================================================================
    """
    def test_high_mutation(self):

        result = calculateProbHBM(0.5, 0.5)
        self.assertAlmostEqual(result, 0.333333, 5)

"""
# =============================================================================

CALCULATE PROBABILITY HOMOLOGOUS K-MERS MATCH

# =============================================================================
"""
class TestCalculateProbHKM(unittest.TestCase):

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the bounds of the function.

    INPUT:
        0: mutation = 0.1, gc = 1.1, k = 21
        1: mutation = 0.1, gc = -0.1, k = 21
        2: mutation = 1.1, gc = 0.5, k = 21
        3: mutation = -0.1, gc = 0.5, k = 21
        4: mutation = 0.1, gc = 0.5, k = 0

    EXPECTED:
        0: RuntimeError
        1: RuntimeError
        2: RuntimeError
        3: RuntimeError
        4: RuntimeError

    # =============================================================================
    """
    def test_bounds(self):

        # 0: mutation = 0.1, gc = 1.1, k = 21 -- gc high
        with self.assertRaises(RuntimeError):
            result = calculateProbHKM(0.01, 1.1, 21)

        # 1: mutation = 0.1, gc = -0.1, k = 21 -- gc low
        with self.assertRaises(RuntimeError):
            result = calculateProbHKM(0.01, -0.1, 21)

        # 2: mutation = 1.1, gc = 0.5, k = 21 -- mutation high
        with self.assertRaises(RuntimeError):
            result = calculateProbHKM(1.1, 0.5, 21)

        # 3: mutation = -0.1, gc = 0.5, k = 21 -- mutation low
        with self.assertRaises(RuntimeError):
            result = calculateProbHKM(-0.1, 0.5, 21)

        # 4: mutation = 0.1, gc = 0.5, k = 0 -- k low
        with self.assertRaises(RuntimeError):
            result = calculateProbHKM(0.01, 0.5, 0)

    """ 
    # =============================================================================

    test_standard

    PURPOSE:
        Tests a standard example.

    INPUT:
        0: mutation = 0.1, gc = 0.5, k = 21 

    EXPECTED:
        0: 0.656128

    # =============================================================================
    """
    def test_standard(self):

        # 0: mutation = 0.1, gc = 0.5, k = 21 
        result = calculateProbHKM(0.01, 0.5, 21)
        self.assertAlmostEqual(result, 0.656128, 5)

    """ 
    # =============================================================================

    test_high_gc

    PURPOSE:
        Tests a high, but acceptable, GC-content.

    INPUT:
        0: mutation = 0.1, gc = 0.9, k = 21 

    EXPECTED:
        0: 0.656591

    # =============================================================================
    """
    def test_high_gc(self):

        # 0: mutation = 0.1, gc = 0.9, k = 21
        result = calculateProbHKM(0.01, 0.9, 21)
        self.assertAlmostEqual(result, 0.656591, 5)

    """ 
    # =============================================================================

    test_high_mutation

    PURPOSE:
        Tests a high, but acceptable, mutation rate.

    INPUT:
        0: mutation = 0.5, gc = 0.5, k = 21

    EXPECTED:
        0: 0.0

    # =============================================================================
    """
    def test_high_mutation(self):

        # 0: mutation = 0.5, gc = 0.5, k = 21
        result = calculateProbHKM(0.5, 0.5, 21)
        self.assertAlmostEqual(result, 0.000000, 5)

    """ 
    # =============================================================================

    test_high_k

    PURPOSE:
        Tests a large k-mer size.

    INPUT:
        0: mutation = 0.01, gc = 0.5, k = 51

    EXPECTED:
        0: 0.359371

    # =============================================================================
    """
    def test_high_k(self):

        # 0: mutation = 0.01, gc = 0.5, k = 51
	    result = calculateProbHKM(0.01, 0.5, 51)
	    self.assertAlmostEqual(result, 0.359371, 5)

    """ 
    # =============================================================================

    test_low_k

    PURPOSE:
        Tests a small k-mer size.

    INPUT:
        0: mutation = 0.01, gc = 0.5, k = 5

    EXPECTED:
        0: 0.904536

    # =============================================================================
    """
    def test_low_k(self):

        # 0: mutation = 0.01, gc = 0.5, k = 5
        result = calculateProbHKM(0.01, 0.5, 5)
        self.assertAlmostEqual(result, 0.904536, 5)

"""
# =============================================================================

ESTIMATE GAP SIZE

# =============================================================================
"""
class TestEstimateGapSize(unittest.TestCase):

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the bounds of the function.

    INPUT:
        confidence = 0.95

        0: mutation = 0.1, gc = 1.1, k = 21
        1: mutation = 0.1, gc = -0.1, k = 21
        2: mutation = 1.1, gc = 0.5, k = 21
        3: mutation = -0.1, gc = 0.5, k = 21
        4: mutation = 0.1, gc = 0.5, k = 0

    EXPECTED:
        0: RuntimeError
        1: RuntimeError
        2: RuntimeError
        3: RuntimeError
        4: RuntimeError

    # =============================================================================
    """
    def test_bounds(self):

        # 0: mutation = 0.1, gc = 1.1, k = 21 -- gc high
        with self.assertRaises(RuntimeError):
            result = estimateGapSize(0.01, 1.1, 21, 0.95)

        # 1: mutation = 0.1, gc = -0.1, k = 21 -- gc low
        with self.assertRaises(RuntimeError):
            result = estimateGapSize(0.01, -0.1, 21, 0.95)

        # 2: mutation = 1.1, gc = 0.5, k = 21 -- mutation high
        with self.assertRaises(RuntimeError):
            result = estimateGapSize(1.1, 0.5, 21, 0.95)

        # 3: mutation = -0.1, gc = 0.5, k = 21 -- mutation low
        with self.assertRaises(RuntimeError):
            result = estimateGapSize(-0.1, 0.5, 21, 0.95)

        # 4: mutation = 0.1, gc = 0.5, k = 0 -- k low
        with self.assertRaises(RuntimeError):
            result = estimateGapSize(0.01, 0.5, 0, 0.95) 

    """ 
    # =============================================================================

    test_trivial

    PURPOSE:
        Tests a trivial example.

    INPUT:
        0: mutation = 0.01, gc = 0.5, k = 5

    EXPECTED:
        0: 11

    # =============================================================================
    """
    def test_trivial(self):

        # 0: mutation = 0.01, gc = 0.5, k = 5
        mutationRate = 0.01
        GC = 0.5
        kmerSize = 5
        confidence = 0.95

        # p = 0.980133; from rate & GC
        # u = 5.31239
        # var^2 = 1.23319
        # std = 1.11049
        # CONFIDENCE = 4.472 # 95% (Chebyshev's Inequality)

        expected = 11
        result = estimateGapSize(mutationRate, GC, kmerSize, confidence)

        self.assertEqual(expected, result)

"""
# =============================================================================

EXTRACT

# =============================================================================
"""
class TestExtract(unittest.TestCase):

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the bounds of the function.

    INPUT:
        0: references = None, (...)
        1: k = 0, (...)
        2: inmers = None, (...)
        3: exmers = None, (...)
        4: size = 0, (...)
        5: gap = 0, (...)
        6: output = None, (...)
        7: normal

    EXPECTED:
        0: RuntimeError
        1: RuntimeError
        2: RuntimeError
        3: RuntimeError
        4: RuntimeError
        5: RuntimeError
        6: RuntimeError
        7: AAAAA

    # =============================================================================
    """
    def test_bounds(self):

        references = {}
        references["1"] = "CCCCCAAAAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["CCA"] = 1
        inmers["CAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        output = io.StringIO()

        # 0: references = None, (...) -- no references
        with self.assertRaises(RuntimeError):
            result = extract(None, k, inmers, exmers, size, gap, output) 

        # k = 0, (...) -- k low
        with self.assertRaises(RuntimeError):
            result = extract(references, 0, inmers, exmers, size, gap, output)

        # inmers = None, (...) -- no inclusion kmers
        with self.assertRaises(RuntimeError):
            result = extract(references, k, None, exmers, size, gap, output)

        # exmers = None, (...) -- no exclusion kmers
        with self.assertRaises(RuntimeError):
            result = extract(references, k, inmers, None, size, gap, output)

        # size = 0, (...) -- size low
        with self.assertRaises(RuntimeError):
            result = extract(references, k, inmers, exmers, 0, gap, output)

        # gap = 0, (...) -- gap low
        with self.assertRaises(RuntimeError):
            result = extract(references, k, inmers, exmers, size, 0, output)

        # output = None, (...) -- no output
        with self.assertRaises(RuntimeError):
            result = extract(references, k, inmers, exmers, size, gap, None)

        # normal
        output = io.StringIO()
        result = extract(references, k, inmers, exmers, size, gap, output)
        lines = output.getvalue().split("\n")

        self.assertEqual(lines[0], ">0 score=0.0000 in=0.0000 ex=0.0000 len=5 ref=1 pos=5")
        self.assertEqual(lines[1], "AAAAA")

        output.close()

    """ 
    # =============================================================================

    test_trivial

    PURPOSE:
        Tests a trivial example.

    INPUT:
        0:

        references["1"] = "CCCCCAAAAACCCCC"

        inmers:
            AAA 1
            AAC 1
            ACC 1
            CAA 1
            CCA 1

        exmers:
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0: "AAAAA"

    # =============================================================================
    """
    def test_trivial(self):

        references = {}
        references["1"] = "CCCCCAAAAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1
        inmers["CAA"] = 1
        inmers["CCA"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        expected = ("AAAAA")
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        signature = output.getvalue().split("\n")[1]
        
        self.assertEqual(signature, expected)

        output.close()

    """ 
    # =============================================================================

    test_small_gap

    PURPOSE:
        Tests extraction with a small gap.

    INPUT:
        0:

        references["1"] = "CCCCCAAAAATAAAAACCCCC"

        inmers:
            AAA 1
            AAC 1
            ACC 1
            CAA 1
            CCA 1

        exmers:
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0: "AAAAATAAAAA"

    # =============================================================================
    """
    def test_small_gap(self):

        references = {}
        references["1"] = "CCCCCAAAAATAAAAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1
        inmers["CAA"] = 1
        inmers["CCA"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        expected = ("AAAAATAAAAA")
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        signature = output.getvalue().split("\n")[1]        

        self.assertEqual(signature, expected)

        output.close()

    """ 
    # =============================================================================

    test_large_gap

    PURPOSE:
        Tests extraction with a large gap. The large gap should break the
        signature into two signatures.

    INPUT:
        0:

        references["1"] = "CCCCCAAAAATATATATAAAAACCCCC"

        inmers:
            AAA 1
            AAC 1
            ACC 1
            CAA 1
            CCA 1

        exmers:
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0:
        
        "AAA"
        "AAA"

    # =============================================================================
    """
    def test_large_gap(self):

        references = {}
        references["1"] = "CCCCCAAAAATATATATAAAAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1
        inmers["CAA"] = 1
        inmers["CCA"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        expected1 = ("AAA")
        expected2 = ("AAA")
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        signature1 = output.getvalue().split("\n")[1]
        signature2 = output.getvalue().split("\n")[3]    

        self.assertEqual(signature1, expected1)
        self.assertEqual(signature2, expected2)

        output.close()

    """ 
    # =============================================================================

    test_start

    PURPOSE:
        Tests that extraction finds a signature at the start of a reference.

    INPUT:
        0:

        references["1"] = "AAAAACCCCC"

        inmers:
            AAA 1
            AAC 1
            ACC 1
            CAA 1
            CCA 1

        exmers:
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0: "AAA"

    # =============================================================================
    """
    def test_start(self):

        references = {}
        references["1"] = "AAAAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1
        inmers["CAA"] = 1
        inmers["CCA"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        expected = ("AAA")
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        signature = output.getvalue().split("\n")[1]
        
        self.assertEqual(signature, expected)

        output.close()

    """ 
    # =============================================================================

    test_end

    PURPOSE:
        Tests the extraction finds a signature at the end of a reference.

    INPUT:
        0:

        references["1"] = "CCCCCAAAAA"

        inmers:
            AAA 1
            AAC 1
            ACC 1
            CAA 1
            CCA 1

        exmers:
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0: "AAA"

    # =============================================================================
    """
    def test_end(self):

        references = {}
        references["1"] = "CCCCCAAAAA"

        k = 3

        inmers = {}
        inmers["AAA"] = 1
        inmers["AAC"] = 1
        inmers["ACC"] = 1
        inmers["CAA"] = 1
        inmers["CCA"] = 1

        exmers = {}
        exmers["CCC"] = 1

        size = 2
        gap = 4

        expected = ("AAA")
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        signature = output.getvalue().split("\n")[1]
        
        self.assertEqual(signature, expected)

        output.close()

    """ 
    # =============================================================================

    test_no_signature

    PURPOSE:
        Tests that extraction doesn't find a signature when there is no
        appropriate signature.

    INPUT:
        0:

        references["1"] = "CCCCCAACCCCC"

        inmers:
            AAA 1

        exmers:
            AAC 1
            ACC 1
            CAA 1
            CCA 1
            CCC 1

        size = 2, gap = 4

    EXPECTED:
        0: NONE

    # =============================================================================
    """
    def test_no_signature(self):

        references = {}
        references["1"] = "CCCCCAACCCCC"

        k = 3

        inmers = {}
        inmers["AAA"] = 1

        exmers = {}
        exmers["CCC"] = 1
        exmers["CCA"] = 1
        exmers["CAA"] = 1
        exmers["AAC"] = 1
        exmers["ACC"] = 1

        size = 2
        gap = 4
        
        output = io.StringIO()
        extract(references, k, inmers, exmers, size, gap, output)
        
        self.assertEqual(output.getvalue(), "")

        output.close()

"""
# =============================================================================

ESTIMATE SIGNATURE SIZE

# =============================================================================
"""
class TestEstimateSignatureSize(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0: k = 51

    EXPECTED:
        0: 204

    # =============================================================================
    """
    def test_simple(self):

        k = 51

        expected = 204
        result = estimateSignatureSize(k)

        self.assertEqual(result, expected)

    """ 
    # =============================================================================

    test_zero

    PURPOSE:
        Tests when k is zero. 

    INPUT:
        0: k = 0

    EXPECTED:
        0: 0

    # =============================================================================
    """
    def test_zero(self):

        k = 0

        expected = 0
        result = estimateSignatureSize(k)

        self.assertEqual(result, expected)

"""
# =============================================================================

ESTIMATE EXCLUSION HITS

# =============================================================================
"""
class TestEstimateExclusionHits(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple use case.

    INPUT:
        0: totalExclusion = 20, rate = 0.01, k = 11

    EXPECTED:
        0: 1

    # =============================================================================
    """
    def test_simple(self):

        totalExclusion = 20
        rate = 0.01
        kmerSize = 11

        expected = 1
        result = estimateExclusionHits(totalExclusion, rate, kmerSize)

        self.assertEqual(result, expected)

"""
# =============================================================================

ESTIMATE INCLUSION HITS

# =============================================================================
"""
class TestEstimateInclusionHits(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        totalInclusion = 100
        mutation = 0.01
        gc = 0.5
        k = 21
        confidence = 0.95

    EXPECTED:
        57

    # =============================================================================
    """
    def test_simple(self):

        totalInclusion = 100
        mutationRate = 0.01
        GC = 0.5
        kmerSize = 21
        confidence = 0.95

        # p = probHKM = 0.656128
        # mean = 64.9567
        # variance = 22.3368
        # stdev = 4.72618

        # estimate = 1 + 64.9567 - 1.645 * 4.72618 = 58.1821339 -> 58
        expected = 58

        result = estimateInclusionHits(
            totalInclusion, mutationRate, GC, kmerSize, confidence)

        self.assertEqual(result, expected)

"""
# =============================================================================

ESTIMATE K

# =============================================================================
"""
class TestEstimateK(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0: "AAA 1 1"

    EXPECTED:
        0: 3

    # =============================================================================
    """
    def test_simple(self):

        buff = io.StringIO()
        buff.write("AAA 1 1\n")
        buff.seek(0)

        result = estimateK(buff)
        expected = 3

        self.assertEqual(result, expected)

        buff.close()

"""
# =============================================================================

BUILD KMERS

# =============================================================================
"""
class TestBuildKMers(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0: 

        AAA 3 4
        ACA 4 3
        CAA 3 4
        CCA 3 3

    EXPECTED:
        0: 

        inmers:
        AAA 3
        ACA 4
        CAA 3
        CCA 3

        exmers:
        AAA 4
        ACA 3
        CAA 4
        CCA 3

    # =============================================================================
    """
    def test_simple(self):

        buff = io.StringIO()
        buff.write("AAA 3 4\n")
        buff.write("ACA 4 3\n")
        buff.write("CAA 3 4\n")
        buff.write("CCA 3 3\n")
        buff.seek(0)

        inmers = {}
        exmers = {}
        inhits = 1
        exhits = 1

        result = buildKMers(buff, inmers, exmers, inhits, exhits)
        
        expected_inmers = {'AAA': 3, 'ACA': 4, 'CAA': 3, 'CCA': 3}
        self.assertEqual(inmers, expected_inmers)

        expected_exmers = {'AAA': 4, 'ACA': 3, 'CAA': 4, 'CCA': 3}
        self.assertEqual(exmers, expected_exmers)

        buff.close()

    """ 
    # =============================================================================

    test_no_inclusion

    PURPOSE:
        Tests when there are no inclusion k-mers.

    INPUT:
        0: 

        AAA 0 4
        ACA 0 3
        CAA 0 4
        CCA 0 3

    EXPECTED:
        0:

        inmers:
        NONE

        exmers:
        AAA 4
        ACA 3
        CAA 4
        CCA 3

    # =============================================================================
    """
    def test_no_inclusion(self):

        buff = io.StringIO()
        buff.write("AAA 0 4\n")
        buff.write("ACA 0 3\n")
        buff.write("CAA 0 4\n")
        buff.write("CCA 0 3\n")
        buff.seek(0)

        inmers = {}
        exmers = {}
        inhits = 1
        exhits = 1

        result = buildKMers(buff, inmers, exmers, inhits, exhits)
        
        expected_inmers = {}
        self.assertEqual(inmers, expected_inmers)

        expected_exmers = {'AAA': 4, 'ACA': 3, 'CAA': 4, 'CCA': 3}
        self.assertEqual(exmers, expected_exmers)

        buff.close()

    """ 
    # =============================================================================

    test_no_exclusion

    PURPOSE:
        Tests when there are no exclusion k-mers.

    INPUT:
        0: 

        AAA 0 4
        ACA 0 3
        CAA 0 4
        CCA 0 3

    EXPECTED:
        0:

        inmers:
        AAA 3
        ACA 4
        CAA 3
        CCA 3

        exmers:
        NONE

    # =============================================================================
    """
    def test_no_exclusion(self):

        buff = io.StringIO()
        buff.write("AAA 3 0\n")
        buff.write("ACA 4 0\n")
        buff.write("CAA 3 0\n")
        buff.write("CCA 3 0\n")
        buff.seek(0)

        inmers = {}
        exmers = {}
        inhits = 1
        exhits = 1

        result = buildKMers(buff, inmers, exmers, inhits, exhits)
        
        expected_inmers = {'AAA': 3, 'ACA': 4, 'CAA': 3, 'CCA': 3}
        self.assertEqual(inmers, expected_inmers)

        expected_exmers = {}
        self.assertEqual(exmers, expected_exmers)

        buff.close()

    """ 
    # =============================================================================

    test_bounds

    PURPOSE:
        Tests the input bounds of the function. This is testing edge cases.

    INPUT:
        0:

        AAA 1 1
        ACA 2 2
        CAA 3 3
        CCA 4 4

        inhits = 3
        exhits = 2

    EXPECTED:
        0:

        inmers:
        CAA 3
        CCA 4

        exmers:
        ACA 2
        CAA 3
        CCA 4

    # =============================================================================
    """
    def test_bounds(self):

        buff = io.StringIO()
        buff.write("AAA 1 1\n")
        buff.write("ACA 2 2\n")
        buff.write("CAA 3 3\n")
        buff.write("CCA 4 4\n")
        buff.seek(0)

        inmers = {}
        exmers = {}
        inhits = 3
        exhits = 2

        result = buildKMers(buff, inmers, exmers, inhits, exhits)
        
        expected_inmers = {'CAA': 3, 'CCA': 4}
        self.assertEqual(inmers, expected_inmers)

        expected_exmers = {'ACA': 2, 'CAA': 3, 'CCA': 4}
        self.assertEqual(exmers, expected_exmers)

        buff.close()

"""
# =============================================================================

REPORT PARAMETERS

# =============================================================================
"""
class TestReportParameters(unittest.TestCase):

    """ 
    # =============================================================================

    test_run

    PURPOSE:
        Tests the report function causes no errors. This test does not contain
        asserts because of the verbosity of the output.

    INPUT:
        0: 

    EXPECTED:
        0: [NO CRASH]

    # =============================================================================
    """
    def test_run(self):

        buff = io.StringIO()

        reportParameters(buff, "A.fasta", 12, 0.01, 
		    3, 3, 2, 2, 3, "some.kmers", 5, 6, 0.5)

        buff.close()

"""
# =============================================================================

MAIN

# =============================================================================
"""
class TestMain(unittest.TestCase):

    """ 
    # =============================================================================

    test_simple

    PURPOSE:
        Tests a simple example.

    INPUT:
        0:

        simple.fasta:
        >0
        ACGTACGTACGT

        alternative.fasta:
        >0
        ATATATATATAT

        simple.kmers:
        ACGTA 4 0
        CGTAC 4 0

    EXPECTED:
        0:

        >0 score=0.0000 in=0.0000 ex=0.0000 len=4 ref=0 pos=4
        ACGT

    # =============================================================================
    """
    def test_simple(self):

        outputLocation = getPath("tests/output/extract/temp.out")

        sys.argv[1:] = [
            REFERENCE_LONG, "tests/data/extract/simple.fasta",
            INCLUSION_LONG, "tests/data/extract/simple.fasta",
            EXCLUSION_LONG, "tests/data/extract/alternative.fasta",
            KMERS_LONG, "tests/data/extract/simple.kmers",
            SIZE_LONG, "4",
            OUTPUT_LONG, outputLocation
            ]

        main()

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ">0 score=0.0000 in=0.0000 ex=0.0000 len=4 ref=0 pos=4\nACGT\n"
            self.assertEqual(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_no_reference

    PURPOSE:
        Tests there is a RuntimeError when the reference file does not exist.

    INPUT:
        0:

        reference = NONE

        simple.fasta:
        >0
        ACGTACGTACGT

        alternative.fasta:
        >0
        ATATATATATAT

        simple.kmers:
        ACGTA 4 0
        CGTAC 4 0

    EXPECTED:
        0: RuntimeError

    # =============================================================================
    """
    def test_no_reference(self):

        outputLocation = getPath("tests/output/extract/temp.out")

        sys.argv[1:] = [
            REFERENCE_LONG, "tests/data/extract/DOES_NOT_EXIST.FAKE",
            INCLUSION_LONG, "tests/data/extract/simple.fasta",
            EXCLUSION_LONG, "tests/data/extract/alternative.fasta",
            KMERS_LONG, "tests/data/extract/simple.kmers",
            OUTPUT_LONG, outputLocation
            ]

        with self.assertRaises(RuntimeError):
            main()

    """ 
    # =============================================================================

    test_

    PURPOSE:
        Tests there is a RuntimeError when the k-mer file does not exist.

    INPUT:
        0:

        simple.fasta:
        >0
        ACGTACGTACGT

        alternative.fasta:
        >0
        ATATATATATAT

        k-mers: NONE

    EXPECTED:
        0: RuntimeError

    # =============================================================================
    """
    def test_no_kmers(self):

        outputLocation = getPath("tests/output/temp.out")

        sys.argv[1:] = [
            REFERENCE_LONG, "tests/data/extract/simple.fasta",
            INCLUSION_LONG, "tests/data/extract/simple.fasta",
            EXCLUSION_LONG, "tests/data/extract/alternative.fasta",
            KMERS_LONG, "tests/data/extract/DOES_NOT_EXIST.FAKE",
            OUTPUT_LONG, outputLocation
            ]

        with self.assertRaises(RuntimeError):
            main()

    """ 
    # =============================================================================

    test_specify_parameters

    PURPOSE:
        Tests specifying many parameters.

    INPUT:
        0:

        simple.fasta:
        >0
        ACGTACGTACGT

        alternative.fasta:
        >0
        ATATATATATAT

        simple.kmers:
        ACGTA 4 0
        CGTAC 4 0

        refSize = 12, rate = 0.01, inhits = 1, exhits = 1, gap = 3, size = 5, gc = 0.5

    EXPECTED:
        0: 

        >0 score=0.0000 in=0.0000 ex=0.0000 len=4 ref=0 pos=4
        ACGT
    # =============================================================================
    """
    def test_specify_parameters(self):

        outputLocation = getPath("tests/output/extract/temp.out")

        sys.argv[1:] = [
            REFERENCE_LONG, "tests/data/extract/simple.fasta",
            INCLUSION_LONG, "tests/data/extract/simple.fasta",
            EXCLUSION_LONG, "tests/data/extract/alternative.fasta",
            KMERS_LONG, "tests/data/extract/simple.kmers",
            OUTPUT_LONG, outputLocation,
            REFERENCE_SIZE_LONG, "12",
            RATE_LONG, "0.01",
            INHITS_LONG, "1",
            EXHITS_LONG, "1",
            GAP_LONG, "3",
            SIZE_LONG, "4",
            GC_LONG, "0.5"
            ]

        main()

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ">0 score=0.0000 in=0.0000 ex=0.0000 len=4 ref=0 pos=4\nACGT\n"
            self.assertEqual(result, expected)

        os.remove(outputLocation)

    """ 
    # =============================================================================

    test_random

    PURPOSE:
        Tests extraction on a randomly-generated sequence.

    INPUT:
        0:

        >random
        GGTTGGTGC TCTAAACTTCAT TTGGTTA

        AAACTTC 1 0
        AAATGAA 1 0
        AACCAAA 1 0
        AACTTCA 1 0
        AAGTTTA 1 0
        AATGAAG 1 0
        ACCAAAT 1 0
        ACCAACC 1 1
        ACTTCAT 1 0
        AGAGCAC 1 0
        AGCACCA 1 0
        AGTTTAG 1 0
        CAAATGA 1 0
        CACCAAC 1 1
        CATTTGG 1 0
        CTCTAAA 1 0
        GAGCACC 1 0
        GCACCAA 1 1
        GCTCTAA 1 0
        GTTTAGA 1 0
        TAACCAA 1 1
        TAGAGCA 1 0

        refSize = 28, rate = 0.01, inhits = 1, exhits = 1, gap = 3, size = 5, gc = 0.5

    EXPECTED:
        0:

        >0 score=0.0000 in=0.0000 ex=0.0000 len=12 ref=random pos=9
        TCTAAACTTCAT

    # =============================================================================
    """
    def test_random(self):

        outputLocation = getPath("tests/output/temp.out")

        sys.argv[1:] = [
            REFERENCE_LONG, "tests/data/random.fasta",
            INCLUSION_LONG, "tests/data/random.fasta",
            EXCLUSION_LONG, "tests/data/alternative.fasta",
            KMERS_LONG, "tests/data/random.kmers",
            OUTPUT_LONG, outputLocation,
            REFERENCE_SIZE_LONG, "28",
            RATE_LONG, "0.01",
            INHITS_LONG, "1",
            EXHITS_LONG, "1",
            GAP_LONG, "3",
            SIZE_LONG, "5",
            GC_LONG, "0.5"
            ]

        main()

        with open (outputLocation, "r") as myfile:

            result = myfile.read()
            expected = ">0 score=0.0000 in=0.0000 ex=0.0000 len=12 ref=random pos=9\nTCTAAACTTCAT\n"
            self.assertEqual(result, expected)

        os.remove(outputLocation)

if __name__ == '__main__':
	unittest.main()
