"""
Unit tests for VCF CallFilter classes.

These tests validate the expected input/output structures for each CallFilter
using real VCF fixture files.

To focus on a single test (similar to fit() in Mocha):

1. Use unittest.skip decorator on other tests:
   @unittest.skip("Temporarily skipping")
   
2. Run specific test from command line:
   python -m unittest publisher.VCF_test.TestTranscriptsCallFilter
   python -m unittest publisher.VCF_test.TestTranscriptsCallFilter.test_transform_output_structure
   
3. Use pytest with -k flag (if pytest is installed):
   pytest publisher/VCF_test.py -k "Transcripts"
   
4. Use environment variable or attribute (demonstrated below with FOCUS_TEST)
"""

import unittest
import os
from pathlib import Path

# Set to True to enable focus mode - only focused tests will run
FOCUS_MODE = os.environ.get('FOCUS_TEST', 'false').lower() == 'true'

def focus(cls):
    """Decorator to mark a test class as focused. Only runs when FOCUS_MODE=true."""
    cls._focused = True
    return cls

def skipUnlessFocused(cls):
    """Decorator to skip test class unless it's focused or FOCUS_MODE is off."""
    if FOCUS_MODE and not getattr(cls, '_focused', False):
        return unittest.skip("Skipping - not focused")(cls)
    return cls

from publisher.VCF_filters import (
    CallFilter,
    GenesCallFilter,
    TranscriptsCallFilter,
    VariantsCallFilter,
    VariantsTranscriptsCallFilter,
    VariantsAnnotationsCallFilter,
    VariantsConsequencesCallFilter,
    SnvsCallFilter,
    MtsCallFilter,
    GenomicBvlFrequenciesCallFilter,
#    GenomicGnomadFrequenciesCallFilter,
    MtBvlFrequenciesCallFilter,
#    MtGnomadFrequenciesCallFilter,
)

# Helper function to get fixture paths
def get_fixture_path(filename: str) -> str:
    """Get the absolute path to a fixture file."""
    return os.path.join(os.path.dirname(__file__), 'fixtures', 'vcf', filename)

@skipUnlessFocused
class TestBaseFilter(unittest.TestCase):
    
    testInstance = None
    class MockFilter(CallFilter):
            
        def load_vcf_file(self, vcf_file_path: str):
            super().load_vcf_file(vcf_file_path)

        def getTableRows(self):
            return []

    def setUp(self):
        self.testInstance = self.MockFilter(get_fixture_path('mock_snv.vcf'))
        
    def test_instance_creation(self):
        self.assertIsInstance(self.testInstance, CallFilter)
        
    def test_has_records(self):
        self.assertEqual(len(self.testInstance.vcf_records), 2)
        
    def test_csq_getter(self):    
        csq_values = self.testInstance.get_csq_values(self.testInstance.vcf_records[1], 'SYMBOL')
        self.assertEqual(csq_values, ['LA16c-60H5.7', 'NBEAP3'])

@skipUnlessFocused
class TestGenesCallFilter(unittest.TestCase):
    """Test GenesCallFilter."""
    
    # Class variables initialized to None
    filter = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.filter = GenesCallFilter(
            get_fixture_path('mock_snv.vcf'),
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("GenesCallFilter.getTableRows() not yet implemented")
        
        # Check result is a list
        self.assertIsInstance(result, list)
        
        # Check each item has expected keys
        if result:
            self.assertIn('short_name', result[0])


@skipUnlessFocused
@focus
class TestTranscriptsCallFilter(unittest.TestCase):
    """Test TranscriptsCallFilter - FOCUSED for demonstration."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = TranscriptsCallFilter(self.vcf_files)
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("TranscriptsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            # Check expected keys
            self.assertIn('transcript_id', result[0])
            self.assertIn('gene', result[0])
            self.assertIn('transcript_type', result[0])
            self.assertIn('tsl', result[0])
            
            # Check transcript type is encoded (E or R)
            self.assertIn(result[0]['transcript_type'], ['.'])
            self.assertIn(result[1]['transcript_type'], ['E'])


@skipUnlessFocused
class TestVariantsCallFilter(unittest.TestCase):
    """Test VariantsCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = VariantsCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant_id', result[0])
            self.assertIn('var_type', result[0])
            
            # Check variant type matches
            self.assertEqual(result[0]['var_type'], 'SNV')


@skipUnlessFocused
class TestVariantsTranscriptsCallFilter(unittest.TestCase):
    """Test VariantsTranscriptsCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = VariantsTranscriptsCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsTranscriptsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('hgvsc', result[0])
            self.assertEqual(result[0]['hgvsc'], 'ENST00000398242.2:n.402G>C')

@skipUnlessFocused
class TestVariantsAnnotationsCallFilter(unittest.TestCase):
    """Test VariantsAnnotationsCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = VariantsAnnotationsCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsAnnotationsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('hgvsp', result[0])
            self.assertIn('sift', result[0])
            self.assertIn('polyphen', result[0])
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])

@skipUnlessFocused
class TestVariantsConsequencesCallFilter(unittest.TestCase):
    """Test VariantsConsequencesCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = VariantsConsequencesCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsConsequencesCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        if result:
            self.assertIn('severity', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('transcript', result[0])
            
            self.assertEqual(result[0]['severity'], 'FIXME - WHY DOESNT THIS FAIL') #fixme

@skipUnlessFocused
class TestSnvsCallFilter(unittest.TestCase):
    """Test SnvsCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = SnvsCallFilter(
            self.vcf_files,
            assembly='GRCh37'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("SnvsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('type', result[0])
            self.assertIn('chr', result[0])
            self.assertIn('pos', result[0])


@skipUnlessFocused
class TestMtsCallFilter(unittest.TestCase):
    """Test MtsCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        self.skipTest("MtsCallFilter.getTableRows() not yet implemented")
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_mt.vcf')
        self.filter = MtsCallFilter(
            self.vcf_files,
            assembly='GRCh37'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("MtsCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('pos', result[0])
            self.assertIn('ref', result[0])
            self.assertIn('alt', result[0])


@skipUnlessFocused
class TestGenomicBvlFrequenciesCallFilter(unittest.TestCase):
    """Test GenomicBvlFrequenciesCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_snv.vcf')
        self.filter = GenomicBvlFrequenciesCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("GenomicBvlFrequenciesCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('af_tot', result[0])
            self.assertIn('ac_tot', result[0])
            self.assertEquals(result[0]['af_tot'], 0.188889)
            self.assertEquals(result[0]['hom_tot'], 34)

@skipUnlessFocused
class TestMtBvlFrequenciesCallFilter(unittest.TestCase):
    """Test MtBvlFrequenciesCallFilter."""
    
    # Class variables initialized to None
    filter = None
    vcf_files = None
    
    def setUp(self):
        self.skipTest("MtBvlFrequenciesCallFilter.getTableRows() not yet implemented")
        """Set up test fixtures."""
        self.vcf_files = get_fixture_path('mock_mt.vcf')
        self.filter = MtBvlFrequenciesCallFilter(
            self.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = self.filter.getTableRows()
        except NotImplementedError:
            self.skipTest("MtBvlFrequenciesCallFilter.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('an', result[0])
            self.assertIn('ac_hom', result[0])
            self.assertIn('ac_het', result[0])

if __name__ == '__main__':
    unittest.main()
