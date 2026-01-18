"""
Unit tests for VCF CallFilter classes.

These tests validate the expected input/output structures for each CallFilter
using real VCF fixture files.

To focus on a single test (similar to fit() in Mocha):

1. Use unittest.skip decorator on other tests:
   @unittest.skip("Temporarily skipping")
   
2. Run specific test from command line:
   python -m unittest publisher.VCF_test_transformers.TestTranscriptsTransformer
   python -m unittest publisher.VCF_test_transformers.TestTranscriptsTransformer.test_transform_output_structure
   
3. Use pytest with -k flag (if pytest is installed):
   pytest publisher/VCF_test_transformers.py -k "Transcripts"
   
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

from publisher.VCF_transformers import (
    GenesTransformer,
    TranscriptsTransformer,
    VariantsTransformer,
    VariantsTranscriptsTransformer,
    VariantsAnnotationsTransformer,
    VariantsConsequencesTransformer,
    SnvsTransformer,
    MtsTransformer,
    GenomicIbvlFrequenciesTransformer,
    GenomicGnomadFrequenciesTransformer,
    MtIbvlFrequenciesTransformer,
    MtGnomadFrequenciesTransformer,
)


# Helper function to get fixture paths
def get_fixture_path(filename: str) -> str:
    """Get the absolute path to a fixture file."""
    return str(Path(__file__).parent / 'fixtures' / 'vcf' / filename)


@skipUnlessFocused
class TestGenesTransformer(unittest.TestCase):
    """Test GenesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestGenesTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestGenesTransformer.transformer = GenesTransformer(TestGenesTransformer.vcf_files)
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestGenesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("GenesTransformer.getTableRows() not yet implemented")
        
        # Check result is a list
        self.assertIsInstance(result, list)
        
        # Check each item has expected keys
        if result:
            self.assertIn('short_name', result[0])


@skipUnlessFocused
@focus
class TestTranscriptsTransformer(unittest.TestCase):
    """Test TranscriptsTransformer - FOCUSED for demonstration."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestTranscriptsTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestTranscriptsTransformer.transformer = TranscriptsTransformer(TestTranscriptsTransformer.vcf_files)
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestTranscriptsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("TranscriptsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            # Check expected keys
            self.assertIn('transcript_id', result[0])
            self.assertIn('gene', result[0])
            self.assertIn('transcript_type', result[0])
            self.assertIn('tsl', result[0])
            
            # Check transcript type is encoded (E or R)
            self.assertIn(result[0]['transcript_type'], ['E', 'R'])


@skipUnlessFocused
class TestVariantsTransformer(unittest.TestCase):
    """Test VariantsTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestVariantsTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestVariantsTransformer.transformer = VariantsTransformer(
            TestVariantsTransformer.vcf_files, 
            variant_type='SNV'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestVariantsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant_id', result[0])
            self.assertIn('var_type', result[0])
            
            # Check variant type matches
            self.assertEqual(result[0]['var_type'], 'SNV')


@skipUnlessFocused
class TestVariantsTranscriptsTransformer(unittest.TestCase):
    """Test VariantsTranscriptsTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestVariantsTranscriptsTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestVariantsTranscriptsTransformer.transformer = VariantsTranscriptsTransformer(
            TestVariantsTranscriptsTransformer.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestVariantsTranscriptsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsTranscriptsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('hgvsc', result[0])


@skipUnlessFocused
class TestVariantsAnnotationsTransformer(unittest.TestCase):
    """Test VariantsAnnotationsTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestVariantsAnnotationsTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestVariantsAnnotationsTransformer.transformer = VariantsAnnotationsTransformer(
            TestVariantsAnnotationsTransformer.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestVariantsAnnotationsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsAnnotationsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('hgvsp', result[0])
            self.assertIn('sift', result[0])
            self.assertIn('polyphen', result[0])
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])


@skipUnlessFocused
class TestVariantsConsequencesTransformer(unittest.TestCase):
    """Test VariantsConsequencesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    severity_table_path = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestVariantsConsequencesTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestVariantsConsequencesTransformer.severity_table_path = get_fixture_path('severity_table.tsv')
        TestVariantsConsequencesTransformer.transformer = VariantsConsequencesTransformer(
            TestVariantsConsequencesTransformer.vcf_files,
            TestVariantsConsequencesTransformer.severity_table_path
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestVariantsConsequencesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("VariantsConsequencesTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('severity', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('transcript', result[0])


@skipUnlessFocused
class TestSnvsTransformer(unittest.TestCase):
    """Test SnvsTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestSnvsTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestSnvsTransformer.transformer = SnvsTransformer(
            TestSnvsTransformer.vcf_files,
            assembly='GRCh37'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestSnvsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("SnvsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('type', result[0])
            self.assertIn('chr', result[0])
            self.assertIn('pos', result[0])


@skipUnlessFocused
class TestMtsTransformer(unittest.TestCase):
    """Test MtsTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    gnomad_file = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestMtsTransformer.vcf_files = [get_fixture_path('mock_mt.vcf')]
        TestMtsTransformer.gnomad_file = get_fixture_path('gnomad_mt.tsv')
        TestMtsTransformer.transformer = MtsTransformer(
            TestMtsTransformer.vcf_files,
            TestMtsTransformer.gnomad_file,
            assembly='GRCh37'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestMtsTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("MtsTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('pos', result[0])
            self.assertIn('ref', result[0])
            self.assertIn('alt', result[0])


@skipUnlessFocused
class TestGenomicIbvlFrequenciesTransformer(unittest.TestCase):
    """Test GenomicIbvlFrequenciesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestGenomicIbvlFrequenciesTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestGenomicIbvlFrequenciesTransformer.transformer = GenomicIbvlFrequenciesTransformer(
            TestGenomicIbvlFrequenciesTransformer.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestGenomicIbvlFrequenciesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("GenomicIbvlFrequenciesTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('af_tot', result[0])
            self.assertIn('ac_tot', result[0])


@skipUnlessFocused
class TestGenomicGnomadFrequenciesTransformer(unittest.TestCase):
    """Test GenomicGnomadFrequenciesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    gnomad_file = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestGenomicGnomadFrequenciesTransformer.vcf_files = [get_fixture_path('mock_snv.vcf')]
        TestGenomicGnomadFrequenciesTransformer.gnomad_file = get_fixture_path('gnomad_snv.tsv')
        TestGenomicGnomadFrequenciesTransformer.transformer = GenomicGnomadFrequenciesTransformer(
            TestGenomicGnomadFrequenciesTransformer.vcf_files,
            TestGenomicGnomadFrequenciesTransformer.gnomad_file,
            assembly='GRCh37'
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestGenomicGnomadFrequenciesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("GenomicGnomadFrequenciesTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('af_tot', result[0])
            self.assertIn('ac_tot', result[0])


@skipUnlessFocused
class TestMtIbvlFrequenciesTransformer(unittest.TestCase):
    """Test MtIbvlFrequenciesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestMtIbvlFrequenciesTransformer.vcf_files = [get_fixture_path('mock_mt.vcf')]
        TestMtIbvlFrequenciesTransformer.transformer = MtIbvlFrequenciesTransformer(
            TestMtIbvlFrequenciesTransformer.vcf_files
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestMtIbvlFrequenciesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("MtIbvlFrequenciesTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('an', result[0])
            self.assertIn('ac_hom', result[0])
            self.assertIn('ac_het', result[0])


@skipUnlessFocused
class TestMtGnomadFrequenciesTransformer(unittest.TestCase):
    """Test MtGnomadFrequenciesTransformer."""
    
    # Class variables initialized to None
    transformer = None
    vcf_files = None
    gnomad_mt_file = None
    
    def setUp(self):
        """Set up test fixtures."""
        TestMtGnomadFrequenciesTransformer.vcf_files = [get_fixture_path('mock_mt.vcf')]
        TestMtGnomadFrequenciesTransformer.gnomad_mt_file = get_fixture_path('gnomad_mt.tsv')
        TestMtGnomadFrequenciesTransformer.transformer = MtGnomadFrequenciesTransformer(
            TestMtGnomadFrequenciesTransformer.vcf_files,
            TestMtGnomadFrequenciesTransformer.gnomad_mt_file
        )
    
    def test_getTableRows_output_structure(self):
        """Test that getTableRows returns correct structure."""
        try:
            result = TestMtGnomadFrequenciesTransformer.transformer.getTableRows()
        except NotImplementedError:
            self.skipTest("MtGnomadFrequenciesTransformer.getTableRows() not yet implemented")
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant', result[0])
            self.assertIn('an', result[0])
            self.assertIn('ac_hom', result[0])
            self.assertIn('ac_het', result[0])


if __name__ == '__main__':
    unittest.main()
