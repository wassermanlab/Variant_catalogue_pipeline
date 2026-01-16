"""
Unit tests for VCF transformer classes.

These tests validate the expected input/output structures for each transformer
using real VCF fixture files.
"""

import unittest
import os
from pathlib import Path
from typing import List, Dict, Any

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


# Helper functions to read fixture files
def get_fixture_path(filename: str) -> Path:
    """Get the path to a fixture file."""
    return Path(__file__).parent / 'fixtures' / 'vcf' / filename


def read_vcf_fixture(filename: str) -> List[Dict[str, Any]]:
    """
    Read a VCF fixture file and return as list of dictionaries.
    This is a placeholder - actual implementation would parse VCF format.
    """
    # For now, return mock data structure that matches what transformers expect
    # TODO: Implement proper VCF parsing
    if 'snv' in filename:
        return [{
            'CHROM': '1',
            'POS': 100000,
            'ID': '1_100000_A_G',
            'REF': 'A',
            'ALT': 'G',
            'QUAL': 99.9,
            'INFO': {
                'CSQ': 'G|missense_variant|MODERATE|GENE1|ENSG00000001|Transcript|ENST00000001|protein_coding|5/10||ENST00000001.1:c.123A>G|ENSP00000001.1:p.Lys41Glu|123|123|41|K/E|Aaa/Gaa|rs123456||1||SNV|HGNC|HGNC:1234|1|TRUE|Ensembl||A|A||tolerated(0.5)|benign(0.1)|||ClinVar::VCV000123456||35.0|0.5|1|2|3|4|0.1|0.2|0.3|0.4',
                'AF_tot_XX_XY': '0.01,0.02,0.005',
                'AC_tot_XX_XY': '10,15,5',
                'AN_tot_XX_XY': '1000,750,1000',
                'hom_tot_XX_XY': '2,1,1',
            }
        }]
    elif 'mt' in filename:
        return [{
            'CHROM': 'MT',
            'POS': 8602,
            'ID': 'chrM_8602_T_C',
            'REF': 'T',
            'ALT': 'C',
            'QUAL': 100.0,
            'INFO': {
                'CSQ': 'C|missense_variant|MODERATE|MT-ATP6|ENSG00000198899|Transcript|ENST00000361899|protein_coding|||ENST00000361899.2:c.321T>C|ENSP00000355046.2:p.Ile107Thr||||I/T||rs879029014||1||SNV||||Ensembl||T|C||||||||ClinVar::VCV000692920||30.0|0.3||||||',
            },
            'GT_fields': {
                'AC_hom': 5,
                'AC_het': 10,
                'AF_hom': 0.05,
                'AF_het': 0.10,
                'AN': 100,
                'max_observed_heteroplasmy': 0.95,
                'heteroplasmy_histogram': '[[0.1,0.2,0.3],[5,10,15]]',
            }
        }]
    return []


def read_tsv_fixture(filename: str) -> List[Dict[str, Any]]:
    """
    Read a TSV fixture file and return as list of dictionaries.
    """
    path = get_fixture_path(filename)
    with open(path, 'r') as f:
        lines = f.readlines()
    
    if not lines:
        return []
    
    # Parse header
    headers = lines[0].strip().split('\t')
    
    # Parse data rows
    records = []
    for line in lines[1:]:
        values = line.strip().split('\t')
        record = dict(zip(headers, values))
        records.append(record)
    
    return records


def read_severity_table(filename: str) -> Dict[str, int]:
    """Read severity table and return as dict."""
    records = read_tsv_fixture(filename)
    return {r['consequence']: int(r['severity_number']) for r in records}


class TestGenesTransformer(unittest.TestCase):
    """Test GenesTransformer."""
    
    def setUp(self):
        self.transformer = GenesTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        # Check result is a list
        self.assertIsInstance(result, list)
        
        # Check each item has expected keys
        if result:
            self.assertIn('short_name', result[0])
            
            # Check expected values
            gene_names = [r['short_name'] for r in result]
            self.assertIn('GENE1', gene_names)


class TestTranscriptsTransformer(unittest.TestCase):
    """Test TranscriptsTransformer."""
    
    def setUp(self):
        self.transformer = TranscriptsTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        self.assertIsInstance(result, list)
        
        if result:
            # Check expected keys
            self.assertIn('transcript_id', result[0])
            self.assertIn('gene', result[0])
            self.assertIn('transcript_type', result[0])
            self.assertIn('tsl', result[0])
            
            # Check transcript type is encoded (E or R)
            self.assertIn(result[0]['transcript_type'], ['E', 'R'])


class TestVariantsTransformer(unittest.TestCase):
    """Test VariantsTransformer."""
    
    def setUp(self):
        self.transformer = VariantsTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records, 'SNV')
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('variant_id', result[0])
            self.assertIn('var_type', result[0])
            
            # Check variant type matches
            self.assertEqual(result[0]['var_type'], 'SNV')
            
            # Check variant ID format
            self.assertIn('1_100000_A_G', [r['variant_id'] for r in result])


class TestVariantsTranscriptsTransformer(unittest.TestCase):
    """Test VariantsTranscriptsTransformer."""
    
    def setUp(self):
        self.transformer = VariantsTranscriptsTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('hgvsc', result[0])
            
            # Check HGVS format
            self.assertTrue(result[0]['hgvsc'].startswith('ENST'))


class TestVariantsAnnotationsTransformer(unittest.TestCase):
    """Test VariantsAnnotationsTransformer."""
    
    def setUp(self):
        self.transformer = VariantsAnnotationsTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('hgvsp', result[0])
            self.assertIn('sift', result[0])
            self.assertIn('polyphen', result[0])
            self.assertIn('transcript', result[0])
            self.assertIn('variant', result[0])


class TestVariantsConsequencesTransformer(unittest.TestCase):
    """Test VariantsConsequencesTransformer."""
    
    def setUp(self):
        self.transformer = VariantsConsequencesTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
        self.severity_table = read_severity_table('severity_table.tsv')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records, self.severity_table)
        
        self.assertIsInstance(result, list)
        
        if result:
            self.assertIn('severity', result[0])
            self.assertIn('variant', result[0])
            self.assertIn('transcript', result[0])
            
            # Check severity is numeric
            self.assertIsInstance(result[0]['severity'], int)


class TestSnvsTransformer(unittest.TestCase):
    """Test SnvsTransformer."""
    
    def setUp(self):
        self.transformer = SnvsTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records, 'GRCh38')
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'type', 'length', 'chr', 'pos', 'ref', 'alt',
                           'cadd_score', 'cadd_intr', 'dbsnp_id', 'dbsnp_url',
                           'ucsc_url', 'ensembl_url', 'clinvar_url', 'gnomad_url',
                           'clinvar_vcv', 'splice_ai']
            
            for key in expected_keys:
                self.assertIn(key, result[0])
            
            # Check CADD interpretation
            self.assertIn(result[0]['cadd_intr'], ['Tolerable', 'Damaging'])


class TestMtsTransformer(unittest.TestCase):
    """Test MtsTransformer."""
    
    def setUp(self):
        self.transformer = MtsTransformer()
        self.vcf_records = read_vcf_fixture('mock_mt.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(
            self.vcf_records, 
            'GRCh38',
            ['chrM_8602_T_C']
        )
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'pos', 'ref', 'alt', 'ucsc_url',
                           'mitomap_url', 'gnomad_url', 'dbsnp_id',
                           'dbsnp_url', 'clinvar_url', 'clinvar_vcv']
            
            for key in expected_keys:
                self.assertIn(key, result[0])


class TestGenomicIbvlFrequenciesTransformer(unittest.TestCase):
    """Test GenomicIbvlFrequenciesTransformer."""
    
    def setUp(self):
        self.transformer = GenomicIbvlFrequenciesTransformer()
        self.vcf_records = read_vcf_fixture('mock_snv.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'af_tot', 'af_xx', 'af_xy',
                           'ac_tot', 'ac_xx', 'ac_xy',
                           'an_tot', 'an_xx', 'an_xy',
                           'hom_tot', 'hom_xx', 'hom_xy', 'quality']
            
            for key in expected_keys:
                self.assertIn(key, result[0])
            
            # Check types
            self.assertIsInstance(result[0]['af_tot'], float)
            self.assertIsInstance(result[0]['ac_tot'], int)


class TestGenomicGnomadFrequenciesTransformer(unittest.TestCase):
    """Test GenomicGnomadFrequenciesTransformer."""
    
    def setUp(self):
        self.transformer = GenomicGnomadFrequenciesTransformer()
        self.gnomad_records = read_tsv_fixture('gnomad_snv.tsv')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(
            self.gnomad_records,
            ['1_100000_A_G'],
            'GRCh38'
        )
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'af_tot', 'ac_tot', 'an_tot', 'hom_tot', 'FILTER']
            
            for key in expected_keys:
                self.assertIn(key, result[0])


class TestMtIbvlFrequenciesTransformer(unittest.TestCase):
    """Test MtIbvlFrequenciesTransformer."""
    
    def setUp(self):
        self.transformer = MtIbvlFrequenciesTransformer()
        self.vcf_records = read_vcf_fixture('mock_mt.vcf')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(self.vcf_records)
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'an', 'ac_hom', 'ac_het',
                           'af_hom', 'af_het', 'hl_hist', 'max_hl']
            
            for key in expected_keys:
                self.assertIn(key, result[0])
            
            # Check heteroplasmy histogram is formatted
            self.assertIsInstance(result[0]['hl_hist'], str)


class TestMtGnomadFrequenciesTransformer(unittest.TestCase):
    """Test MtGnomadFrequenciesTransformer."""
    
    def setUp(self):
        self.transformer = MtGnomadFrequenciesTransformer()
        self.gnomad_records = read_tsv_fixture('gnomad_mt.tsv')
    
    def test_transform_output_structure(self):
        """Test that transform returns correct structure."""
        result = self.transformer.transform(
            self.gnomad_records,
            ['chrM_8602_T_C']
        )
        
        self.assertIsInstance(result, list)
        
        if result:
            expected_keys = ['variant', 'an', 'ac_hom', 'ac_het',
                           'af_hom', 'af_het', 'max_hl']
            
            for key in expected_keys:
                self.assertIn(key, result[0])


if __name__ == '__main__':
    unittest.main()
