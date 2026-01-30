"""
CallFilter base class and global VCF options.
"""


# Global VCF options (imported from constants)
from constants import HYPEN_VARIANT_NOTATION, CHR_NOTATION, NA

import vcfpy
import os
import logging

import gzip
import io
from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)

class CallFilter(ABC):
    vcf_records: List[vcfpy.Record]
    
    def __init__(self, vcf_file_path: str):
        self.vcf_records = []
        self.csq_fields = []
        self.csq_index_map = {}
        self.severity_map = {}
        
        child_class_name = self.__class__.__name__
        logger.info("booting up %s", child_class_name)
        
        #read severity table file
        severity_table_path = os.path.join(os.path.dirname(__file__), "../severities.tsv")
        try:
            with open(severity_table_path, "r") as f:
                for line in f.readlines()[1:]:
                    parts = line.strip().split("\t")
                    if len(parts) == 3:
                        severity, consequence = parts
                        self.severity_map[consequence] = int(severity)
        except FileNotFoundError:
            logger.warning("Severity table file not found: %s", severity_table_path)
        
        self.load_vcf_file(vcf_file_path)
    
    def load_vcf_file(self, file, type = "SNV"):
        
        def read_vcf(reader: vcfpy.Reader):
                
            if type == "SNV":
            
                csq = reader.header.get_info_field_info("CSQ")
                csq_elements = csq.description.split("Format: ")[1]
                self.csq_fields = csq_elements.split("|")
                self.csq_index_map = {field: index for index, field in enumerate(self.csq_fields)}
                
                for record in reader:
                    try:
                        self.vcf_records.append(record)
                    except Exception as e:
                        logger.error(f"Error processing record {l} in file {file}: {e}")
                        
            elif type == "MT":
                # ????? TBA
                
                self.csq_index_map = {field: index for index, field in enumerate(self.csq_fields)}

                for record in reader:    
                    self.vcf_records.append(record)
                
        if file.endswith('.gz'):
            with gzip.open(file, 'rb') as gz:
                with io.TextIOWrapper(gz, encoding='utf-8', errors='replace') as f:  # or errors='ignore'
                    vcf_reader = vcfpy.Reader(stream=f)
                    read_vcf(vcf_reader)
                        
        else:
            with io.TextIOWrapper(open(file, 'rb'), encoding='utf-8', errors='replace') as f:
                vcf_reader = vcfpy.Reader(stream=f)
                read_vcf(vcf_reader)
        
    
    def get_csq_values(self, record: vcfpy.Record, field_name: str) -> List[str]:
        """
        Helper method to extract a specific CSQ field value from a VCF record.
        
        Args:
            record: VCF record object
            field_name: Name of the CSQ field to extract
        """
        index = self.csq_index_map.get(field_name)
        values = []
        if index is None:
            return []
        csq_list = record.INFO.get("CSQ", [])
        if not csq_list:
            return []
        for list in csq_list:
            csq_parts = list.split("|")
            if index >= len(csq_parts):
                return []
            else:
                values.append(csq_parts[index])
        return values
    
    def get_info_value(self, record: vcfpy.Record, field_name: str, fallback = None) -> str:
        """
        Helper method to extract a specific INFO field value from a VCF record.
        
        Args:
            record: VCF record object
            field_name: Name of the INFO field to extract
        """
        return record.INFO.get(field_name, fallback)

    def make_variant_id(self, record: vcfpy.Record) -> str:
        """
        Helper method to construct a variant ID from VCF record fields.
        
        Args:
            record: VCF record object
        """
        if HYPEN_VARIANT_NOTATION:
            if CHR_NOTATION:
                chrom = record.CHROM
            else:
                chrom = record.CHROM.replace("chr", "")
            pos = record.POS
            ref = record.REF
            alt = record.ALT[0].value  # assuming single ALT allele
            variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        else:
            variant_id = record.ID[0]
        return variant_id
    @abstractmethod
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract table rows from the loaded VCF data.
        
        Returns:
            List of dictionaries where each dict represents a row in the table.
        """
        pass