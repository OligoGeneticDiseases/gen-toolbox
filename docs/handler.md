# CommandHandler and hail_methods.py Overview

This document provides an overview of the `CommandHandler` class and the `hail_methods.py` module, which together streamline the processing of VCF files in the genomics toolbox.

## CommandHandler

The `CommandHandler` class is responsible for handling the different commands used in the genomics toolbox. The main function of this class is to process VCF files more efficiently by following a logical workflow.

### Methods

- `__init__(self, args)`: Initializes the CommandHandler object with the provided command-line arguments.
- `handle_find_type_command(self)`: Handles the `find_type` command.
- `handle_read_vcfs_command(self)`: Handles the `read_vcfs` command. This method has been updated to follow the new logic for processing VCF files more efficiently.
- `handle_load_db_command(self)`: Handles the `load_db` command. 

## hail_methods.py

The `hail_methods.py` module contains a set of functions that implement the core logic for processing VCF files using the Hail library. These functions are used by the `CommandHandler` class to perform various tasks in the VCF processing pipeline.

### Functions

- `import_and_annotate_vcf(vcf_path)`: Imports a VCF file and annotates it using `hail.import_vcf()` and `hail.VEP()`. Takes a VCF file path as input and returns an annotated MatrixTable.
- `merge_matrix_tables(matrix_tables)`: Merges a list of MatrixTables into a single MatrixTable with the same header and the same set of samples. Takes a list of MatrixTables as input and returns the merged MatrixTable.
- `reduce_to_2d_table(mt)`: Reduces a MatrixTable to a 2D MatrixTable with gene and frequency as keys. Takes a MatrixTable as input and returns the reduced MatrixTable.
- `create_frequency_bins(mt, num_bins=16)`: Creates a frequency table with the specified number of bins (default: 16). Takes a MatrixTable and an optional number of bins as input and returns the frequency table as a MatrixTable.

By using the `CommandHandler` class in conjunction with the functions in the `hail_methods.py` module, the genomics toolbox can efficiently process VCF files and perform various tasks, such as merging MatrixTables and creating frequency tables, with improved performance and maintainability.
