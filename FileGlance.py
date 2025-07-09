# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# FileGlance.py
# class to glance at the file you are working with

import gzip

class FileGlance:
    def __init__(self, filepath, numlines=10):
        self.filepath = filepath
        self.numlines = numlines
        self.openfunc = self.getopenfunc()
        self.glance()

    def getopenfunc(self):
        if self.filepath.endswith('.gz'):
            return gzip.open
        else:
            return open

    def glance(self):
        """Preview the file contents."""
        try:
            with self.openfunc(self.filepath, 'rt') as f:
                print(f"Previewing first {self.numlines} lines of: {self.filepath}")
                for i, line in enumerate(f):
                    if i >= self.numlines:
                        break
                    print(f"Line {i + 1}: {line.strip()}")
        except Exception as e:
            print(f"Error reading file: {e}")




