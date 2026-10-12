"""Shared paths for scripts; makes src/ importable."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
DATA = os.path.join(ROOT, "data")
NOTES = os.path.join(ROOT, "notes")

if SRC not in sys.path:
    sys.path.insert(0, SRC)
