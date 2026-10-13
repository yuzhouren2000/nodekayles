"""Check the product-pairing theorem against enumeration and simulation.

Run: python3 tests/test_pairing_classification.py
The exhaustive branch bypasses the theorem; simulate() checks every first
move with a separately implemented neighbourhood calculation.
"""

import sys
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from pairing_certificates import construct_product_certificate, involutions, search
from verify_certificates import simulate


def main():
    boards = 0
    for m in range(1, 7):
        for n in range(m, 7):
            expected = construct_product_certificate(m, n)
            enumerated = search(m, n, exhaustive=True)
            simulated = {
                (tau, sigma)
                for tau, sigma in product(involutions(m), involutions(n))
                if simulate(m, n, tau, sigma)
            }
            assert set(enumerated) == simulated, (m, n)
            assert bool(enumerated) == (expected is not None), (m, n)
            boards += 1
    print(f"classification agrees with exhaustive search on {boards} boards")

    # Paths are the obstruction in the proof, including losing paths such as
    # P8 which nevertheless have no certificate of this form.
    for n in range(1, 11):
        assert bool(search(1, n, exhaustive=True)) == (n == 4), n
    print("path lemma agrees with enumeration for n=1..10")

    checked = 0
    for m in range(1, 15):
        for n in range(1, 15):
            cert = construct_product_certificate(m, n)
            if cert is not None:
                assert simulate(m, n, *cert), (m, n, cert)
                checked += 1
    # Long boards check the construction independently of the SG solver's
    # bitboard size and the finite range of the original enumeration.
    for m, n in [(4, 101), (101, 4), (6, 100)]:
        cert = construct_product_certificate(m, n)
        assert cert is not None and simulate(m, n, *cert), (m, n)
        checked += 1
    print(f"{checked} constructed strategies pass independent simulation")

    for m, n in [(5, 8), (5, 12), (7, 8), (9, 12), (3, 22), (3, 26)]:
        assert search(m, n, limit=1) == [], (m, n)
    print("mixed-parity boards outside the classified families have no certificate")


if __name__ == "__main__":
    main()
