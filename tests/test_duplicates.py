import pytest

from rdkit import Chem
from rdkit.Chem import RegistrationHash

from common_code import load_templates


@pytest.mark.xfail(reason='we know we have some duplicates')
def test_check_duplicates():
    """
    Check if there are any duplicate templates.

    This test is slow, since it has to look at all the templates and
    generate hashes, so we want to run it only if we know all templates
    pass the basic requisites.
    """
    all_templates = {}
    duplicates = 0
    for i, _, cxsmiles in load_templates():
        mol = Chem.MolFromSmiles(cxsmiles)
        mol_layers = RegistrationHash.GetMolLayers(mol,
                                                   enable_tautomer_hash_v2=True)
        mol_hash = RegistrationHash.GetMolHash(mol_layers)
        if (seen_idx := all_templates.get(mol_hash, None)) is not None:
            print(f'Template #{i} is a duplicate of line #{seen_idx}!')
            duplicates += 1
        else:
            all_templates[mol_hash] = i
    assert duplicates == 0
