import os
import git
import json
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

img_width = 400
img_height = 400


def get_new_templates(base_commit_hash):
    """
    Finds the new templates that have been added since the base_commit
    """
    template_file = 'templates.smi'

    repo = git.Repo('.')
    base_commit = repo.commit(base_commit_hash)
    tpl_file_diff = repo.git.diff(base_commit, repo.head, template_file)

    # Make sure we only generate images for things that have been added,
    # not those we moved around or changed by adding a new line at the end
    added = []
    removed = []

    # We want to skip the diff header, as we don't want to hit the line
    # that says '+++ b/templates.smi'
    changes = tpl_file_diff.split('\n')
    for line in changes[4:]:
        if line.startswith('-') and (cxsmiles := line[1:].strip()):
            removed.append(cxsmiles)
        elif line.startswith('+') and (cxsmiles := line[1:].strip()):
            added.append(cxsmiles)

    change_id = 0
    for cxsmiles in added:
        if cxsmiles not in removed:
            yield change_id, cxsmiles
            change_id += 1


def draw_mol(cxsmiles, idx, output_dir='.'):
    smiles = cxsmiles.split()[0]
    legend = f"#{idx} {smiles}"
    mol = Chem.MolFromSmiles(cxsmiles)
    if not mol.GetNumConformers():
        raise ValueError("SMILES must include coordinates")

    print(f'Creating PNG image #{idx} for SMILES "{smiles}"')

    fname = os.path.join(output_dir, f'{idx}.png')
    png = MolToImage(mol, size=(img_width, img_height), legend=legend)
    png.save(fname)
    return fname


def export_generated_imgs(fnames):
    if gh_output := os.environ.get('GITHUB_OUTPUT', ''):
        json_list = json.dumps(fnames)
        with open(gh_output, 'a') as f:
            f.write(f'imgs={json_list}')


def main(args):
    if len(args) < 2:
        print(f'Usage {__file__} [base commit hash]')
        exit(1)

    print(f'Differencing with sha {args[1]}')

    count = 0
    exported_files = []
    tmpdir = tempfile.mkdtemp()
    for idx, new_template in get_new_templates(args[1]):
        fname = draw_mol(new_template, idx, tmpdir)
        exported_files.append(fname)
        count += 1

    print(f'Wrote {count} images to path {tmpdir}')
    export_generated_imgs(exported_files)


if __name__ == '__main__':
    main(sys.argv)
