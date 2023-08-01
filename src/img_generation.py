import os
import tempfile

import git
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from imgur_upload import upload_img

img_width = 400
img_height = 400
templates_file = 'templates.smi'


def get_new_templates(base_commit_hash):
    """
    Finds the new templates that have been added since the base_commit
    """

    with open(templates_file) as f:
        line_nums = {line.strip(): i for i, line in enumerate(f, 1)}

    repo = git.Repo('.')
    base_commit = repo.commit(base_commit_hash)
    tpl_file_diff = repo.git.diff(base_commit, repo.head, templates_file)

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

    for cxsmiles in added:
        if cxsmiles not in removed:
            yield line_nums[cxsmiles], cxsmiles


def draw_mol(cxsmiles, idx, output_dir='.'):
    smiles = cxsmiles.split()[0]
    legend = f"#{idx} {smiles}"
    mol = Chem.MolFromSmiles(cxsmiles)
    if not mol.GetNumConformers():
        raise ValueError("SMILES must include coordinates")

    print(f'Creating PNG image #{idx} for SMILES "{smiles}"')

    fpath = os.path.join(output_dir, f'{idx}.png')
    png = MolToImage(mol, size=(img_width, img_height), legend=legend)
    png.save(fpath)
    return fpath, legend


def export_image_urls(template_imgs):
    if gh_output := os.environ.get('GITHUB_OUTPUT', ''):
        markdown_imgs = [f'![{title}]({url})' for url, title in template_imgs]
        with open(gh_output, 'a') as f:
            f.write(f"template_imgs=\"{''.join(markdown_imgs)}\"")


def main():
    base_sha = os.environ.get('GH_PR_BASE', None)
    if not base_sha:
        print(
            "The base commit hash must be provided via the GH_PR_BASE env var")
        exit(1)
    print(f'Differencing with sha {base_sha}')

    template_imgs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, cxsmiles in get_new_templates(base_sha):
            fpath, title = draw_mol(cxsmiles, idx, tmpdir)
            img_url, title = upload_img(fpath, title)
            template_imgs.append((img_url, title))

    export_image_urls(template_imgs)
    print(f'{len(template_imgs)} images generated and uploaded')


if __name__ == '__main__':
    main()
