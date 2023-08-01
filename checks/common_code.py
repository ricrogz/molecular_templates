def load_templates():
    with open("templates.smi") as f:
        for i, line in enumerate(f, 1):
            cxsmiles = line.strip()
            if not cxsmiles:
                continue
            smiles = cxsmiles.split('|', 1)[0]
            yield i, smiles, cxsmiles
